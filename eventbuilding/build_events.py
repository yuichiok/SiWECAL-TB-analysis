#!/usr/bin/env python
from __future__ import print_function

import collections

import numpy as np
import ROOT

from bcid_handling import BCIDHandler
from ecal_config import EcalConfig
from parse_config import create_cli_from_default_config
from util import aligned_path, EventBuildingException, speed_warning_if_python2
from util import get_tree_spills, get_tree_spills_no_progress_info


def get_hits_per_event(entry, bcid_handler, ecal_config, zero_suppress=True):
    event = {}
    bcid_merge_end = collections.defaultdict(int)
    n_chips = ecal_config._n_chips
    n_channels = ecal_config._n_channels
    n_scas = ecal_config._n_scas
    channel_arange = np.arange(n_channels)

    # It is faster to fill the arrays once outside the loop.
    all_gain_hit = np.array(entry.hitbit_high)
    all_adc_high = np.array(entry.adc_high)
    all_adc_low = np.array(entry.adc_low)

    # assert bcid_handler.merged_bcid.shape == np.array(event.bcid).shape
    for i_chip, scas, bcid, bcid_before_merge in bcid_handler.yield_from_spill():
        if len(scas) == 1:
            bcid_channel_id = (i_chip * n_scas + scas[0]) * n_channels + channel_arange
        else:
            all_scas_bcid_channel_id = (
                np.repeat(i_chip * n_scas + np.array(scas), n_channels)
                * n_channels
            ).reshape(-1, n_channels) + channel_arange
            hitbit_high = all_gain_hit[all_scas_bcid_channel_id]
            # Pick the first trigger of each channel within the BCID window.
            bcid_channel_id = all_scas_bcid_channel_id[
                np.argmax(hitbit_high, axis=0), channel_arange
            ]
        hitbit_high = all_gain_hit[bcid_channel_id]
        if np.all(hitbit_high < 0):
            raise EventBuildingException(
                "Should have been caught in BCIDHandler",
                i_chip, scas, bcid, bcid_before_merge, hitbit_high
            )
        # hitbit_high == hitbit_low (up to tiny errors), confirmed by Stephane Callier.
        is_hit = np.array(hitbit_high > 0)
        if zero_suppress:
            if is_hit.sum() == 0:
                raise EventBuildingException(
                    "Should have been caught in BCIDHandler",
                    i_chip, scas, bcid, bcid_before_merge, hitbit_high
                )

            def ext(name, array):
                return event[bcid][name].extend(array[is_hit])
        else:
            def ext(name, array):
                return event[bcid][name].extend(array)

        if bcid not in event:
            event[bcid] = collections.defaultdict(list)
        slab = entry.slboard_id[(i_chip // n_chips)]
        slab_id = ecal_config._slabs.index(slab)
        chip = entry.chipid[i_chip]
        sca = (bcid_channel_id // n_channels)  % n_scas
        per_sca_picker = (slab_id, chip, channel_arange, sca)

        ext("isHit", is_hit)
        ext("isMasked", np.array(ecal_config.masked_map[slab_id, chip], dtype=int))
        ext("isCommissioned", np.array(ecal_config.is_commissioned_map[per_sca_picker], dtype=int))
        ext("x", np.array(ecal_config.x[slab_id, chip]))
        ext("y", np.array(ecal_config.y[slab_id, chip]))
        ext("z", np.array(ecal_config.z[slab_id, chip]))

        adc_high = all_adc_high[bcid_channel_id]
        energy = adc_high.astype(float)
        energy -= ecal_config.pedestal_map[per_sca_picker]
        energy /= ecal_config.mip_map[slab_id, chip]
        ext("adc_high", adc_high)
        ext("energy", energy)
        adc_low = all_adc_low[bcid_channel_id]
        if not ecal_config.no_lg:
            energy_lg = adc_low.astype(float)
            energy_lg -= ecal_config.pedestal_lg_map[per_sca_picker]
            energy_lg /= ecal_config.mip_lg_map[slab_id, chip]
            ext("energy_lg", energy_lg)
        ext("adc_low", adc_low)

        n_in_batch = n_channels
        ext("slab", np.full(n_in_batch, slab, dtype=int))
        ext("chip", np.full(n_in_batch, chip, dtype=int))
        ext("sca", np.full(n_in_batch, sca, dtype=int))
        ext("chan", np.arange(n_channels, dtype=int))
        n_scas_filled = entry.nColumns[i_chip]
        ext("n_scas_filled", np.full(n_in_batch, n_scas_filled, dtype=int))

        bcid_merge_end[bcid] = max(bcid_before_merge, bcid_merge_end[bcid])

    for bcid, event_arrays in event.items():
        for key in event_arrays:
            event_arrays[key] = np.array(event_arrays[key])
        event[bcid] = event_arrays
    return event, bcid_merge_end


def _get_hit_branches(out_arrays):
    """Get the names of the branches that are filled per-hit (through EcalHit).

    At the same time, this function checks that EcalHit and the branches
    in BuildEvents actually agree on what these branches in should be.
    This would not be necessary as a runtime check, but is fast and should be
    a good hint in case of erroneous code changes.
    """
    hit_branches = {br[4:] for br in out_arrays if br.startswith("hit_")}
    return hit_branches


class BuildEvents:
    _in_tree_name = "siwecaldecoded"
    _out_tree_name = "ecal"

    _branch_tags = [
        # "event info":
            "event/I",
            "spill/I",
            "cycle/I",
            "bcid/I",
            "bcid_first_sca_full/I",
            "bcid_merge_end/I",
            "id_run/I",
            "id_dat/I",
        # "hit summary":
            "nhit_slab/I",
            "nhit_chip/I",
            "nhit_chan/I",
            "nhit_len/I",
            "sum_energy/F",
            "sum_energy_lg/F",
        # 64: the number of channels on a chip.
        # "hit id":
            "hit_slab[nhit_len]/I",
            "hit_chip[nhit_len]/I",
            "hit_chan[nhit_len]/I",
            "hit_sca[nhit_len]/I",
        # "hit coord":
            "hit_x[nhit_len]/F",
            "hit_y[nhit_len]/F",
            "hit_z[nhit_len]/F",
        # "hit readout":
            "hit_adc_high[nhit_len]/I",
            "hit_adc_low[nhit_len]/I",
            "hit_energy[nhit_len]/F",
            "hit_energy_lg[nhit_len]/F",
            "hit_n_scas_filled[nhit_len]/I",
        # "hit booleans":
            "hit_isHit[nhit_len]/I",
            "hit_isMasked[nhit_len]/I",
            "hit_isCommissioned[nhit_len]/I",
    ]

    _branch_tags_lg_only = ["sum_energy_lg/F", "hit_energy_lg[nhit_len]/F"]
    # Fast checks to avoid introducing errors.
    assert set(_branch_tags_lg_only).issubset(_branch_tags)


    def __init__(
        self,
        config_parser,
    ):
        self._config_parser = config_parser
        eb_config = config_parser["eventbuilding"]
        self.file_name = eb_config["converted_path"]
        self.w_config = eb_config["w_config"]
        self.max_entries = int(eb_config["max_entries"])
        self.out_file_name = eb_config["build_path"]
        self.min_slabs_hit = int(eb_config["min_slabs_hit"])
        self._zero_suppress = eb_config.getboolean("zero_suppress")
        self._merge_within_chip = eb_config.getboolean("merge_within_chip")
        if self._merge_within_chip and not self._zero_suppress:
            raise EventBuildingException(
                "`merge_within_chip` requires `zero_suppress`."
            )

        self._previous_cycle = -1
        self.event_counter = 0

        self.redo_config = eb_config.getboolean("redo_config")
        self.in_tree = self._get_tree(self.file_name)
        slabs = self._get_slabs(self.in_tree)

        self._no_progress_info = eb_config.getboolean("no_progress_info")
        self._id_dat = int(eb_config["id_dat"])
        self._id_run = int(eb_config["id_run"])
        speed_warning_if_python2()
        self.ecal_config = EcalConfig(
            calibration_files=eb_config,
            slabs=slabs,
            asu_versions=eb_config["asu_versions"],
            geometry_config = self._config_parser["geometry"],
            commissioning_config = self._config_parser["commissioning"],
            no_lg=eb_config.getboolean("no_lg"),
        )

    def _get_tree(self, file_name):
        self.in_file = ROOT.TFile(file_name,"read")
        if self.redo_config:
            tree_name = self._out_tree_name
        else:
            tree_name = self._in_tree_name
        self.tree = self.in_file.Get(tree_name)
        if not self.tree:
            trees_available = [k.GetName() for k in self.in_file.GetListOfKeys()]
            print("Found tree names:", trees_available)
            ex_txt = "Tree %s not found in %s" %(tree_name, file_name)
            raise EventBuildingException(ex_txt)
        return self.tree


    def _get_slabs(self, tree):
        if self.redo_config:
            return self._get_slabs_from_buildfile()
        tree.Draw("slboard_id >> slot_hist", "", "goff")
        hist = ROOT.gDirectory.Get("slot_hist")
        slabs = []
        for i in range(1, hist.GetNbinsX() + 1):  # 0 is underflow bin.
            if hist.GetBinContent(i) > 0:
                slabs.append(int(np.ceil(hist.GetBinLowEdge(i))))
        if -1 in slabs:
            slabs.remove(-1)  # Used to indicate dummy entries.
        assert len(set(slabs)) == len(slabs)
        return slabs


    def _get_slabs_from_buildfile(self):
        hist = self.in_file.Get("config").Get("w_in_front")
        slabs = []
        for i in range(1, hist.GetNbinsX() + 1):  # 0 is underflow bin.
            slabs.append(int(np.ceil(hist.GetBinLowEdge(i))))
        return slabs


    def _add_branch(self, tag):
        name, branch_type = tag.split("/")
        branch_type = {"I": "i", "F": "f"}[branch_type]
        array_indicator = "[nhit_len]"
        if array_indicator in name:
            name = name.replace(array_indicator, "")
            starting_value = 10000 * [0]
        else:
            starting_value = [0]
        self.out_arrays[name] = np.array(starting_value, dtype=branch_type)
        self.out_tree.Branch(name, self.out_arrays[name], tag)


    def _create_out_tree(self, out_file_name, file_name):
        if out_file_name is None:
            if file_name.endswith("_converted.root"):
                out_file_name = file_name[:-len("_converted.root")] + "_build.root"
            elif file_name.endswith(".root"):
                out_file_name = file_name[:-len(".root")] + "_build.root"
            else:
                raise EventBuildingException("Unexpected file extension: %s" %file_name)
        print(aligned_path("# Creating ecal tree in file ", out_file_name))
        self.out_file = ROOT.TFile(out_file_name,"recreate")
        self.out_tree = ROOT.TTree(self._out_tree_name, "Build ecal events")

        for branch_tag in self._branch_tags:
            if self.ecal_config.no_lg:
                if branch_tag in self._branch_tags_lg_only:
                    continue
            self._add_branch(branch_tag)
        self._hit_branches = _get_hit_branches(self.out_arrays)
        self._other_branches = sorted(set(self.out_arrays.keys()) - set(["hit_" + x for x in self._hit_branches]))
        return self.out_tree


    def _write_and_close(self):
        self.out_tree.Write()
        print("# Created tree with %i events." % self.out_tree.GetEntries())
        self.out_file.Close()
        self.in_file.Close()


    def build_events(self, file_name=None, max_entries=None, out_file_name=None):
        if file_name is None:
            file_name = self.file_name
        if out_file_name is None:
            out_file_name = self.out_file_name
        self.in_tree = self._get_tree(file_name)
        self.out_arrays = {}
        self.out_tree = self._create_out_tree(out_file_name, file_name)
        self._create_config_hists()

        if max_entries is None:
            max_entries = self.max_entries
        if max_entries < 0:
            max_entries = self.in_tree.GetEntries()

        if self._no_progress_info:
            _get_tree_spills = get_tree_spills_no_progress_info
        else:
            _get_tree_spills = get_tree_spills

        if self.redo_config:
            for _, event in _get_tree_spills(self.in_tree, max_entries):
                self._redo_fill_event(event)
        else:
            for i_spill, entry in _get_tree_spills(self.in_tree, max_entries):
                self._fill_spill(i_spill, entry)
        self._write_and_close()


    def _create_config_hists(self):
        self.out_file.mkdir("config")
        self.out_file.GetDirectory("config").cd()
        self._fill_w_config_hist()
        self._fill_config_maps()
        self.out_file.cd()


    def _fill_w_config_hist(self):
        slabs = self.ecal_config._slabs
        bin_centers = np.array(slabs, dtype="float64")
        bin_edge_candidates = np.concatenate([bin_centers - 0.5, bin_centers + 0.5])
        bin_edges = np.sort(np.unique(bin_edge_candidates))

        w_hist = ROOT.TH1F("w_in_front","w_in_front",len(bin_edges) - 1, bin_edges)
        if "," in self.w_config:
            abs_thick = np.array(list(map(float, self.w_config.split(","))))
        else:
            abs_thick = np.full_like(slabs, float(self.w_config))
        assert len(slabs) == len(abs_thick)
        for i in range(len(slabs)):
            w_hist.Fill(slabs[i], abs_thick[i])
        w_hist.Write()
        print("# W thickness [mm]:", abs_thick)
        # w_x0 = 1 / 3.5  # 0.56 #Xo per mm of W.
        # pos_x0 = np.cumsum(abs_thick) * w_x0
        # print("# -> X0:", pos_x0)

    def _fill_config_maps(self):
        def is_config_map(name):
            return not name.startswith("_") and "_map" in name

        for name in dir(self.ecal_config):
            if is_config_map(name):
                config_map = getattr(self.ecal_config, name)
                self._fill_config_map(config_map, name)


    def _fill_config_map(self, config_map, name):
        # There might be a faster implementation than this for-loop.
        # But at it is only run once (per config map), this should be ok.
        if len(config_map.shape) == 4:
            # There is no TH4. Instead, build a TH3 for each sca.
            assert config_map.shape[3] == self.ecal_config._n_scas
            for i in range(config_map.shape[3]):
                new_name = name + "_sca{:02d}".format(i)
                self._fill_config_map(config_map[:,:,:, i], new_name)
            return
        if len(config_map.shape) != 3:
            print("Warning: Storing of this config map is not implemented:", name)
        assert config_map.shape[0] == self.ecal_config._n_slabs
        assert config_map.shape[1] == self.ecal_config._n_chips
        assert config_map.shape[2] == self.ecal_config._n_channels
        if config_map.dtype == int or config_map.dtype == bool:
            hist_fct = ROOT.TH3I
        else:
            hist_fct = ROOT.TH3F
        hist = hist_fct(
            name, name,
            self.ecal_config._n_slabs, np.arange(-0.5, self.ecal_config._n_slabs),
            self.ecal_config._n_chips, np.arange(-0.5, self.ecal_config._n_chips),
            self.ecal_config._n_channels, np.arange(-0.5, self.ecal_config._n_channels),
        )
        for i_slab in range(self.ecal_config._n_slabs):
            for i_chip in range(self.ecal_config._n_chips):
                for i_channel in range(self.ecal_config._n_channels):
                    hist.Fill(i_slab, i_chip, i_channel, config_map[i_slab, i_chip, i_channel])
        hist.Write()


    def _redo_fill_event(self, event):
        if event.nhit_slab < self.min_slabs_hit:
            return
        for branch_name in self._hit_branches:
            branch_name = "hit_"+ branch_name
            if branch_name == "hit_energy":
                mip = self.ecal_config.mip_map[event.hit_slab, event.hit_chip, event.hit_chan]
                pedestal = self.ecal_config.pedestal_map[event.hit_slab, event.hit_chip, event.hit_chan, event.hit_sca]
                arr = (event.hit_adc_high - pedestal) / mip
            elif branch_name == "hit_energy_lg":
                mip = self.ecal_config.mip_lg_map[event.hit_slab, event.hit_chip, event.hit_chan]
                pedestal = self.ecal_config.pedestal_lg_map[event.hit_slab, event.hit_chip, event.hit_chan, event.hit_sca]
                arr = (event.hit_adc_low - pedestal) / mip
            elif branch_name == "hit_isCommissioned":
                arr = self.ecal_config.is_commissioned_map[event.hit_slab, event.hit_chip, event.hit_chan, event.hit_sca]
            elif branch_name == "hit_isMasked":
                arr = self.ecal_config.masked_map[event.hit_slab, event.hit_chip, event.hit_chan]
            elif branch_name == "hit_x":
                arr = self.ecal_config.x[event.hit_slab, event.hit_chip, event.hit_chan]
            elif branch_name == "hit_y":
                arr = self.ecal_config.y[event.hit_slab, event.hit_chip, event.hit_chan]
            elif branch_name == "hit_z":
                arr = self.ecal_config.z[event.hit_slab, event.hit_chip, event.hit_chan]
            else:
                arr = getattr(event, branch_name)
            self.out_arrays[branch_name][:len(arr)] = arr
        for branch_name in self._other_branches:
            if branch_name == "sum_energy":
                val = np.sum(self.out_arrays["hit_energy"][:event.nhit_len][np.array(event.hit_isHit, dtype=bool)])
            elif branch_name == "sum_energy_lg":
                val = np.sum(self.out_arrays["hit_energy_lg"][:event.nhit_len][np.array(event.hit_isHit, dtype=bool)])
            elif branch_name == "id_dat" and self._id_dat != -1:
                val = self._id_dat
            elif branch_name == "id_run" and self._id_run != -1:
                val = self._id_run
            else:
                val = getattr(event, branch_name)
            self.out_arrays[branch_name][0] = val
        self.out_tree.Fill()


    def _fill_spill(self, spill, entry):
        b = self.out_arrays
        b["spill"][0] = spill
        bcid_handler = BCIDHandler(
            self._config_parser["bcid"],
            self._config_parser["geometry"],
            self.min_slabs_hit,
            self._merge_within_chip,
        )
        bcid_handler.load_spill(entry)
        hits_per_event, bcid_merge_end = get_hits_per_event(
            entry, bcid_handler, self.ecal_config, self._zero_suppress
        )

        max_hits = int(self._config_parser["eventbuilding"]["max_hits_per_event"])
        for bcid in sorted(hits_per_event):
            hits = hits_per_event[bcid]
            if len(hits) == 0:
                raise EventBuildingException("Event with 0 hits. Why?")
            tmp_for_unique = hits["slab"]
            nhit_slab = np.unique(tmp_for_unique).size
            if nhit_slab < self.min_slabs_hit:
                # This check is not redundant: here we filter out cases were all
                # channels in a chip/slab had isHit == False.
                continue
            cycle = entry.acqNumber
            if self._previous_cycle < cycle:
                self.event_counter = 0
                self._previous_cycle = cycle
            else:
                self.event_counter += 1
            b["cycle"][0] = cycle
            b["event"][0] = self.event_counter
            b["bcid"][0] = bcid
            b["bcid_first_sca_full"][0] = bcid_handler.bcid_first_sca_full
            b["bcid_merge_end"][0] = bcid_merge_end[bcid]
            b["id_run"][0] = self._id_run
            b["id_dat"][0] = self._id_dat

            # count hits per slab/chan/chip
            b["nhit_slab"][0] = nhit_slab
            tmp_for_unique = tmp_for_unique * self.ecal_config._n_chips + hits["chip"]
            b["nhit_chip"][0] = np.unique(tmp_for_unique).size
            b["nhit_len"][0] = hits["isHit"].size
            tmp_for_unique = tmp_for_unique * self.ecal_config._n_channels + hits["chan"]
            b["nhit_chan"][0] = np.unique(tmp_for_unique[hits["isHit"]]).size
            b["sum_energy"][0] = hits["energy"][hits["isHit"]].sum()
            if "energy_lg" in hits:
                b["sum_energy_lg"][0] = hits["energy_lg"][hits["isHit"]].sum()

            if max_hits > 0 and len(hits) > max_hits:
                txt = "Suspicious number of hits! %i " %len(hits)
                txt += "for bcid %i." %b["bcid"][0]
                txt += "Skipping event %i." % b["event"][0]
                print(txt)
                continue

            for hit_field in self._hit_branches:
                ev_val = hits[hit_field]
                b["hit_" + hit_field][:ev_val.size] = ev_val
            self.out_tree.Fill()


if __name__ == "__main__":
    parser = create_cli_from_default_config()
    BuildEvents(parser).build_events()
