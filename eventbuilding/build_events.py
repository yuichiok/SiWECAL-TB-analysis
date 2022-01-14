#!/usr/bin/env python
from __future__ import print_function

import argparse
import collections
import sys

import numpy as np
import ROOT as rt
from help_tools import *


try:
    from tqdm.autonotebook import tqdm

    def get_tree_spills(tree, max_entries):
        for i, spill in enumerate(tqdm(tree, desc="# Build events", total=max_entries, unit=" spills")):
            if i > max_entries:
                break
            yield i, spill
except ImportError:
    from datetime import datetime

    def get_tree_spills(tree, max_entries):
        print("# Going to analyze %i entries..." %max_entries)
        print("# For better progress information: `pip install tqdm`.")
        print("# Start time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        progress_bar = ""
        for i, spill in enumerate(tree):
            if i > max_entries:
                break
            if i%10 == 0 or i == max_entries:
                progress_bar = "#" * int(30 * i / max_entries)
                print("# Build events [{}] Spill {}/{}".format(
                        progress_bar.ljust(30),
                        str(i).rjust(len(str(max_entries))),
                        max_entries,
                    ),
                    end="\r",
                )
                sys.stdout.flush()
            yield i, spill
        print("\n# Final time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


class BCIDHandler:
    def __init__(self, bcid_numbers, min_slabs_hit=1):
        self._N = bcid_numbers
        self._min_slabs_hit = min_slabs_hit


    def _apply_overflow_correction(self, arr):
        """"Within a chips memory, the BCID must increase. A drop indicates overflow."""
        # assert np.max(arr) < self._N.bcid_overflow
        n_memory_cycles = np.cumsum(arr[:,1:] - arr[:,:-1] < 0, axis=1)
        is_filled = arr[:,1:] != self._N.bcid_bad_value
        arr[:,1:] = arr[:,1:] + self._N.bcid_overflow * n_memory_cycles * is_filled
        return arr


    def _is_retrigger(self, arr):
        """Within each chip's memory, look for consecutive BCIDs."""
        delta = arr[:,1:] - arr[:,:-1]
        return np.logical_and(delta > 0, delta <= self._N.bcid_drop_retrigger_delta)


    def _find_bcid_filling_last_sca(self, arr):
        bcid_in_last_sca = arr[:,-1]
        bcid_in_last_sca = bcid_in_last_sca[bcid_in_last_sca != self._N.bcid_bad_value]
        if len(bcid_in_last_sca):
            bcid_in_last_sca = bcid_in_last_sca.min()
        else:
            bcid_in_last_sca = self._N.bcid_bad_value
        return bcid_in_last_sca


    def _get_bcids(self, raw_bcid_array):
        raw_bcids = np.array(list(raw_bcid_array)).reshape(-1, self._N.n_scas)
        chip_ids = np.nonzero(
            np.any(raw_bcids != self._N.bcid_bad_value, axis=1)
        )[0]
        bcids = np.copy(raw_bcids[chip_ids])
        bcids = self._apply_overflow_correction(bcids)
        bcids[:,1:][self._is_retrigger(bcids)] = self._N.bcid_bad_value
        bcids[bcids < self._N.bcid_skip_noisy_acquisition_start] = self._N.bcid_bad_value
        for drop_bcid in self._N.bcid_drop:
            bcids[bcids == drop_bcid] = self._N.bcid_bad_value
        return chip_ids, bcids


    def _get_bcid_start_stop(self, arr):
        all_bcids = np.unique(arr)
        all_bcids = all_bcids[all_bcids != self._N.bcid_bad_value]
        if len(all_bcids) == 0:
            return dict()
        previous_bcid = all_bcids[0]
        bcid_start_stop = {previous_bcid: previous_bcid}
        for current_bcid in all_bcids[1:]:
            if current_bcid - bcid_start_stop[previous_bcid] < self._N.bcid_merge_delta:
                bcid_start_stop[previous_bcid] = current_bcid
            else:
                bcid_start_stop[current_bcid] = current_bcid
                previous_bcid = current_bcid
        return bcid_start_stop


    def _get_array_positions_per_bcid(self, chip_ids, bcids):
        bcid_start_stop = self._get_bcid_start_stop(bcids)
        _pos_per_bcid = {k: [] for k in bcid_start_stop}
        bcid_before_merge = {k: [] for k in bcid_start_stop}
        for bcid_start in bcid_start_stop:
            i_chip_local, i_sca = np.nonzero(np.logical_and(
                bcids >= bcid_start, bcids <= bcid_start_stop[bcid_start]
            ))
            bcid_before_merge[bcid_start].extend(bcids[i_chip_local, i_sca])
            i_chip = chip_ids[i_chip_local]
            _pos_per_bcid[bcid_start].extend(i_chip * self._N.n_scas + i_sca)

        pos_per_bcid = {}
        for k, v in _pos_per_bcid.items():
            pos = np.array(v)
            slabs_hit = pos // (self._N.n_chips * self._N.n_scas)
            if len(np.unique(slabs_hit)) >= self._min_slabs_hit:
                pos_per_bcid[k] = pos
        return pos_per_bcid, bcid_before_merge


    def load_spill(self, spill_entry):
        """False if the current spill is empty."""
        chip_ids, bcids = self._get_bcids(spill_entry.bcid)
        self.bcid_first_sca_full = self._find_bcid_filling_last_sca(bcids)

        self.pos_per_bcid, self._bcid_before_merge = \
            self._get_array_positions_per_bcid(chip_ids, bcids)
        self.spill_bcids = sorted(self.pos_per_bcid.keys())
        return len(self.pos_per_bcid) > 0


    def yield_from_spill(self):
        for bcid in self.spill_bcids:
            before_merge = self._bcid_before_merge[bcid]
            for i, pos in enumerate(self.pos_per_bcid[bcid]):
                yield pos, bcid, before_merge[i]


    def previous_bcid(self, bcid):
        i_bcid = self.spill_bcids.index(bcid)
        if i_bcid == 0:
            return -1
        else:
            return self.spill_bcids[i_bcid - 1]

    def next_bcid(self, bcid):
        i_bcid = self.spill_bcids.index(bcid)
        if i_bcid == len(self.spill_bcids) - 1:
            return -1
        else:
            return self.spill_bcids[i_bcid + 1]


def get_hits_per_event(entry, bcid_handler, ecal_config):
    event = {}
    bcid_merge_end = collections.defaultdict(int)
    n_chips = ecal_config._N.n_chips
    n_channels = ecal_config._N.n_channels
    n_scas = ecal_config._N.n_scas

    # It is faster to fill the arrays once outside the loop.
    all_gain_hit = np.array(entry.gain_hit_high)
    all_hg = np.array(entry.charge_hiGain)
    all_lg = np.array(entry.charge_lowGain)

    # assert bcid_handler.merged_bcid.shape == np.array(event.bcid).shape
    for i, bcid, bcid_before_merge in bcid_handler.yield_from_spill():
        bcid_channel_id = slice(i * n_channels, (i + 1) * n_channels)
        gain_hit_high = all_gain_hit[bcid_channel_id]
        if np.all(gain_hit_high < 0):
            continue
        # gain_hit_high == gain_hit_low (up to tiny errors), confirmed by Stephane Callier.
        is_hit = np.array(gain_hit_high > 0)
        if ecal_config.zero_suppress:
            if is_hit.sum() == 0:
                # It is not clear to me why chips with no hit-bit for any of their channels
                # are ever written to memory, but the DAQ writes such lines.
                continue

            def ext(name, array):
                return event[bcid][name].extend(array[is_hit])
        else:
            def ext(name, array):
                return event[bcid][name].extend(array)

        if bcid not in event:
            event[bcid] = collections.defaultdict(list)
        slab = entry.slot[(i // n_scas // n_chips)]
        slab_id = ecal_config._N.slabs.index(slab)
        chip = entry.chipid[i // n_scas]
        sca = i % n_scas

        ext("isHit", is_hit)
        ext("isMasked", np.array(ecal_config.masked_map[slab_id, chip], dtype=int))
        ext("isCommissioned", np.array(ecal_config.is_commissioned_map[slab_id, chip, :, sca], dtype=int))
        ext("x", np.array(ecal_config.x[slab_id, chip]))
        ext("y", np.array(ecal_config.y[slab_id, chip]))
        ext("z", np.array(ecal_config.z[slab_id, chip]))

        hg = all_hg[bcid_channel_id]
        energy = hg.astype(float)
        energy -= ecal_config.pedestal_map[slab_id, chip, :, sca]
        energy /= ecal_config.mip_map[slab_id, chip]
        ext("hg", hg)
        ext("energy", energy)
        lg = all_lg[bcid_channel_id]
        if not ecal_config.no_lg:
            energy_lg = hg.astype(float)
            energy_lg -= ecal_config.pedestal_lg_map[slab_id, chip, :, sca]
            energy_lg /= ecal_config.mip_lg_map[slab_id, chip]
            ext("energy_lg", energy_lg)
        ext("lg", lg)

        n_in_batch = n_channels
        ext("slab", np.full(n_in_batch, slab, dtype=int))
        ext("chip", np.full(n_in_batch, chip, dtype=int))
        ext("sca", np.full(n_in_batch, sca, dtype=int))
        ext("chan", np.arange(n_channels, dtype=int))
        n_scas_filled = entry.nColumns[i // n_scas]
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
            "bcid/I",
            "bcid_first_sca_full/I",
            "bcid_merge_end/I",
            "bcid_prev/I",
            "bcid_next/I",
        # "hit summary":
            "nhit_slab/I",
            "nhit_chip/I",
            "nhit_chan/I",
            "nhit_len/I",
            "sum_hg/F",
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
            "hit_hg[nhit_len]/I",
            "hit_lg[nhit_len]/I",
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
        file_name,
        w_config=-1,
        max_entries=-1,
        out_file_name=None,
        commissioning_folder=None,
        min_slabs_hit=4,
        cob_positions_string="",
        ecal_numbers=None,  # Not provided in CLI. Mainly useful for debugging/changing.
        no_zero_suppress=False,
        no_lg=False,
        **config_file_kws
    ):
        self.file_name = file_name
        self.w_config = w_config
        self.max_entries = max_entries
        self.out_file_name = out_file_name
        self.min_slabs_hit = min_slabs_hit
        self.event_counter = 0

        self.in_tree = self._get_tree(file_name)
        slabs = self._get_slabs(self.in_tree)
        cob_slabs = set(map(int, filter(None, cob_positions_string.split(" "))))
        if ecal_numbers is None:
           ecal_numbers = EcalNumbers(slabs=slabs, cob_slabs=cob_slabs)
        else:
            assert ecal_numbers.slabs == slabs
            assert ecal_numbers.cob_slabs == cob_slabs

        speed_warning_if_python2()
        self.ecal_config = EcalConfig(
            commissioning_folder=commissioning_folder,
            numbers=ecal_numbers,
            no_lg=no_lg,
            zero_suppress=not bool(no_zero_suppress),
            **config_file_kws
        )

    def _get_tree(self, file_name):
        self.in_file = rt.TFile(file_name,"read")
        self.tree = self.in_file.Get(self._in_tree_name)
        if not self.tree:
            trees_available = [k.GetName() for k in self.in_file.GetListOfKeys()]
            print("Found tree names:", trees_available)
            ex_txt = "Tree %s not found in %s" %(self._in_tree_name, file_name)
            raise EventBuildingException(ex_txt)
        return self.tree


    def _get_slabs(self, tree):
        tree.Draw("slot >> slot_hist", "", "goff")
        hist = rt.gDirectory.Get("slot_hist")
        slabs = []
        for i in range(1, hist.GetNbinsX() + 1):  # 0 is undeflow bin.
            if hist.GetBinContent(i) > 0:
                slabs.append(int(np.ceil(hist.GetBinLowEdge(i))))
        if -1 in slabs:
            slabs.remove(-1)  # Used to indicate dummy entries.
        assert len(set(slabs)) == len(slabs)
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
        self.out_file = rt.TFile(out_file_name,"recreate")
        self.out_tree = rt.TTree(self._out_tree_name, "Build ecal events")

        for branch_tag in self._branch_tags:
            if self.ecal_config.no_lg:
                if branch_tag in self._branch_tags_lg_only:
                    continue
            self._add_branch(branch_tag)
        self._hit_branches = _get_hit_branches(self.out_arrays)
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

        for i_spill, entry in get_tree_spills(self.in_tree, max_entries):
            self._fill_spill(i_spill, entry)
        self._write_and_close()


    def _create_config_hists(self):
        self.out_file.mkdir("config")
        self.out_file.GetDirectory("config").cd()
        self._fill_w_config_hist()
        self._fill_config_maps()
        self.out_file.cd()


    def _fill_w_config_hist(self):
        slabs = self.ecal_config._N.slabs
        bin_centers = np.array(slabs, dtype="float64")
        bin_edge_candidates = np.concatenate([bin_centers - 0.5, bin_centers + 0.5])
        bin_edges = np.sort(np.unique(bin_edge_candidates))

        w_hist = rt.TH1F("w_in_front","w_in_front",len(bin_edges) - 1, bin_edges)
        if self.w_config in self.ecal_config._N.w_config_options.keys():
            abs_thick = self.ecal_config._N.w_config_options[self.w_config]
        elif self.w_config == 0:
            abs_thick = np.zeros_like(slabs)
        else:
            raise EventBuildingException("Not a valid W config:", self.w_config)
        assert len(slabs) == len(abs_thick)
        for i in range(len(slabs)):
            w_hist.Fill(slabs[i], abs_thick[i])
        # w_x0 = 1 / 3.5  # 0.56 #Xo per mm of W.
        # pos_x0 = np.cumsum(abs_thick) * w_x0
        w_hist.Write()
        print("W config %i used." %self.w_config, end=" ")
        print("Absolute thickness:", abs_thick)


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
            assert config_map.shape[3] == self.ecal_config._N.n_scas
            for i in range(config_map.shape[3]):
                new_name = name + "_sca{:02d}".format(i)
                self._fill_config_map(config_map[:,:,:, i], new_name)
            return
        if len(config_map.shape) != 3:
            print("Warning: Storing of this config map is not implemented:", name)
        assert config_map.shape[0] == self.ecal_config._N.n_slabs
        assert config_map.shape[1] == self.ecal_config._N.n_chips
        assert config_map.shape[2] == self.ecal_config._N.n_channels
        if config_map.dtype == int or config_map.dtype == bool:
            hist_fct = rt.TH3I
        else:
            hist_fct = rt.TH3F
        hist = hist_fct(
            name, name,
            self.ecal_config._N.n_slabs, np.arange(-0.5, self.ecal_config._N.n_slabs),
            self.ecal_config._N.n_chips, np.arange(-0.5, self.ecal_config._N.n_chips),
            self.ecal_config._N.n_channels, np.arange(-0.5, self.ecal_config._N.n_channels),
        )
        for i_slab in range(self.ecal_config._N.n_slabs):
            for i_chip in range(self.ecal_config._N.n_chips):
                for i_channel in range(self.ecal_config._N.n_channels):
                    hist.Fill(i_slab, i_chip, i_channel, config_map[i_slab, i_chip, i_channel])
        hist.Write()


    def _fill_spill(self, spill, entry):
        b = self.out_arrays
        b["spill"][0] = spill
        bcid_handler = BCIDHandler(self.ecal_config._N, self.min_slabs_hit)
        bcid_handler.load_spill(entry)
        hits_per_event, bcid_merge_end = get_hits_per_event(
            entry, bcid_handler, self.ecal_config
        )

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

            self.event_counter += 1
            b["event"][0] = self.event_counter
            b["bcid"][0] = bcid
            b["bcid_first_sca_full"][0] = bcid_handler.bcid_first_sca_full
            b["bcid_merge_end"][0] = bcid_merge_end[bcid]
            b["bcid_prev"][0] = bcid_handler.previous_bcid(bcid)
            b["bcid_next"][0] = bcid_handler.next_bcid(bcid)

            # count hits per slab/chan/chip
            b["nhit_slab"][0] = nhit_slab
            tmp_for_unique = tmp_for_unique * self.ecal_config._N.n_chips + hits["chip"]
            b["nhit_chip"][0] = np.unique(tmp_for_unique).size
            b["nhit_len"][0] = hits["isHit"].size
            tmp_for_unique = tmp_for_unique * self.ecal_config._N.n_channels + hits["chan"]
            b["nhit_chan"][0] = np.unique(tmp_for_unique[hits["isHit"]]).size
            b["sum_hg"][0] = hits["hg"][hits["isHit"]].sum()
            b["sum_energy"][0] = hits["energy"][hits["isHit"]].sum()
            if "energy_lg" in hits:
                b["sum_energy_lg"][0] = hits["energy_lg"][hits["isHit"]].sum()

            if len(hits) > self.ecal_config._N.bcid_too_many_hits:
                txt = "Suspicious number of hits! %i " %len(hits)
                txt += "for bcid %i and previous bcid %i" %(b["bcid"][0], b["prev_bcid"][0])
                txt += "\nSkipping event %i." % b["event"][0]
                print(txt)
                continue

            for hit_field in self._hit_branches:
                ev_val = hits[hit_field]
                b["hit_" + hit_field][:ev_val.size] = ev_val
            self.out_tree.Fill()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build an event-level rootfile (smaller) from the raw rootfile.",
    )
    parser.add_argument("file_name", help="The raw rootfile from converter_SLB")
    parser.add_argument("-n", "--max_entries", default=-1, type=int)
    parser.add_argument("-w", "--w_config", default=-1, type=int)
    parser.add_argument("-o", "--out_file_name", default=None)
    parser.add_argument("-c", "--commissioning_folder", default=None)
    parser.add_argument("-s", "--min_slabs_hit", default=4, type=int)
    parser.add_argument("--cob_positions_string", default="")
    parser.add_argument("--no_zero_suppress", action="store_true", help="Store all channels on hit chip.")
    parser.add_argument("--no_lg", action="store_true", help="Ignore low gain")
    # Run ./build_events.py --help to see all options.
    for config_option, config_value in dummy_config.items():
        parser.add_argument("--" + config_option, default=config_value)
    BuildEvents(**vars(parser.parse_args())).build_events()
