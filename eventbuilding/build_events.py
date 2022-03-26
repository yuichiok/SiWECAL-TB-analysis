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


def get_tree_spills_no_progress_info(tree, max_entries):
    for i, spill in enumerate(tree):
        if i > max_entries:
            break
        yield i, spill


class BCIDHandler:
    def __init__(self, bcid_numbers, min_slabs_hit=1):
        self._N = bcid_numbers
        self._min_slabs_hit = min_slabs_hit


    def _calculate_overflow_correction(self, arr):
        """"Within a chips memory, the BCID must increase. A drop indicates overflow."""
        # assert np.max(arr) < self._N.bcid_overflow
        n_memory_cycles = np.cumsum(arr[:,1:] - arr[:,:-1] < 0, axis=1)
        is_filled = arr[:,1:] != self._N.bcid_bad_value
        overflown_bcids = np.copy(arr)
        overflown_bcids[:,1:] = arr[:,1:] + self._N.bcid_overflow * n_memory_cycles * is_filled
        bcid_to_overflow_adjusted = {}
        unique_overflow = np.unique(overflown_bcids)
        unique_overflow = unique_overflow[unique_overflow >= self._N.bcid_overflow]
        for of_bcid in np.unique(overflown_bcids):
            bcid_to_overflow_adjusted[of_bcid % self._N.bcid_overflow] = of_bcid
        return overflown_bcids, bcid_to_overflow_adjusted


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
        return pos_per_bcid, bcid_before_merge, bcid_start_stop


    def load_spill(self, spill_entry):
        """False if the current spill is empty."""
        chip_ids, bcids = self._get_bcids(spill_entry.bcid)
        overflown_bcids, self.bcid_to_overflow_adjusted = \
            self._calculate_overflow_correction(bcids)
        self.bcid_first_sca_full = self._find_bcid_filling_last_sca(overflown_bcids)
        self.pos_per_bcid, self._bcid_before_merge, bcid_start_stop = \
            self._get_array_positions_per_bcid(chip_ids, bcids)
        self.spill_bcids_overflown = []
        for start_bcid in self.pos_per_bcid.keys():
            n_memory_cycles = max([
                self.bcid_to_overflow_adjusted.get(b, b)
                for b in range(start_bcid, bcid_start_stop[start_bcid] + 1)
            ]) // self._N.bcid_overflow
            _overflown = start_bcid + n_memory_cycles * self._N.bcid_overflow
            self.spill_bcids_overflown.append(_overflown)
        self.spill_bcids_overflown = sorted(self.spill_bcids_overflown)
        return len(self.pos_per_bcid) > 0


    def yield_from_spill(self):
        for overflown_bcid in self.spill_bcids_overflown:
            bcid = overflown_bcid % self._N.bcid_overflow
            before_merge = self._bcid_before_merge[bcid]
            plus_overflow = (overflown_bcid // self._N.bcid_overflow) * self._N.bcid_overflow
            for i, pos in enumerate(self.pos_per_bcid[bcid]):
                yield pos, bcid + plus_overflow, before_merge[i] + plus_overflow


    def previous_bcid(self, overflown_bcid):
        i_bcid = self.spill_bcids_overflown.index(overflown_bcid)
        if i_bcid == 0:
            return -1
        else:
            return self.spill_bcids_overflown[i_bcid - 1]

    def next_bcid(self, overflown_bcid):
        i_bcid = self.spill_bcids_overflown.index(overflown_bcid)
        if i_bcid == len(self.spill_bcids_overflown) - 1:
            return -1
        else:
            return self.spill_bcids_overflown[i_bcid + 1]


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
        slab = entry.slboard_id[(i // n_scas // n_chips)]
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
            "cycle/I",
            "bcid/I",
            "bcid_first_sca_full/I",
            "bcid_merge_end/I",
            "bcid_prev/I",
            "bcid_next/I",
            "id_run/I",
            "id_dat/I",
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
        asu_version="",
        ecal_numbers=None,  # Not provided in CLI. Mainly useful for debugging/changing.
        no_zero_suppress=False,
        no_lg=False,
        redo_config=False,
        no_progress_info=False,
        id_dat=-1,
        id_run=-1,
        **config_file_kws
    ):
        self.file_name = file_name
        self.w_config = w_config
        self.max_entries = max_entries
        self.out_file_name = out_file_name
        self.min_slabs_hit = min_slabs_hit
        self._previous_cycle = 0
        self.event_counter = 0

        self.redo_config = redo_config
        self.in_tree = self._get_tree(file_name)
        slabs = self._get_slabs(self.in_tree)
        if ecal_numbers is None:
           ecal_numbers = EcalNumbers(slabs=slabs, asu_version=asu_version)
        else:
            assert ecal_numbers.slabs == slabs
            assert ecal_numbers.asu_version == asu_version

        self._no_progress_info = no_progress_info
        self._id_dat = id_dat
        self._id_run = id_run
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
        hist = rt.gDirectory.Get("slot_hist")
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
        self.out_file = rt.TFile(out_file_name,"recreate")
        self.out_tree = rt.TTree(self._out_tree_name, "Build ecal events")

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


    def _redo_fill_event(self, event):
        if event.nhit_slab < self.min_slabs_hit:
            return
        for branch_name in self._hit_branches:
            branch_name = "hit_"+ branch_name
            if branch_name == "hit_energy":
                mip = self.ecal_config.mip_map[event.hit_slab, event.hit_chip, event.hit_chan]
                pedestal = self.ecal_config.pedestal_map[event.hit_slab, event.hit_chip, event.hit_chan, event.hit_sca]
                arr = (event.hit_hg - pedestal) + mip
            elif branch_name == "hit_energy_lg":
                mip = self.ecal_config.mip_lg_map[event.hit_slab, event.hit_chip, event.hit_chan]
                pedestal = self.ecal_config.pedestal_lg_map[event.hit_slab, event.hit_chip, event.hit_chan, event.hit_sca]
                arr = (event.hit_lg - pedestal) + mip
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
            b["bcid_prev"][0] = bcid_handler.previous_bcid(bcid)
            b["bcid_next"][0] = bcid_handler.next_bcid(bcid)
            b["id_run"][0] = self._id_run
            b["id_dat"][0] = self._id_dat

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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file_name", help="The raw rootfile from converter_SLB.")
    parser.add_argument("-n", "--max_entries", default=-1, type=int)
    parser.add_argument("-w", "--w_config", default=-1, type=int)
    parser.add_argument("-o", "--out_file_name", default=None)
    parser.add_argument("-c", "--commissioning_folder", default=None)
    parser.add_argument("-s", "--min_slabs_hit", default=4, type=int)
    _help = "For each layer, one of 10,11,12,13,COB. Comma seperated list. "
    _help += "When left empty, all layers are assumed to be 12 (FEV12)."
    parser.add_argument("-a", "--asu_version", default="", help=_help)
    parser.add_argument("--no_zero_suppress", action="store_true", help="Store all channels on hit chip.")
    parser.add_argument("--no_lg", action="store_true", help="Ignore low gain.")
    _help ="Do not (re)-build the events but only change the configuration options. Then file_name should be a build.root file already."
    parser.add_argument("--redo_config", action="store_true", help=_help)
    _help = "Less verbose output, especially for batch processing."
    parser.add_argument("--no_progress_info", action="store_true", help=_help)
    parser.add_argument("--id_run", default=-1, type=int, help="Integer ID of the run within the testbeam campaign.")
    parser.add_argument("--id_dat", default=-1, type=int, help="Integer ID for piece-by-piece eventbuilding within a run.")
    # Run ./build_events.py --help to see all options.
    for config_option, config_value in dummy_config.items():
        parser.add_argument("--" + config_option, default=config_value)
    BuildEvents(**vars(parser.parse_args())).build_events()
