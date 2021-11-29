#!/usr/bin/env python
from __future__ import print_function

import argparse
import collections
import numpy as np
import ROOT as rt
from array import array
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
        print(datetime.now())
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
            yield i, spill


class BCIDHandler:
    def __init__(self, bcid_skip_noisy_acquisition_start, delta_merge):
        self.bcid_skip_noisy_acquisition_start = bcid_skip_noisy_acquisition_start
        self.delta_merge = delta_merge
        self.bad_bcid_value = -999


    def _get_corrected_bcids(self, bcids):
        """Using corrected_bcid, these corrections should already be in place."""
        bcids = np.array(list(bcids))
        bcids[bcids < self.bcid_skip_noisy_acquisition_start] = self.bad_bcid_value
        return bcids


    def _get_good_bcids(self, bcids, bad_bcids):
        """Returns Counter for #events within the spill with same good BCID."""
        good_bcid_counts = collections.Counter(bcids)
        # Apparently np.arrray(array.array) does not work, as numpy thinks the
        # array/buffer is longer, and includes the following memory-nonsense to.
        bad_bcid_values = set(np.array(list(bcids))[np.array(list(bad_bcids)) != 0])
        for bad_bcid in bad_bcid_values:
            good_bcid_counts.pop(bad_bcid)
        return good_bcid_counts


    def _choose_main_bcid_for_merger(self, bcid_group, good_bcid_counts):
        if len(bcid_group) == 1:
            main_bcid = bcid_group[0]
        else:
            group_counts = np.array([good_bcid_counts[x] for x in bcid_group])
            is_main_bcid_candidate = group_counts == max(group_counts)
            if sum(is_main_bcid_candidate) == 1:
                main_bcid = bcid_group[np.argmax(is_main_bcid_candidate)]
            else:
                bcid_group = np.array(bcid_group)
                weighted_mean = bcid_group * group_counts / len(bcid_group)
                weighted_mean -= 0.01  # Choose the smaller BCID when in between.
                main_bcid = bcid_group[np.argmin(np.abs(bcid_group - weighted_mean))]
        return main_bcid


    def merge_bcids(self, bcids, bad_bcids):
        bcids = self._get_corrected_bcids(bcids)
        good_bcid_counts = self._get_good_bcids(bcids, bad_bcids)
        good_bcids = np.array(sorted(good_bcid_counts))
        delta_good_bcids = good_bcids[1:] - good_bcids[:-1]
        merge_with_following = delta_good_bcids <= self.delta_merge

        if len(good_bcids) == 0: 
            return np.full_like(bcids, self.bad_bcid_value)
        good_bcid_groups = [[good_bcids[0]]]
        for i, do_merge in enumerate(merge_with_following, start=1):
            if do_merge:
                good_bcid_groups[-1].append(good_bcids[i])
            else:
                good_bcid_groups.append([good_bcids[i]])

        map_to_main_bcid_in_group = collections.defaultdict(lambda: self.bad_bcid_value)
        for group in good_bcid_groups:
            main_bcid = self._choose_main_bcid_for_merger(group, good_bcid_counts)
            for bcid_in_group in group:
                map_to_main_bcid_in_group[bcid_in_group] = main_bcid

        merged_bcids = np.array([map_to_main_bcid_in_group[b] for b in bcids])
        return merged_bcids


    def load_spill(self, spill_entry):
        self.merged_bcid = self.merge_bcids(spill_entry.corrected_bcid, spill_entry.badbcid)
        self.spill_bcids = sorted(set(self.merged_bcid) - {self.bad_bcid_value})

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
    n_chips = ecal_config._N.n_chips
    n_channels = ecal_config._N.n_channels
    n_scas = ecal_config._N.n_scas
    
    # It is faster to fill the arrays once outside the loop.
    all_gain_hit = np.array(entry.gain_hit_high)
    all_hg = np.array(entry.charge_hiGain)
    all_lg = np.array(entry.charge_lowGain)

    # assert bcid_handler.merged_bcid.shape == np.array(event.bcid).shape
    for i, bcid in enumerate(bcid_handler.merged_bcid):
        if bcid == bcid_handler.bad_bcid_value:
            continue
        bcid_channel_id = slice(i * n_channels, (i + 1) * n_channels)
        gain_hit_high = all_gain_hit[bcid_channel_id]
        if np.all(gain_hit_high < 0):
            continue
        if bcid not in event:
            event[bcid] = collections.defaultdict(list)

        slab = entry.slot[(i // n_scas // n_chips)]
        slab_id = ecal_config._N.slabs.index(slab)
        chip = entry.chipid[i // n_scas]
        sca = i % n_scas

        def ext(name, array):
            return event[bcid][name].extend(array)

        ext("isHit", gain_hit_high > 0)
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
        ext("lg", lg)

        n_in_batch = n_channels
        ext("slab", np.full(n_in_batch, slab, dtype=int))
        ext("chip", np.full(n_in_batch, chip, dtype=int))
        ext("sca", np.full(n_in_batch, sca, dtype=int))
        ext("chan", np.arange(n_channels, dtype=int))
        n_scas_filled = entry.nColumns[i // n_scas]
        ext("n_scas_filled", np.full(n_in_batch, n_scas_filled, dtype=int))

    for bcid, event_arrays in event.items():
        for key in event_arrays:
            event_arrays[key] = np.array(event_arrays[key])
        event[bcid] = event_arrays
    return event


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

    _branch_tags = {
        "event info": [
            "event/I",
            "spill/I",
            "bcid/I",
            "prev_bcid/I",
            "next_bcid/I",
        ],
        "hit summary": [
            "nhit_slab/I",
            "nhit_chip/I",
            "nhit_chan/I",
            "nhit_len/I",
            "sum_hg/F",
            "sum_energy/F",
        ],
        # 64: the number of channels on a chip.
        "hit id": [
            "hit_slab[nhit_len]/I",
            "hit_chip[nhit_len]/I",
            "hit_chan[nhit_len]/I",
            "hit_sca[nhit_len]/I",
        ],
        "hit coord": [
            "hit_x[nhit_len]/F",
            "hit_y[nhit_len]/F",
            "hit_z[nhit_len]/F",
        ],
        "hit readout": [
            "hit_hg[nhit_len]/I",
            "hit_lg[nhit_len]/I",
            "hit_energy[nhit_len]/F",
            "hit_n_scas_filled[nhit_len]/I",
        ],
        "hit booleans": [
            "hit_isHit[nhit_len]/I",
            "hit_isMasked[nhit_len]/I",
            "hit_isCommissioned[nhit_len]/I",
        ],
    }

    def __init__(
        self,
        file_name,
        w_config=-1,
        max_entries=-1,
        out_file_name=None,
        commissioning_folder=None,
        cob_positions_string="",
        ecal_numbers=None,  # Not provided in CLI. Mainly useful for debugging/changing.
        **config_file_kws,
    ):
        self.file_name = file_name
        self.w_config = w_config
        self.max_entries = max_entries
        self.out_file_name = out_file_name
        self.event_counter = 0

        self.in_tree = self._get_tree(file_name)
        slabs = self._get_slabs(self.in_tree)
        cob_slabs = set(map(int, filter(None, cob_positions_string.split(" "))))
        if ecal_numbers is None:
            ecal_numbers = EcalNumbers(slabs=slabs, cob_slabs=cob_slabs)
        else:
            assert ecal_numbers.slabs == slabs
            assert ecal_numbers.cob_slabs == cob_slabs

        self.ecal_config = EcalConfig(
            commissioning_folder=commissioning_folder,
            numbers=ecal_numbers,
            **config_file_kws,
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
        print("# Creating ecal tree in file %s" %out_file_name)
        self.out_file = rt.TFile(out_file_name,"recreate")
        self.out_tree = rt.TTree(self._out_tree_name, "Build ecal events")

        for branch_tag in [bt for bts in self._branch_tags.values() for bt in bts]:
            self._add_branch(branch_tag)
        self._hit_branches = _get_hit_branches(self.out_arrays)
        return self.out_tree


    def _write_and_close(self):
        self.out_tree.Write()
        print("\n# Created tree with %i events." % self.out_tree.GetEntries())
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
        self._fill_w_config_hist()

        if max_entries is None:
            max_entries = self.max_entries
        if max_entries < 0:
            max_entries = self.in_tree.GetEntries()

        for i_spill, entry in get_tree_spills(self.in_tree, max_entries):
            self._fill_spill(i_spill, entry)
        self._write_and_close()


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


    def _fill_spill(self, spill, entry):
        b = self.out_arrays
        b["spill"][0] = spill
        bcid_handler = BCIDHandler(
            self.ecal_config._N.bcid_skip_noisy_acquisition_start,
            self.ecal_config._N.bcid_merge_delta,
        )
        bcid_handler.load_spill(entry)
        hits_per_event = get_hits_per_event(entry, bcid_handler, self.ecal_config)

        for bcid in sorted(hits_per_event):
            hits = hits_per_event[bcid]
            if len(hits) == 0:
                raise EventBuildingException("Event with 0 hits. Why?")
            self.event_counter += 1
            b["event"][0] = self.event_counter
            b["bcid"][0] = bcid
            b["prev_bcid"][0] = bcid_handler.previous_bcid(bcid)
            b["next_bcid"][0] = bcid_handler.next_bcid(bcid)

            # count hits per slab/chan/chip
            tmp_for_unique = hits["slab"]
            b["nhit_slab"][0] = np.unique(tmp_for_unique).size
            tmp_for_unique = tmp_for_unique * self.ecal_config._N.n_chips + hits["chip"]
            b["nhit_chip"][0] = np.unique(tmp_for_unique).size
            b["nhit_len"][0] = hits["isHit"].size
            tmp_for_unique = tmp_for_unique * self.ecal_config._N.n_channels + hits["chan"]
            b["nhit_chan"][0] = np.unique(tmp_for_unique[hits["isHit"]]).size
            b["sum_hg"][0] = hits["hg"][hits["isHit"]].sum()
            b["sum_energy"][0] = hits["energy"][hits["isHit"]].sum()

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
    parser.add_argument("--cob_positions_string", default="")
    # Run ./build_events.py --help to see all options.
    for config_option, config_value in dummy_config.items():
        parser.add_argument("--" + config_option, default=config_value)
    BuildEvents(**vars(parser.parse_args())).build_events()
