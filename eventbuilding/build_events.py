#!/usr/bin/env python
from __future__ import print_function

import argparse
import collections
import sys
import numpy as np
import ROOT as rt
from array import array
from help_tools import *


try:
    from tqdm.autonotebook import tqdm

    def get_tree_events(tree, max_entries):
        for i, event in enumerate(tqdm(tree, desc="# Build events: ", total=max_entries)):
            if i > max_entries:
                break
            yield i, event
except ImportError:
    def get_tree_events(tree, max_entries):
        print("# Going to analyze %i entries..." %max_entries)
        print("# For better progress information: `pip install tqdm`.")
        for i, event in enumerate(tree):
            if i > max_entries:
                break
            if i%100 == 0:
                print("Entry %i" %i)
            yield i, event


class BCIDHandler:
    def __init__(self, val_evt, delta_merge):
        self.val_evt = val_evt
        self.delta_merge = delta_merge
        self.bad_bcid_value = -999


    def _get_corrected_bcid(self, bcid):
        if bcid < self.val_evt:
            # TODO: I understand why bcid = -999 is the filler value from converter_SLB. But why is bcid < 50 not ok?
            bcid = self.bad_bcid_value
        elif bcid >= 4096:  # 4096 = 2^12
            raise EventBuildingException("BCID Overflow:", bcid)
        return bcid


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
        bcids = [self._get_corrected_bcid(b) for b in bcids]
        good_bcid_counts = self._get_good_bcids(bcids, bad_bcids)
        good_bcids = np.array(sorted(good_bcid_counts))
        delta_good_bcids = good_bcids[1:] - good_bcids[:-1]
        merge_with_following = delta_good_bcids <= self.delta_merge

        good_bcid_groups = [[good_bcids[0]]]
        for i, do_merge in enumerate(merge_with_following, start=1):
            if do_merge:
                good_bcid_groups[-1].append(good_bcids[i])
            else:
                good_bcid_groups.append([good_bcids[i]])

        map_to_main_bcid_in_group = collections.defaultdict(lambda:self.bad_bcid_value)
        for group in good_bcid_groups:
            main_bcid = self._choose_main_bcid_for_merger(group, good_bcid_counts)
            for bcid_in_group in group:
                map_to_main_bcid_in_group[bcid_in_group] = main_bcid

        merged_bcids = np.array([map_to_main_bcid_in_group[b] for b in bcids])
        return merged_bcids


    def load_spill(self, spill_entry):
        # TODO: Should corrected_bcid be used instead of bcid?
        self.merged_bcid = self.merge_bcids(spill_entry.bcid, spill_entry.badbcid)
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
    event = collections.defaultdict(list)
    gain_hit_high = entry.gain_hit_high
    charge_hiGain = entry.charge_hiGain
    charge_lowGain = entry.charge_lowGain

    n_slabs = len(ecal_config._N.slab_map)
    n_chips = ecal_config._N.n_chips
    n_channels = ecal_config._N.n_channels
    n_scas = ecal_config._N.n_scas
    for i_slab in range(n_slabs):
        for i_chip in range(n_chips):
            for i_sca in range(n_scas):
                index_sca = (i_slab * n_chips + i_chip) * n_scas + i_sca
                bcid = bcid_handler.merged_bcid[index_sca]
                if bcid == bcid_handler.bad_bcid_value:
                    continue
                ## energies
                for i_channel in range(n_channels):
                    index_channel = index_sca * n_channels + i_channel
                    hit = EcalHit(
                        i_slab,
                        i_chip,
                        i_channel,
                        i_sca,
                        charge_hiGain[index_channel],
                        charge_lowGain[index_channel],
                        gain_hit_high[index_channel],
                        ecal_config,
                    )
                    event[bcid].append(hit)
    return event


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
            "sum_hg/F",
            "sum_energy/F",
        ],
        "hit id": [
            "hit_slab[nhit_chan]/I",
            "hit_chip[nhit_chan]/I",
            "hit_chan[nhit_chan]/I",
            "hit_sca[nhit_chan]/I",
        ],
        "hit coord": [
            "hit_x[nhit_chan]/F",
            "hit_y[nhit_chan]/F",
            "hit_z[nhit_chan]/F",
            "hit_x0[nhit_chan]/F",
        ],
        "hit readout": [
            "hit_hg[nhit_chan]/F",
            "hit_lg[nhit_chan]/F",
            "hit_energy[nhit_chan]/F",
        ],
        "hit booleans": [
            "hit_isHit[nhit_chan]/I",
            "hit_isMasked[nhit_chan]/I",
        ],
    }

    def __init__(
        self,
        file_name,
        w_config=-1,
        max_entries=-1,
        out_file_name=None,
        commissioning_folder=None,
        ecal_numbers=None,  # Not provided in CLI. Mainly useful for debugging/changing.
    ):
        self.ecal_config = EcalConfig(
            w_config=w_config,
            commissioning_folder=commissioning_folder,
            numbers=ecal_numbers,
        )
        self.file_name = file_name
        self.max_entries = max_entries
        self.out_file_name = out_file_name
        self.event_counter = 0


    def _get_tree(self, file_name):
        self.in_file = rt.TFile(file_name,"read")
        self.tree = self.in_file.Get(self._in_tree_name)
        if not self.tree:
            trees_available = [k.GetName() for k in self.in_file.GetListOfKeys()]
            print("Found tree names:", trees_available)
            ex_txt = "Tree %s not found in %s" %(self._in_tree_name, file_name)
            raise EventBuildingException(ex_txt)
        return self.tree


    def _add_branch(self, tag):
        name, branch_type = tag.split("/")
        branch_type = {"I": "i", "F": "f"}[branch_type]
        array_indicator = "[nhit_chan]"
        if array_indicator in name:
            name = name.replace(array_indicator, "")
            starting_value = 10000 * [0]
        else:
            starting_value = [0]
        self.out_arrays[name] = array(branch_type, starting_value)
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
        self._hit_branches = [br[4:] for br in self.out_arrays if br.startswith("hit_")]
        return self.out_tree


    def _write_and_close(self):
        self.out_tree.Write()
        # self.out_tree.Print()
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

        if max_entries is None:
            max_entries = self.max_entries
        if max_entries < 0:
            max_entries = self.in_tree.GetEntries()

        for i_spill, entry in get_tree_events(self.in_tree, max_entries):
            self._fill_spill(i_spill, entry)
        self._write_and_close()


    def _fill_spill(self, spill, entry):
        b = self.out_arrays
        b["spill"][0] = spill
        bcid_handler = BCIDHandler(
            self.ecal_config._N.bcid_val_event,
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
            b["nhit_slab"][0] = len(set([hit.slab for hit in hits]))
            b["nhit_chip"][0] = len(set([(hit.slab, hit.chip) for hit in hits]))
            b["nhit_chan"][0] = len(set([(hit.slab, hit.chip, hit.chan) for hit in hits]))
            b["sum_hg"][0] = sum([hit.hg for hit in hits])
            b["sum_energy"][0] = sum([hit.energy for hit in hits])

            if len(hits) > self.ecal_config._N.bcid_too_many_hits:
                txt = "Suspicious number of hits! %i " %len(hits)
                txt += "for bcid %i and previous bcid %i" %(b["bcid"][0], b["prev_bcid"][0])
                txt += "\nSkipping event %i." % b["event"][0]
                print(txt)
                continue

            for i,hit in enumerate(hits):
                for hit_attr in self._hit_branches:
                    b["hit_" + hit_attr] = getattr(hit, hit_attr)
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
    BuildEvents(**vars(parser.parse_args())).build_events()
