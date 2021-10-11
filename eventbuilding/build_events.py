#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
import numpy as np
import ROOT as rt
from array import array
from help_tools import *

BCID_VALEVT = 50  # TODO: Put this into EcalNumbers?


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


def get_corr_bcid(bcid):
    #if bcid > 10: return bcid
    if bcid > BCID_VALEVT: return bcid
    else: return -9999
    #if bcid < 0: return 0
    #else: return bcid + 4096

def merge_bcids(bcid_cnts):
    ## Set of BCIDs present in this entry
    bcids_unique = set(bcid_cnts.keys())

    #bcids_cnts = bcids
    ## format: bcid: (counts,corresponding bcid)
    ## initialize corresponding with itself
    new_bcid_cnts = bcid_cnts#{bcid:(bcid_cnts[bcid],bcid) for bcid in bcid_cnts}
    bcid_map = {bcid:bcid for bcid in bcid_cnts}

    ## Do BCID matching
    for i,bcid in enumerate(bcids_unique):

        for bcid_close in [bcid-1,bcid+1,bcid-2,bcid+2,bcid-3,bcid+3]:
            if bcid_close in new_bcid_cnts:
                #print("Found nearby", bcid, bcid_close)
                #if bcids_cnts[bcid_close] < 1: continue

                # found bcid nearby
                # merge bcids based on "occupancy" counter
                if new_bcid_cnts[bcid_close] >= new_bcid_cnts[bcid]:
                    # nearby bcid has more counts
                    # -> assign this bcid to nearby bcid
                    new_bcid_cnts[bcid_close] += new_bcid_cnts[bcid]
                    new_bcid_cnts[bcid] -= new_bcid_cnts[bcid]

                    bcid_map[bcid] = bcid_close
                break

    return bcid_map

def get_good_bcids(entry):

    all_bcids = {}
    entry_badbcid = entry.badbcid
    entry_nhits = entry.nhits

    for i,bcid in enumerate(entry.bcid):
        if bcid < 0: continue

        bcid_flag = 0 #0 is OK!

        if entry_badbcid[i] > 0 or entry_badbcid[i] < 0: bcid_flag = 1
        #if entry_nhits[i] > 20: bcid_flag = 1

        bcid = get_corr_bcid(bcid)

        if bcid in all_bcids: all_bcids[bcid].append(bcid_flag)
        else: all_bcids[bcid] = [bcid_flag]

    ## make counter
    good_bcids = {}
    for bcid,flags in all_bcids.items():
        if sum(flags) == 0:
            good_bcids[bcid] = len(flags)

    return good_bcids

def get_hits(entry, bcid_map, ecal_config):
    ## Collect hits in bcid containe
    #print(bcid_map)
    event = {bcid:[] for bcid in bcid_map if bcid_map[bcid] > 0} # bcid : hits
    entry_bcids = entry.bcid
    gain_hit_low = entry.gain_hit_low
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
                if index_sca >= len(entry_bcids):
                    continue
                bcid = get_corr_bcid(entry_bcids[index_sca])
                # filter bad bcids
                if bcid not in bcid_map:
                    continue
                # get assigned bcid
                bcid = bcid_map[bcid]
                if bcid not in event:
                    continue
                ## energies
                for i_channel in range(n_channels):
                    index_channel = index_sca * n_channels + i_channel
                    #if not entry.gain_hit_low[index_channel]: continue
                    isHit = gain_hit_high[index_channel]
                    #if not isHit: continue
                    hg_ene = charge_hiGain[index_channel]
                    lg_ene = charge_lowGain[index_channel]
                    hit = EcalHit(i_slab, i_chip, i_channel, i_sca, hg_ene, lg_ene, isHit, ecal_config)
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
        ecal_numbers=None,  # Not provided in CLI. Mainly useful for debugging/changing.
    ):
        self.ecal_config = EcalConfig(w_config=w_config, numbers=ecal_numbers)
        self.file_name = file_name
        self.max_entries = max_entries
        self.out_file_name = out_file_name


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

        print(max_entries)
        for i_spill, entry in get_tree_events(self.in_tree, max_entries):
            self._fill_spill(i_spill, entry)
        self._write_and_close()


    def _fill_spill(self, spill, entry):
        b = self.out_arrays
        b["spill"][0] = spill

        ## BCID
        bcids = get_good_bcids(entry)
        bcid_map = merge_bcids(bcids)

        ## Collect hits in bcid container
        ev_hits = get_hits(entry, bcid_map, self.ecal_config)

        #for bcid,hits in ev_hits.items():
        for ibc,bcid in enumerate(sorted(ev_hits)):

            hits = ev_hits[bcid]

            if len(hits) == 0: continue

            ## each bcid -- single event
            corr_bcid = get_corr_bcid(bcid)
            global event_counter
            event_counter = event_counter+1
            b["event"][0] = event_counter#int(spill[0]*10000 + corr_bcid)
            b["bcid"][0] = corr_bcid

            ## store distance to previous bcid
            if ibc > 0:
                prev_bcid = sorted(ev_hits)[ibc -1]
                b["prev_bcid"][0] = get_corr_bcid(prev_bcid)
            else:
                b["prev_bcid"][0] = -1

            if ibc + 1 < len(ev_hits):
                next_bcid = sorted(ev_hits)[ibc +1]
                b["next_bcid"][0] = get_corr_bcid(next_bcid)
                #print("ibc=%i length=%i bcid=%i spill=%i"%(ibc,len(ev_hits),b["next_bcid"][0],spill[0]))
            else:
                b["next_bcid"][0] = -1

            # count hits per slab/chan/chip
            b["nhit_slab"][0] = len(set([hit.slab for hit in hits]))
            b["nhit_chip"][0] = len(set([(hit.slab*self.ecal_config._N.n_chips + hit.chip) for hit in hits]))
            # b["nhit_chan"][0] = len(set([(hit.slab*self.ecal_config._N.n_chips + hit.chip)*self.ecal_config._N.n_channels + hit.chan for hit in hits]))
            b["nhit_chan"][0] = len(hits)
            b["sum_hg"][0] = sum([hit.hg for hit in hits])
            b["sum_energy"][0] = sum([hit.energy for hit in hits])

            if len(hits) > 8000:
                print("Suspicious number of hits! %i for bcid %i and previous bcid %i" %(len(hits),b["bcid"][0],b["prev_bcid"][0]))
                print("Skipping event %i" % b["event"][0] )
                continue

            for i,hit in enumerate(hits):
                b["hit_slab"][i] = hit.slab; b["hit_chip"][i] = hit.chip; b["hit_chan"][i] = hit.chan; b["hit_sca"][i] = hit.sca
                b["hit_x"][i] = hit.x; b["hit_y"][i] = hit.y; b["hit_z"][i] = hit.z; b["hit_x0"][i] = hit.x0
                b["hit_hg"][i] = hit.hg; b["hit_lg"][i] = hit.lg
                b["hit_isHit"][i] = hit.isHit; b["hit_isMasked"][i] = hit.isMasked

                b["hit_energy"][i] = hit.energy

            self.out_tree.Fill()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build an event-level rootfile (smaller) from the raw rootfile.",
    )
    parser.add_argument("file_name", help="The raw rootfile from converter_SLB")
    parser.add_argument("-n", "--max_entries", default=-1, type=int)
    parser.add_argument("-w", "--w_config", default=-1, type=int)
    parser.add_argument("-o", "--out_file_name", default=None)
    BuildEvents(**vars(parser.parse_args())).build_events()
