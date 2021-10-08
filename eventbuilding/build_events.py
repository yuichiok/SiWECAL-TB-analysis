#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
import numpy as np
import ROOT as rt
from array import array
from help_tools import *

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

def get_hits(entry,bcid_map):
    ## Collect hits in bcid containe
    #print(bcid_map)
    event = {bcid:[] for bcid in bcid_map if bcid_map[bcid] > 0} # bcid : hits
    entry_bcids = entry.bcid
    gain_hit_low = entry.gain_hit_low
    gain_hit_high = entry.gain_hit_high
    charge_hiGain = entry.charge_hiGain
    charge_lowGain = entry.charge_lowGain
    tdc = entry.tdc

    for slab in range(NSLAB):
        for chip in range(NCHIP):
            for sca in range(NSCA):

                sca_indx = (slab * NCHIP + chip) * NSCA + sca
                #print("%i eee"%sca_indx)
                #print("%i %i %i %i %i"%(slab,NCHIP,chip,NSCA,sca))

                if sca_indx >= len(entry_bcids): continue
                #print(entry_bcids[sca_indx])
                bcid = get_corr_bcid(entry_bcids[sca_indx])
                #print(bcid)
                #print(bcid_map)

                # filter bad bcids
                if bcid not in bcid_map: continue
                # get assigned bcid
                bcid = bcid_map[bcid]

                if bcid not in event: continue

                ## energies
                for chan in range(NCHAN):
                    chan_indx = sca_indx * NCHAN + chan

                    #if not entry.gain_hit_low[chan_indx]: continue
                    isHit = gain_hit_high[chan_indx]
                    #if not isHit: continue

                    hg_ene = charge_hiGain[chan_indx]
                    lg_ene = charge_lowGain[chan_indx]
                    tdc_cp = tdc[chan_indx]

                    hit = EcalHit(slab,chip,chan,sca,hg_ene,lg_ene,tdc_cp,isHit)
                    event[bcid].append(hit)

    return event

def build_events(file_name, max_entries=-1, w_config=-1, out_file_name=None):

    ## Build tungsten config
    build_w_config(w_config)
    ## Read channel mapping
    read_mapping()
    ## Read channel mapping, cob
    read_mapping_cob()
    ## Read masked channels
    read_masked()
    ## Read pedestals
    read_pedestals()
    ## Read mip MPV values
    read_mip_values()

    # Get ttree
    tfile = rt.TFile(file_name,"read")
    tfile = rt.TFile(filename,"read")
    treename = "ecal_raw"
    tree = tfile.Get(treename)
    if not tree:
        print("No tree found in %s" %file_name)
        print(tree)
        exit(0)

    ##### TREE #####
    if out_file_name is None:
        if file_name.endswith("_converted.root"):
            out_file_name = file_name[:-len("_converted.root")] + "_build.root"
        elif file_name.endswith(".root"):
            out_file_name = file_name[:-len(".root")] + "_build.root"
        else:
            raise Exception("Unexpected file extension:", file_name)
    outf = rt.TFile(out_file_name,"recreate")
    outtree = rt.TTree("ecal","Build ecal events")

    print("# Creating ecal tree in file %s" %out_file_name)

    #### BRANCHES
    # event info
    event = array('i', [0]); outtree.Branch( 'event', event, 'event/I' )
    spill = array('i', [0]); outtree.Branch( 'spill', spill, 'spill/I' )
    bcid_b = array('i', [0]); outtree.Branch( 'bcid', bcid_b, 'bcid/I' )
    prev_bcid_b = array('i', [0]); outtree.Branch( 'prev_bcid', prev_bcid_b, 'prev_bcid/I' )
    next_bcid_b = array('i', [0]); outtree.Branch( 'next_bcid', next_bcid_b, 'next_bcid/I' )

    # occupancy/hit info
    nhit_slab = array('i', [0]); outtree.Branch( 'nhit_slab', nhit_slab, 'nhit_slab/I' )
    nhit_chip = array('i', [0]); outtree.Branch( 'nhit_chip', nhit_chip, 'nhit_chip/I' )
    nhit_chan = array('i', [0]); outtree.Branch( 'nhit_chan', nhit_chan, 'nhit_chan/I' )
    sum_hg = array('f', [0]); outtree.Branch( 'sum_hg', sum_hg, 'sum_hg/F' )
    sum_energy = array('f', [0]); outtree.Branch( 'sum_energy', sum_energy, 'sum_energy/F' )

    ## hit information
    # detid
    hit_slab = array('i', 10000*[0]); outtree.Branch( 'hit_slab', hit_slab, 'hit_slab[nhit_chan]/I' )
    hit_chip = array('i', 10000*[0]); outtree.Branch( 'hit_chip', hit_chip, 'hit_chip[nhit_chan]/I' )
    hit_chan = array('i', 10000*[0]); outtree.Branch( 'hit_chan', hit_chan, 'hit_chan[nhit_chan]/I' )
    hit_sca = array('i', 10000*[0]); outtree.Branch( 'hit_sca', hit_sca, 'hit_sca[nhit_chan]/I' )
    # coord
    hit_x = array('f', 10000*[0]); outtree.Branch( 'hit_x', hit_x, 'hit_x[nhit_chan]/F' )
    hit_y = array('f', 10000*[0]); outtree.Branch( 'hit_y', hit_y, 'hit_y[nhit_chan]/F' )
    hit_z = array('f', 10000*[0]); outtree.Branch( 'hit_z', hit_z, 'hit_z[nhit_chan]/F' )
    hit_x0 = array('f', 10000*[0]); outtree.Branch( 'hit_x0', hit_x0, 'hit_x0[nhit_chan]/F' )
    # energy
    hit_hg = array('f', 10000*[0]); outtree.Branch( 'hit_hg', hit_hg, 'hit_hg[nhit_chan]/F' )
    hit_lg = array('f', 10000*[0]); outtree.Branch( 'hit_lg', hit_lg, 'hit_lg[nhit_chan]/F' )
    hit_tdc = array('f', 10000*[0]); outtree.Branch( 'hit_tdc', hit_tdc, 'hit_tdc[nhit_chan]/F' )
    hit_energy = array('f', 10000*[0]); outtree.Branch( 'hit_energy', hit_energy, 'hit_energy[nhit_chan]/F' )
    # boolean
    hit_isHit = array('i', 10000*[0]); outtree.Branch( 'hit_isHit', hit_isHit, 'hit_isHit[nhit_chan]/I' )
    hit_isMasked = array('i', 10000*[0]); outtree.Branch( 'hit_isMasked', hit_isMasked, 'hit_isMasked[nhit_chan]/I' )

    if max_entries < 0: max_entries = tree.GetEntries()
    #else: max_entries = 1000

    spill_cnt = 0

    print("# Going to analyze %i entries..." %max_entries )
    for ientry,entry in enumerate(tree):#.GetEntries():

        if ientry > max_entries: break
        #if ientry != 9: continue

        if ientry%100 == 0: print("Entry %i" %ientry)

        ## BCID
        bcids = get_good_bcids(entry)
        bcid_map = merge_bcids(bcids)

        spill[0] = spill_cnt
        spill_cnt += 1

        ## Collect hits in bcid container
        ev_hits = get_hits(entry,bcid_map)

        #for bcid,hits in ev_hits.items():
        for ibc,bcid in enumerate(sorted(ev_hits)):

            hits = ev_hits[bcid]

            if len(hits) == 0: continue

            ## each bcid -- single event
            corr_bcid = get_corr_bcid(bcid)
            global event_counter
            event_counter = event_counter+1
            event[0] = event_counter#int(spill[0]*10000 + corr_bcid)
            bcid_b[0] = corr_bcid

            ## store distance to previous bcid
            if ibc > 0:
                prev_bcid = sorted(ev_hits)[ibc -1]
                prev_bcid_b[0] = get_corr_bcid(prev_bcid)
            else:
                prev_bcid_b[0] = -1

            if ibc + 1 < len(ev_hits):
                next_bcid = sorted(ev_hits)[ibc +1]
                next_bcid_b[0] = get_corr_bcid(next_bcid)
                #print("ibc=%i length=%i bcid=%i spill=%i"%(ibc,len(ev_hits),next_bcid_b[0],spill[0]))
            else:
                next_bcid_b[0] = -1

            # count hits per slab/chan/chip
            nhit_slab[0] = len(set([hit.slab for hit in hits]))
            nhit_chip[0] = len(set([(hit.slab*NCHIP + hit.chip) for hit in hits]))
            #nhit_chan[0] = len(set([(hit.slab*NCHIP + hit.chip)*NCHAN + hit.chan for hit in hits]))
            nhit_chan[0] = len(hits)
            sum_hg[0] = sum([hit.hg for hit in hits])
            sum_energy[0] = sum([hit.energy for hit in hits])

            if len(hits) > 8000:
                print("Suspicious number of hits! %i for bcid %i and previous bcid %i" %(len(hits),bcid_b[0],prev_bcid_b[0]))
                print("Skipping event %i" % event[0] )
                continue

            for i,hit in enumerate(hits):
                hit_slab[i] = hit.slab; hit_chip[i] = hit.chip; hit_chan[i] = hit.chan; hit_sca[i] = hit.sca
                hit_x[i] = hit.x; hit_y[i] = hit.y; hit_z[i] = hit.z; hit_x0[i] = hit.x0
                hit_hg[i] = hit.hg; hit_lg[i] = hit.lg; hit_tdc[i]=hit.tdc
                hit_isHit[i] = hit.isHit; hit_isMasked[i] = hit.isMasked

                hit_energy[i] = hit.energy

            outtree.Fill()

    outtree.Write()
    #outtree.Print()
    print("# Created tree with %i events" % outtree.GetEntries())
    outf.Close()

    tfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build an event-level rootfile (smaller) from the raw rootfile.",
    )
    parser.add_argument("file_name", help="The raw rootfile from converter_SLB")
    parser.add_argument("-n", "--max_entries", default=-1, type=int)
    parser.add_argument("-w", "--w_config", default=-1, type=int)
    parser.add_argument("-o", "--out_file_name", default=None)
    build_events(**vars(parser.parse_args()))
