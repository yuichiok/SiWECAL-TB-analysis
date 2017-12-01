#!/usr/bin/env python
import sys
import numpy as np
import ROOT as rt
from array import array
from help_tools import *

def get_corr_bcid(bcid):
    if bcid < 0: return bcid
    if bcid > BCID_VALEVT: return bcid
    else: return bcid + 4096

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

        for bcid_close in [bcid-1,bcid+1,bcid-2,bcid+2]:
            if bcid_close in new_bcid_cnts:
                #print "Found nearby", bcid, bcid_close
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

        if entry_badbcid[i] > 1 or entry_badbcid[i] < 0: bcid_flag = 1
        #if entry_nhits[i] > 20: bcid_flag = 1

        bcid = get_corr_bcid(bcid)

        if bcid in all_bcids: all_bcids[bcid].append(bcid_flag)
        else: all_bcids[bcid] = [bcid_flag]

    ## make counter
    good_bcids = {}
    for bcid,flags in all_bcids.iteritems():
        if sum(flags) == 0:
            good_bcids[bcid] = len(flags)

    return good_bcids

def get_hits(entry,bcid_map):
    ## Collect hits in bcid containe

    event = {bcid:[] for bcid in bcid_map if bcid_map[bcid] > 0} # bcid : hits
    entry_bcids = entry.bcid
    gain_hit_low = entry.gain_hit_low
    gain_hit_high = entry.gain_hit_high
    charge_hiGain = entry.charge_hiGain
    charge_lowGain = entry.charge_lowGain

    for slab in xrange(NSLAB):
        for chip in xrange(NCHIP):
            for sca in xrange(NSCA):

                sca_indx = (slab * NCHIP + chip) * NSCA + sca
                bcid = get_corr_bcid(entry_bcids[sca_indx])
                
                # filter bad bcids
                if bcid not in bcid_map: continue
                # get assigned bcid
                bcid = bcid_map[bcid]

                if bcid not in event: continue

                ## energies
                for chan in xrange(NCHAN):
                    chan_indx = sca_indx * NCHAN + chan
                    
                    #if not entry.gain_hit_low[chan_indx]: continue
                    isHit = gain_hit_high[chan_indx]
                    #if not isHit: continue
                    
                    hg_ene = charge_hiGain[chan_indx]
                    lg_ene = charge_lowGain[chan_indx]
                    
                    hit = EcalHit(slab,chip,chan,sca,hg_ene,lg_ene,isHit)
                    event[bcid].append(hit)

    return event

def build_events(filename, maxEntries = -1, w_config = -1):

    ## Build tungsten config
    build_w_config(w_config)
    ## Read channel mapping
    read_mapping()
    ## Read masked channels
    read_masked()
    ## Read pedestals
    read_pedestals()
    ## Read mip MPV values
    read_mip_values()

    # Get ttree
    tfile = rt.TFile(filename,"read")
    treename = "fev10"
    tree = tfile.Get(treename)
    if not tree:
        print("No tree found in ")
        print(tree)
        exit(0)

    ##### TREE #####
    outfname = filename.replace("merge","build")
    outf = rt.TFile(outfname,"recreate")
    outtree = rt.TTree("ecal","Build ecal events")

    print("# Creating ecal tree in file %s" %outfname)

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
    hit_energy = array('f', 10000*[0]); outtree.Branch( 'hit_energy', hit_energy, 'hit_energy[nhit_chan]/F' )
    # boolean
    hit_isHit = array('i', 10000*[0]); outtree.Branch( 'hit_isHit', hit_isHit, 'hit_isHit[nhit_chan]/I' )
    hit_isMasked = array('i', 10000*[0]); outtree.Branch( 'hit_isMasked', hit_isMasked, 'hit_isMasked[nhit_chan]/I' )

    if maxEntries < 0: maxEntries = tree.GetEntries()
    #else: maxEntries = 1000

    spill_cnt = 0

    print("# Going to analyze %i entries..." %maxEntries )
    for ientry,entry in enumerate(tree):#.GetEntries():

        if ientry > maxEntries: break
        #if ientry != 9: continue

        if ientry%100 == 0: print("Entry %i" %ientry)

        ## BCID
        bcids = get_good_bcids(entry)
        bcid_map = merge_bcids(bcids)

        spill[0] = spill_cnt
        spill_cnt += 1

        ## Collect hits in bcid container
        ev_hits = get_hits(entry,bcid_map)

        #for bcid,hits in ev_hits.iteritems():
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
                hit_hg[i] = hit.hg; hit_lg[i] = hit.lg
                hit_isHit[i] = hit.isHit; hit_isMasked[i] = hit.isMasked

                hit_energy[i] = hit.energy

            outtree.Fill()

    outtree.Write()
    #outtree.Print()
    print("# Created tree with %i events" % outtree.GetEntries())
    outf.Close()

    tfile.Close()


if __name__ == "__main__":


    filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9__merge.root"
    maxEntries = -1
    w_config = 0

    if len(sys.argv) < 3:
        filename = sys.argv[1]
    elif len(sys.argv) < 4:
        filename = sys.argv[1]
        maxEntries = int(sys.argv[2])
    elif len(sys.argv) < 5:
        filename = sys.argv[1]
        maxEntries = int(sys.argv[2])
        w_config = int(sys.argv[3])

    print("# Input file is %s" % filename)
    print("# maxEntries is %i" % maxEntries)
    print("# W-config is %i" % w_config)


    if os.path.exists(filename):
        build_events(filename,maxEntries,w_config)
    else:
        print("The file does not exist!")
