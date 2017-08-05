#!/usr/bin/env python
import sys
import numpy as np
import ROOT as rt
from array import array
from help_tools import *

def get_corr_bcid(bcid):
    return bcid if bcid > BCID_VALEVT else bcid + 4096
    #return bcid

def merge_bcids(bcids):
    ## Set of BCIDs present in this entry
    #entry_bcids = [bcid for bcid in entry.bcid if bcid > -1]

    entry_bcids = bcids

    entry_bcids_unique = set(entry_bcids)
    #entry_bcids_cnt = [entry_bcids.count(bcid) for bcid in entry_bcids_unique]
    entry_bcids_cnts = {bcid:entry_bcids.count(bcid) for bcid in entry_bcids_unique}

    ## Do BCID matching
    for i,bcid in enumerate(entry_bcids):
        #if bcid < 0: continue

        for bcid_close in [bcid-1,bcid+1,bcid-2,bcid+2]:
            if bcid_close in entry_bcids_cnts:
                #print "Found nearby", bcid, bcid_close
                #if entry_bcids_cnts[bcid_close] < 1: continue

                # found bcid nearby
                # merge bcids based on "occupancy" counter
                if entry_bcids_cnts[bcid_close] >= entry_bcids_cnts[bcid]:
                    # nearby bcid has more counts
                    # -> assign this bcid to nearby bcid
                    entry_bcids_cnts[bcid_close] += 1
                    entry_bcids_cnts[bcid] -= 1
                break

    return entry_bcids_cnts

def get_good_bcids(entry):

    bcids = []
    entry_badbcid = entry.badbcid
    entry_nhits = entry.nhits

    for i,bcid in enumerate(entry.bcid):

        if bcid < 0: continue
        if entry_badbcid[i] > 1 or entry_badbcid[i] < 0: continue
        if entry_nhits[i] > 20: continue

        bcids.append(get_corr_bcid(bcid))

    return bcids

def get_hits(entry,bcids):
    ## Collect hits in bcid containe

    event = {bcid:[] for bcid in bcids if bcids[bcid] > 0} # bcid : hits
    entry_bcids = entry.bcid
    gain_hit_low = entry.gain_hit_low
    charge_hiGain = entry.charge_hiGain
    charge_lowGain = entry.charge_lowGain

    for slab in xrange(NSLAB):
        for chip in xrange(NCHIP):
            for sca in xrange(NSCA):

                sca_indx = (slab * NCHIP + chip) * NSCA + sca
                bcid = get_corr_bcid(entry_bcids[sca_indx])
                
                # filter bad bcids
                if bcid not in bcids: continue
                # filter merged bcids
                #if bcids[bcid] == 0: continue
                
                ## if merged bcid, find closeby bcid
                if bcids[bcid] == 0:
                    for bcid_close in [bcid-1,bcid+1,bcid-2,bcid+2]:
                        if bcid_close in bcids:
                            if bcids[bcid_close] > 0: bcid = bcid_close
                            
                            ## energies
                for chan in xrange(NCHAN):
                    chan_indx = sca_indx * NCHAN + chan
                    
                    #if not entry.gain_hit_low[chan_indx]: continue
                    isHit = gain_hit_low[chan_indx]
                    if not isHit: continue
                    
                    hg_ene = charge_hiGain[chan_indx]
                    lg_ene = charge_lowGain[chan_indx]
                    
                    hit = EcalHit(slab,chip,chan,sca,hg_ene,lg_ene,isHit)
                    event[bcid].append(hit)

    return event

def build_events(filename, maxEntries = -1, w_config = 0):

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
        bcid_cnts = merge_bcids(bcids)

        ## reset counters
        #nhit_slab[0] = nhit_chip[0] = nhit_chan[0] = nhit_sca[0] = 0

        spill[0] = spill_cnt
        spill_cnt += 1

        ## Collect hits in bcid container
        #ev_hits = get_hits(entry,bcid_cnts)
        ev_hits = get_hits(entry,bcid_cnts)

        #for bcid,hits in ev_hits.iteritems():
        for ibc,bcid in enumerate(sorted(ev_hits)):
            hits = ev_hits[bcid]

            if bcid_cnts[bcid] < 1:
                print "Hits with empty bcid! %i ! This should not happen." %bcid
                continue #skip emptied bcids

            ## each bcid -- single event
            corr_bcid = get_corr_bcid(bcid)
            event[0] = int(spill[0]*5000 + corr_bcid)
            bcid_b[0] = corr_bcid

            ## store distance to previous bcid
            if ibc > 0:
                prev_bcid = sorted(ev_hits)[ibc -1]
                prev_bcid_b[0] = get_corr_bcid(prev_bcid)
            else:
                prev_bcid_b[0] = -1

            # count hits per slab/chan/chip
            nhit_slab[0] = len(set([hit.slab for hit in hits]))
            nhit_chip[0] = len(set([(hit.slab*NCHIP + hit.chip) for hit in hits]))
            #nhit_chan[0] = len(set([(hit.slab*NCHIP + hit.chip)*NCHAN + hit.chan for hit in hits]))
            nhit_chan[0] = len(hits)
            sum_hg[0] = sum([hit.hg for hit in hits])
            sum_energy[0] = sum([hit.energy for hit in hits])

            if len(hits) > 1000:
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

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else: 
        filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9__merge.root"
        maxEntries = -1
        w_config = 0

    if len(sys.argv) > 2:
        maxEntries = int(sys.argv[2])
    else:
        maxEntries = -1
        w_config = 0

    if len(sys.argv) > 3:
        print(sys.argv[3])
        w_config = int(sys.argv[3])
    else:
        w_config = 0

    print("# Input file is %s" % filename)
    print("# maxEntries is %i" % maxEntries)
    print("# W-config is %i" % w_config)


    if os.path.exists(filename):
        build_events(filename,maxEntries)
    else:
        print("The file does not exist!")
