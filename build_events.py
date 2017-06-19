#!/usr/bin/env python
import sys
import numpy as np
import ROOT as rt
from array import array

NSLAB = 7
NCHIP = 16
NSCA = 15
NCHAN = 64

BCID_VALEVT = 1245

chan_map = {}

#pos_xzero = [2,2,2,2,4,4,6]#*0.56
pos_xzero = [2,4,6,8,12,16,22]#*0.56


class EcalHit:
    def __init__(self,slab,chip,chan,sca,hg,lg,isHit):
        self.slab = slab
        self.chip = chip
        self.chan = chan
        self.sca = sca
        self.hg = hg
        self.lg = lg
        self.isHit = isHit

        #global pos_xzero
        ## get x-y coordinates
        self.z = pos_xzero[slab] * 0.56
        (self.x,self.y) = chan_map[(chip,chan)]

class EcalEvent:

    def __init__(self,hits):
        self.hits = hits

def read_mapping(fname = "fev10_chip_channel_x_y_mapping.txt"):

    global chan_map# = {}

    with open(fname) as fmap:
        for i,line in enumerate(fmap.readlines()):
            if i == 0: continue

            # items: chip x0 y0 channel x y
            items = line.split()

            chip = int(items[0]); chan = int(items[3])
            x = float(items[4]); y = float(items[5])

            chan_map[(chip,chan)] = (x,y)

    return chan_map

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

            else:
                # found no bcids closeby
                # -> bcid counter untouched
                pass

    return entry_bcids_cnts

def get_good_bcids(entry):

    bcids = []

    for i,bcid in enumerate(entry.bcid):

        #chip = i % NCHIP
        #if chip not in [3,5,10,12]: continue
        #if chip != 5: continue

        if bcid < 0: continue
        if entry.badbcid[i] != 0: continue
        if entry.nhits[i] > 20: continue

        #if i%7 == 6: print "HERE"

        bcids.append(bcid)

    '''
    #bcids = [bcid for bcid in entry.bcid if bcid > -1]
    for slab in range(NSLAB):
        #if slab != 1: continue
        for chip in range(NCHIP):

            #if chip in [3,5,10,12]: continue
            if chip != 12: continue
            for sca in range(NSCA):

                bcid_indx = slab * NCHIP * NSCA + chip * NSCA + sca
                bcid = entry.bcid[bcid_indx]
                if bcid < 0: continue

                if entry.badbcid[bcid_indx] != 0: continue
                if entry.nhits[bcid_indx] > 20: continue

                bcids.append(bcid)
    '''
    return bcids

def get_hits(entry,bcids):
    ## Collect hits in bcid containe

    event = {bcid:[] for bcid in bcids if bcids[bcid] > 0} # bcid : hits

    slab_hit_cnts = NSLAB*[0]

    for slab in range(NSLAB):
        for chip in range(NCHIP):
            for sca in range(NSCA):

                sca_indx = (slab * NCHIP + chip) * NSCA + sca
                bcid = entry.bcid[sca_indx]

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
                for chan in range(NCHAN):
                    chan_indx = sca_indx * NCHAN + chan

                    #if not entry.gain_hit_low[chan_indx]: continue
                    isHit = entry.gain_hit_low[chan_indx]
                    if not isHit: continue

                    hg_ene = entry.charge_hiGain[chan_indx]
                    lg_ene = entry.charge_lowGain[chan_indx]

                    hit = EcalHit(slab,chip,chan,sca,hg_ene,lg_ene,isHit)
                    event[bcid].append(hit)

    return event

def build_events(filename, maxEntries = -1):

    ## Read channel mapping
    #chan_map = read_mapping()
    read_mapping()

    # Get ttree
    tfile = rt.TFile(filename,"read")
    treename = "fev10"
    tree = tfile.Get(treename)

    ##### TREE #####
    outfname = filename.replace("merge","build")
    #outf = rt.TFile("event_tree.root","recreate")
    outf = rt.TFile(outfname,"recreate")
    outtree = rt.TTree("ecal","Build ecal events")

    #### BRANCHES
    # event info
    event = array('i', [0]); outtree.Branch( 'event', event, 'event/I' )
    spill = array('i', [0]); outtree.Branch( 'spill', spill, 'spill/I' )
    bcid_b = array('i', [0]); outtree.Branch( 'bcid', bcid_b, 'bcid/I' )

    # occupancy/hit info
    nhit_slab = array('i', [0]); outtree.Branch( 'nhit_slab', nhit_slab, 'nhit_slab/I' )
    nhit_chip = array('i', [0]); outtree.Branch( 'nhit_chip', nhit_chip, 'nhit_chip/I' )
    nhit_chan = array('i', [0]); outtree.Branch( 'nhit_chan', nhit_chan, 'nhit_chan/I' )

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
    # energy
    hit_hg = array('f', 10000*[0]); outtree.Branch( 'hit_hg', hit_hg, 'hit_hg[nhit_chan]/F' )
    hit_lg = array('f', 10000*[0]); outtree.Branch( 'hit_lg', hit_lg, 'hit_lg[nhit_chan]/F' )
    #
    hit_isHit = array('i', 10000*[0]); outtree.Branch( 'hit_isHit', hit_isHit, 'hit_isHit[nhit_chan]/I' )

    if maxEntries == -1: maxEntries = tree.GetEntries()
    #else: maxEntries = 1000

    for ientry,entry in enumerate(tree):#.GetEntries():

        if ientry > maxEntries: break
        #if ientry != 9: continue

        if ientry%100 == 0: print("Entry %i" %ientry)

        ## BCID
        bcids = get_good_bcids(entry)
        bcid_cnts = merge_bcids(bcids)

        ## reset counters
        #nhit_slab[0] = nhit_chip[0] = nhit_chan[0] = nhit_sca[0] = 0
        spill[0] = entry.acqNumber

        ## Collect hits in bcid container
        ev_hits = get_hits(entry,bcid_cnts)
        for bcid,hits in ev_hits.iteritems():
            if bcid_cnts[bcid] < 1: continue #skip emptied bcids

            ## each bcid -- single event
            corr_bcid = bcid if bcid > 1245 else bcid + 4096
            event[0] = entry.acqNumber*10000 + corr_bcid
            bcid_b[0] = corr_bcid

            # count hits per slab/chan/chip
            nhit_slab[0] = len(set([hit.slab for hit in hits]))
            nhit_chip[0] = len(set([(hit.slab*NCHIP + hit.chip) for hit in hits]))
            #nhit_chan[0] = len(set([(hit.slab*NCHIP + hit.chip)*NCHAN + hit.chan for hit in hits]))
            nhit_chan[0] = len(hits)

            for i,hit in enumerate(hits):
                hit_slab[i] = hit.slab; hit_chip[i] = hit.chip; hit_chan[i] = hit.chan; hit_sca[i] = hit.sca
                hit_x[i] = hit.x; hit_y[i] = hit.y; hit_z[i] = hit.z
                hit_hg[i] = hit.hg; hit_lg[i] = hit.lg
                hit_isHit[i] = hit.isHit

            outtree.Fill()

    outtree.Write()
    outtree.Print()
    outf.Close()

    tfile.Close()


if __name__ == "__main__":

    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        #filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9_dif_1_1_1.raw.root"
        filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9__merge.root"
        #filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_10_all_difs_merge.root"
        #filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_10__merge.root"
    print("# Input file is %s" % filename)

    maxEntries = -1
    build_events(filename,maxEntries)
