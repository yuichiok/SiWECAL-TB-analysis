#!/usr/bin/env python
import numpy as np
import ROOT as rt

NSLAB = 7
NCHIP = 16
NSCA = 15
NCHAN = 64

BCID_VALEVT = 1245

chan_map = {}

class EcalHit:
    def __init__(self,slab,chip,chan,hg,lg):
        self.slab = slab
        self.chip = chip
        self.chan = chan
        self.hg = hg
        self.lg = lg

        ## get x-y coordinates
        self.z = slab#*10.
        #print chan_map
        (self.x,self.y) = chan_map[(chip,chan)]
        #if (chip,chan) in chan_map:
        #    print chan_map[(chip,chan)]
        #else:
        #    print (chip,chan), " missing"

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

    #print bcids

    ## Set of BCIDs present in this entry
    #entry_bcids = [bcid for bcid in entry.bcid if bcid > -1]

    #print entry_bcids
    #print len(entry.bcid), len(entry.badbcid)
    '''
    entry_bcids = []#[bcid for i,bcid in enumerate(entry.bcid) if (bcid > -1 and entry.badbcid[i] == 0)]
    for i,bcid in enumerate(entry.bcid):
        if bcid < 0: continue
        if entry.badbcid[i] != 0: continue
        entry_bcids.append(bcid)
    '''

    entry_bcids = bcids

    entry_bcids_unique = set(entry_bcids)
    #entry_bcids_cnt = [entry_bcids.count(bcid) for bcid in entry_bcids_unique]
    entry_bcids_cnts = {bcid:entry_bcids.count(bcid) for bcid in entry_bcids_unique}

    #print entry_bcids
    #print sorted(entry_bcids_unique)
    #print entry_bcids_cnts

    ## Do BCID matching
    #for i,bcid in enumerate(entry.bcid):
    for i,bcid in enumerate(entry_bcids):

        #if bcid < 0: continue
        #if entry.badbcid[i] != 0: continue

        #if bcid not in entry_bcids_cnts: print "HELP"
        #print "Orig", bcid, entry_bcids_cnts[bcid]

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
                '''
                if entry_bcids_cnts[bcid_close] < entry_bcids_cnts[bcid]:
                    continue
                    # nearby bcid has less counts
                    # delete nearby bcid, increment this bcid counter
                    entry_bcids_cnts[bcid] += 1
                    entry_bcids_cnts[bcid_close] -= 1
                else:
                    # nearby bcid has more counts
                    # -> assign this bcid to nearby bcid
                    entry_bcids_cnts[bcid_close] += 1
                    entry_bcids_cnts[bcid] -= 1
                '''
                break

            else:
                # found no bcids closeby
                # ->
                pass
                #print "Found no nearby", bcid


        '''
        if entry_bcids_cnts[bcid] == 1:
            # only this chip fired this bcid
            # -> search for nearby bcids
            for bcid_close in range(bcid-2,bcid+3):
                if bcid_close in entry_bcids_cnts:
                    # found bcid nearby
                    if entry_bcids_cnts[bcid_close] == 1:
                        # nearby bcid has only 1 count
                        # delete nearby bcid, increment this bcid
                        entry_bcids_cnts[bcid] += 1
                        del entry_bcids_cnts[bcid_close]
                    else:
                        # nearby bcid has more counts
                        # -> assign this bcid to nearby bcid
                        entry_bcids_cnts[bcid_close] += 1
                        del entry_bcids_cnts[bcid]
                else:
                    # found no bcids closeby
                    # ->
                    pass
        else:
            # multiple chips fired this bcid
            pass

        '''
    '''
    for bcid,cnt in entry_bcids_cnts.iteritems():
        hist.Fill(bcid,cnt)
        #hist.Fill(cnt,bcid)
        #if cnt > 1: print bcid,":", cnt,",",
    #print
    '''

    '''

    for slab in range(NSLAB):
        for chip in range(NCHIP):
            for sca in range(NSCA):

                bcid_indx = slab * NCHIP * NSCA + chip * NSCA + sca
                bcid = entry.bcid[bcid_indx]

    for slab in range(NSLAB):
        for chip in range(NCHIP):
            for sca in range(NSCA):
                for chan in range(NCHAN):

                    # example
                    glob_indx = slab * NSLAB * NCHIP * NSCA * NCHAN
                    a = entry.bcid[ glob_indx ]

    '''

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
                if bcids[bcid] == 0: continue

                ## energies
                for chan in range(NCHAN):
                    chan_indx = sca_indx * NCHAN + chan

                    if not entry.gain_hit_low[chan_indx]: continue

                    #slab_hit_cnts[slab] += 1

                    hg_ene = entry.charge_hiGain[chan_indx]
                    lg_ene = entry.charge_lowGain[chan_indx]

                    hit = EcalHit(slab,chip,chan,hg_ene,lg_ene)
                    event[bcid].append(hit)

    return event

if __name__ == "__main__":

    ## Read channel mapping
    #chan_map = read_mapping()
    read_mapping()

    # Get ttree
    #filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9_dif_1_1_1.raw.root"
    filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9__merge.root"
    #filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/test/run_10_all_difs_merge.root"

    tfile = rt.TFile(filename,"read")
    treename = "fev10"
    tree = tfile.Get(treename)

    #tree.Print()

    glob_event_counter = 0

    #bcid_cnts = {}
    hBCID = rt.TH2F("hBCID","hist",4096,0,4095,20,0,20)
    #hist = rt.TH2F("hist","hist",50,0,50,4096,0,4095)
    hXY = rt.TH2F("hXY","XY",32,-88,88,32,-88,88)
    #hXY = rt.TH2F("hXY","XY",100,-100,100,100,-100,100)
    hXYZ = rt.TH3F("hXYZ","XY",7,0,7,32,-88,88,32,-88,88)

    hXZ = rt.TH2F("hXZ","XY",7,0,7,32,-88,88)

    all_hits = []

    for ientry,entry in enumerate(tree):#.GetEntries():

        if ientry > 50: break
        #if ientry != 9: continue

        if ientry%100 == 0: print("Entry %i" %ientry)

        #if len([1 for nhits in entry.nhits if nhits > 30]) > 20: continue
        #print sum([1 for nhits in entry.nhits if nhits > 20])

        ## CHIP-wise vars
        bcids = get_good_bcids(entry)
        bcid_cnts = merge_bcids(bcids)
        for bcid,cnt in bcid_cnts.iteritems():
            if cnt > 0: hBCID.Fill(bcid,cnt)

        ## Collect hits in bcid containe
        ev_hits = get_hits(entry,bcid_cnts)
        for bcid,hits in ev_hits.iteritems():
            if bcid_cnts[bcid] < 5: continue
            for hit in hits:
                hXY.Fill(hit.x,hit.y)#,hit.hg)
                hXZ.Fill(hit.z,hit.x,hit.hg)
                hXYZ.Fill(hit.slab,hit.x,hit.y)#,hit.hg)

                all_hits.append(hit)

    hBCID.Draw("colz")
    #hXY.Draw("colz")
    hXYZ.Draw("lego")
    #hXZ.Draw("colz")

    gr = rt.TGraph2D()
    gr.SetMarkerStyle(20)
    gr.SetName("gr_event")
    for i, hit in enumerate(all_hits):
        if hit.hg > 420:
            gr.SetPoint(i,hit.z,hit.x,hit.y)
    gr.Draw("pcol")

    q = raw_input("Exit")

    tfile.Close()
