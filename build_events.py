#!/usr/bin/env python
import numpy as np
import ROOT as rt

from array import array

NSLAB = 7
NCHIP = 16
NSCA = 15
NCHAN = 64

BCID_VALEVT = 1245

chan_map = {}

class EcalHit:
    def __init__(self,slab,chip,chan,sca,hg,lg):
        self.slab = slab
        self.chip = chip
        self.chan = chan
        self.sca = sca
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
                #if bcids[bcid] == 0: continue

                ## if merged bcid, find closeby bcid
                if bcids[bcid] == 0:
                    for bcid_close in [bcid-1,bcid+1,bcid-2,bcid+2]:
                        if bcid_close in bcids:
                            if bcids[bcid_close] > 0: bcid = bcid_close

                ## energies
                for chan in range(NCHAN):
                    chan_indx = sca_indx * NCHAN + chan

                    if not entry.gain_hit_low[chan_indx]: continue

                    #slab_hit_cnts[slab] += 1

                    hg_ene = entry.charge_hiGain[chan_indx]
                    lg_ene = entry.charge_lowGain[chan_indx]

                    hit = EcalHit(slab,chip,chan,sca,hg_ene,lg_ene)
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

    ##### TREE #####
    outf = rt.TFile("event_tree.root","recreate")
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
    #nhit_sca = array('i', [0]); outtree.Branch( 'nhit_sca', nhit_sca, 'nhit_sca/I' )

    ## hit information
    # eleid
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

    for ientry,entry in enumerate(tree):#.GetEntries():

        if ientry > 100: break
        #if ientry != 9: continue

        if ientry%100 == 0: print("Entry %i" %ientry)

        #if len([1 for nhits in entry.nhits if nhits > 30]) > 20: continue
        #print sum([1 for nhits in entry.nhits if nhits > 20])

        ## event counter
        #acqNumber = entry.acqNumber

        ## BCID
        bcids = get_good_bcids(entry)
        bcid_cnts = merge_bcids(bcids)

        #> drawing
        for bcid,cnt in bcid_cnts.iteritems():
            if cnt > 0: hBCID.Fill(bcid,cnt)

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
            '''
            if nhit_chip[0]!= bcid_cnts[bcid] and nhit_chip[0] < 4:
                print bcid,nhit_chip[0], bcid_cnts[bcid]
                print bcid_cnts
                break
            '''
            for i,hit in enumerate(hits):
                hit_slab[i] = hit.slab; hit_chip[i] = hit.chip; hit_chan[i] = hit.chan; hit_sca[i] = hit.sca
                hit_x[i] = hit.x; hit_y[i] = hit.y; hit_z[i] = hit.z
                hit_hg[i] = hit.hg; hit_lg[i] = hit.lg

            outtree.Fill()


            '''
            if bcid_cnts[bcid] < 5: continue
            for hit in hits:
                hXY.Fill(hit.x,hit.y)#,hit.hg)
                hXZ.Fill(hit.z,hit.x,hit.hg)
                hXYZ.Fill(hit.slab,hit.x,hit.y)#,hit.hg)

                all_hits.append(hit)
            '''

    outtree.Write()
    outtree.Print()
    outf.Close()

    '''

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
    '''

    tfile.Close()
