#!/usr/bin/env python
import numpy as np
import ROOT as rt

NSLAB = 7
NCHIP = 16
NSCA = 15
NCHAN = 64

BCID_VALEVT = 1245

def read_mapping(fname = "fev10_chip_channel_x_y_mapping.txt"):

    chan_map = {}

    with open(fname) as fmap:
        for i,line in enumerate(fmap.readlines()):
            if i == 0: continue

            # items: chip x0 y0 channel x y
            items = line.split()

            chip = items[0]; chan = items[3]
            x = items[4]; y = items[5]

            chan_map[(chip,chan)] = (x,y)

    return chan_map

if __name__ == "__main__":

    ## Read channel mapping
    #chan_map = read_mapping()

    # Get ttree
    #filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9_dif_1_1_1.raw.root"
    #filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/run_9__merge.root"
    filename = "/Users/artur/cernbox/CALICE/TB2017/data/Jun_2017_TB/BT2017/findbeam/test/run_10_all_difs_merge.root"

    tfile = rt.TFile(filename,"read")
    treename = "fev10"
    tree = tfile.Get(treename)

    #tree.Print()

    glob_event_counter = 0

    hist = rt.TH2F("hist","hist",4096,0,4095,50,0,50)
    #hist = rt.TH2F("hist","hist",50,0,50,4096,0,4095)

    for ientry,entry in enumerate(tree):#.GetEntries():

        if ientry > 1000: break
        #if ientry != 9: continue

        if ientry%100 == 0: print("Entry %i" %ientry)

        #if len([1 for nhits in entry.nhits if nhits > 30]) > 20: continue
        #print sum([1 for nhits in entry.nhits if nhits > 20])

        ## CHIP-wise vars
        #for i,bcid in enumerate(entry.bcid):
        #    a = i,bcid,entry.badbcid[i]
        #exit(0)

        ## Set of BCIDs present in this entry
        entry_bcids = [bcid for bcid in entry.bcid if bcid > -1]
        #print len(entry.bcid), len(entry.badbcid)
        entry_bcids = []#[bcid for i,bcid in enumerate(entry.bcid) if (bcid > -1 and entry.badbcid[i] == 0)]
        for i,bcid in enumerate(entry.bcid):
            if bcid < 0: continue
            if entry.badbcid[i] != 0: continue
            entry_bcids.append(bcid)

        entry_bcids_unique = set(entry_bcids)
        #entry_bcids_cnt = [entry_bcids.count(bcid) for bcid in entry_bcids_unique]
        entry_bcids_cnts = {bcid:entry_bcids.count(bcid) for bcid in entry_bcids_unique}

        #print entry_bcids
        #print sorted(entry_bcids_unique)
        #print entry_bcids_cnts

        ## Do BCID matching
        for i,bcid in enumerate(entry.bcid):

            if bcid < 0: continue
            if entry.badbcid[i] != 0: continue

            #if bcid not in entry_bcids_cnts: print "HELP"
            #print "Orig", bcid, entry_bcids_cnts[bcid]

            for bcid_close in [bcid-1,bcid+1]:#,bcid-2,bcid+2]:
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

        #print entry_bcids_cnts

        for bcid,cnt in entry_bcids_cnts.iteritems():
            hist.Fill(bcid,cnt)
            #hist.Fill(cnt,bcid)
            #if cnt > 1: print bcid,":", cnt,",",
        #print
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

    hist.Draw("colz")

    q = raw_input("Exit")

    tfile.Close()
