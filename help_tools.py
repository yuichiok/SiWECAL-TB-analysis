#!/usr/bin/env python
import os
import numpy as np

NSLAB = 7
NCHIP = 16
NSCA = 15
NCHAN = 64

BCID_VALEVT = 1245

chan_map = {}
#ped_map = read_pedestals()
ped_map = {}

# SLAB positions
pos_z = [0,1,2,3,4,5,9] * 15#mm gap

## Tungsten / W configuration
# Config 1
abs_thick = [2,2,2,2,4,4,6]
# Config 2
#abs_thick = [4,2,2,4,4,6,6]
# Config 3
#abs_thick = [6,2,4,4,6,6,6]

## sum up thickness
w_xzero = 0.56#Xo per mm of W
pos_xzero = [sum(abs_thick[:i+1])*w_xzero for i in range(len(abs_thick))]
## Print
print("W config used:")
print(abs_thick, pos_xzero)

class EcalHit:
    def __init__(self,slab,chip,chan,sca,hg,lg,isHit):
        self.slab = slab
        self.chip = chip
        self.chan = chan
        self.sca = sca
        self.hg = hg
        self.lg = lg
        self.isHit = isHit

        ## get x-y coordinates
        self.x0 = pos_xzero[slab]
        self.z = pos_z[slab]
        (self.x,self.y) = chan_map[(chip,chan)]

        # do pedestal subtraction
        self.hg -= ped_map[self.slab][self.chip][self.chan][self.sca]

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

def read_pedestals(indir_prefix = "./pedestals/"):

    global ped_map

    slab_map = {
        0: '_dif_1_1_1.txt',
        1: '_dif_1_1_2.txt',
        2: '_dif_1_1_3.txt',
        3: '_dif_1_1_4.txt',
        4: '_dif_1_1_5.txt',
        5: '_dif_1_2_1.txt',
        6: '_dif_1_2_2.txt'
    }

    ## pedestal map (n-dim numpy array)
    pedestal_map = np.zeros((NSLAB,NCHIP,NCHAN,NSCA))

    for slab in slab_map:
        fname = indir_prefix + "Pedestal" + slab_map[slab]
        print("Reading pedestals for %s from %s" %(slab,fname))
        if not os.path.exists(fname):
            print fname, " does not exist"
            continue

        with open(fname) as fmap:
            for i,line in enumerate(fmap.readlines()):
                if '#' in line: continue

                items = [float(item) for item in line.split()]

                chip,chan = int(items[0]),int(items[1])
                peds = items[2::2]
                peds_err = items[3::2]
                pedestal_map[slab][chip][chan] = peds

    ped_map = pedestal_map
    return pedestal_map

if __name__ == "__main__":

    print read_pedestals()
