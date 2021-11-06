# SiWECAL-TB-analysis

> _A. Irles, 1st August 2019._ irles NOT SPAM @lal.in2p3.fr

Scripts for SiW ECAL test beam analysis

## Single slab analysis

Analysis tools are provided in

cd singleslab/

--> pedestal studies
--> mip studies.
(see README in the folder)

## Build Events

See [eventbuilding README](eventbuilding/README.md).

### Steering the build process through a Makefile

To keep your changes in variables consistent throughout the event building
for cosmics, a Makefile can be used.
Most likely you will find the variables that you have to change 
in the first rows of the [`cosmics.mk`](cosmics.mk) Makefile.
Usage example:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
make cosmics.mk build -j8 \
  WOLFRAM_CONFIG=0 \
  RAW_DATA_DIR=/data_ilc/flc/ECAL/cosmics/TB2021-11/run_050016_10192021_21h49min_Ascii
```

It is recommended to split the make run into parts:

1. `make cosmics.mk converted`
2. `make cosmics.mk masked`
3. `make cosmics.mk pedestals`
4. `make cosmics.mk mip_calib`
5. `make cosmics.mk run`

Then, check the results of each step before proceeding to the next.


## Analyze build events

In progress,

folder proto
