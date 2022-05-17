# Event building

> _A. Irles, 1st August 2019. J. Kunath 2021._

This is the last step in the (central) data preparation.
The commissioning has to be performed before.


## Prerequisites

- The DAQ's data converted into root format
- A masked channels file
- A pedestal file
- A MIP calibration file

## Usage

Either change the default values (mainly in [`default_eventbuilding.cfg`](./default_eventbuilding.cfg))
or provide the changes through the command line interface.
For an exhaustive list of steering options, see the program's help:

```python
python buildevents.py --help
```

You will need a version of CERN root with python bindings. E.g.:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
```

## Features

- Build events from the raw data using calibration information.
- Apply a BCID merging of nearby, good BCIDs.
- Change the configuration at any time through the `--redo-config` mode.
- Compare two `build.root` files:
  - Tree comparison (e.g. the `ecal` tree) with
    `root -l -q -b test/AssertTreeEquality.C\(\"previous.root\",\"new.root\"\)`
  - (Calibration) histogram comparison with
    `./test/compare_histograms.py previous.root new.root --root_directory config`

## Dictionary

- `cycle`: `acqNumber` in converted.root files. Counts the acquisitions done by the DAQ.
- `event`: Event counter within a cycle.
  (_new_. Previously, the counter was not reset for the next cycle)
- `id_run`: E.g. 50123 for all events in a file. Counter that increases by 1 for each new run that we start.
- `id_dat`: Within a run, the DAQ will write the data into partial files.
   This tracks which of the partial files the event is from. 
   When we run the event building in parallel, or as part of the monitoring, the order might be shuffled.
   In most situations, this branch is not interesting (more for debugging).
- `nhit_len`: Is more an implementation detail (needed for the ROOT TBranches that take an array (e.g. hit_adc_high). 
   It tells ROOT the length of a specific leaf's array). 
   `nhit_len != nhit_chan` if you choose to also record the non-triggered channels (`zero_suppress = False`).
- `nhit_slab`: Number of slabs where at least one chip recorded a hit.
   E.g. for a MIP you would hope that this is close to 15. 
   In a 3GeV shower, 10 might be good (as the shower depletes).
- `spill`: Only counts those DAQ cycles where data was written.
   If the event building is performed on parts of the run
   (multiprocessing in [SiWECAL-TB-monitoring](https://github.com/SiWECAL-TestBeam/SiWECAL-TB-monitoring)),
   then this counter will restart on each new run.
