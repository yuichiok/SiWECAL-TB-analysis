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
- `spill`: Only counts those DAQ cycles where data was written.
  If the event building is performed on parts of the run
  (multiprocessing in [SiWECAL-TB-monitoring](https://github.com/SiWECAL-TestBeam/SiWECAL-TB-monitoring)),
  then this counter will restart on each new run.
- `event`: Event counter within a cycle.
  (_new_. Previously, the counter was not reset for the next cycle)
