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

Either change the default values (mainly in [`help_tools.py`](./help_tools.py))
or provide the changed through the command line interface.
For an exhaustive list of steering options, see the program's help:

```python
python buildevents.py --help
```

The top-level [`cosmics.mk`](../cosmics.mk) Makefile demonstrates the whole data preparation chain,
starting with the raw `.dat` files from the DAQ.

Note that this Makefile is not intended to be used automatically from raw to build.
Instead, use the intermediate steps that are provided and confirm their results
before continuing to the next step.

## Features

- Build events from the raw data using calibration information.
- Apply a BCID merging of nearby, good BCIDs.

Merges all 9 files in one (adding one extra dimension to the arrays).
It does also the inversion between high gain and tdc if this is the case.
It also makes the bcid-offset correction between fev13s and slbs

## (Old) ToDos

1. What about the auto gain? This bit is not included in the standard RAW2ROOT output, 
  so I did not used it here.
