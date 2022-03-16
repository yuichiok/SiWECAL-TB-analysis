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
or provide the changes through the command line interface.
For an exhaustive list of steering options, see the program's help:

```python
python buildevents.py --help
```

The top-level [`cosmics.mk`](../cosmics.mk) Makefile demonstrates the whole data preparation chain,
starting with the raw `.datXXXX` files from the DAQ.

Note that this Makefile is not intended to be used automatically from raw to build.
Instead, use the intermediate steps that are provided and confirm their results
before continuing to the next step.

## Features

- Build events from the raw data using calibration information.
- Apply a BCID merging of nearby, good BCIDs.
- Change the configuration at any time through the `--redo-config` mode.
