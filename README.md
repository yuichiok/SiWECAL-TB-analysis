# SiWECAL-TB-analysis
Scripts for SiW ECAL test beam analysis

## Merge single dif/slab files
Use `mergeRootFiles.C` to merge the root files from different slabs:

#### python version
`$ ./mergeRootFiles.py /path/to/root/files/prefix_ (before dif_1)..)`

#### original C macro version
```
root[] .L mergeRootFiles.C
root[] m = new mergeRootFiles()
root[] m->Merge("path/to/dir/with/files/run_prefix_")
```

## Build events (python builder)
Use build_events.py to build BCID-based events with hits from the merged file.
Run as:
`$ ./build_events.py merged_filename.root`
One can modify the `maxEntries` to process.

The pedestal file directory has to be updated in the `help_tools.py` file.
The tungsten/W configuration is hard-coded in the `build_events.py`.

## Analyze build events
For now simple analysis is possible directly from the ROOT prompt with `ecal->Draw()`.
A python analyzer will be added soon.

# TODO:
- per-channel MIP calibration
- python event analyzer
