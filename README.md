# SiWECAL-TB-analysis
Scripts for SiW ECAL test beam analysis

###################################
## Single slab analysis
RAW2ROOT.cc is the macro that reads the raw files and makes covnersion to root files, for single slabs.
Example of how to run it is in ConvertDirectory.cc:
root -l -q /home/irles/cernbox/TB2017/SiWECAL-TB-analysis-dev/ConvertDirectory.cc\(\"/Folder_path/"/"\",\""dif_1_1_1.raw"\",bcid_threshold\)
The bcid threshold is used to refine retriggers (retriggers are consecutive triggers with bcid[sca+1}-bcid[sca] < bcid_trheshold

Analysis tools are provided in

cd singleslab/

--> pedestal studies
--> mip studies.
(see README in the folder)


##################################3
## Built Events

## Batch conversion and event building
One can use the `batch_convert.py` script to automatically convert (raw2root), merge (single slab files) and build events.
The usage is:
```
./batch_convert.py path/to/dir/with/runs [options]
```
The options can be seen by executing without arguments.

# Single steps
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
- python event analyzer
