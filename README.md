# SiWECAL-TB-analysis
Scripts for SiW ECAL test beam analysis

## Merge single dif/slab files
Use `mergeRootFiles.C` to merge the root files from different slabs:

### python version
`$ ./mergeRootFiles.py /path/to/root/files/prefix_ (before dif_1)..)`

### original C macro version
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
