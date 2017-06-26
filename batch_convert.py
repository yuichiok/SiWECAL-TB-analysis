#!/usr/bin/env python
import sys,os,glob,fnmatch
import ROOT as rt
rt.gROOT.LoadMacro("mergeRootFiles.C+")
from build_events import *

maxEntries = -1#300
dif_pattern = "_dif_1_1_2.raw.root"
force_merge = True
force_build = True

if __name__ == "__main__":

    if len(sys.argv) > 1:
        indir = sys.argv[1]
        print("# Input dir is" + indir)
    else:
        print("No input directory given!")
        exit(0)

    ## search input dir for files matching run*_dif_1_1_1.raw.root
    pattern = "run_1*" + dif_pattern

    fnames = []
    for root, dirs, files in os.walk(indir):
        #print root,dirs
        for items in fnmatch.filter(files, pattern):
            fnames.append(root + "/" + items)

    print("# Found %i runs" % len(fnames))

    for fname in fnames:
        # analyze run
        run_name = os.path.basename(fname).replace(dif_pattern,"")
        run_dir =  os.path.dirname(fname)

        print 80*"#"
        print("## Analyzing run " + run_name)
        merge = rt.mergeRootFiles()
        run_pattern = "%s/%s_" %(run_dir,run_name)
        merge_fname = run_pattern + "_merge.root"

        ## Merge single slab files
        if not os.path.exists(merge_fname) or force_merge:
            print("### Going to merge run pattern %s" %run_pattern)
            merge.Merge(run_pattern)
        else:
            print("### Already merged!")

        ## Build events
        build_fname = run_pattern + "_build.root"
        if not os.path.exists(build_fname) or force_build:
            print("### Building events from merged file")
            build_events(merge_fname,maxEntries)
        else:
            print("### Already built!")

