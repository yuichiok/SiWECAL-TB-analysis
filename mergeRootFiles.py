#!/usr/bin/env python
import sys
import ROOT as rt
rt.gROOT.LoadMacro("mergeRootFiles.cc+")

if __name__ == "__main__":

    if len(sys.argv) > 1:
        indir = sys.argv[1]
        print("# Input dir (and prefix) is " + indir)
    else:
        print("No input directory given!")
        exit(0)

    if len(sys.argv) == 3:
        outname = sys.argv[2]
        print("# Output name is " + outname)
    else:
        outdir = "default"

    merge = rt.mergeRootFiles()
    merge.Merge(indir,outname)
