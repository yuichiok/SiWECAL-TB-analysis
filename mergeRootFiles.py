#!/usr/bin/env python
import sys
import ROOT as rt
rt.gSystem.Load("mergeRootFiles_C.so")

if __name__ == "__main__":

    if len(sys.argv) > 1:
        indir = sys.argv[1]
        print("# Input dir (and prefix) is" + indir)
    else:
        print("No input directory given!")
        exit(0)

    merge = rt.mergeRootFiles()
    merge.Merge(indir)
