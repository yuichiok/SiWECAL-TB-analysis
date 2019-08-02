#!/usr/bin/env python
import sys
import ROOT as rt
rt.gROOT.LoadMacro("mergeRootFiles.cc+")

if __name__ == "__main__":

    if len(sys.argv) > 5:
        run = sys.argv[1]
        indir1 = sys.argv[2]
        indir2 = sys.argv[3]
        outdir = sys.argv[4]
        mode   = sys.argv[5]
        print("# Run name " + run)
        print("# Input dir for DIFs is " + indir1)
        print("# Input dir for SLBs is " + indir2)
        print("# Output dir is " + outdir)       
        print("# Data Mode for DIF slabs is " + mode)
    else:
        print("You need to parse at least five arguments, 5 strings!")
        print("          1: runname")
        print("          2: Input dir for DIFs")
        print("          3: Input dir for SLBs")
        print("          4: Output dir")
        print("          5: Data Mode for DIF slabs ")
        exit(0)

    if mode != "TDC" and mode != "HighLow":
        print("# Wrong mode " + mode + "  ONLY accepted modes are TDC or HighLow")
        exit(0)


    merge = rt.mergeRootFiles()
    merge.Merge(run,indir1,indir2,outdir,mode)
