#!/usr/bin/env python
import sys,os,glob,fnmatch,argparse
import ROOT as rt
rt.gROOT.LoadMacro("mergeRootFiles.C+")
from build_events import *

def convert_dir(indir,opts):

    ## search input dir for files matching run*_dif_1_1_1.raw.root
    pattern = "run_*" + opts.dif_suffix

    fnames = []
    if os.path.exists(indir):
        for root, dirs, files in os.walk(indir):
            for items in fnmatch.filter(files, pattern):
                fnames.append(root + "/" + items)
    else:
        # if input pattern already given
        fnames = glob.glob(indir + opts.dif_suffix)

    print("# Found %i runs" % len(fnames))

    for fname in fnames:
        # analyze run
        run_name = os.path.basename(fname).replace(opts.dif_suffix,"")
        run_dir =  os.path.dirname(fname)

        print 80*"#"
        print("## Analyzing run " + run_name)
        merge = rt.mergeRootFiles()
        run_pattern = "%s/%s_" %(run_dir,run_name)
        merge_fname = run_pattern + "_merge.root"

        ## Merge single slab files
        if not os.path.exists(merge_fname) or opts.force_merge:
            print("### Going to merge run pattern %s" %run_pattern)
            merge.Merge(run_pattern)
        else:
            print("### Already merged!")

        ## Build events
        build_fname = run_pattern + "_build.root"
        if not os.path.exists(build_fname) or opts.force_build:
            print("### Building events from merged file")
            build_events(merge_fname,opts.maxEntries,opts.w_config)
        else:
            print("### Already built!")



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Batch convert (merge&build) ecal runs.')
    parser.add_argument('indirs', nargs='+', help='input run dirs')
    parser.add_argument('--force-merge', dest='force_merge', action='store_true', help='Force merge files (default: False)')
    parser.add_argument('--force-build', dest='force_build', action='store_true', help='Force build files (default: False)')
    parser.add_argument('-n','--maxEntries', dest='maxEntries', type=int, default = -1, help='maximum number of entries to process (default: all)')
    parser.add_argument('--dif-suffix', dest='dif_suffix', default = "_dif_1_1_2.raw.root", help='single dif/slab suffix pattern (default: _dif_1_1_2.raw.root)')
    parser.add_argument('-c','--w_config', dest='w_config', type=int, default = 1, help=' (default: 1)')

    args = parser.parse_args()

    indirs = args.indirs
    print("Input dirs are ")
    print(indirs)

    for indir in indirs:
        convert_dir(indir,args)
