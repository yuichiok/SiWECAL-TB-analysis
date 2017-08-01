#!/usr/bin/env python
import sys,os,glob,fnmatch,argparse
import ROOT as rt
rt.gROOT.LoadMacro("RAW2ROOT.cc+")
rt.gROOT.LoadMacro("mergeRootFiles.cc+")
from build_events import *

def new_path(old_name,opts):

    if opts.outdir == "default":
        return old_name

    ## original filename: old_name
    ## we want to create the same filename in another directory
    new_name = opts.outdir + "/"
    dirname = os.path.dirname(old_name)
    if opts.cutoff == "default":
        cutoff_dirs = dirname
        new_name += cutoff_dirs + "/"
    else:
        cutoff_dirs = dirname[dirname.find(opts.cutoff):]
        new_name += cutoff_dirs + "/"

    ## new to create new dir
    if not os.path.exists(new_name):
        print("# Creating output directory %s" %new_name )
        os.makedirs(new_name)

    new_name += os.path.basename(old_name)

    print("# New path: %s" % new_name)

    return new_name


def convert_run_raw(run_pattern,opts):

    fnames_raw = glob.glob(run_pattern + "*.raw")
    fnames_root = glob.glob(run_pattern + "*.raw.root")

    raw2root = rt.RAW2ROOT()
    for fname_raw in fnames_raw:
        if (fname_raw + ".root" not in fnames_root) or opts.force_convert:
            print("# Need to convert " + fname_raw)
            if opts.outdir == "default":
                raw2root.ReadFile(fname_raw,optsforce_convert)
            else:
                outfname = new_path(fname_raw,opts) + ".root"
                raw2root.ReadFile(fname_raw,opts.force_convert,outfname)


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

    ## loop through files and convert, merge and build them
    for fname in fnames:
        # analyze run
        run_name = os.path.basename(fname).replace(opts.dif_suffix,"")
        run_dir =  os.path.dirname(fname)

        print 80*"#"
        print("## Analyzing run " + run_name)
        run_pattern = "%s/%s_" %(run_dir,run_name)

        ## convert missing files
        convert_run_raw(run_pattern,opts)

        merge = rt.mergeRootFiles()
        merge_fname = run_pattern + "merge.root"
        # new path
        #merge_fname = new_path(merge_fname,opts)

        ## Merge single slab files
        if opts.force_merge:
            merge_fname = new_path(merge_fname,opts)

        if not os.path.exists(merge_fname) or opts.force_merge:
            print("### Going to merge run pattern %s" %run_pattern)
            merge.Merge(run_pattern,merge_fname)
        else:
            print("### Already merged!")

        ## Build events
        build_fname = merge_fname.replace("merge","build")
        if opts.force_build:
            build_fname = new_path(build_fname,opts)

        if not os.path.exists(build_fname) or opts.force_build:
            print("### Building events from merged file")
            build_events(merge_fname,opts.maxEntries,opts.w_config)
        else:
            print("### Already built!")



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Batch convert (merge&build) ecal runs.')
    parser.add_argument('indirs', nargs='+', help='input run dirs')
    parser.add_argument('--force-r2r','--force-convert', dest='force_convert', action='store_true', help='Force convert files (default: False)')
    parser.add_argument('--force-merge', dest='force_merge', action='store_true', help='Force merge files (default: False)')
    parser.add_argument('--force-build', dest='force_build', action='store_true', help='Force build files (default: False)')
    parser.add_argument('-n','--maxEntries', dest='maxEntries', type=int, default = -1, help='maximum number of entries to process (default: all)')
    parser.add_argument('--dif-suffix', dest='dif_suffix', default = "_dif_1_1_2.raw", help='single dif/slab suffix pattern (default: _dif_1_1_2.raw)')
    parser.add_argument('-c','--w_config', dest='w_config', type=int, default = 1, help=' (default: 1)')
    parser.add_argument('--outdir', dest='outdir', default = "default", help="output directory")
    parser.add_argument('--cutoff', dest='cutoff', default = "default", help="specify from which dir in the indir path to preserve the dir structure")

    args = parser.parse_args()

    indirs = args.indirs
    print("Input dirs are :")
    print(indirs)

    if args.force_convert: args.force_merge = True
    if args.force_merge:   args.force_build = True

    '''
    ##### fiddle with outdir if it is not default
    if args.outdir != "default":
        if not os.path.exists(args.outdir):
            print("# Creating output directory")
            os.makedirs(args.outdir)
    '''

    for indir in indirs:
        convert_dir(indir,args)
