#!/usr/bin/env python
from __future__ import print_function

import ROOT


def main(file1, file2, root_directory, silent=False):
    """Use Kolmogorov test score to conveniently check for equality-"""
    f1 = ROOT.TFile(file1)
    f2 = ROOT.TFile(file2)
    f1.cd(root_directory)
    scores = []
    n_failed = 0
    for k in ROOT.gDirectory.GetListOfKeys():
        if not k.GetClassName().startswith("TH"):
            continue
        name = root_directory + "/" + k.GetName()
        h1 = f1.Get(name)
        h2 = f2.Get(name)
        for hist, file in [(h1, file1), (h2, file2)]:
            if not bool(hist):
                print("{0} does not exist in {1}!".format(name, file))
                n_failed += 1
                break
        else:
            scores.append((h1.KolmogorovTest(h2, "UO"), name))
    if not silent:
        for score, name in sorted(scores)[::-1]:
            print("{0}: {1:.2f}".format(name, score))
    n_smaller_one = len([s for s in scores if s[0] < 1])
    print("Kolmogorov test: {0}/{1} tests failed.".format(n_smaller_one, len(scores)))
    return n_smaller_one + n_failed


if __name__ == "__main__":
    import argparse
    import sys
    _help = "Compare histograms through Kolmogorov test. "
    _help += "Equal histograms should return exactly 1."
    arg_parser = argparse.ArgumentParser(description=_help)
    arg_parser.add_argument("file1")
    arg_parser.add_argument("file2")
    _help = "E.g. 'config' to compare (only) the calibration of two build files."
    arg_parser.add_argument("-d", "--root_directory", default="config")
    arg_parser.add_argument("--silent", action="store_true")
    args = arg_parser.parse_args()
    sys.exit(main(**vars(args)))
