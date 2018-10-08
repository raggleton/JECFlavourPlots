#!/usr/bin/env python

"""Extract only the respone hists from a ROOT file, save to new ROOT file"""


from __future__ import print_function
import argparse
import os
import sys
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


def get_list_of_objects(tobj):
    return [x.GetName() for x in tobj.GetListOfKeys()]


def extract_response_hists(input_filename, output_filename):
    inf = ROOT.TFile(input_filename)
    outf = ROOT.TFile(output_filename, "RECREATE")
    print("Writing output file", output_filename)
    dirs = get_list_of_objects(inf)
    name_starts = ("ud_", "b_", "c_", "s_", "g_", "RelRsp")
    for dirname in dirs:
        print(dirname)
        tdir = inf.Get(dirname)
        outf.mkdir(dirname)
        outdir = outf.Get(dirname)
        outdir.cd()
        hists = get_list_of_objects(tdir)
        hists = [h for h in hists if "RelRsp_JetEta" in h and h.startswith(name_starts)]
        counter = 0
        for hname in hists:
            # print(hname)
            h = tdir.Get(hname)
            if h.GetEntries() > -1:
                h.Write()
                counter += 1
        print("Copying", counter, "histograms")


def main(in_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("input",
                        help="Input ROOT files produced by jet_response_analyser_x", 
                        nargs='*')
    args = parser.parse_args(in_args)

    for filename in args.input:
        if not os.path.isfile(filename):
            raise IOError("Cannot find input file %s" % filename)
        parts = os.path.splitext(filename)
        output_filename = parts[0]+"_responseOnly"+parts[1]
        extract_response_hists(filename, output_filename)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))