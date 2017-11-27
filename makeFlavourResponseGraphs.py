#!/Usr/bin/env python

"""Make lots of plots for flavour-binned response.

Uses output ROOT files from
"""

# split by:
# eta binned
# pt binned
# inclusive over pt/eta
#
# want graphs, and comparison fits/hists
#
# do same for resolution?
#

import os
import sys
import argparse
from copy import deepcopy
import re

import ROOT
from MyStyle import My_Style

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
My_Style.cd()



def check_file_exists(filepath):
    """Check if file exists. Can do absolute or relative file paths."""
    return os.path.isfile(os.path.realpath(filepath))


def check_dir_exists_create(filepath):
    """Check if directory exists. If not, create it."""
    if not os.path.isdir(filepath):
        os.makedirs(os.path.realpath(filepath))

#
# ROOT specific fns, like opening files safely
#
def open_root_file(filename, mode="READ"):
    """Safe way to open ROOT file. Could be improved."""
    if mode in ["READ", "UPDATE"]:
        if not check_file_exists(filename):
            raise IOError("No such file %s" % filename)
    f = ROOT.TFile(filename, mode)
    if f.IsZombie() or not f:
        raise IOError("Can't open TFile %s" % filename)
    return f


def exists_in_file(tfile, obj_name):
    """Check if object exists in TFile.

    Also handles directory structure, e.g. histograms/central/pt_1
    """
    parts = obj_name.split("/")
    current_obj = tfile
    for p in parts:
        if current_obj.GetListOfKeys().Contains(p):
            current_obj = current_obj.Get(p)
        else:
            return False
    return True


def get_from_tfile(tfile, obj_name, info=False):
    """Get some object from ROOT TFile with checks."""
    if info:
        print "Getting %s" % obj_name
    if not exists_in_file(tfile, obj_name):
        raise IOError("No object named %s in %s" % (obj_name, tfile.GetName()))
    else:
        return tfile.Get(obj_name)


def grab_obj_from_file(file_name, obj_name):
    """Get object names obj_name from ROOT file file_name"""
    input_file = open_root_file(file_name)
    obj = get_from_tfile(input_file, obj_name)
    if isinstance(obj, (ROOT.TH1)):
        obj.SetDirectory(0)  # Ownership kludge
        input_file.Close()
        return obj.Clone(ROOT.TUUID().AsString())
    else:
        return obj


def get_list_of_element_names(thing):
    """Get list of key names for given thing in a ROOT file"""
    return [x.GetName() for x in thing.GetListOfKeys()]


def get_list_of_objects_in_dir(filename, dirname):
    """Get list of elements in a TDirectory"""
    f = open_root_file(filename)  # need to keep this in memory otherwise the get_list... segfaults
    d = get_from_tfile(f, dirname)
    return get_list_of_element_names(d)


def get_common_eta_bins(obj_list):
    """Get list of common eta bins from list of object names"""
    eta_bins = []
    for x in obj_list:
        m = re.search(r'JetEta[0-9.]+to[0-9.]+', x)
        if m:
            eta_bins.append(m.group(0))
    return list(set(eta_bins))


def get_common_pt_bins(obj_list):
    """Get list of common pt bins from list of object names"""
    pt_bins = []
    for x in obj_list:
        m = re.search(r'RefPt[0-9.]+to[0-9.]+', x)
        if m:
            pt_bins.append(m.group(0))
    return list(set(pt_bins))


def do_comparison_graph(entries, output_filename, title="", xtitle="", ytitle="", other_elements=None, logx=False, logy=False, do_line=True):
    """Draw several graphs on one canvas and save to file

    Parameters
    ----------
    entries : [dict]
        List of entries to plot. Each is represented by a dict, with the graph,
        label, and various other style options to be applied.
    output_filename : str
        Name of output plot file
    title : str
        Overall plot title
    xtitle : str
        X axis label
    ytitle : str
        y axis label
    other_elements : [TObject], optional
        Other elements to Draw on the canvas
    logx : bool, optional
        Log x axis
    logy : bool, optional
        Log y axis
    do_line : bool, optional
        Do horizontal line at 1

    """
    mg = ROOT.TMultiGraph()
    mg.SetTitle(";".join([title, xtitle, ytitle]))
    leg = ROOT.TLegend(0.7, 0.55, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    for entry in entries:
        entry['graph'].SetLineColor(entry.get('line_color', ROOT.kBlack))
        entry['graph'].SetLineStyle(entry.get('line_style', 1))
        entry['graph'].SetLineWidth(entry.get('line_width', 1))

        entry['graph'].SetMarkerColor(entry.get('marker_color', ROOT.kBlack))
        entry['graph'].SetMarkerStyle(entry.get('marker_style', 1))
        entry['graph'].SetMarkerSize(entry.get('marker_size', 1))

        mg.Add(entry['graph'])
        leg.AddEntry(entry['graph'], entry.get('label', entry['graph'].GetName()), "LP")

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetTicks(1, 1)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()
    mg.Draw("ALP")
    leg.Draw()
    if do_line:
        x_min, x_max = mg.GetXaxis().GetXmin(), mg.GetXaxis().GetXmax()
        line = ROOT.TLine(x_min, 1, x_max, 1)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw()
    if other_elements:
        for ele in other_elements:
            ele.Draw()
    canv.SaveAs(output_filename)


def main(in_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--outputDir", help="Output directory for plots", default=os.getcwd())
    parser.add_argument("--inputFits", help="Input ROOT file with response histograms & fits (from jet_response_fitter_x)")
    parser.add_argument("--inputGraphs", help="Input ROOT file with response & resolution graphs (from jet_response_and_resolution_x)")
    args = parser.parse_args(in_args)

    check_dir_exists_create(args.outputDir)

    if args.inputGraphs:

        lw = 2
        entry_dicts = [
            {"flav": "RelRspVsRefPt", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 2, "line_width": lw, "marker_size": 1.2},
            {"flav": "ud_RspVsRefPt_RelRsp", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 1.2},
            {"flav": "s_RRspVsRefPt_RelRsp", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.2},
            {"flav": "c_RRspVsRefPt_RelRsp", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.2},
            {"flav": "b_RRspVsRefPt_RelRsp", "label": "b", "colour": ROOT.kMagenta, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.6},
            {"flav": "g_RRspVsRefPt_RelRsp", "label": "g", "colour": ROOT.kAzure+1, "marker_style": ROOT.kFullCrossX, "line_style": 1, "line_width": lw, "marker_size": 1.6},
        ]

        dirs = get_list_of_element_names(open_root_file(args.inputGraphs))
        for mydir in dirs:
            # mydir = "ak4pfchsl1l2l3"
            plot_dir = os.path.join(args.outputDir, mydir)
            check_dir_exists_create(plot_dir)

            obj_list = get_list_of_objects_in_dir(args.inputGraphs, mydir)

            # Do all flavs for given eta bins
            common_eta_bins = get_common_eta_bins(obj_list)
            for eta_bin in common_eta_bins:
                entries = []
                for fdict in entry_dicts:
                    entry = deepcopy(fdict)
                    entry["graph"] = grab_obj_from_file(args.inputGraphs, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                    entry["line_color"] = fdict['colour']
                    entry["marker_color"] = fdict['colour']
                    entries.append(entry)
                title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
                do_comparison_graph(entries, title=title,
                                    xtitle="p_{T}^{Gen} [GeV]", ytitle="Response", logx=True,
                                    output_filename=os.path.join(plot_dir, "rsp_vs_pt_%s.pdf" % (eta_bin)))

            # Do all flavs for given pt bins
            common_pt_bins = get_common_pt_bins(obj_list)
            for pt_bin in common_pt_bins:
                entries = []
                for fdict in entry_dicts:
                    entry = deepcopy(fdict)
                    entry["graph"] = grab_obj_from_file(args.inputGraphs, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin))
                    entry["line_color"] = fdict['colour']
                    entry["marker_color"] = fdict['colour']
                    entries.append(entry)
                title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
                do_comparison_graph(entries, title=title + " GeV", xtitle="|#eta|", ytitle="Response",
                                    output_filename=os.path.join(plot_dir, "rsp_vs_eta_%s.pdf" % (pt_bin)))

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
