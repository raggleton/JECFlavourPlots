#!/usr/bin/env python

"""Plot response & correciton graphs using ROOT file from output of jet_l2_correction_x"""


from __future__ import print_function
import os
import sys
import argparse
from copy import deepcopy
import re
from array import array
import numpy as np

import ROOT
from MyStyle import My_Style
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
My_Style.cd()

OPT_FIT_STYLE = 1111
ROOT.gStyle.SetOptFit(OPT_FIT_STYLE)

FONT_SIZE = 0.032


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


def construct_inverse_graph(graph):
    """Construct graph with y values = 1/y of input graph"""
    n = graph.GetN()
    x, y = cu.get_xy(graph)
    # y = np.ndarray(n, 'd', graph.GetY())
    new_y = array('d', [1/old_y if old_y != 0 else 0 for old_y in y])
    ex = array('d', [0] * n)
    ey = array('d', [0] * n)
    gr = ROOT.TGraphErrors(n, x, new_y, ex, ey)
    return gr


def do_comparison_graph(entries, output_filename, bin_title="", xtitle="", ytitle="",
                        other_elements=None, logx=False, logy=False,
                        do_line=True, xlimits=None, ylimits=None,
                        y_limit_protection=None, draw_fits=True):
    """Draw several graphs on one canvas and save to file

    Parameters
    ----------
    entries : [dict]
        List of entries to plot. Each is represented by a dict, with the graph,
        label, and various other style options to be applied.
    output_filename : str
        Name of output plot file
    bin_title : str
        Bin title e.g. x < pT < y
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
    xlimits : (min, max), optional
        Set hard x limits
    ylimits : (min, max), optional
        Set hard y limits
    y_limit_protection : (y_min, y_max), optional
        Set minimum and maximum y values in the event of a huge stat error or weird point
    draw_fits : bool, optional
        Draw fitted functions or not

    """
    mg = ROOT.TMultiGraph()
    mg.SetTitle(";".join(["", xtitle, ytitle]))
    delta = 0.12
    middle = 0.77
    leg = ROOT.TLegend(middle-delta, 0.75, middle+delta, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)
    leg.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignCenter)
    
    # VITAL to stop ROOT destroying all the graphs inside the TMultiGraph as well
    ROOT.SetOwnership(mg, False)

    if len(entries) == 1:
        ROOT.gStyle.SetOptFit(OPT_FIT_STYLE)
    else:
        ROOT.gStyle.SetOptFit(0)

    for entry in entries:
        default_colour = ROOT.kBlack
        graph = entry['graph']
        graph.SetLineColor(entry.get('line_color', default_colour))
        graph.SetLineStyle(entry.get('line_style', 1))
        graph.SetLineWidth(entry.get('line_width', 1))

        graph.SetMarkerColor(entry.get('marker_color', default_colour))
        graph.SetMarkerStyle(entry.get('marker_style', 1))
        graph.SetMarkerSize(entry.get('marker_size', 1))

        if graph.GetListOfFunctions().GetSize() > 0:
            func = graph.GetListOfFunctions().Last()
            if draw_fits:
                func.SetLineColor(entry.get('line_color', default_colour))
                func.SetLineStyle(4)
                func.SetLineWidth(2)
            else:
                # Broken for now
                graph.GetListOfFunctions().Remove(func)

        c_dummy = ROOT.TCanvas("dummy"+ROOT.TUUID().AsString(), "", 800, 800)
        graph.Draw("ALP")
        ROOT.gPad.Modified()
        ROOT.gPad.Update()
        fit_stats = graph.FindObject("stats")
        fit_stats.SetBorderSize(0)
        fit_stats.SetFillStyle(0)
        fit_stats.SetFillColor(ROOT.kWhite)
        fit_stats.SetX1NDC(0.55)
        fit_stats.SetX2NDC(0.88)
        fit_stats.SetY1NDC(0.55)
        fit_stats.SetY2NDC(0.75)
        if len(entries) > 1:
            # This is an awful hack to not show the fit stats boxes,
            # but I don't know how to do it otherwise
            fit_stats.SetTextColorAlpha(ROOT.kWhite, 0)
        fit_stats.SetOptStat(0)  # Doesn't work
        entry['stats'] = fit_stats
        
        mg.Add(graph)
        leg.AddEntry(graph, entry.get('label', graph.GetName()), "LP")
    
    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetTicks(1, 1)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()

    mg.Draw("ALP")

    # Little extra breathing room
    mg.GetHistogram().SetMaximum(mg.GetYaxis().GetXmax() * 1.05)

    # Protection in case y limits are dominated by large stat error
    if y_limit_protection and len(y_limit_protection) == 2:
        y_min, y_max = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
        y_lim_lower, y_lim_upper = y_limit_protection
        if y_max > y_lim_upper:
            mg.GetHistogram().SetMaximum(y_lim_upper)
        if y_min < y_lim_lower:
            mg.GetHistogram().SetMinimum(y_lim_lower)

    if xlimits and len(xlimits) == 2:
        mg.GetXaxis().SetLimits(*xlimits)

    if ylimits and len(ylimits) == 2:
        mg.GetHistogram().SetMaximum(ylimits[1])
        mg.GetHistogram().SetMinimum(ylimits[0])

    leg.Draw()

    if do_line:
        y_min, y_max = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
        if y_min < 1 and y_max > 1:
            x_min, x_max = mg.GetXaxis().GetXmin(), mg.GetXaxis().GetXmax()
            line = ROOT.TLine(x_min, 1, x_max, 1)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray+2)
            line.Draw()

    cms_text = ROOT.TPaveText(0.17, 0.84, 0.2, 0.85, "NDC")
    cms_text.AddText("CMS Simulation")
    cms_text.SetTextFont(62)
    cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_text.SetTextSize(FONT_SIZE)
    cms_text.SetBorderSize(0)
    cms_text.SetFillStyle(0)
    cms_text.Draw()

    bin_text = ROOT.TPaveText(0.17, 0.8, 0.2, 0.81, "NDC")
    bin_text.AddText(bin_title)
    bin_text.SetTextFont(42)
    bin_text.SetTextSize(FONT_SIZE)
    bin_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    bin_text.SetBorderSize(0)
    bin_text.SetFillStyle(0)
    bin_text.Draw()

    if other_elements:
        for ele in other_elements:
            ele.Draw()

    canv.SaveAs(output_filename)


def main(in_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input ROOT file with response & resolution graphs (from jet_response_and_resolution_x)", required=True)
    parser.add_argument("--outputDir", help="Output directory for plots", default=os.getcwd())
    parser.add_argument("--title", help="Title string for plots", default="")
    parser.add_argument("--sampleName", help="Sample name string for plots", default="")
    args = parser.parse_args(in_args)

    cu.check_dir_exists_create(args.outputDir)

    lw = 1
    entry_dicts = [
        {"flav": "AbsCorVsJetPt", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 1, "line_width": lw, "marker_size": 1.2},
        {"flav": "ud_AbsCorVsJetPt", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 1.2},
        {"flav": "s_AbsCorVsJetPt", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "c_AbsCorVsJetPt", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "b_AbsCorVsJetPt", "label": "b", "colour": ROOT.kOrange-3, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.6},
        {"flav": "g_AbsCorVsJetPt", "label": "g", "colour": ROOT.kAzure+1, "marker_style": 29, "line_style": 1, "line_width": lw, "marker_size": 1.8},
    ]

    # Loop through all different ak4pfchs, etc
    dirs = cu.get_list_of_element_names(cu.open_root_file(args.input))
    for mydir in dirs[:]:
        jec_text = ROOT.TPaveText(0.14, 0.91, 0.2, 0.92, "NDC")
        jec_text.AddText(args.title)
        jec_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        jec_text.SetTextFont(42)
        jec_text.SetTextSize(FONT_SIZE)
        jec_text.SetBorderSize(0)
        jec_text.SetFillStyle(0)

        dir_text = ROOT.TPaveText(0.17, 0.75, 0.2, 0.77, "NDC")
        dir_label = mydir.upper().replace("PFCHS", " PF CHS").replace("PUPPI", " PUPPI").replace("L1L2L3", " + L1L2L3")
        dir_text.AddText(dir_label)
        dir_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        dir_text.SetTextFont(42)
        dir_text.SetTextSize(FONT_SIZE)
        dir_text.SetBorderSize(0)
        dir_text.SetFillStyle(0)

        sample_text = ROOT.TPaveText(0.89, 0.91, 0.9, 0.92, "NDC")
        # sample_text.AddText("Flat QCD 13 TeV")
        sample_text.AddText(args.sampleName + " 13 TeV")
        sample_text.SetTextFont(42)
        sample_text.SetTextSize(FONT_SIZE)
        sample_text.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        sample_text.SetBorderSize(0)
        sample_text.SetFillStyle(0)
        sample_text.Draw()

        other_elements = [jec_text, dir_text, sample_text]

        plot_dir = os.path.join(args.outputDir, mydir)
        cu.check_dir_exists_create(plot_dir)

        obj_list = cu.get_list_of_objects_in_dir(args.input, mydir)

        # Do all flavs rsp vs pt for given eta bin
        common_eta_bins = get_common_eta_bins(obj_list)
        for eta_bin in common_eta_bins:
            entries = []
            
            # repsonse
            # for fdict in entry_dicts:
            #     entry = deepcopy(fdict)
            #     entry["graph"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
            #     entry["line_color"] = fdict['colour']
            #     entry["marker_color"] = fdict['colour']
            #     entries.append(entry)
            # bin_title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
            # do_comparison_graph(entries, bin_title=bin_title,
            #                     xtitle="p_{T}^{Gen} [GeV]", ytitle="Response", logx=True,
            #                     xlimits=(10, 5000), y_limit_protection=(0.5, 1.5),
            #                     other_elements=other_elements,
            #                     output_filename=os.path.join(plot_dir, "rsp_vs_pt_%s.pdf" % (eta_bin)))

            # correction
            entries = []
            bin_title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
            for fdict in entry_dicts:
                entry = deepcopy(fdict)
                entry["graph"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                entry["line_color"] = fdict['colour']
                entry["marker_color"] = fdict['colour']
                entries.append(entry)
                do_comparison_graph([entry], bin_title=bin_title,
                                    xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                    y_limit_protection=(0.7, 1.8),
                                    other_elements=other_elements,
                                    output_filename=os.path.join(plot_dir, "%s_%s_logX.pdf" % (fdict['flav'], eta_bin)))

                do_comparison_graph([entry], bin_title=bin_title,
                                    xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=False,
                                    y_limit_protection=(0.7, 1.8),
                                    other_elements=other_elements,
                                    output_filename=os.path.join(plot_dir, "%s_%s_linX.pdf" % (fdict['flav'], eta_bin)))

            do_comparison_graph(entries, bin_title=bin_title,
                                xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                y_limit_protection=(0.7, 1.8),
                                other_elements=other_elements,
                                output_filename=os.path.join(plot_dir, "corr_vs_pt_%s_logX.pdf" % (eta_bin)))

            do_comparison_graph(entries, bin_title=bin_title,
                                xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=False,
                                y_limit_protection=(0.7, 1.8),
                                other_elements=other_elements,
                                output_filename=os.path.join(plot_dir, "corr_vs_pt_%s_linX.pdf" % (eta_bin)))

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
