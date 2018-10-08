#!/usr/bin/env python

"""Compare + and - eta correction graphs as made by jet_l2_correction_x"""


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


FONT_SIZE = 0.032


def make_pos_neg_eta_pairs(common_eta_bins):
    pos_bins = []
    for entry in common_eta_bins:
        if "-" not in entry:
            pos_bins.append(entry)

    pairs = []
    for pos_entry in pos_bins:
        m = re.search(r'JetEta([0-9.]+)to([0-9.]+)', pos_entry)
        if m:
            eta_low = m.group(1)
            eta_high = m.group(2)
            eta_high_str = "-" + eta_high
            eta_low_str = eta_low
            if eta_low != "0":
                eta_low_str = "-" + eta_low
            neg_entry = "JetEta%sto%s" % (eta_high_str, eta_low_str)
            if neg_entry in common_eta_bins:
                pairs.append((pos_entry, neg_entry))
            else:
                raise RuntimeError("Cannot find corresponding bin %s", neg_entry)
    return pairs


def do_comparison_graph(entries, output_filename, title="", xtitle="", ytitle="",
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
    delta = 0.21
    middle = 0.69
    leg = ROOT.TLegend(middle-delta, 0.75, middle+delta, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if len(entries) > 4:
        leg.SetNColumns(2)
    leg.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignCenter)

    # VITAL to stop ROOT destroying all the graphs inside the TMultiGraph as well
    ROOT.SetOwnership(mg, False)

    for ind, entry in enumerate(entries):
        default_colour = ROOT.kBlack
        graph = entry['graph']
        graph.SetLineColor(entry.get('line_color', default_colour))
        graph.SetLineStyle(entry.get('line_style', 1))
        graph.SetLineWidth(entry.get('line_width', 1))

        graph.SetMarkerColor(entry.get('marker_color', default_colour))
        graph.SetMarkerStyle(entry.get('marker_style', 1))
        graph.SetMarkerSize(entry.get('marker_size', 1)*0.5)

        graph.SetFillColorAlpha(entry.get('fill_color', default_colour), entry.get('fill_alpha', 1))
        graph.SetFillStyle(entry.get('fill_style', 1001))

        if graph.GetListOfFunctions().GetSize() > 0:
            func = graph.GetListOfFunctions().Last()
            if draw_fits:
                func.SetLineColor(entry.get('line_color', default_colour))
                func.SetLineStyle(entry.get('line_style', 1))
                # func.SetLineStyle(3)
                func.SetLineWidth(2)
            else:
                # Broken for now
                graph.GetListOfFunctions().Remove(func)

        mg.Add(graph)
        leg.AddEntry(graph, entry.get('label', graph.GetName()), "LFP")

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetTicks(1, 1)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()
    mg.Draw("AP2")

    # Little extra breathing room
    mg.GetHistogram().SetMaximum(mg.GetYaxis().GetXmax() * 1.1)

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

    if title:
        bin_text = ROOT.TPaveText(0.17, 0.8, 0.2, 0.81, "NDC")
        bin_text.AddText(title)
        bin_text.SetTextFont(42)
        bin_text.SetTextSize(FONT_SIZE)
        bin_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        bin_text.SetBorderSize(0)
        bin_text.SetFillStyle(0)
        bin_text.Draw()

    # sample_text = ROOT.TPaveText(0.65, 0.91, 0.67, 0.92, "NDC")
    # sample_text.AddText("QCD 13 TeV")
    # sample_text.SetTextFont(42)
    # sample_text.SetTextSize(FONT_SIZE)
    # sample_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    # sample_text.SetBorderSize(0)
    # sample_text.SetFillStyle(0)
    # sample_text.Draw()

    if other_elements:
        for ele in other_elements:
            ele.Draw()

    canv.SaveAs(output_filename)


def get_open_marker(marker):
    opposites = {
        20: 24,
        21: 25,
        22: 26,
        23: 32,
        33: 27,
        29: 30,
        34: 28,
    }
    return opposites[marker]


def main(in_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input ROOT file with correction graphs", required=True)
    parser.add_argument("--outputDir", help="Output directory for plots", default=os.getcwd())
    parser.add_argument("--title", help="Title string for plots", default="")
    parser.add_argument("--sampleName", help="Sample name string for plots", default="")
    args = parser.parse_args(in_args)

    cu.check_dir_exists_create(args.outputDir)

    lw = 1
    entry_dicts = [
        {"flav": "AbsCorVsJetPt", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 2, "line_width": lw, "marker_size": 1.2},
        {"flav": "ud_AbsCorVsJetPt", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 1.2},
        {"flav": "s_AbsCorVsJetPt", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "c_AbsCorVsJetPt", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "b_AbsCorVsJetPt", "label": "b", "colour": ROOT.kOrange-3, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.6},
        {"flav": "g_AbsCorVsJetPt", "label": "g", "colour": ROOT.kAzure+1, "marker_style": 29, "line_style": 1, "line_width": lw, "marker_size": 1.8},
    ]


    all_dirs = cu.get_list_of_element_names(cu.open_root_file(args.input))
    dirs = sorted(list(all_dirs))[:1]
    print("Doing: ", dirs)
    # Loop through all different ak4pfchs, etc
    for mydir in dirs:
        jec_text = ROOT.TPaveText(0.15, 0.91, 0.2, 0.92, "NDC")
        jec_text.AddText(args.title)
        jec_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        jec_text.SetTextFont(42)
        jec_text.SetTextSize(FONT_SIZE)
        jec_text.SetBorderSize(0)
        jec_text.SetFillStyle(0)

        dir_text = ROOT.TPaveText(0.17, 0.76, 0.2, 0.77, "NDC")
        dir_label = mydir.upper().replace("PFCHS", " PF CHS").replace("PUPPI", " PUPPI").replace("L1L2L3", " + L1L2L3")
        dir_text.AddText(dir_label)
        dir_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        dir_text.SetTextFont(42)
        dir_text.SetTextSize(FONT_SIZE)
        dir_text.SetBorderSize(0)
        dir_text.SetFillStyle(0)

        sample_text = ROOT.TPaveText(0.93, 0.91, 0.94, 0.92, "NDC")
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

        # ylimits = (float(args.ylim[0]), float(args.ylim[1])) if args.ylim else None
        ylimits = None

        X_MIN, X_MAX = 8, 5000
        # For limit protection:
        Y_MIN, Y_MAX = 0.8, 1.6

        # Do all flavs corr vs pt for given eta bin
        common_eta_bins = cu.sort_human(cu.get_common_eta_bins(obj_list))
        # Find the matching + and - eta bins
        eta_bin_pairs = make_pos_neg_eta_pairs(common_eta_bins)
        print(eta_bin_pairs)

        for eta_bin_pair in eta_bin_pairs:
            all_entries = []
            
            title = eta_bin_pair[0].replace("to", " < |#eta| < ").replace("JetEta", "")
            
            # Do a per-flavour comparison plot
            for fdict in entry_dicts:
                entries = []
                for ind, (eta_bin, label) in enumerate(zip(eta_bin_pair, ['+#eta', '-#eta'])):
                    entry = deepcopy(fdict)
                    entry["graph"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                    entry['label'] += " [%s]" % label
                    entry["line_color"] = fdict['colour']+ind
                    entry["marker_color"] = fdict['colour']+ind
                    entry["fill_color"] = fdict['colour']+ind
                    entry["fill_style"] = 1001
                    entry["fill_alpha"] = 0.7
                    if ind == 1:
                        entry["marker_style"] = get_open_marker(entry['marker_style'])
                        entry["line_style"] += 1
                    if ind == 2:
                        entry["line_style"] += 1
                    if ind == 3:
                        entry["line_style"] += 1
                        entry["marker_style"] = get_open_marker(entry['marker_style'])
                    entries.append(entry)
                    all_entries.append(entry)

                do_comparison_graph(entries, title=title,
                                    xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                    xlimits=(X_MIN, X_MAX), y_limit_protection=(Y_MIN, Y_MAX), 
                                    ylimits=ylimits,
                                    other_elements=other_elements,
                                    output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_pos_neg_%s.pdf" % (eta_bin_pair[0], fdict['label'])))

            # Do a plot with all flavours
            do_comparison_graph(all_entries, title=title,
                                xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                xlimits=(X_MIN, X_MAX), y_limit_protection=(Y_MIN, Y_MAX),
                                other_elements=other_elements,
                                output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_pos_neg_allFlavs.pdf" % (eta_bin_pair[0])))
            
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
