#!/usr/bin/env python

"""Compare correction graphs as made by jet_l2_correction_x"""


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
from comparator import Contribution, Plot


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
My_Style.cd()


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


def construct_difference_graph(graph, other_graph):
    x, y = cu.get_xy(graph)
    x_other, y_other = cu.get_xy(other_graph)
    n = graph.GetN()
    if len(y) != len(y_other):
        # If different # points, just use the smaller
        if len(x) < len(x_other):
            n = len(x)
            x_other = x_other[:len(x)]
            y_other = y_other[:len(y)]
        elif len(x) > len(x_other):
            n = len(x_other)
            x = x[:len(x_other)]
            y = y[:len(y_other)]
    if x != x_other:
        raise RuntimeError("x values different")
    diff_y = [y1-y2 for y1, y2 in zip(y, y_other)]
    ex = [0]*n
    ey = [0]*n
    gr = ROOT.TGraphErrors(n, array('d', x), array('d', diff_y), array('d', ex), array('d', ey))
    return gr


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
        graph.SetMarkerSize(entry.get('marker_size', 1))

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
        leg.AddEntry(graph, entry.get('label', graph.GetName()), "LP")

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetTicks(1, 1)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()
    mg.Draw("AP")

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
    cms_text.AddText("CMS")
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

    sample_text = ROOT.TPaveText(0.65, 0.91, 0.67, 0.92, "NDC")
    sample_text.AddText("QCD 13 TeV")
    sample_text.SetTextFont(42)
    sample_text.SetTextSize(FONT_SIZE)
    sample_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    sample_text.SetBorderSize(0)
    sample_text.SetFillStyle(0)
    sample_text.Draw()

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
    parser.add_argument("--input", help="Input ROOT file with correction graphs", action='append')
    parser.add_argument("--label", help="Label for input file", action='append')
    parser.add_argument("--outputDir", help="Output directory for plots", default=os.getcwd())
    parser.add_argument("--title", help="Title string for plots")
    parser.add_argument("--chi2", help="Do Chi2 comparison plot", action='store_true')
    parser.add_argument("--ylim", help="Set explicit y limits", nargs=2)
    args = parser.parse_args(in_args)

    cu.check_dir_exists_create(args.outputDir)

    if len(args.input) != len(args.label):
        raise RuntimeError("Need a --label for each --input")

    if args.input:

        lw = 1
        entry_dicts = [
            {"flav": "AbsCorVsJetPt", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 2, "line_width": lw, "marker_size": 1.2},
            {"flav": "ud_AbsCorVsJetPt", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 1.2},
            {"flav": "s_AbsCorVsJetPt", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.4},
            {"flav": "c_AbsCorVsJetPt", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.4},
            {"flav": "b_AbsCorVsJetPt", "label": "b", "colour": ROOT.kOrange-3, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.6},
            {"flav": "g_AbsCorVsJetPt", "label": "g", "colour": ROOT.kAzure+1, "marker_style": 29, "line_style": 1, "line_width": lw, "marker_size": 1.8},
        ]


        all_dirs = [set(cu.get_list_of_element_names(cu.open_root_file(infile))) for infile in args.input]
        dirs = all_dirs[0]
        for d in all_dirs[1:]:
            dirs = dirs & d
        dirs = sorted(list(dirs))[:1]
        print("Doing: ", dirs)
        # Loop through all different ak4pfchs, etc
        for mydir in dirs:
            jec_text = ROOT.TPaveText(0.17, 0.91, 0.2, 0.92, "NDC")
            # jec_label = "Without JEC"
            jec_label = "With JEC"
            # jec_label = "Summer16_23Sep2016V4"
            # jec_label = "Summer16_03Feb2017_V8"
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

            other_elements = [jec_text, dir_text]

            plot_dir = os.path.join(args.outputDir, mydir)
            cu.check_dir_exists_create(plot_dir)

            obj_list = cu.get_list_of_objects_in_dir(args.input[0], mydir)

            ylimits = (float(args.ylim[0]), float(args.ylim[1])) if args.ylim else None

            X_MIN, X_MAX = 8, 5000
            # For limit protection:
            Y_MIN, Y_MAX = 0.8, 1.6

            # Do all flavs corr vs pt for given eta bin
            common_eta_bins = cu.sort_human(get_common_eta_bins(obj_list))
            for eta_bin in common_eta_bins:
                # Do a per-flavour comparison plot
                for fdict in entry_dicts:
                    entries = []
                    chi2entries = []
                    for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
                        entry = deepcopy(fdict)
                        entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                        entry['label'] += " [%s]" % label
                        entry["line_color"] = fdict['colour']
                        entry["marker_color"] = fdict['colour']
                        if ind == 1:
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                            entry["line_style"] += 1
                        if ind == 2:
                            entry["line_style"] += 1
                        if ind == 3:
                            entry["line_style"] += 1
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                        entries.append(entry)

                        if args.chi2:
                            chi2entry = deepcopy(entry)
                            chi2entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("corrVs", "Chi2NDoFVs"), eta_bin))
                            chi2entries.append(chi2entry)

                    title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
                    do_comparison_graph(entries, title=title,
                                        xtitle="p_{T}^{Gen} [GeV]", ytitle="Correction", logx=True,
                                        xlimits=(X_MIN, X_MAX), y_limit_protection=(Y_MIN, Y_MAX), 
                                        ylimits=ylimits,
                                        other_elements=other_elements,
                                        output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s.pdf" % (eta_bin, fdict['label'])))
                    if args.chi2:
                        do_comparison_graph(chi2entries, title=title,
                                            xtitle="p_{T}^{Gen} [GeV]", ytitle="#chi^{2}/N_{DoF}", logx=True,
                                            xlimits=(X_MIN, X_MAX),
                                            other_elements=other_elements,
                                            output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_chi2.pdf" % (eta_bin, fdict['label'])))

                # Do a plot with all flavours
                entries = []
                ud_entries, g_entries = [], []
                for fdict in entry_dicts:
                    for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
                        entry = deepcopy(fdict)
                        entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                        entry['label'] += " [%s]" % label
                        entry["line_color"] = fdict['colour']
                        entry["marker_color"] = fdict['colour']
                        if ind == 1:
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                            entry["line_style"] += 1
                        if ind == 2:
                            entry["line_style"] += 1
                        if ind == 3:
                            entry["line_style"] += 1
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                        entries.append(entry)
                        if fdict['label'] == "ud":
                            ud_entries.append(entry)
                        elif fdict['label'] == "g":
                            g_entries.append(entry)

                title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
                do_comparison_graph(entries, title=title,
                                    xtitle="p_{T}^{Gen} [GeV]", ytitle="Correction", logx=True,
                                    xlimits=(X_MIN, X_MAX), y_limit_protection=(Y_MIN, Y_MAX),
                                    other_elements=other_elements,
                                    output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_allFlavs.pdf" % (eta_bin)))
                
                # diff_entries = []
                # for ind, (ud, g, label) in enumerate(zip(ud_entries, g_entries, args.label)):
                #     diff_graph = construct_difference_graph(ud['graph'], g['graph'])
                #     diff = deepcopy(ud)
                #     diff['graph'] = diff_graph
                #     diff['label'] = label
                #     diff_entries.append(diff)
                # do_comparison_graph(diff_entries, title=title,
                #                     xtitle="p_{T}^{Gen} [GeV]", ytitle="Correction (ud) - Correction (g)", logx=True,
                #                     xlimits=(X_MIN, X_MAX), ylimits=(0, 0.08),
                #                     other_elements=other_elements,
                #                     output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_ud_g_diff.pdf" % (eta_bin)))
            
            return    

            # Do all flavs corr vs eta for given pt bin
            common_pt_bins = cu.sort_human(get_common_pt_bins(obj_list))
            for pt_bin in common_pt_bins:
                # Do a per-flavour comparison plot
                for fdict in entry_dicts:
                    entries = []
                    chi2entries = []
                    for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
                        entry = deepcopy(fdict)
                        entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin))
                        entry['label'] += " [%s]" % label
                        entry["line_color"] = fdict['colour']
                        entry["marker_color"] = fdict['colour']
                        if ind == 1:
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                        if ind == 2:
                            entry["line_style"] += 1
                        if ind == 3:
                            entry["line_style"] += 1
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                        entries.append(entry)
                        if args.chi2:
                            chi2entry = deepcopy(entry)
                            chi2entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta").replace("corrVs", "Chi2NDoFVs"), pt_bin))
                            chi2entries.append(chi2entry)
                    
                    title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
                    do_comparison_graph(entries, title=title + " GeV",
                                        xtitle="|#eta|", ytitle="Correction",
                                        xlimits=(0, 5.2), y_limit_protection=(Y_MIN, Y_MAX),
                                        other_elements=other_elements, ylimits=ylimits,
                                        output_filename=os.path.join(plot_dir, "compare_corr_vs_eta_%s_%s.pdf" % (pt_bin, fdict['label'])))
                    if args.chi2:
                        do_comparison_graph(entries, title=title,
                                            xtitle="|#eta|", ytitle="#chi^{2}/N_{DoF}",
                                            xlimits=(0, 5.2),
                                            other_elements=other_elements,
                                            output_filename=os.path.join(plot_dir, "compare_corr_vs_eta_%s_%s_chi2.pdf" % (pt_bin, fdict['label'])))

                # Do a plot with all flavours
                entries = []
                for fdict in entry_dicts:
                    for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
                        entry = deepcopy(fdict)
                        entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin))
                        entry['label'] += " [%s]" % label
                        entry["line_color"] = fdict['colour']
                        entry["marker_color"] = fdict['colour']
                        if ind == 1:
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                        if ind == 2:
                            entry["line_style"] += 1
                        if ind == 3:
                            entry["line_style"] += 1
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                        entries.append(entry)
                title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
                do_comparison_graph(entries, title=title + " GeV",
                                    xtitle="|#eta|", ytitle="Correction",
                                    xlimits=(0, 5.2), y_limit_protection=(Y_MIN, Y_MAX),
                                    other_elements=other_elements, ylimits=ylimits,
                                    output_filename=os.path.join(plot_dir, "compare_corr_vs_eta_%s_allFlavs.pdf" % (pt_bin)))


    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
