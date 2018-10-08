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
from math import sqrt, fabs

import ROOT
from MyStyle import My_Style
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
My_Style.cd()


FONT_SIZE = 0.032


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


def construct_graph_func_diff_graph(graph, func):
    """Construct a graph of function - graph for each (x, y) point in graph"""
    x, y = cu.get_xy(graph)
    n = len(x)
    diff = array('d', [func.Eval(this_x) - this_y for this_x, this_y in zip(x, y)])
    ex = array('d', [0] * n)
    ey = array('d', [0] * n)
    gr = ROOT.TGraphErrors(n, x, diff, ex, ey)
    return gr


def construct_graph_func_ratio_graph(graph, func):
    """Construct a graph of function / graph for each (x, y) point in graph"""
    x, y = cu.get_xy(graph)
    n = len(x)
    diff = array('d', [func.Eval(this_x) / this_y for this_x, this_y in zip(x, y)])
    ex = array('d', [0] * n)
    ey = array('d', [0] * n)
    gr = ROOT.TGraphErrors(n, x, diff, ex, ey)
    return gr


def construct_graph_ratio(graph, ref_graph):
    """Contruct a TGraph that is the ratio of the graph : ref_graph

    Interpolation will take place since x values may not align.
    The x values from ref_graph will be used as the reference values.
    """
    x, y = cu.get_xy(graph)
    ex, ey = cu.get_exey(graph)
    x_ref, y_ref = cu.get_xy(ref_graph)
    ex_ref, ey_ref = cu.get_exey(ref_graph)
    n = graph.GetN()
    # handle different length arrays
    if len(y) != len(y_ref):
        # If different # points, just use the smaller
        if len(x) < len(x_ref):
            n = len(x)
            x_ref = x_ref[:len(x)]
            y_ref = y_ref[:len(y)]
            ex_ref = ex_ref[:len(y)]
            ey_ref = ey_ref[:len(y)]
        elif len(x) > len(x_ref):
            n = len(x_ref)
            x = x[:len(x_ref)]
            y = y[:len(y_ref)]
            ex = ex[:len(x_ref)]
            ey = ey[:len(y_ref)]
    # if x != x_ref:
        # raise RuntimeError("x values different")
    # Interpolate
    y = [graph.Eval(xi) for xi in x_ref]
    ratio_y = [y1/y2 for y1, y2 in zip(y, y_ref)]
    ex = [0]*n
    ey = [(iy/iy_ref)*sqrt(pow(iey/iy, 2) + pow(iey_ref/iy_ref, 2)) for iy, iey, iy_ref, iey_ref in zip(y, ey, y_ref, ey_ref)]
    gr = ROOT.TGraphErrors(n, array('d', x_ref), array('d', ratio_y), array('d', ex), array('d', ey))
    return gr


def construct_func_ratio(func, ref_func):
    """Contruct a TF1 that is the ratio of the func : ref_func TF1s"""
    func_formula = func.GetExpFormula("P")
    ref_func_formula = ref_func.GetExpFormula("P")
    xmin = max([func.GetXmin(), ref_func.GetXmin()])
    xmax = min([func.GetXmax(), ref_func.GetXmax()])
    new_func = ROOT.TF1("ratio_%s_%s" % (func.GetName(), ref_func.GetName()),
                        "(%s)/(%s)" % (func_formula, ref_func_formula),
                        xmin, xmax)
    return new_func


def rescale_plot_labels(container, factor):
    # What a pile of poop, why does ROOT scale all these sizes?
    container.GetXaxis().SetLabelSize(container.GetXaxis().GetLabelSize()/factor)
    container.GetXaxis().SetTitleSize(container.GetXaxis().GetTitleSize()/factor)
    container.GetXaxis().SetTitleOffset(container.GetXaxis().GetTitleOffset()*factor)  # doesn't seem to work?
    container.GetXaxis().SetTickLength(container.GetXaxis().GetTickLength()/factor)

    container.GetYaxis().SetLabelSize(container.GetYaxis().GetLabelSize()/factor)
    container.GetYaxis().SetTitleSize(container.GetYaxis().GetTitleSize()/factor)
    container.GetYaxis().SetTitleOffset(container.GetYaxis().GetTitleOffset()*factor)
    # container.GetYaxis().SetTickLength(0.03/factor)


def do_comparison_graph(entries, output_filename, title="", xtitle="", ytitle="",
                        other_elements=None, logx=False, logy=False,
                        do_line=False, xlimits=None, ylimits=None,
                        y_limit_protection=None, draw_fits=True,
                        vertical_line=None,
                        do_fit_graph_ratio=False, 
                        do_ratio_plots=False, ratio_title=""):
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
    vertical_line : int, optional
        Draw a vertical dashed line at this x value
    do_fit_graph_ratio : bool, optional
        Draw subplot of fit to graph ratio for each entry
    do_ratio_plots : bool, optional
        Draw subplot of ratio of graph (& its fit) to some reference object
        that should be defined in each entry of `entries` under `ratio` key
        Cannot use do_fit_graph_ratio and do_ratio_plots simultaneously
    ratio_title : str, optional
        Y axis title for ratio subplot. Only used if do_ratio_plots is True
    
    Raises
    ------
    RuntimeError
        If you specify both do_fit_graph_ratio and do_ratio_plots to be True
    
    """
    if do_fit_graph_ratio and do_ratio_plots:
        raise RuntimeError("Cannot use both do_fit_graph_ratio and do_ratio_plots - choose only one")

    mg = ROOT.TMultiGraph()
    mg.SetTitle(";".join(["", xtitle, ytitle]))
    delta = 0.21
    middle = 0.77
    leg = ROOT.TLegend(middle-delta, 0.7, middle+delta, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if len(entries) > 4:
        leg.SetNColumns(2)
    leg.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignCenter)

    # VITAL to stop ROOT destroying all the graphs inside the TMultiGraph as well
    ROOT.SetOwnership(mg, False)

    subplot_entries = []
    func_mins, func_maxs = [], []
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
            # func = graph.FindObject("fit")
            # if type(func) != ROOT.TF1:
            #     raise RuntimeError("func has type " + str(type(func)))

            func = graph.GetListOfFunctions().Last()
            if draw_fits:
                func.SetLineColor(entry.get('line_color', default_colour))
                func.SetLineStyle(entry.get('line_style', 1))
                # func.SetLineStyle(3)
                func.SetLineWidth(2)
                # func_mins.append(func.GetMinimum())
                # func_maxs.append(func.GetMaximum())
            else:
                # Broken for now
                graph.GetListOfFunctions().Remove(func)

            if do_fit_graph_ratio:
                diff_graph = construct_graph_func_ratio_graph(graph, func)
                diff_graph.SetLineColor(entry.get('line_color', default_colour))
                diff_graph.SetLineStyle(entry.get('line_style', 1))
                diff_graph.SetLineWidth(entry.get('line_width', 1))
                diff_graph.SetMarkerColor(entry.get('marker_color', default_colour))
                diff_graph.SetMarkerStyle(entry.get('marker_style', 1))
                diff_graph.SetMarkerSize(entry.get('marker_size', 1))
                diff_entry = deepcopy(entry)
                diff_entry['graph'] = diff_graph
                subplot_entries.append(diff_entry)
            elif do_ratio_plots:
                ref_graph = entry.get('ratio', None)
                if ref_graph:
                    ratio_graph = construct_graph_ratio(graph, ref_graph)
                    ratio_graph.SetLineColor(entry.get('line_color', default_colour))
                    ratio_graph.SetLineStyle(entry.get('line_style', 1))
                    ratio_graph.SetLineWidth(entry.get('line_width', 1))
                    ratio_graph.SetMarkerColor(entry.get('marker_color', default_colour))
                    ratio_graph.SetMarkerStyle(entry.get('marker_style', 1))
                    ratio_graph.SetMarkerSize(entry.get('marker_size', 1))
                    
                    if ref_graph.GetListOfFunctions().GetSize() > 0:
                        ref_func = ref_graph.GetListOfFunctions().Last()
                        ratio_func = construct_func_ratio(func, ref_func)
                        ratio_func.SetLineColor(entry.get('line_color', default_colour))
                        ratio_func.SetLineStyle(entry.get('line_style', 1))
                        ratio_func.SetLineWidth(entry.get('line_width', 2))
                        ratio_graph.GetListOfFunctions().AddLast(ratio_func)
                    ratio_entry = deepcopy(entry)
                    ratio_entry['graph'] = ratio_graph
                    subplot_entries.append(ratio_entry)

        mg.Add(graph)
        leg.AddEntry(graph, entry.get('label', graph.GetName()), "LP")

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    # canv.SetRightMargin(0.03)
    right_margin = 0.03
    subplot_pad_height = 0.32
    subplot_pad_fudge = 0.01
    main_pad, subplot_pad = None, None
    if do_fit_graph_ratio or do_ratio_plots:
        main_pad = ROOT.TPad("main_pad", "", 0, subplot_pad_height+subplot_pad_fudge, 1, 1)
        ROOT.SetOwnership(main_pad, False)
        main_pad.SetTicks(1, 1)
        main_pad.SetBottomMargin(2*subplot_pad_fudge)
        main_pad.SetTopMargin(0.12)
        main_pad.SetRightMargin(right_margin)
        canv.cd()
        main_pad.Draw()
        subplot_pad = ROOT.TPad("subplot_pad", "", 0, 0, 1, subplot_pad_height-subplot_pad_fudge)
        ROOT.SetOwnership(subplot_pad, False)
        subplot_pad.SetTicks(1, 1)
        subplot_pad.SetFillColor(0)
        subplot_pad.SetFillStyle(0)
        subplot_pad.SetTopMargin(5*subplot_pad_fudge)
        subplot_pad.SetRightMargin(right_margin)
        subplot_pad.SetBottomMargin(0.35)
        canv.cd()
        subplot_pad.Draw()
    else:
        main_pad = ROOT.TPad("main_pad", "", 0, 0, 1, 1)
        ROOT.SetOwnership(main_pad, False)
        main_pad.SetRightMargin(right_margin)
        main_pad.SetTicks(1, 1)
        main_pad.Draw()
        # main_pad.SetBottomMargin(2*subplot_pad_fudge)
        main_pad.SetTopMargin(0.12 * (1-subplot_pad_height))
        main_pad.SetRightMargin(right_margin)

    main_pad.cd()
    if logx:
        main_pad.SetLogx()
    if logy:
        main_pad.SetLogy()

    mg.Draw("AP")

    if do_fit_graph_ratio or do_ratio_plots:
        rescale_plot_labels(mg, 1-subplot_pad_height)
        # global FONT_SIZE
        # FONT_SIZE = FONT_SIZE / 1-subplot_pad_height

    # Little extra breathing room
    mg.GetHistogram().SetMaximum(mg.GetYaxis().GetXmax() * 1.1)
    mg.GetHistogram().SetMinimum(mg.GetYaxis().GetXmin() / 1.05)

    # Protection in case y limits are dominated by large stat error
    if y_limit_protection and len(y_limit_protection) == 2:
        y_min, y_max = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
        y_lim_lower, y_lim_upper = y_limit_protection
        if y_max > y_lim_upper:
            mg.GetHistogram().SetMaximum(y_lim_upper)
        if y_min < y_lim_lower:
            mg.GetHistogram().SetMinimum(y_lim_lower)

    xlimits = list(xlimits)
    if xlimits and len(xlimits) == 2:
        if xlimits[0] is None:
            xlimits[0] = mg.GetXaxis().GetXmin() / 1.05
        if xlimits[1] is None:
            xlimits[1] = mg.GetXaxis().GetXmax() * 1.2
        mg.GetXaxis().SetLimits(*xlimits)

    ylimits = list(ylimits)
    if ylimits and len(ylimits) == 2:
        if ylimits[0] is None:
            ylimits[0] = mg.GetYaxis().GetXmin() / 1.05
        if ylimits[1] is None:
            ylimits[1] = mg.GetYaxis().GetXmax() * 1.2

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

    canv.cd()
    cms_text = ROOT.TPaveText(0.17, 0.86, 0.2, 0.88, "NDC")
    cms_text.AddText("CMS Simulation Preliminary")
    cms_text.SetTextFont(62)
    cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_text.SetTextSize(FONT_SIZE)
    cms_text.SetBorderSize(0)
    cms_text.SetFillStyle(0)
    cms_text.Draw()

    if title:
        bin_text = ROOT.TPaveText(0.17, 0.82, 0.2, 0.83, "NDC")
        bin_text.AddText(title)
        bin_text.SetTextFont(42)
        bin_text.SetTextSize(FONT_SIZE)
        bin_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        bin_text.SetBorderSize(0)
        bin_text.SetFillStyle(0)
        bin_text.Draw()

    sample_text = ROOT.TPaveText(0.78, 0.93, 0.8, 0.94, "NDC")
    sample_text.AddText("QCD 13 TeV")
    sample_text.SetTextFont(42)
    sample_text.SetTextSize(FONT_SIZE)
    sample_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    sample_text.SetBorderSize(0)
    sample_text.SetFillStyle(0)
    sample_text.Draw()

    canv.Update()

    main_pad.cd()
    if vertical_line:
        y_min, y_max = ROOT.gPad.GetUymin(), ROOT.gPad.GetUymax()
        vert_line = ROOT.TLine(vertical_line, y_min, vertical_line, y_max)
        vert_line.SetLineStyle(2)
        vert_line.SetLineWidth(2)
        vert_line.SetLineColor(ROOT.kGray+2)
        vert_line.Draw()

    canv.cd()
    if other_elements:
        for ele in other_elements:
            ele.Draw()

    if do_fit_graph_ratio or do_ratio_plots:
        # Remove x axis labels & title
        mg.GetHistogram().GetXaxis().SetLabelSize(0)
        mg.GetXaxis().SetLabelSize(0)
        mg.GetHistogram().GetXaxis().SetTitleSize(0)
        mg.GetXaxis().SetTitleSize(0)

        subplot_pad.cd()
        if logx:
            subplot_pad.SetLogx()

        mg_sub = ROOT.TMultiGraph()
        if do_fit_graph_ratio:
            mg_sub.SetTitle(";".join(["", xtitle, "Fit / graph"]))
        elif do_ratio_plots:
            mg_sub.SetTitle(";".join(["", xtitle, ratio_title]))
        ROOT.SetOwnership(mg_sub, False)
        for entry in subplot_entries:
            mg_sub.Add(entry['graph'])
        if do_fit_graph_ratio:
            mg_sub.Draw("ALP")
        else:
            mg_sub.Draw("AP")

        # Set x axis range
        mg_sub.GetXaxis().SetLimits(mg.GetXaxis().GetXmin(), mg.GetXaxis().GetXmax())

        # Set y axis range
        if do_fit_graph_ratio:
            mg_sub.GetHistogram().SetMaximum(1.02)
            mg_sub.GetHistogram().SetMinimum(0.98)
        else:
            mg_sub.GetHistogram().SetMaximum(1.04)
            mg_sub.GetHistogram().SetMinimum(0.96)

        # Make it look sensible
        mg_sub.GetXaxis().SetTitleOffset(mg_sub.GetXaxis().GetTitleOffset()*2.8)
        mg_sub.GetYaxis().SetNdivisions(505)
        # mg_sub.GetYaxis().CenterTitle()

        # Draw a line at 1
        y_min, y_max = mg_sub.GetYaxis().GetXmin(), mg_sub.GetYaxis().GetXmax()
        if y_min < 1 and y_max > 1:
            x_min, x_max = mg_sub.GetXaxis().GetXmin(), mg_sub.GetXaxis().GetXmax()
            line = ROOT.TLine(x_min, 1, x_max, 1)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray+2)
            line.Draw()

        if vertical_line:
            y_min, y_max = mg_sub.GetHistogram().GetMinimum(), mg_sub.GetHistogram().GetMaximum()
            sub_vert_line = ROOT.TLine(vertical_line, y_min, vertical_line, y_max)
            sub_vert_line.SetLineStyle(2)
            sub_vert_line.SetLineWidth(2)
            sub_vert_line.SetLineColor(ROOT.kGray+2)
            sub_vert_line.Draw()

        rescale_plot_labels(mg_sub, subplot_pad_height)

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
    parser.add_argument("--slides", help="Dump JSON snippet to file for use with beamer-plot-slides", default=None)
    parser.add_argument("--vertLine", help="Add a vertical dashed line at this x value", type=float)
    args = parser.parse_args(in_args)
    print(args)

    cu.check_dir_exists_create(args.outputDir)

    if len(args.input) != len(args.label):
        raise RuntimeError("Need a --label for each --input")

    do_slides = args.slides not in (None, "")

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
        slides_contents = []
        for mydir in dirs:
            jec_text = ROOT.TPaveText(0.15, 0.93, 0.2, 0.94, "NDC")
            # jec_label = "Without JEC"
            # jec_label = "With JEC"
            # jec_label = "Summer16_23Sep2016V4"
            # jec_label = "Summer16_03Feb2017_V8"
            jec_text.AddText(args.title)
            jec_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
            jec_text.SetTextFont(42)
            jec_text.SetTextSize(FONT_SIZE)
            jec_text.SetBorderSize(0)
            jec_text.SetFillStyle(0)

            dir_text = ROOT.TPaveText(0.17, 0.78, 0.2, 0.79, "NDC")
            dir_label = mydir.upper().replace("PFCHS", " PF CHS").replace("PUPPI", " PUPPI").replace("L1L2L3", " + L1L2L3").replace("L1", " + L1")
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

            X_MIN, X_MAX = 4, None
            # For limit protection:
            Y_MIN, Y_MAX = 0.8, 1.6

            # Do all flavs corr vs pt for given eta bin
            common_eta_bins = cu.sort_human(cu.get_common_eta_bins(obj_list))
            for eta_bin in common_eta_bins:
                title = eta_bin.replace("to", " < #eta < ").replace("JetEta", "")
                # Do a per-flavour comparison plot
                outputs = []
                for fdict in entry_dicts:
                    entries = []
                    chi2entries = []
                    for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
                        entry = deepcopy(fdict)
                        entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                        entry['label'] += " [%s]" % label
                        entry["line_color"] = fdict['colour']+ind
                        entry["marker_color"] = fdict['colour']+ind
                        entry["fill_color"] = fdict['colour']+ind
                        entry["fill_style"] = 1001
                        entry["fill_alpha"] = 0.7
                        entry["ratio"] = None
                        if ind == 1:
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                            entry["line_style"] += 1
                        if ind == 2:
                            entry["line_style"] += 1
                        if ind == 3:
                            entry["line_style"] += 1
                            entry["marker_style"] = get_open_marker(entry['marker_style'])
                        entries.append(entry)
                        if ind != 0:
                            entries[-1]['ratio'] = entries[0]['graph']


                        if args.chi2:
                            chi2entry = deepcopy(entry)
                            chi2entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("corrVs", "Chi2NDoFVs"), eta_bin))
                            chi2entries.append(chi2entry)

                    output_filename = os.path.abspath(os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_funcGraphRatio.pdf" % (eta_bin, fdict['label'])))
                    outputs.append(output_filename)
                    do_comparison_graph(entries, title=title,
                                        xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                        xlimits=(X_MIN, X_MAX),
                                        # xlimits=None,
                                        y_limit_protection=(Y_MIN, Y_MAX),
                                        ylimits=ylimits,
                                        other_elements=other_elements,
                                        output_filename=output_filename,
                                        vertical_line=args.vertLine,
                                        do_fit_graph_ratio=True)

                    output_filename = os.path.abspath(os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_pyHerwigRatio.pdf" % (eta_bin, fdict['label'])))
                    do_comparison_graph(entries, title=title,
                                        xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                        xlimits=(X_MIN, X_MAX),
                                        # xlimits=None,
                                        y_limit_protection=(Y_MIN, Y_MAX),
                                        ylimits=ylimits,
                                        other_elements=other_elements,
                                        output_filename=output_filename,
                                        vertical_line=args.vertLine,
                                        do_ratio_plots=True, ratio_title="H++ / Py8")
                    if args.chi2:
                        do_comparison_graph(chi2entries, title=title,
                                            xtitle="p_{T}^{Reco} [GeV]", ytitle="#chi^{2}/N_{DoF}", logx=True,
                                            xlimits=(X_MIN, X_MAX),
                                            other_elements=other_elements,
                                            output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_chi2.pdf" % (eta_bin, fdict['label'])))

                # make 1 slide with all flavours for this eta bin
                if do_slides:
                    plots = ",\n".join(['            ["{fname}", "{title}"]'.format(fname=fname, title="") for fname in outputs])
                    text = """    {{
        "title": "${thistitle}$",
        "plots": [
{plots}
        ]
    }}""".format(thistitle=title, plots=plots)
                    slides_contents.append(text)

                # Do a plot with all flavours
                entries = []
                ud_entries, g_entries = [], []
                for fdict in entry_dicts:
                    for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
                        entry = deepcopy(fdict)
                        entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                        entry['label'] += " [%s]" % label
                        entry["line_color"] = fdict['colour']+ind
                        entry["marker_color"] = fdict['colour']+ind
                        entry["fill_color"] = fdict['colour']+ind
                        entry["fill_style"] = 1001
                        entry["fill_alpha"] = 0.8
                        entry['ratio'] = None
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
                        if ind > 0:
                            entries[-1]['ratio'] = entries[-2]['graph']

                title = eta_bin.replace("to", " < #eta < ").replace("JetEta", "")
                do_comparison_graph(entries, title=title,
                                    xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                    xlimits=(X_MIN, X_MAX),
                                    y_limit_protection=(Y_MIN, Y_MAX), ylimits=ylimits,
                                    other_elements=other_elements,
                                    output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_allFlavs_funcGraphRatio.pdf" % (eta_bin)),
                                    vertical_line=args.vertLine,
                                    do_fit_graph_ratio=True)

                do_comparison_graph(entries, title=title,
                                    xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction", logx=True,
                                    xlimits=(X_MIN, X_MAX),
                                    y_limit_protection=(Y_MIN, Y_MAX), ylimits=ylimits,
                                    other_elements=other_elements,
                                    output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_allFlavs_pyHerwigRatio.pdf" % (eta_bin)),
                                    vertical_line=args.vertLine,
                                    do_ratio_plots=True, ratio_title="H++ / Py8")

                # diff_entries = []
                # for ind, (ud, g, label) in enumerate(zip(ud_entries, g_entries, args.label)):
                #     diff_graph = construct_difference_graph(ud['graph'], g['graph'])
                #     diff = deepcopy(ud)
                #     diff['graph'] = diff_graph
                #     diff['label'] = label
                #     diff_entries.append(diff)
                # do_comparison_graph(diff_entries, title=title,
                #                     xtitle="p_{T}^{Reco} [GeV]", ytitle="Correction (ud) - Correction (g)", logx=True,
                #                     xlimits=(X_MIN, X_MAX), ylimits=(0, 0.08),
                #                     other_elements=other_elements,
                #                     output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_ud_g_diff.pdf" % (eta_bin)))

            # # Do all flavs corr vs eta for given pt bin
            # common_pt_bins = cu.sort_human(cu.get_common_pt_bins(obj_list))
            # for pt_bin in common_pt_bins:
            #     # Do a per-flavour comparison plot
            #     for fdict in entry_dicts:
            #         entries = []
            #         chi2entries = []
            #         for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
            #             entry = deepcopy(fdict)
            #             entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin))
            #             entry['label'] += " [%s]" % label
            #             entry["line_color"] = fdict['colour']
            #             entry["marker_color"] = fdict['colour']
            #             if ind == 1:
            #                 entry["marker_style"] = get_open_marker(entry['marker_style'])
            #             if ind == 2:
            #                 entry["line_style"] += 1
            #             if ind == 3:
            #                 entry["line_style"] += 1
            #                 entry["marker_style"] = get_open_marker(entry['marker_style'])
            #             entries.append(entry)
            #             if args.chi2:
            #                 chi2entry = deepcopy(entry)
            #                 chi2entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta").replace("corrVs", "Chi2NDoFVs"), pt_bin))
            #                 chi2entries.append(chi2entry)

            #         title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
            #         do_comparison_graph(entries, title=title + " GeV",
            #                             xtitle="|#eta|", ytitle="Correction",
            #                             xlimits=(0, 5.2), y_limit_protection=(Y_MIN, Y_MAX),
            #                             other_elements=other_elements, ylimits=ylimits,
            #                             output_filename=os.path.join(plot_dir, "compare_corr_vs_eta_%s_%s.pdf" % (pt_bin, fdict['label'])))
            #         if args.chi2:
            #             do_comparison_graph(entries, title=title,
            #                                 xtitle="|#eta|", ytitle="#chi^{2}/N_{DoF}",
            #                                 xlimits=(0, 5.2),
            #                                 other_elements=other_elements,
            #                                 output_filename=os.path.join(plot_dir, "compare_corr_vs_eta_%s_%s_chi2.pdf" % (pt_bin, fdict['label'])))

            #     # Do a plot with all flavours
            #     entries = []
            #     for fdict in entry_dicts:
            #         for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
            #             entry = deepcopy(fdict)
            #             entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin))
            #             entry['label'] += " [%s]" % label
            #             entry["line_color"] = fdict['colour']
            #             entry["marker_color"] = fdict['colour']
            #             if ind == 1:
            #                 entry["marker_style"] = get_open_marker(entry['marker_style'])
            #             if ind == 2:
            #                 entry["line_style"] += 1
            #             if ind == 3:
            #                 entry["line_style"] += 1
            #                 entry["marker_style"] = get_open_marker(entry['marker_style'])
            #             entries.append(entry)
            #     title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
            #     do_comparison_graph(entries, title=title + " GeV",
            #                         xtitle="|#eta|", ytitle="Correction",
            #                         xlimits=(0, 5.2), y_limit_protection=(Y_MIN, Y_MAX),
            #                         other_elements=other_elements, ylimits=ylimits,
            #                         output_filename=os.path.join(plot_dir, "compare_corr_vs_eta_%s_allFlavs.pdf" % (pt_bin)))

    if do_slides:
        print("Writing JSON to", args.slides)
        with open(args.slides, "w") as fout:
            fout.write(",\n".join(slides_contents))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
