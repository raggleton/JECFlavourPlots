#!/usr/bin/env python

"""Plot response & resolution graphs using ROOT file from output of jet_response_and_resolution_x"""


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

FONT_SIZE = 0.03


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
                        y_limit_protection=None, draw_fits=True, 
                        do_ratio=True, ratio_limits=None):
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
    do_ratio : bool, optional
        Add ratio subplot
    
    """

    plot = Plot(entries, what='graph', 
                title=None, xtitle=xtitle, ytitle=ytitle,
                xlim=xlimits, ylim=ylimits,
                legend=True,
                subplot=entries[0] if do_ratio else None, subplot_type="ratio",
                subplot_title="Ratio to All")

    # replace legend with our own for now
    delta = 0.12
    middle = 0.77
    plot.legend = ROOT.TLegend(middle-delta, 0.75, middle+delta, 0.88)
    plot.legend.SetBorderSize(0)
    plot.legend.SetFillStyle(0)
    plot.legend.SetNColumns(2)
    plot.legend.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignCenter)

    plot.plot()
    if logx:
        plot.set_logx()
    if logy:
        plot.set_logy()

    plot.main_pad.cd()
    if not ylimits:
        plot.container.GetHistogram().SetMaximum(plot.container.GetYaxis().GetXmax() * 1.05)

    # Protection in case y limits are dominated by large stat error
    if y_limit_protection and len(y_limit_protection) == 2:
        y_min, y_max = plot.container.GetYaxis().GetXmin(), plot.container.GetYaxis().GetXmax()
        y_lim_lower, y_lim_upper = y_limit_protection
        if y_max > y_lim_upper:
            plot.container.GetHistogram().SetMaximum(y_lim_upper)
        if y_min < y_lim_lower:
            plot.container.GetHistogram().SetMinimum(y_lim_lower)

    # add protection for subplot
    if not ratio_limits:
        low_lim = 0.8
        upper_lim = 1.2
        y_min, y_max = plot.subplot_container.GetYaxis().GetXmin(), plot.subplot_container.GetYaxis().GetXmax()
        plot.subplot_container.GetYaxis().SetRangeUser(max(low_lim, y_min), min(y_max, upper_lim))  # set limits doesn't work for y axis
    elif len(ratio_limits) == 2:
        plot.subplot_container.GetYaxis().SetRangeUser(*ratio_limits)  # set limits doesn't work for y axis

    plot.subplot_pad.Update()
    plot.canvas.Update()

    # if do_line:
    #     y_min, y_max = plot.container.GetYaxis().GetXmin(), plot.container.GetYaxis().GetXmax()
    #     if y_min < 1 and y_max > 1:
    #         x_min, x_max = plot.container.GetXaxis().GetXmin(), plot.container.GetXaxis().GetXmax()
    #         line = ROOT.TLine(x_min, 1, x_max, 1)
    #         line.SetLineStyle(2)
    #         line.SetLineColor(ROOT.kGray+2)
    #         line.Draw()

    plot.canvas.cd()

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

    plot.save(output_filename)


def main(in_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input ROOT file with response & resolution graphs (from jet_response_and_resolution_x)", required=True)
    parser.add_argument("--outputDir", help="Output directory for plots", default=os.getcwd())
    parser.add_argument("--title", help="Title string for plots", default="")
    parser.add_argument("--sampleName", help="Sample name string for plots", default="")
    args = parser.parse_args(in_args)

    cu.check_dir_exists_create(args.outputDir)

    lw = 2
    entry_dicts = [
        {"flav": "RelRspVsRefPt", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 2, "line_width": lw, "marker_size": 1.2},
        {"flav": "ud_RspVsRefPt_RelRsp", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 1.2},
        {"flav": "s_RRspVsRefPt_RelRsp", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "c_RRspVsRefPt_RelRsp", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "b_RRspVsRefPt_RelRsp", "label": "b", "colour": ROOT.kOrange-3, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.6},
        {"flav": "g_RRspVsRefPt_RelRsp", "label": "g", "colour": ROOT.kAzure+1, "marker_style": 29, "line_style": 1, "line_width": lw, "marker_size": 1.8},
    ]

    # no JEC
    eta_bin_limits = {
        'JetEta0to0.783': (0.9, 1.15),
        'JetEta0.783to1.305': None,
        'JetEta1.305to1.93': (0.85, 1.05),
        'JetEta1.93to2.5': None,
        'JetEta2.5to2.964': None,
        'JetEta2.964to3.139': None,
        'JetEta3.139to5.191': (0.98, 1.8),
    }
    eta_bin_ratio_limits = {
        'JetEta0to0.783': (0.96, 1.07),
        'JetEta0.783to1.305': None,
        'JetEta1.305to1.93': (0.95, 1.07),
        'JetEta1.93to2.5': None,
        'JetEta2.5to2.964': None,
        'JetEta2.964to3.139': None,
        'JetEta3.139to5.191': (0.93, 1.1),
    }
    
    # with L1 JEC
    # eta_bin_limits = {
    #     'JetEta0to0.783': (0.76, 1.05),
    #     'JetEta0.783to1.305': (0.77, 1),
    #     'JetEta1.305to1.93': (0.7, 1.),
    #     'JetEta1.93to2.5': (0.68, 1.1),
    #     'JetEta2.5to2.964': (0.55, 1.1),
    #     'JetEta2.964to3.139': (0.6, 1.1),
    #     'JetEta3.139to5.191': (0.75, 1.1),
    # }
    # eta_bin_ratio_limits = {
    #     'JetEta0to0.783': (0.96, 1.08),
    #     'JetEta0.783to1.305': (0.96, 1.08),
    #     'JetEta1.305to1.93': (0.95, 1.09),
    #     'JetEta1.93to2.5': (0.95, 1.13),
    #     'JetEta2.5to2.964': (0.93, 1.15),
    #     'JetEta2.964to3.139': (0.9, 1.12),
    #     'JetEta3.139to5.191': (0.93, 1.1),
    # }

    # default
    # eta_bin_limits = {
    #     'JetEta0to0.783': None,
    #     'JetEta0.783to1.305': None,
    #     'JetEta1.305to1.93': None,
    #     'JetEta1.93to2.5': None,
    #     'JetEta2.5to2.964': None,
    #     'JetEta2.964to3.139': None,
    #     'JetEta3.139to5.191': None,
    # }
    # eta_bin_ratio_limits = {
    #     'JetEta0to0.783': None,
    #     'JetEta0.783to1.305': None,
    #     'JetEta1.305to1.93': None,
    #     'JetEta1.93to2.5': None,
    #     'JetEta2.5to2.964': None,
    #     'JetEta2.964to3.139': None,
    #     'JetEta3.139to5.191': None,
    # }


    # Loop through all different ak4pfchs, etc
    dirs = cu.get_list_of_element_names(cu.open_root_file(args.input))
    for mydir in dirs[:]:
        jec_text = ROOT.TPaveText(0.14, 0.93, 0.2, 0.96, "NDC")
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

        sample_text = ROOT.TPaveText(0.94, 0.93, 0.95, 0.96, "NDC")
        # sample_text.AddText("Flat QCD 13 TeV")
        sample_text.AddText(args.sampleName) # + " 13 TeV")
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
        common_eta_bins = cu.sort_human(cu.get_common_eta_bins(obj_list))
        for eta_bin in common_eta_bins[:]:
            entries = []
            contributions = []
            for fdict in entry_dicts:
                entry = deepcopy(fdict)
                # entry["graph"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                entry["obj"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                entry["line_color"] = fdict['colour']
                entry["marker_color"] = fdict['colour']
                # entries.append(entry)

                del entry['flav']
                del entry['colour']
                this_contrib = Contribution(**entry)
                contributions.append(this_contrib)

            bin_title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
            do_comparison_graph(contributions, bin_title=bin_title,
                                xtitle="p_{T}^{Gen} [GeV]", ytitle="Response", logx=True,
                                xlimits=(10, 5000), 
                                ylimits=eta_bin_limits.get(eta_bin, None),
                                y_limit_protection=(0.5, 1.7),
                                other_elements=other_elements,
                                output_filename=os.path.join(plot_dir, "rsp_vs_pt_%s.pdf" % (eta_bin)),
                                ratio_limits=eta_bin_ratio_limits.get(eta_bin, None)
                                )

            # Inverse response ie ~ correction
            # entries = []
            # for fdict in entry_dicts:
            #     entry = deepcopy(fdict)
            #     entry["graph"] = construct_inverse_graph(cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin)))
            #     entry["line_color"] = fdict['colour']
            #     entry["marker_color"] = fdict['colour']
            #     entries.append(entry)
            # bin_title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
            # do_comparison_graph(entries, bin_title=bin_title,
            #                     xtitle="p_{T}^{Gen} [GeV]", ytitle="1/Response", logx=True,
            #                     y_limit_protection=(0.8, 1.2),
            #                     other_elements=other_elements,
            #                     output_filename=os.path.join(plot_dir, "inv_rsp_vs_pt_%s.pdf" % (eta_bin)))

        # Do all flavs rsp vs eta for given pt bin
        # common_pt_bins = cu.get_common_pt_bins(obj_list)
        # for pt_bin in common_pt_bins:
        #     entries = []
        #     contributions = []
        #     for fdict in entry_dicts:
        #         entry = deepcopy(fdict)
        #         # entry["graph"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin))
        #         entry["obj"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin))
        #         entry["line_color"] = fdict['colour']
        #         entry["marker_color"] = fdict['colour']
        #         # entries.append(entry)

        #         del entry['flav']
        #         del entry['colour']
        #         this_contrib = Contribution(**entry)
        #         contributions.append(this_contrib)
            
        #     bin_title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
        #     do_comparison_graph(contributions, bin_title=bin_title + " GeV",
        #                         xtitle="|#eta|", ytitle="Response",
        #                         xlimits=(0, 5.2), y_limit_protection=(0.5, 1.5),
        #                         other_elements=other_elements,
        #                         output_filename=os.path.join(plot_dir, "rsp_vs_eta_%s.pdf" % (pt_bin)))

            # Inverse response ie ~ correction
            # entries = []
            # for fdict in entry_dicts:
            #     entry = deepcopy(fdict)
            #     entry["graph"] = construct_inverse_graph(cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta"), pt_bin)))
            #     entry["line_color"] = fdict['colour']
            #     entry["marker_color"] = fdict['colour']
            #     entries.append(entry)
            # bin_title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
            # do_comparison_graph(entries, bin_title=bin_title + " GeV",
            #                     xtitle="|#eta|", ytitle="1/Response",
            #                     y_limit_protection=(0.8, 1.6),
            #                     other_elements=other_elements,
            #                     output_filename=os.path.join(plot_dir, "inv_rsp_vs_eta_%s.pdf" % (pt_bin)))


        # Do all flavs resolution vs pt for given eta bin
        # for eta_bin in common_eta_bins:
        #     entries = []
        #     for fdict in entry_dicts:
        #         entry = deepcopy(fdict)
        #         entry["graph"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'].replace("Rsp", "Res", 1), eta_bin))
        #         entry["line_color"] = fdict['colour']
        #         entry["marker_color"] = fdict['colour']
        #         entries.append(entry)
        #     bin_title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
        #     do_comparison_graph(entries, bin_title=bin_title,
        #                         xtitle="p_{T}^{Gen} [GeV]", ytitle="Relative resolution", logx=True,
        #                         y_limit_protection=(0, 0.5), draw_fits=False,
        #                         ylimits=(0, 0.5),
        #                         xlimits=(10, 5000), other_elements=other_elements,
        #                         output_filename=os.path.join(plot_dir, "res_vs_pt_%s.pdf" % (eta_bin)))

        # Do all flavs resolution plots vs eta for given pt bin
        # for pt_bin in common_pt_bins:
        #     entries = []
        #     for fdict in entry_dicts:
        #         entry = deepcopy(fdict)
        #         entry["graph"] = cu.grab_obj_from_file(args.input, "%s/%s_%s" % (mydir, fdict['flav'].replace("Rsp", "Res", 1).replace("RefPt", "JetEta"), pt_bin))
        #         entry["line_color"] = fdict['colour']
        #         entry["marker_color"] = fdict['colour']
        #         entries.append(entry)
        #     bin_title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
        #     do_comparison_graph(entries, bin_title=bin_title + " GeV",
        #                         xtitle="|#eta|", ytitle="Relative resolution", other_elements=other_elements,
        #                         y_limit_protection=(0, 0.5), draw_fits=True, xlimits=(0, 5.2),
        #                         output_filename=os.path.join(plot_dir, "res_vs_eta_%s.pdf" % (pt_bin)))

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
