#!/usr/bin/env python

"""Plot response hists using ROOT file from output of jet_response_fitter"""


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


def do_comparison_hist(entries, output_filename, title="", xtitle="", ytitle="",
                       other_elements=None, logx=False, logy=False,
                       do_line=True, xlimits=None, ylimits=None,
                       y_limit_protection=None, draw_fits=True, normalise=True):
    """Draw several hists on one canvas and save to file
    
    Parameters
    ----------
    entries : [dict]
        List of entries to plot. Each is represented by a dict, with the hist,
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
    normalise : bool, optional
        Normalise each hist so that integral = 1
    
    """
    hst = ROOT.THStack(ROOT.TUUID().AsString(), ";".join(["", xtitle, ytitle]))
    delta = 0.12
    middle = 0.77
    leg = ROOT.TLegend(middle-delta, 0.75, middle+delta, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)
    leg.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignCenter)

    funcs = []  # have to save them here, THStack doesn't draw them by default for some stupid reason
    for entry in entries:
        default_colour = ROOT.kBlack
        hist = entry['hist']
        hist.SetLineColor(entry.get('line_color', default_colour))
        hist.SetLineStyle(entry.get('line_style', 1))
        hist.SetLineWidth(entry.get('line_width', 1))

        hist.SetMarkerColor(entry.get('marker_color', default_colour))
        hist.SetMarkerStyle(entry.get('marker_style', 1))
        hist.SetMarkerSize(entry.get('marker_size', 1))

        if hist.GetListOfFunctions().GetSize() > 0:
            func = hist.GetListOfFunctions().Last()
            if draw_fits:
                func.SetLineColor(entry.get('line_color', default_colour))
                # func.SetLineStyle(3)
                func.SetLineWidth(1)
                funcs.append(func)
            else:
                # Broken for now
                hist.GetListOfFunctions().Remove(func)

        hst.Add(hist)
        this_label = entry.get('label', hist.GetName())
        for i, substr in enumerate(this_label.split("\n")):
            if i == 0:
                leg.AddEntry(hist, substr, "LP")
            else:
                leg.AddEntry(0, substr, "")

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetTicks(1, 1)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()
    hst.Draw("NOSTACK HISTE")

    # Little extra breathing room
    hst.SetMaximum(hst.GetMaximum("NOSTACK") * 1.1)

    if draw_fits:
        for f in funcs:
            f.Draw("SAME")

    # Protection in case y limits are dominated by large stat error
    # if y_limit_protection and len(y_limit_protection) == 2:
    #     y_min, y_max = hst.GetYaxis().GetXmin(), hst.GetYaxis().GetXmax()
    #     y_lim_lower, y_lim_upper = y_limit_protection
    #     if y_max > y_lim_upper:
    #         hst.GetHistogram().SetMaximum(y_lim_upper)
    #     if y_min < y_lim_lower:
    #         hst.GetHistogram().SetMinimum(y_lim_lower)

    if xlimits and len(xlimits) == 2:
        hst.GetXaxis().SetLimits(*xlimits)

    if ylimits and len(ylimits) == 2:
        hst.SetMaximum(ylimits[1])
        hst.SetMinimum(ylimits[0])

    leg.Draw()

    if do_line:
        x_min, x_max = hst.GetXaxis().GetXmin(), hst.GetXaxis().GetXmax()
        if x_min < 1 and x_max > 1:
            y_min, y_max = hst.GetYaxis().GetXmin(), hst.GetHistogram().GetMaximum()
            line = ROOT.TLine(1, y_min, 1, y_max)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray+2)
            line.Draw()

    cms_text = ROOT.TPaveText(0.17, 0.84, 0.2, 0.85, "NDC")
    cms_text.AddText("CMS")
    cms_text.SetTextFont(62)
    cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_text.SetTextSize(0.035)
    cms_text.SetBorderSize(0)
    cms_text.SetFillStyle(0)
    cms_text.Draw()

    n = title.count("\n")
    # HERE BE MAGIC. Seriously though, it appears this works, so dont touch it
    y1 = 0.79 - (n*0.04)
    y2 = y1 + ((n+1)*0.04)
    bin_text = ROOT.TPaveText(0.17, y1, 0.2, y2, "NDC")
    for i, substr in enumerate(title.split("\n")):
        bin_text.AddText(substr)
    bin_text.SetTextFont(42)
    bin_text.SetTextSize(0.035)
    bin_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    bin_text.SetBorderSize(0)
    bin_text.SetFillStyle(0)
    bin_text.Draw()

    sample_text = ROOT.TPaveText(0.65, 0.91, 0.67, 0.92, "NDC")
    sample_text.AddText("Flat QCD 13 TeV")
    sample_text.SetTextFont(42)
    sample_text.SetTextSize(0.035)
    sample_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    sample_text.SetBorderSize(0)
    sample_text.SetFillStyle(0)
    sample_text.Draw()

    if other_elements:
        for ele in other_elements:
            ele.Draw()

    canv.SaveAs(output_filename)


def main(in_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputHists", help="Input ROOT file with response & resolution hists (from jet_response_fitter_x)")
    parser.add_argument("--outputDir", help="Output directory for plots", default=os.getcwd())
    parser.add_argument("--title", help="Title string for plots", default="")
    args = parser.parse_args(in_args)

    cu.check_dir_exists_create(args.outputDir)

    if args.inputHists:

        lw = 2
        entry_dicts = [
            {"flav": "RelRsp", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 2, "line_width": lw, "marker_size": 1.2},
            {"flav": "ud_RelRsp", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 1.2},
            {"flav": "s_RelRsp", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.4},
            {"flav": "c_RelRsp", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.4},
            {"flav": "b_RelRsp", "label": "b", "colour": ROOT.kOrange-3, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.6},
            {"flav": "g_RelRsp", "label": "g", "colour": ROOT.kAzure+1, "marker_style": 29, "line_style": 1, "line_width": lw, "marker_size": 1.8},
        ]


        # Loop through all different ak4pfchs, etc
        dirs = cu.get_list_of_element_names(cu.open_root_file(args.inputHists))
        for mydir in dirs[:]:
            jec_text = ROOT.TPaveText(0.17, 0.91, 0.2, 0.92, "NDC")
            jec_text.AddText(args.title)
            jec_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
            jec_text.SetTextFont(42)
            jec_text.SetTextSize(0.035)
            jec_text.SetBorderSize(0)
            jec_text.SetFillStyle(0)

            dir_text = ROOT.TPaveText(0.17, 0.72, 0.2, 0.73, "NDC")
            dir_label = mydir.upper().replace("PFCHS", " PF CHS").replace("PUPPI", " PUPPI").replace("L1L2L3", " + L1L2L3")
            dir_text.AddText(dir_label)
            dir_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
            dir_text.SetTextFont(42)
            dir_text.SetTextSize(0.035)
            dir_text.SetBorderSize(0)
            dir_text.SetFillStyle(0)

            other_elements = [jec_text, dir_text]

            plot_dir = os.path.join(args.outputDir, mydir)
            cu.check_dir_exists_create(plot_dir)

            obj_list = cu.get_list_of_objects_in_dir(args.inputHists, mydir)

            # Do separate dir for each eta bin, then separate plot for each pt bin
            common_eta_bins = get_common_eta_bins(obj_list)
            common_pt_bins = get_common_pt_bins(obj_list)
            for eta_bin in common_eta_bins:
                this_plot_dir = os.path.join(plot_dir, eta_bin)
                cu.check_dir_exists_create(this_plot_dir)

                for pt_bin in common_pt_bins:
                    entries = []
                    for fdict in entry_dicts:
                        entry = deepcopy(fdict)
                        entry["hist"] = cu.grab_obj_from_file(args.inputHists, "%s/%s_%s_%s" % (mydir, fdict['flav'], eta_bin, pt_bin))
                        entry["line_color"] = fdict['colour']
                        entry["marker_color"] = fdict['colour']
                        entries.append(entry)
                    title = eta_bin.replace("to", " < |#eta| < ").replace("JetEta", "")
                    title += "\n"
                    title += pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
                    title += " GeV"
                    do_comparison_hist(entries, title=title,
                                        xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", ytitle="N",
                                        other_elements=other_elements,
                                        output_filename=os.path.join(this_plot_dir, "rsp_vs_pt_%s.pdf" % (pt_bin)))


    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
