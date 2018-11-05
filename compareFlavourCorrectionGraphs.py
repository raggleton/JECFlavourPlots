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


def get_eta_mid(s):
    # find numbers
    parts = cu.alphanum_key(s)
    nums = [i for i in parts if isinstance(i, float)]
    if len(nums) > 2:
        raise RuntimeError("no idea how to get midpoint of bin %s" % s)
    return 0.5*(nums[0] + nums[1])


def get_eta_half_width(s):
    # find half width between bin edges
    parts = cu.alphanum_key(s)
    nums = [i for i in parts if isinstance(i, float)]
    if len(nums) > 2:
        raise RuntimeError("no idea how to get width of bin %s" % s)
    return fabs(0.5*(nums[0] - nums[1]))


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

    if not args.input or len(args.input) == 0:
        return 1

    cu.check_dir_exists_create(args.outputDir)

    if len(args.input) != len(args.label):
        raise RuntimeError("Need a --label for each --input")

    do_slides = args.slides not in (None, "")

    lw = 1
    entry_dicts = [
        # {"flav": "AbsCorVsJetPt", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 2, "line_width": lw, "marker_size": 1.2},
        {"flav": "ud_AbsCorVsJetPt", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 0.8},
        {"flav": "g_AbsCorVsJetPt", "label": "g", "colour": ROOT.kAzure+1, "marker_style": 29, "line_style": 1, "line_width": lw, "marker_size": 1.},
        {"flav": "s_AbsCorVsJetPt", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.},
        {"flav": "c_AbsCorVsJetPt", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.},
        {"flav": "b_AbsCorVsJetPt", "label": "b", "colour": ROOT.kOrange-3, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.2},
    ][:]

    sample_str = "QCD MC"

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

        X_MIN, X_MAX = 7, None
        # For limit protection:
        Y_MIN, Y_MAX = 0.8, 1.6

        py_herwig_ratio_title = "Parton / Hadron"
        py_herwig_ratio_title = "H++ / Py8"

        # Store chi2 for fits over pT
        # Each entry is for a flavour,
        # and is a dict of file label : chi2/dof values
        fit_chi2entries = {}
        for fdict in entry_dicts:
            fit_chi2entries[fdict['label']] = {}
            for label in args.label:
                fit_chi2entries[fdict['label']][label] = []

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
                    new_colour = cu.get_alternate_colour(fdict['colour'], ind)
                    entry["line_color"] = new_colour
                    entry["marker_color"] = new_colour
                    entry["fill_color"] = new_colour
                    entry["fill_style"] = 1001
                    entry["fill_alpha"] = 0.7
                    entry["ratio"] = None
                    if ind == 1:
                        entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
                        entry["line_style"] += 1
                    if ind == 2:
                        entry["line_style"] += 1
                    if ind == 3:
                        entry["line_style"] += 1
                        entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
                    entries.append(entry)
                    if ind != 0:
                        entries[-1]['ratio'] = entries[0]['graph']


                    if args.chi2:
                        # chi2entry = deepcopy(entry)
                        # chi2entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("corrVs", "Chi2NDoFVs"), eta_bin))
                        # chi2entries.append(chi2entry)
                        if entry['graph'].GetListOfFunctions().GetSize() > 0:
                            fit = entry['graph'].GetListOfFunctions().Last()
                            # chi2 = entry['graph'].Chisquare(fit)  # can't do this as we don't want to use the plateau in chi2
                            chi2 = fit.GetChisquare()
                            ndof = fit.GetNDF()
                            fit_chi2entries[fdict['label']][label].append(chi2 / ndof)

                output_filename = os.path.abspath(os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_funcGraphRatio.pdf" % (eta_bin, fdict['label'])))
                cu.do_comparison_graph(entries, title=title, sample_title=sample_str,
                                       xtitle="p_{T}^{Reco} [GeV]",
                                       ytitle="Correction",
                                       logx=True,
                                       xlimits=(X_MIN, X_MAX),
                                       # xlimits=None,
                                       y_limit_protection=(Y_MIN, Y_MAX),
                                       ylimits=ylimits,
                                       other_elements=other_elements,
                                       output_filename=output_filename,
                                       vertical_line=args.vertLine,
                                       do_fit_graph_ratio=True)
                outputs.append(output_filename)

                output_filename = os.path.abspath(os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_pyHerwigRatio.pdf" % (eta_bin, fdict['label'])))
                cu.do_comparison_graph(entries, title=title, sample_title=sample_str,
                                       xtitle="p_{T}^{Reco} [GeV]",
                                       ytitle="Correction",
                                       logx=True,
                                       xlimits=(X_MIN, X_MAX),
                                       # xlimits=None,
                                       y_limit_protection=(Y_MIN, Y_MAX),
                                       ylimits=ylimits,
                                       other_elements=other_elements,
                                       output_filename=output_filename,
                                       vertical_line=args.vertLine,
                                       draw_fits=True,
                                       do_ratio_plots=True,
                                       ratio_title=py_herwig_ratio_title)
 
                output_filename = os.path.abspath(os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_pull.pdf" % (eta_bin, fdict['label'])))
                cu.do_comparison_graph(entries, title=title, sample_title=sample_str,
                                       xtitle="p_{T}^{Reco} [GeV]",
                                       ytitle="Correction",
                                       logx=True,
                                       xlimits=(X_MIN, X_MAX),
                                       # xlimits=None,
                                       y_limit_protection=(Y_MIN, Y_MAX),
                                       ylimits=ylimits,
                                       other_elements=other_elements,
                                       output_filename=output_filename,
                                       vertical_line=args.vertLine,
                                       do_fit_graph_pull=True)
                # if args.chi2:
                #     cu.do_comparison_graph(chi2entries, title=title, sample_title=sample_str,
                #                         xtitle="p_{T}^{Reco} [GeV]", ytitle="#chi^{2}/N_{DoF}", logx=True,
                #                         xlimits=(X_MIN, X_MAX),
                #                         other_elements=other_elements,
                #                         output_filename=os.path.join(plot_dir, "compare_corr_vs_pt_%s_%s_chi2.pdf" % (eta_bin, fdict['label'])))

            # Do a plot with all flavours
            entries = []
            ud_entries, g_entries = [], []
            for fdict in entry_dicts[:]:
                for ind, (input_filename, label) in enumerate(zip(args.input, args.label)):
                    entry = deepcopy(fdict)
                    entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'], eta_bin))
                    entry['label'] += " [%s]" % label
                    new_colour = cu.get_alternate_colour(fdict['colour'], ind)
                    entry["line_color"] = new_colour
                    entry["marker_color"] = new_colour
                    entry["fill_color"] = new_colour
                    entry["fill_style"] = 1001
                    entry["fill_alpha"] = 0.8
                    entry['ratio'] = None
                    if ind == 1:
                        entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
                        entry["line_style"] += 1
                    if ind == 2:
                        entry["line_style"] += 1
                    if ind == 3:
                        entry["line_style"] += 1
                        entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
                    entries.append(entry)
                    if fdict['label'] == "ud":
                        ud_entries.append(entry)
                    elif fdict['label'] == "g":
                        g_entries.append(entry)
                    if ind > 0:
                        entries[-1]['ratio'] = entries[-2]['graph']

            title = eta_bin.replace("to", " < #eta < ").replace("JetEta", "")
            output_filename = os.path.abspath(os.path.join(plot_dir, "compare_corr_vs_pt_%s_allFlavs_funcGraphRatio.pdf" % (eta_bin)))
            cu.do_comparison_graph(entries, title=title, sample_title=sample_str,
                                   xtitle="p_{T}^{Reco} [GeV]",
                                   ytitle="Correction",
                                   logx=True,
                                   xlimits=(X_MIN, X_MAX),
                                   y_limit_protection=(Y_MIN, Y_MAX), ylimits=ylimits,
                                   other_elements=other_elements,
                                   output_filename=output_filename,
                                   vertical_line=args.vertLine,
                                   do_fit_graph_ratio=True)
            outputs.insert(0, output_filename)

            output_filename = os.path.abspath(os.path.join(plot_dir, "compare_corr_vs_pt_%s_allFlavs_pyHerwigRatio.pdf" % (eta_bin)))
            cu.do_comparison_graph(entries, title=title, sample_title=sample_str,
                                   xtitle="p_{T}^{Reco} [GeV]",
                                   ytitle="Correction",
                                   logx=True,
                                   xlimits=(X_MIN, X_MAX),
                                   y_limit_protection=(Y_MIN, Y_MAX), ylimits=ylimits,
                                   other_elements=other_elements,
                                   output_filename=output_filename,
                                   vertical_line=args.vertLine,
                                   draw_fits=True,
                                   do_ratio_plots=True, 
                                   ratio_title=py_herwig_ratio_title)
            
            # make 1 slide with all flavours for this eta bin
            if do_slides:
                plots = ",\n".join(['            ["{{{fname}}}{ext}", "{title}"]'.format(
                                        fname=os.path.splitext(fname)[0],
                                        ext=os.path.splitext(fname)[1],
                                        title="")
                                    for fname in outputs])
                text = """    {{
    "title": "${thistitle}$",
    "plots": [
{plots}
        ]
    }}""".format(thistitle=title.replace("#", r"\\"), plots=plots)
                
                slides_contents.append(text)


        if args.chi2:
            # for each flav, do a plot of chi2 of total fit vs eta for all files
            all_entries = []
            all_ymax = 0
            for fdict in entry_dicts:
                entries = []
                file_labels = list(fit_chi2entries[fdict['label']].keys())
                ymax = 0
                for ind, flabel in enumerate(file_labels):
                    eta_bins = [get_eta_mid(x) for x in common_eta_bins]
                    y = fit_chi2entries[fdict['label']][flabel]
                    ymax = max(ymax, max(y))
                    ey = [0] * len(y)
                    ex = [get_eta_half_width(x) for x in common_eta_bins]
                    gr = ROOT.TGraphErrors(len(y), array('d', eta_bins), array('d', y), array('d', ex), array('d', ey))
                    entry = deepcopy(fdict)
                    entry["graph"] = gr
                    entry['label'] = "%s [%s]" % (fdict['label'], flabel)
                    new_colour = cu.get_alternate_colour(fdict['colour'], ind)
                    entry["line_color"] = new_colour
                    entry["marker_color"] = new_colour
                    entry["fill_color"] = new_colour
                    entry["fill_style"] = 1001
                    entry["fill_alpha"] = 0.8
                    if ind == 1:
                        entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
                        entry["line_style"] += 1
                    if ind == 2:
                        entry["line_style"] += 1
                    if ind == 3:
                        entry["line_style"] += 1
                        entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
                    entries.append(entry)
                all_ymax = max(all_ymax, ymax)
                ylimits = [0, ymax*1.3]
                # if ylimits[1] > 20:
                #     ylimits[1] = 20
                output_filename = os.path.abspath(os.path.join(plot_dir, "compare_chi2_vs_eta_%s.pdf" % (fdict['label'])))
                cu.do_comparison_graph(entries, title="", sample_title=sample_str,
                                       xtitle="#eta^{Reco}", ytitle="#chi^{2} / N_{dof}",
                                       draw_opt="ALP",
                                       ylimits=ylimits,
                                       other_elements=other_elements,
                                       output_filename=output_filename)
                all_entries.extend(entries)

            # do a single plot with all flavs
            ylimits = [0, all_ymax*1.3]
            if ylimits[1] > 20:
                ylimits[1] = 20
            output_filename = os.path.abspath(os.path.join(plot_dir, "compare_chi2_vs_eta_allFlavs.pdf"))
            cu.do_comparison_graph(all_entries, title="", sample_title=sample_str,
                                   xtitle="#eta^{Reco}", ytitle="#chi^{2} / N_{dof}",
                                   draw_opt="ALP",
                                   ylimits=ylimits,
                                   other_elements=other_elements,
                                   output_filename=output_filename)


            # diff_entries = []
            # for ind, (ud, g, label) in enumerate(zip(ud_entries, g_entries, args.label)):
            #     diff_graph = construct_difference_graph(ud['graph'], g['graph'])
            #     diff = deepcopy(ud)
            #     diff['graph'] = diff_graph
            #     diff['label'] = label
            #     diff_entries.append(diff)
            # cu.do_comparison_graph(diff_entries, title=title, sample_title=sample_str,
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
        #                 entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
        #             if ind == 2:
        #                 entry["line_style"] += 1
        #             if ind == 3:
        #                 entry["line_style"] += 1
        #                 entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
        #             entries.append(entry)
        #             if args.chi2:
        #                 chi2entry = deepcopy(entry)
        #                 chi2entry["graph"] = cu.grab_obj_from_file(input_filename, "%s/%s_%s" % (mydir, fdict['flav'].replace("RefPt", "JetEta").replace("corrVs", "Chi2NDoFVs"), pt_bin))
        #                 chi2entries.append(chi2entry)

        #         title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
        #         cu.do_comparison_graph(entries, title=title + " GeV", sample_title=sample_str,
        #                             xtitle="|#eta|", ytitle="Correction",
        #                             xlimits=(0, 5.2), y_limit_protection=(Y_MIN, Y_MAX),
        #                             other_elements=other_elements, ylimits=ylimits,
        #                             output_filename=os.path.join(plot_dir, "compare_corr_vs_eta_%s_%s.pdf" % (pt_bin, fdict['label'])))
        #         if args.chi2:
        #             cu.do_comparison_graph(entries, title=title, sample_title=sample_str,
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
        #                 entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
        #             if ind == 2:
        #                 entry["line_style"] += 1
        #             if ind == 3:
        #                 entry["line_style"] += 1
        #                 entry["marker_style"] = cu.get_open_marker(entry['marker_style'])
        #             entries.append(entry)
        #     title = pt_bin.replace("to", " < p_{T} < ").replace("RefPt", "")
        #     cu.do_comparison_graph(entries, title=title + " GeV", sample_title=sample_str,
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
