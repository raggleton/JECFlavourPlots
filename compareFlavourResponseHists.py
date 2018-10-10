#!/usr/bin/env python

"""Plot response hists using ROOT file from output of jet_response_analyzer_x,
also plot flavour fraction graphs"""


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
# Move exponent label as it overlaps title box
ROOT.TGaxis.SetExponentOffset(-0.06, 0, "y")

FONT_SIZE = 0.03


def do_comparison_hist(entries, output_filename, bin_title="", xtitle="", ytitle="",
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
    num_entries = 0
    for entry in entries:
        default_colour = ROOT.kBlack
        hist_name = entry['hist'].GetName()
        hist = entry['hist'].Clone(hist_name+"Clone")
        hist.SetLineColor(entry.get('line_color', default_colour))
        hist.SetLineStyle(entry.get('line_style', 1))
        hist.SetLineWidth(entry.get('line_width', 1))

        hist.SetMarkerColor(entry.get('marker_color', default_colour))
        hist.SetMarkerStyle(entry.get('marker_style', 1))
        hist.SetMarkerSize(entry.get('marker_size', 1))

        # hist.Rebin(2)

        if hist.Integral() == 0:
            continue

        scale_factor = 1./hist.Integral()
        if normalise:
            hist.Scale(scale_factor)

        if hist.GetListOfFunctions().GetSize() > 0:
            func = hist.GetListOfFunctions().Last()
            if draw_fits:
                func.SetLineColor(entry.get('line_color', default_colour))
                # func.SetLineStyle(3)
                func.SetLineWidth(1)
                if normalise:
                    func.SetParameter(0, func.GetParameter(0)*scale_factor)
                funcs.append(func)
            else:
                # Broken for now
                hist.GetListOfFunctions().Remove(func)

        hst.Add(hist)
        num_entries += 1
        this_label = entry.get('label', hist_name)
        for i, substr in enumerate(this_label.split("\n")):
            if i == 0:
                leg.AddEntry(hist, substr, "LP")
            else:
                leg.AddEntry(0, substr, "")

    # if hst.GetHists().GetSize() == 0:  # segfaults, completely unobviously
    if num_entries == 0:
        return

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetTicks(1, 1)
    canv.SetRightMargin(0.03)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()
    hst.Draw("NOSTACK HISTE")

    # Little extra breathing room
    hst.SetMaximum(hst.GetMaximum("NOSTACK") * 1.4)

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
    cms_text.AddText("CMS Simulation Preliminary")
    cms_text.SetTextFont(62)
    cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_text.SetTextSize(FONT_SIZE)
    cms_text.SetBorderSize(0)
    cms_text.SetFillStyle(0)
    cms_text.Draw()

    n = bin_title.count("\n")
    # HERE BE MAGIC. Seriously though, it appears this works, so dont touch it
    y1 = 0.79 - (n*0.04)
    y2 = y1 + ((n+1)*0.04)
    bin_text = ROOT.TPaveText(0.17, y1, 0.2, y2, "NDC")
    for i, substr in enumerate(bin_title.split("\n")):
        bin_text.AddText(substr)
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


def construct_flavour_fractions(entries, bin_names, add_unknown=True):
    """Construct dict of flavour fractions from histograms.
    
    Parameters
    ----------
    entries : list[list[dict]]
        List of list of entry dicts which must have a 'hist' key
        Assumes that each top-level entry corresponds to a given pt bin
        Then each entry in that top-level entry corresponds to a different flavour
        God this is overly complicated
    bin_names : list[(float, float)]
        List of bin edges, should correspond to order of entries
    add_unknown : bool, optional
        Include an unknown fraction
    
    Returns
    -------
    dict{dict}
        Dict of { bin center : {flav : fraction } }

    """
    # use a tuple so can be used as dict keys
    bin_edges = [tuple([int(x) for x in re.search(r'([0-9.]+)to([0-9.]+)', bn).groups()]) 
                 for bn in bin_names]
    x = [(i[0]+i[1])*0.5 for i in sorted(bin_edges)]
    # entries may come in unsorted. this will help us sort them based on x value
    matched_entries = {(i[0]+i[1])*0.5: ent for i, ent in zip(bin_edges, entries)}
    # for each x bin, store dict of {flav: fraction}
    # key is midpoint of pt bin to get ordering
    # so overall looks like:
    # { pt : { flav : fraction } }
    all_entries = {}
    for x_value in sorted(matched_entries.keys()):
        values = {}
        total = 0
        for entry in matched_entries[x_value]:
            # flav_name = entry['label'].split()[0].strip()
            flav_name = entry['label']
            if "All" in flav_name:
                total = entry['hist'].Integral()
            else:
                values[flav_name] = entry['hist'].Integral()
        if total == 0:
            break
        else:
            for k in values:
                values[k] = values[k]/total

            if add_unknown:
                values['Unknown'] = 1 - sum(values.values())

        all_entries[x_value] = values

    return all_entries


def do_flavour_fraction_graph(entries, bin_names, output_filename, add_unknown=True,
                             title="", xtitle="", ytitle="", other_elements=None,
                             logx=False, logy=False,):
    """Do plot of flavour fraction vs something
    
    Parameters
    ----------
    entries : list[list[dict]]
        List of list of entry dicts which must have a 'hist' key
        Assumes that each top-level entry corresponds to a given pt bin
        Then each entry in that top-level entry corresponds to a different flavour
        Order is thus:
    bin_names : TYPE
        Description
    output_filename : TYPE
        Description
    add_unknown : bool, optional
        Description
    title : str, optional
        Description
    xtitle : str, optional
        Description
    ytitle : str, optional
        Description
    other_elements : None, optional
        Description
    logx : bool, optional
        Description
    logy : bool, optional
        Description
    """
    all_entries = []

    for ind, (label, these_entries) in enumerate(entries.items()):
        these_flav_fractions = construct_flavour_fractions(these_entries, bin_names, add_unknown)
        flavs = [x for x in next(iter(these_flav_fractions.values())).keys() if "All" not in x]

        for find, flav in enumerate(flavs):
            # Turn flavour fraction values into graphs
            pt_vals = sorted(these_flav_fractions.keys())
            y = [these_flav_fractions[mid_val].get(flav, 0) for mid_val in pt_vals]
            gr = ROOT.TGraph(len(y), array('d', pt_vals), array('d', y))

            default_colour = 13
            if "Unknown" in flav:
                new_entry = {
                    'graph': gr,
                    'line_color': cu.get_alternate_colour(default_colour, ind),
                    'line_style': 1,
                    'line_width': 2,
                    
                    'marker_color': cu.get_alternate_colour(default_colour, ind),
                    'marker_style': 34,
                    'marker_size': 1,
                    'label': "Unknown [%s]" % label,
                }
                if ind == 1:
                    new_entry["marker_style"] = cu.get_open_marker(new_entry['marker_style'])
                    new_entry["line_style"] += 1
                if ind > 0:
                    new_entry['ratio'] = all_entries[len(flavs)*(ind-1) + find]['graph']
                all_entries.append(new_entry)
            else:
                this_entry = [e for e in these_entries[0] if e['label'] == flav][0]
                new_entry = deepcopy(this_entry)
                new_entry['graph'] = gr
                if ind > 0:
                    new_entry['ratio'] = all_entries[len(flavs)*(ind-1) + find]['graph']
                all_entries.append(new_entry)

    ylim = (0.001, 10) if logy else (0, 1.2)
    cu.do_comparison_graph(all_entries, 
                           title=title,
                           sample_title="QCD MC",
                           xtitle=xtitle,
                           ytitle=ytitle,
                           logx=logx,
                           logy=logy,
                           xlimits=(10, None),
                           ylimits=ylim,
                           draw_opt="ALP",
                           other_elements=other_elements,
                           output_filename=output_filename,
                           do_ratio_plots=True,
                           ratio_title="H++ / Py8",
                           ratio_limits=(0.5, 1.5),
                           ratio_draw_opt="ALP"
                           )


def do_plot_N_graph(entries, bin_names, output_filename, which,
                    title="", xtitle="", other_elements=None,
                    logx=False, logy=False,):
    # dict of methodname : y-axis label
    opts = {'GetN': "N", "GetEffectiveEntries": "N_{effective}", "Integral": "Integral"}
    if which not in opts:
        raise RuntimeError("`which` arg of `do_plot_N_graph` must be one of %s" % list(opts.keys()))

    # use a tuple so can be used as dict keys
    bin_edges = [tuple([int(x) for x in re.search(r'([0-9.]+)to([0-9.]+)', bn).groups()]) for bn in bin_names]
    x = [(i[0]+i[1])*0.5 for i in sorted(bin_edges)]

    if len(entries) != len(bin_names):
        raise RuntimeError ("Collections need to be same length)")

    # entries come in unsorted. this will help us sort them based on x value
    matched_entries = {(i[0]+i[1])*0.5: ent for i, ent in zip(bin_edges, entries)}

    # for each x bin, store dict of {flav: N}
    # key is midpoint of pt bin to get ordering
    all_entries = {}
    for x_value in sorted(matched_entries.keys()):
        values = {}
        for entry in matched_entries[x_value]:
            num = getattr(entry['hist'], which)()
            values[entry['label']] = num
        all_entries[x_value] = values

    delta = 0.12
    middle = 0.77
    leg = ROOT.TLegend(middle-delta, 0.75, middle+delta, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)
    leg.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignCenter)

    # flavs = all_entries.values()[0].keys()  # screws up ordering
    flavs = [e['label'] for e in entries[0] ]
    mg = ROOT.TMultiGraph()
    mg.SetTitle(";".join(["", xtitle, opts[which]]))
    graphs = []
    y_max, y_min = 0, 99999
    for flav in flavs:
        y = [all_entries[mid_val].get(flav, 0) for mid_val in sorted(all_entries.keys())]
        if max(y) > y_max:
            y_max = max(y)
        if min(y) < y_min:
            y_min = min(y)
        gr = ROOT.TGraph(len(all_entries), array('d', x), array('d', y))

        default_colour = ROOT.kBlack
        this_entry = [e for e in entries[0] if e['label'] == flav][0]
        gr.SetLineColor(this_entry.get('line_color', default_colour))
        gr.SetLineStyle(this_entry.get('line_style', 1))
        gr.SetLineWidth(this_entry.get('line_width', 1))

        gr.SetMarkerColor(this_entry.get('marker_color', default_colour))
        gr.SetMarkerStyle(this_entry.get('marker_style', 1))
        gr.SetMarkerSize(this_entry.get('marker_size', 1))

        graphs.append(gr)
        mg.Add(gr)
        leg.AddEntry(gr, flav, "LP")

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetRightMargin(0.03)
    canv.SetTicks(1, 1)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()

    mg.Draw("ALP")
    print(y_min, y_max)
    if logy:
        mg.GetHistogram().SetMaximum(y_max*30)
        # need little offset to avoid setting to 0 in log scale
        if y_min == 0:
            y_min = 0.1
        mg.GetHistogram().SetMinimum(y_min / 1.1)
    else:
        mg.GetHistogram().SetMaximum(y_max)
        mg.GetHistogram().SetMinimum(y_min)


    x_max = mg.GetXaxis().GetXmax()
    mg.GetXaxis().SetLimits(10, x_max)
    leg.Draw()

    cms_text = ROOT.TPaveText(0.17, 0.84, 0.2, 0.85, "NDC")
    cms_text.AddText("CMS Simulation")
    cms_text.SetTextFont(62)
    cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_text.SetTextSize(FONT_SIZE)
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
    bin_text.SetTextSize(FONT_SIZE)
    bin_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    bin_text.SetBorderSize(0)
    bin_text.SetFillStyle(0)
    bin_text.Draw()

    if other_elements:
        for ele in other_elements:
            ele.Draw()

    canv.SaveAs(output_filename)


def do_plot_RMS_graph(entries, bin_names, output_filename,
                      title="", xtitle="", other_elements=None,
                      logx=False, logy=False,):
    # use a tuple so can be used as dict keys
    bin_edges = [tuple([int(x) for x in re.search(r'([0-9.]+)to([0-9.]+)', bn).groups()]) for bn in bin_names]
    x = [(i[0]+i[1])*0.5 for i in sorted(bin_edges)]

    if len(entries) != len(bin_names):
        raise RuntimeError ("Collections need to be same length)")

    # entries come in unsorted. this will help us sort them based on x value
    matched_entries = {(i[0]+i[1])*0.5: ent for i, ent in zip(bin_edges, entries)}

    # for each x bin, store dict of {flav: N}
    # key is midpoint of pt bin to get ordering
    all_entries = {}
    for x_value in sorted(matched_entries.keys()):
        values = {}
        for entry in matched_entries[x_value]:
            values[entry['label']] = entry['hist'].GetRMS()
        all_entries[x_value] = values

    delta = 0.21
    middle = 0.77
    leg = ROOT.TLegend(middle-delta, 0.65, middle+delta, 0.86)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if len(entries) > 4:
        leg.SetNColumns(2)
    leg.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignCenter)

    # flavs = all_entries.values()[0].keys()  # screws up ordering
    flavs = [e['label'] for e in entries[0] ]
    mg = ROOT.TMultiGraph()
    mg.SetTitle(";".join(["", xtitle, "RMS"]))
    graphs = []
    for flav in flavs:
        y = [all_entries[mid_val].get(flav, 0) for mid_val in sorted(all_entries.keys())]
        gr = ROOT.TGraph(len(all_entries), array('d', x), array('d', y))

        default_colour = ROOT.kBlack
        this_entry = [e for e in entries[0] if e['label'] == flav][0]
        gr.SetLineColor(this_entry.get('line_color', default_colour))
        gr.SetLineStyle(this_entry.get('line_style', 1))
        gr.SetLineWidth(this_entry.get('line_width', 1))

        gr.SetMarkerColor(this_entry.get('marker_color', default_colour))
        gr.SetMarkerStyle(this_entry.get('marker_style', 1))
        gr.SetMarkerSize(this_entry.get('marker_size', 1))

        graphs.append(gr)
        mg.Add(gr)
        leg.AddEntry(gr, flav, "LP")

    canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canv.SetRightMargin(0.03)
    canv.SetTicks(1, 1)
    if logx:
        canv.SetLogx()
    if logy:
        canv.SetLogy()

    mg.Draw("ALP")
    if logy:
        mg.GetHistogram().SetMaximum(2)
        mg.GetHistogram().SetMinimum(0.01)
    # else:
    #     mg.GetHistogram().SetMaximum(1.2)
    #     mg.GetHistogram().SetMinimum(0)
    x_max = mg.GetXaxis().GetXmax()
    mg.GetXaxis().SetLimits(10, x_max)
    leg.Draw()

    cms_text = ROOT.TPaveText(0.17, 0.84, 0.2, 0.85, "NDC")
    cms_text.AddText("CMS Simulation")
    cms_text.SetTextFont(62)
    cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_text.SetTextSize(FONT_SIZE)
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
    parser.add_argument("--input", help="Input ROOT file with response hists", action='append')
    parser.add_argument("--label", help="Label for input file", action='append')
    parser.add_argument("--outputDir", help="Output directory for plots", default=os.getcwd())
    parser.add_argument("--title", help="Title string for plots")
    parser.add_argument("--sampleName", help="Sample name string for plots", default="")
    parser.add_argument("--slides", help="Dump JSON snippet to file for use with beamer-plot-slides", default=None)
    parser.add_argument("--plotResponseHists", help="Plot individual response hists", action='store_true')
    parser.add_argument("--plotN", help="Plot graph of N vs pT that goes into the median", action='store_true')
    parser.add_argument("--plotRMS", help="Plot graph of raw RMS vs pT that goes into memdian", action='store_true')
    args = parser.parse_args(in_args)

    if not args.input or len(args.input) == 0:
        return 1
    
    cu.check_dir_exists_create(args.outputDir)

    if len(args.input) != len(args.label):
        raise RuntimeError("Need a --label for each --input")

    lw = 2
    entry_dicts = [
        {"flav": "RelRsp", "label": "All", "colour": ROOT.kBlack, "marker_style": ROOT.kFullCircle, "line_style": 1, "line_width": lw, "marker_size": 1.2},
        {"flav": "ud_RelRsp", "label": "ud", "colour": ROOT.kRed, "marker_style": ROOT.kFullSquare, "line_style": 1, "line_width": lw, "marker_size": 1.2},
        {"flav": "g_RelRsp", "label": "g", "colour": ROOT.kAzure+1, "marker_style": 29, "line_style": 1, "line_width": lw, "marker_size": 1.8},
        {"flav": "s_RelRsp", "label": "s", "colour": ROOT.kBlue, "marker_style": ROOT.kFullTriangleUp, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "c_RelRsp", "label": "c", "colour": ROOT.kGreen+2, "marker_style": ROOT.kFullTriangleDown, "line_style": 1, "line_width": lw, "marker_size": 1.4},
        {"flav": "b_RelRsp", "label": "b", "colour": ROOT.kOrange-3, "marker_style": ROOT.kFullDiamond, "line_style": 1, "line_width": lw, "marker_size": 1.6},
    ]

    # Loop through all different ak4pfchs, etc
    all_dirs = [set(cu.get_list_of_element_names(cu.open_root_file(infile))) for infile in args.input]
    dirs = all_dirs[0]
    for d in all_dirs[1:]:
        dirs = dirs & d
    dirs = sorted(list(dirs))[:1]
    print("Doing: ", dirs)
    for mydir in dirs[:]:
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

        # Do separate dir for each eta bin, then separate plot for each pt bin
        common_eta_bins = cu.sort_human(cu.get_common_eta_bins(obj_list))
        common_pt_bins = cu.sort_human(cu.get_common_pt_bins(obj_list))
        print("Doing eta bins", common_eta_bins)
        print("Doing pt bins", common_pt_bins)

        open_tfiles = []
        for input_filename in args.input:
            open_tfiles.append(cu.open_root_file(input_filename))

        for eta_bin in common_eta_bins[:]:
            this_plot_dir = os.path.join(plot_dir, eta_bin)
            cu.check_dir_exists_create(this_plot_dir)

            # Dict of all entries
            # key is input label
            # value is list of list of dicts
            all_entries = {}
            for label in args.label:
                all_entries[label] = []

            for pt_bin in common_pt_bins[:]:
                entries = []
                for ind, label in enumerate(args.label):
                    this_input_entries = []
                    for find, fdict in enumerate(entry_dicts):
                        entry = deepcopy(fdict)
                        try:
                            entry["hist"] = cu.get_from_tfile(open_tfiles[ind], "%s/%s_%s_%s" % (mydir, fdict['flav'], eta_bin, pt_bin))
                            entry["hist"].Rebin(int(entry["hist"].GetNbinsX()/200))
                            entry['label'] += " [%s]" % label
                            new_colour = cu.get_alternate_colour(fdict['colour'], ind)
                            entry["line_color"] = new_colour
                            entry["marker_color"] = new_colour
                            entry["fill_color"] = new_colour
                            entry["fill_style"] = 1001
                            entry["fill_alpha"] = 0.7
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
                            this_input_entries.append(entry)
                        # if ind > 0:
                        #     entries[-1]['ratio'] = entries[-2]['hist']
                            # all_entries[label].append(deepcopy(entry))
                        except IOError:
                            pass
                    all_entries[label].append(this_input_entries)

                if len(entries) == 0:
                    continue

                bin_title = eta_bin.replace("to", " < #eta < ").replace("JetEta", "")
                bin_title += "\n"
                bin_title += pt_bin.replace("to", " < p^{Gen}_{T} < ").replace("RefPt", "")
                bin_title += " GeV"
                norm_entries = deepcopy(entries)
                rsp_str = "Response (p_{T}^{Reco} / p_{T}^{Gen})"
                if args.plotResponseHists:
                    do_comparison_hist(entries,
                                       bin_title=bin_title,
                                       xlimits=(0, 2),
                                       xtitle=rsp_str,
                                       ytitle="N",
                                       other_elements=other_elements,
                                       normalise=False,
                                       output_filename=os.path.join(this_plot_dir, "rsp_vs_pt_%s.pdf" % (pt_bin)))

                    do_comparison_hist(norm_entries,
                                       bin_title=bin_title,
                                       xlimits=(0, 2),
                                       xtitle=rsp_str,
                                       ytitle="p.d.f",
                                       other_elements=other_elements,
                                       normalise=True,
                                       draw_fits=False,
                                       output_filename=os.path.join(this_plot_dir, "rsp_vs_pt_%s_normed.pdf" % (pt_bin)))

            this_dir_text = dir_text.Clone()
            # this_dir_text.SetY1(0.76)
            # this_dir_text.SetY2(0.77)
            bin_title = eta_bin.replace("to", " < #eta < ").replace("JetEta", "")
            do_flavour_fraction_graph(all_entries, common_pt_bins, title=bin_title,
                                      xtitle="p^{Gen}_{T} [GeV]", ytitle="Flavour fraction",
                                      other_elements=other_elements,
                                      add_unknown=True, logx=True, logy=True,
                                      output_filename=os.path.join(plot_dir, "flav_frac_vs_pt_%s.pdf" % (eta_bin)))

            # if args.plotN:
            #     do_plot_N_graph(all_pt_entries, common_pt_bins, title=bin_title,
            #                     which="GetEffectiveEntries",
            #                     xtitle="p^{Gen}_{T} [GeV]",
            #                     other_elements=[jec_text, this_dir_text],
            #                     logx=True, logy=True,
            #                     output_filename=os.path.join(plot_dir, "Neff_vs_pt_%s.pdf" % (eta_bin)))

            # if args.plotRMS:
            #     do_plot_RMS_graph(all_pt_entries, common_pt_bins, title=bin_title,
            #                       xtitle="p^{Gen}_{T} [GeV]",
            #                       other_elements=[jec_text, this_dir_text],
            #                       logx=True, logy=True,
            #                       output_filename=os.path.join(plot_dir, "RMS_vs_pt_%s.pdf" % (eta_bin)))


    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
