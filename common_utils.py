"""Set of common functions that are used in loads of scripts."""


from __future__ import print_function
import ROOT
import os
from subprocess import call
from sys import platform as _platform
import numpy as np
import argparse
from array import array
import re
from math import fabs, log10, frexp, hypot
from copy import deepcopy


from MyStyle import My_Style
My_Style.cd()


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)

FONT_SIZE = 0.032


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Make argparse respect space formatting (newlines, etc) AND show defaults"""
    pass


def open_pdf(pdf_filename):
    """Open a PDF file using system's default PDF viewer."""
    if _platform.startswith("linux"):
        # linux
        call(["xdg-open", pdf_filename])
    elif _platform == "darwin":
        # OS X
        call(["open", pdf_filename])
    elif _platform == "win32":
        # Windows
        call(["start", pdf_filename])


#
# Filepath/directory fns
#
def cleanup_filepath(filepath):
    """Resolve any env vars, ~, etc, and return absolute path."""
    return os.path.abspath(os.path.expandvars(os.path.expanduser(filepath)))


def get_full_path(filepath):
    """Return absolute directory of filepath.
    Resolve any environment vars, ~, sym links(?)"""
    return os.path.dirname(cleanup_filepath(filepath))


def check_file_exists(filepath):
    """Check if file exists. Can do absolute or relative file paths."""
    return os.path.isfile(cleanup_filepath(filepath))


def check_dir_exists(filepath):
    """Check if directory exists."""
    return os.path.isdir(cleanup_filepath(filepath))


def check_dir_exists_create(filepath):
    """Check if directory exists. If not, create it."""
    if not check_dir_exists(filepath):
        os.makedirs(cleanup_filepath(filepath))


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
        print("Getting %s" % obj_name)
    obj = tfile.Get(obj_name) 
    if obj == None:
        raise IOError("No object named %s in %s" % (obj_name, tfile.GetName()))
    else:
        return obj


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


def check_exp(n):
    """
    Checks if number has stupidly larger exponent

    Can occur is using buffers - it just fills unused bins with crap
    """

    m, e = frexp(n)
    return fabs(log10(pow(2, e))) < 10


def get_xy(graph):
    """
    Return lists of x, y points from a graph, because it's such a PITA
    """
    # xarr = list(np.ndarray(graph.GetN(), 'd', graph.GetX()))
    # yarr = list(np.ndarray(graph.GetN(), 'd', graph.GetY()))
    xarr = array('d', graph.GetX())
    yarr = array('d', graph.GetY())
    return xarr, yarr


def get_exey(graph):
    """
    Return lists of errors on x, y points from a graph, because it's such a PITA
    """
    # xarr = list(np.ndarray(graph.GetN(), 'd', graph.GetEX()))
    # yarr = list(np.ndarray(graph.GetN(), 'd', graph.GetEY()))
    xarr = array('d', graph.GetEX())
    yarr = array('d', graph.GetEY())
    return xarr, yarr


def th2_to_arr(h):
    """Convert TH2 to 2D numpy array"""
    arr = np.zeros((h.GetNbinsX(), h.GetNbinsY()))
    for x_ind in xrange(1, h.GetNbinsX() + 1):
        for y_ind in xrange(1, h.GetNbinsY() + 1):
            arr[x_ind-1][y_ind-1] = h.GetBinContent(x_ind, y_ind)
    return arr


def make_normalised_TH2(hist, norm_axis, recolour=True):
    norm_axis = norm_axis.upper()
    possible_opts = ['X', 'Y']
    if norm_axis not in possible_opts:
        raise RuntimeError("norm_axis must be one of %s" % possible_opts)
    norm_axis = norm_axis.upper()

    # easiest way to cope with x or y is to just get a 2D matrix of values,
    # can then do transpose if necessary
    arr = th2_to_arr(hist)

    if norm_axis == 'Y':
        arr = arr.T
    if recolour:
        # can set so the maximum in each bin is the same,
        # scale other bins accordingly
        # this retain the colour scheme for each set of bins
        for ind, xbin in enumerate(arr):
            if xbin.max() > 0:
                arr[ind] = xbin / xbin.max()
    else:
        # alternatively, can rescale so sum over bins = 1
        for ind, xbin in enumerate(arr):
            if xbin.sum() != 0:
                arr[ind] = xbin / xbin.sum()

    if norm_axis == 'Y':
        arr = arr.T

    # Create new hist object - MUST do it this way to get Z range correct
    new_histname = hist.GetName() + "_norm" + norm_axis
    # hnew = ROOT.TH2F(hist)  # TODO: determine if TH2F, TH2D...
    if type(hist) == ROOT.TH2F:
        hnew = ROOT.TH2F(hist)  # TODO: determine if TH2F, TH2D...
    elif type(hist) == ROOT.TH2D:
        hnew = ROOT.TH2D(hist)  # TODO: determine if TH2F, TH2D...
    else:
        raise RuntimeError("Unknown 2D hist type")
    hnew.SetName(new_histname)
    for x_ind, x_arr in enumerate(arr, 1):
        for y_ind, val in enumerate(x_arr, 1):
            hnew.SetBinContent(x_ind, y_ind, val)
#     hnew.SetAxisRange(0.5, 1., 'Z')
    return hnew


def th1_to_arr(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)])


def tryfloat(s):
    try:
        return float(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryfloat(c) for c in re.split('([0-9.-]+)', s) ]


def sort_human(l):
    """ Sort the given list in the way that humans expect.
    """
    return sorted(l, key=alphanum_key)


def get_common_eta_bins(obj_list):
    """Get list of common eta bins from list of object names"""
    eta_bins = []
    for x in obj_list:
        m = re.search(r'JetEta[0-9.-]+to[0-9.-]+', x)
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


def get_alternate_colour(colour, modifier):
    """Get alternate colour, but keeps white and black invariant.

    Parameters
    ----------
    colour : int
        Original colour
    modifier : int
        Int to modify colour by

    Returns
    -------
    int
        New colour
    """
    if colour < 2:
        return colour
    else:
        return colour + int(modifier)


def construct_difference_graph(graph, other_graph):
    x, y = get_xy(graph)
    x_other, y_other = get_xy(other_graph)
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
    x, y = get_xy(graph)
    n = len(x)
    diff = array('d', [func.Eval(this_x) - this_y for this_x, this_y in zip(x, y)])
    ex = array('d', [0] * n)
    ey = array('d', [0] * n)
    gr = ROOT.TGraphErrors(n, x, diff, ex, ey)
    return gr


def construct_graph_func_ratio_graph(graph, func):
    """Construct a graph of function / graph for each (x, y) point in graph"""
    x, y = get_xy(graph)
    ex, ey = get_exey(graph)
    n = len(x)
    ratio = array('d', [func.Eval(this_x) / this_y for this_x, this_y in zip(x, y)])
    if len(ratio) != n:
        raise IndexError("Length of ratio not correct")
    # ey = array('d', [0] * n)
    new_ey = array('d', [f*e for f, e in zip(ratio, ey)])
    gr = ROOT.TGraphErrors(n, x, ratio, ex, new_ey)
    return gr


def construct_graph_ratio(graph, ref_graph):
    """Contruct a TGraph that is the ratio of the graph : ref_graph

    Interpolation will take place since x values may not align.
    The x values from ref_graph will be used as the reference values.
    """
    x, y = get_xy(graph)
    ex, ey = get_exey(graph)
    x_ref, y_ref = get_xy(ref_graph)
    ex_ref, ey_ref = get_exey(ref_graph)
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
    ratio_y = []
    for y1, y2 in zip(y, y_ref):
        if y2 == 0:
            ratio_y.append(0)
        else:
            ratio_y.append(y1/y2)
    # ex = [0]*n
    # ey = [0]*n
    ratio_ey = []
    for iy, iey, iy_ref, iey_ref in zip(y, ey, y_ref, ey_ref):
        if iy_ref == 0 or iy == 0:
            err = 0
        else:
            err = (iy/iy_ref)*hypot(iey/iy, iey_ref/iy_ref)
        ratio_ey.append(err)
    gr = ROOT.TGraphErrors(n, array('d', x_ref), array('d', ratio_y), array('d', ex), array('d', ratio_ey))
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


def do_comparison_graph(entries,
                        output_filename,
                        title="",
                        sample_title="",
                        xtitle="",
                        ytitle="",
                        other_elements=None,
                        logx=False,
                        logy=False,
                        xlimits=None,
                        ylimits=None,
                        y_limit_protection=None,
                        draw_fits=True,
                        do_horizontal_line=False,
                        vertical_line=None,
                        draw_opt="AP",
                        do_fit_graph_ratio=False,
                        do_ratio_plots=False,
                        ratio_limits=None,
                        ratio_title="",
                        ratio_draw_opt=None):
    """Draw several graphs on one canvas and save to file
    
    Parameters
    ----------
    entries : [dict]
        List of entries to plot. Each is represented by a dict, with the graph,
        label, and various other style options to be applied.
    output_filename : str
        Name of output plot file
    title : str
        Overall plot title, goes inside main axes in top left
    sample_title : str, optional
        Text to put with "13 TeV", e.g. QCD MC, goes top right
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
    xlimits : (min, max), optional
        Set hard x limits. Specifying either upper or lower limit as None lets
        the code make a sensible guess
    ylimits : (min, max), optional
        Set hard y limits. Specifying either upper or lower limit as None lets
        the code make a sensible guess
    do_horizontal_line : bool, optional
        Do horizontal line at 1
    y_limit_protection : (y_min, y_max), optional
        Set minimum and maximum y values in the event of a huge stat error or weird point
    draw_fits : bool, optional
        Draw fitted functions or not
    vertical_line : int, optional
        Draw a vertical dashed line at this x value
    draw_opt : str, optional
        Option to give to TMultiGraph.Draw()
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
    xdelta = 0.2
    xmiddle = 0.76
    ydelta = 0.09
    ymiddle = 0.77
    # if len(entries) > 4:
    #     xdelta *= 1.1
    #     # ymiddle -= 0.02
    #     ydelta *= 1.3
    leg = ROOT.TLegend(xmiddle-xdelta, ymiddle-ydelta, xmiddle+xdelta, ymiddle+ydelta)
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
        graph.SetMarkerSize(entry.get('marker_size', 1))

        graph.SetFillColorAlpha(entry.get('fill_color', default_colour), entry.get('fill_alpha', 1))
        graph.SetFillStyle(entry.get('fill_style', 1001))

        if graph.GetListOfFunctions().GetSize() > 0:
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
        
        if do_ratio_plots:
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
                    ratio_func.SetLineWidth(entry.get('line_width', 1))
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

    mg.Draw(draw_opt)

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

    # Handle x & y axis limits
    # None for either edge means let it decide "intelligently"
    if xlimits:
        xlimits = list(xlimits)
        if len(xlimits) == 2:
            if xlimits[0] is None:
                xlimits[0] = mg.GetXaxis().GetXmin() / 1.05
            if xlimits[1] is None:
                xlimits[1] = mg.GetXaxis().GetXmax() * 1.2
            mg.GetXaxis().SetLimits(*xlimits)

    if ylimits:
        ylimits = list(ylimits)
        if len(ylimits) == 2:
            if ylimits[0] is None:
                ylimits[0] = mg.GetYaxis().GetXmin() / 1.05
            if ylimits[1] is None:
                ylimits[1] = mg.GetYaxis().GetXmax() * 1.2

            mg.GetHistogram().SetMaximum(ylimits[1])
            mg.GetHistogram().SetMinimum(ylimits[0])

    leg.Draw()

    if do_horizontal_line:
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
        n = title.count("\n")
        # HERE BE MAGIC. Seriously though, it appears this works, so dont touch it
        gap = 0.055 / (1-subplot_pad_height)
        y1 = 0.79 - (n*gap)
        y2 = y1 + ((n+1)*gap)
        bin_text = ROOT.TPaveText(0.17, y1, 0.2, y2, "NDC")
        for i, substr in enumerate(title.split("\n")):
            bin_text.AddText(substr)
        bin_text.SetTextFont(42)
        bin_text.SetTextSize(FONT_SIZE)
        bin_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        bin_text.SetBorderSize(0)
        bin_text.SetFillStyle(0)
        bin_text.Draw()

    sample_text = ROOT.TPaveText(0.95, 0.93, 0.96, 0.94, "NDC")
    sample_text.AddText(sample_title+" 13 TeV")
    sample_text.SetTextFont(42)
    sample_text.SetTextSize(FONT_SIZE)
    sample_text.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
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

        if ratio_draw_opt:
            mg_sub.Draw(ratio_draw_opt)
        else:
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

        if ratio_limits and len(ratio_limits) == 2:
            if ratio_limits[1] is not None:
                mg_sub.GetHistogram().SetMaximum(ratio_limits[1])
            if ratio_limits[0] is not None:
                mg_sub.GetHistogram().SetMinimum(ratio_limits[0])

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
