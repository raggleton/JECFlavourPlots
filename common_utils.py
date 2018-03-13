"""Set of common functions that are used in loads of scripts."""


from __future__ import print_function
import ROOT
import os
from subprocess import call
from sys import platform as _platform
import numpy as np
import math
import argparse
from array import array
import re


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.TH1.SetDefaultSumw2(True)


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


def check_exp(n):
    """
    Checks if number has stupidly larger exponent

    Can occur is using buffers - it just fills unused bins with crap
    """

    from math import fabs, log10, frexp
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


def tryint(s):
    try:
        return int(s)
    except:
        return s
     

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]


def sort_human(l):
    """ Sort the given list in the way that humans expect.
    """
    return sorted(l, key=alphanum_key)
