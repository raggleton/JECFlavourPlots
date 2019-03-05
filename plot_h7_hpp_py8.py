#!/usr/bin/env python

"""Do H7 vs H++ vs Py8 plots"""


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
sys.path.append("..")
from MyStyle import My_Style
import common_utils as cu
from comparator import Contribution, Plot


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
My_Style.cd()


FONT_SIZE = 0.032

lw = 2

def label(generator, flav, jec):
    generator_dict = {
    "h7": "Herwig7",
    "hpp": "Herwig++",
    "py16": "Pythia8 (CUETP8M1)",
    "py18": "Pythia8 (CP5)",
    } 
    jec_dict = {
        "": "No JEC/NoPU",
        "L1": "L1FastJet",
        "L1L2L3": "L1FastJetL2L3",
    }
    return ", ".join([generator_dict[generator], flav, jec_dict[jec]])


def do_1d_plot(entries, histname, output_filename, plot_kwargs, logx=False, logy=False):
    conts = []
    for ent in entries:
        hist = cu.get_from_tfile(ent['file'], "ak4pfchs/"+ent.get('histname', histname))
        if hist.GetEntries() == 0:
            ent['used'] = False
            continue
        contrib = Contribution(hist, label=ent['label'],
                               line_width=lw, line_color=ent['colour'], line_style=ent.get('line_style', 1),
                               marker_size=ent.get('marker_size', 0), marker_color=ent['colour'], marker_style=ent.get('marker_style', 1),
                               normalise_hist=True, rebin_hist=ent.get('rebin', None))
        conts.append(contrib)
        ent['used'] = True

    entries = [e for e in entries if e['used']]
    for ind, ent in enumerate(entries):
        if not ent['used']:
            continue
        if ent.get('subplot_ind', -1) >= 0:
            conts[ind].subplot = conts[ent['subplot_ind']].obj

    plot = Plot(conts, what="hist", has_data=False, **plot_kwargs)
    # if is_pdgid_plot and conts[0].obj.GetNbinsX() > 10:
    #     plot.default_canvas_size = (800, 800)
    # else:
    plot.default_canvas_size = (450, 600)
    plot.legend.SetX1(0.45)
    plot.legend.SetX2(0.97)
    if len(entries) > 4:
        plot.legend.SetY1(0.7)
        plot.legend.SetNColumns(2)
    else:
        plot.legend.SetY1(0.75)
    plot.legend.SetY2(0.87)
    plot.plot("HISTE NOSTACK")
    if logx:
        plot.set_logx()
    if logy:
        plot.set_logy()
    plot.save(output_filename)



def make_projection_plot(entries, output_filename, plot_kwargs, projection_axis='x', start_val=None, end_val=None, is_pdgid_plot=False, logy=False, histname=""):
    """Make a plot from entries. Each one is a projection of a 2D plot.
    
    Parameters
    ----------
    entries : list[dict]
        List of plot entries. Each is represented by a dict
    output_filename : str
        Filename for output plot
    plot_kwargs : dict
        kwargs for Plot constructor
    projection_axis : str, optional
        Axis to use for projection. If None, must be specified in entry.
    start_val : None, optional
        If set, use this to make 1D projection plot. Cuts on X axis. 
        Otherwise should be set per entry.
    end_val : None, optional
        If set, use this to make 1D projection plot. Cuts on X axis. 
        Otherwise should be set per entry.
    is_pdgid_plot : bool, optional
        True if for PDGIDs (sets special size & binning)
    logy : bool, optional
        Description
    
    Raises
    ------
    RuntimeError
        Description
    """
    conts = []
    if not is_pdgid_plot:
        for ent in entries:
            # h2d = cu.get_from_tfile(ent['file'], "ak4pfchs/"+ent.get('histname', histname)+"_JetEta0to0.783")
            h2d = cu.get_from_tfile(ent['file'], "ak4pfchs/"+ent.get('histname', histname)+"_JetEta2.5to2.964")
            start_val = ent.get('start_val', start_val)
            end_val = ent.get('end_val', end_val)
            if start_val is None or end_val is None:
                raise RuntimeError("Expected start_val and end_val")
            hist = cu.get_projection_plot(h2d, start_val, end_val, ent.get('projection_axis', projection_axis))
            if hist.GetEntries() == 0:
                ent['used'] = False
                continue
            contrib = Contribution(hist, label=ent['label'],
                                   line_width=lw, line_color=ent['colour'], line_style=ent.get('line_style', 1),
                                   marker_size=ent.get('marker_size', 0), marker_color=ent['colour'], marker_style=ent.get('marker_style', 1),
                                   normalise_hist=True, rebin_hist=ent.get('rebin', None))
            conts.append(contrib)
            ent['used'] = True

    else:
        # For PDGID plots, we want a different canvas aspect ratio, and only to include bins
        # with non-0 contents. We also relabel bins.
        custom_bins = []
        for ent in entries:
            h2d = cu.get_from_tfile(ent['file'], "ak4pfchs/"+ent['histname']+"_JetEta0to0.783")
            start_val = ent.get('start_val', start_val)
            end_val = ent.get('end_val', end_val)
            if not start_val or not end_val:
                raise RuntimeError("Expected start_val and end_val")
            hist = cu.get_projection_plot(h2d, start_val, end_val, ent.get('projection_axis', projection_axis))
            ax = hist.GetXaxis()
            bins = dict()
            for i in range(1, hist.GetNbinsX()+1):
                value = ax.GetBinLowEdge(i)
                cont = hist.GetBinContent(i)
                err = hist.GetBinError(i)
                if cont > 0:
                    bins[value] = [cont, err]
                    # print(custom_bins[-1])

            custom_bins.append(bins)

        nbins = max([len(d) for d in custom_bins])
        first_keys = set(custom_bins[0].keys())
        for cb in custom_bins[1:]:
            first_keys = first_keys.union(set(cb.keys()))
        all_keys = sorted(list(first_keys))
        print(all_keys)

        for cbin, ent in zip(custom_bins, entries):
            h = ROOT.TH1D(cu.get_unique_str(), ";PDGID;N", nbins, 0, nbins)
            for ind, k in enumerate(all_keys, 1):
                content, error = cbin.get(k, [0, 0])
                h.SetBinContent(ind, content)
                h.SetBinError(ind, error)
            ax = h.GetXaxis()
            for i in range(1, nbins+1):
                # ax.SetBinLabel(i, str(int(all_keys[i-1])))
                pdgid = int(all_keys[i-1])
                ax.SetBinLabel(i, PDGID_STR.get(pdgid, str(pdgid)))
            h.LabelsOption("v")
            contrib = Contribution(h, label=ent['label'],
                                   line_width=lw, line_color=ent['colour'], line_style=ent.get('line_style', 1),
                                   marker_size=ent.get('marker_size', 0), marker_color=ent['colour'], marker_style=ent.get('marker_style', 1),
                                   normalise_hist=True, rebin_hist=ent.get('rebin', None))
            conts.append(contrib)
            ent['used'] = True

    entries = [e for e in entries if e['used']]
    for ind, ent in enumerate(entries):
        if not ent['used']:
            continue
        if ent.get('subplot_ind', -1) >= 0:
            conts[ind].subplot = conts[ent['subplot_ind']].obj

    plot = Plot(conts, what="hist", ytitle="p.d.f.", has_data=False, **plot_kwargs)
    if is_pdgid_plot and conts[0].obj.GetNbinsX() > 10:
        plot.default_canvas_size = (800, 800)
    else:
        plot.default_canvas_size = (450, 600)
    plot.legend.SetX1(0.5)
    plot.legend.SetX2(0.97)
    if len(entries) > 4:
        plot.legend.SetX1(0.35)
        plot.legend.SetY1(0.75)
        plot.legend.SetNColumns(2)
    else:
        plot.legend.SetY1(0.75)
    plot.legend.SetY2(0.9)
    plot.plot("HISTE NOSTACK")
    if logy:
        plot.set_logy()
    plot.save(output_filename)


if __name__ == "__main__":

    # P8M1 2016 PYTHIA
    PYTHIA_16_L1_DIR = "QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJet"
    pythia_16_all_L1 = cu.open_root_file(os.path.join(PYTHIA_16_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    pythia_16_ud_L1 = cu.open_root_file(os.path.join(PYTHIA_16_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    pythia_16_g_L1 = cu.open_root_file(os.path.join(PYTHIA_16_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    pythia_16_s_L1 = cu.open_root_file(os.path.join(PYTHIA_16_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    pythia_16_c_L1 = cu.open_root_file(os.path.join(PYTHIA_16_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    pythia_16_b_L1 = cu.open_root_file(os.path.join(PYTHIA_16_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))
    
    PYTHIA_16_L1L2L3_DIR = "QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3"
    pythia_16_all_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_16_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    pythia_16_ud_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_16_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    pythia_16_g_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_16_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    pythia_16_s_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_16_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    pythia_16_c_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_16_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    pythia_16_b_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_16_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))

    # pythia_16_all_NoPU = cu.open_root_file("QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3/Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root")
    # pythia_16_ud_NoPU = cu.open_root_file("QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3/Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root")
    # pythia_16_g_NoPU = cu.open_root_file("QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3/Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root")
    # pythia_16_s_NoPU = cu.open_root_file("QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3/Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root")
    # pythia_16_c_NoPU = cu.open_root_file("QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3/Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root")
    # pythia_16_b_NoPU = cu.open_root_file("QCD_Pt_16_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3/Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root")

    # CP5 2018 PYTHIA8
    PYTHIA_18_L1_DIR = "QCD_Pt_18_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Autumn18_V1_MC_L1FastJet"
    pythia_18_all_L1 = cu.open_root_file(os.path.join(PYTHIA_18_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    pythia_18_ud_L1 = cu.open_root_file(os.path.join(PYTHIA_18_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    pythia_18_g_L1 = cu.open_root_file(os.path.join(PYTHIA_18_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    pythia_18_s_L1 = cu.open_root_file(os.path.join(PYTHIA_18_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    pythia_18_c_L1 = cu.open_root_file(os.path.join(PYTHIA_18_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    pythia_18_b_L1 = cu.open_root_file(os.path.join(PYTHIA_18_L1_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))
    
    PYTHIA_18_L1L2L3_DIR = "QCD_Pt_18_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Autumn18_V1_MC_L1FastJetL2L3"
    pythia_18_all_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_18_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    pythia_18_ud_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_18_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    pythia_18_g_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_18_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    pythia_18_s_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_18_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    pythia_18_c_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_18_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    pythia_18_b_L1L2L3 = cu.open_root_file(os.path.join(PYTHIA_18_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15toInf_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))

    # HERWIG++ 2016
    HPP_L1_DIR = "QCD_Pt_Herwigpp_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJet"
    herwigpp_all_L1 = cu.open_root_file(os.path.join(HPP_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    herwigpp_ud_L1 = cu.open_root_file(os.path.join(HPP_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    herwigpp_g_L1 = cu.open_root_file(os.path.join(HPP_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    herwigpp_s_L1 = cu.open_root_file(os.path.join(HPP_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    herwigpp_c_L1 = cu.open_root_file(os.path.join(HPP_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    herwigpp_b_L1 = cu.open_root_file(os.path.join(HPP_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))

    HPP_L1L2L3_DIR = "QCD_Pt_Herwigpp_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Summer16_07Aug2017_V20_MC_L1FastJetL2L3"
    herwigpp_all_L1L2L3 = cu.open_root_file(os.path.join(HPP_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    herwigpp_ud_L1L2L3 = cu.open_root_file(os.path.join(HPP_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    herwigpp_g_L1L2L3 = cu.open_root_file(os.path.join(HPP_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    herwigpp_s_L1L2L3 = cu.open_root_file(os.path.join(HPP_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    herwigpp_c_L1L2L3 = cu.open_root_file(os.path.join(HPP_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    herwigpp_b_L1L2L3 = cu.open_root_file(os.path.join(HPP_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_ext_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))

    HPP_NoPU_DIR = "QCD_Pt_Herwigpp_NoJEC_NoPU_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k"
    herwigpp_all_NoPU = cu.open_root_file(os.path.join(HPP_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    herwigpp_ud_NoPU = cu.open_root_file(os.path.join(HPP_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    herwigpp_g_NoPU = cu.open_root_file(os.path.join(HPP_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    herwigpp_s_NoPU = cu.open_root_file(os.path.join(HPP_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    herwigpp_c_NoPU = cu.open_root_file(os.path.join(HPP_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    herwigpp_b_NoPU = cu.open_root_file(os.path.join(HPP_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwigpp_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))


    # 2018 HERWIG7
    H7_L1_DIR = "QCD_Pt_Herwig7_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Autumn18_V1_MC_L1FastJet"
    herwig7_all_L1 = cu.open_root_file(os.path.join(H7_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    herwig7_ud_L1 = cu.open_root_file(os.path.join(H7_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    herwig7_g_L1 = cu.open_root_file(os.path.join(H7_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    herwig7_s_L1 = cu.open_root_file(os.path.join(H7_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    herwig7_c_L1 = cu.open_root_file(os.path.join(H7_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    herwig7_b_L1 = cu.open_root_file(os.path.join(H7_L1_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))

    H7_L1L2L3_DIR = "QCD_Pt_Herwig7_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_Autumn18_V1_MC_L1FastJetL2L3"
    herwig7_all_L1L2L3 = cu.open_root_file(os.path.join(H7_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    herwig7_ud_L1L2L3 = cu.open_root_file(os.path.join(H7_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    herwig7_g_L1L2L3 = cu.open_root_file(os.path.join(H7_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    herwig7_s_L1L2L3 = cu.open_root_file(os.path.join(H7_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    herwig7_c_L1L2L3 = cu.open_root_file(os.path.join(H7_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    herwig7_b_L1L2L3 = cu.open_root_file(os.path.join(H7_L1L2L3_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))

    H7_NoPU_DIR = "QCD_Pt_Herwig7_NoJEC_NoPU_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k"
    herwig7_all_NoPU = cu.open_root_file(os.path.join(H7_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_all.root"))
    herwig7_ud_NoPU = cu.open_root_file(os.path.join(H7_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_ud.root"))
    herwig7_g_NoPU = cu.open_root_file(os.path.join(H7_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_g.root"))
    herwig7_s_NoPU = cu.open_root_file(os.path.join(H7_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_s.root"))
    herwig7_c_NoPU = cu.open_root_file(os.path.join(H7_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_c.root"))
    herwig7_b_NoPU = cu.open_root_file(os.path.join(H7_NoPU_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_NoPU_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin10_nrefmax0_drmin0p8_ptGenOverlap10_b.root"))

    OUTPUT_DIR = "h7_hpp_py8_studies"
    cu.check_dir_exists_create(OUTPUT_DIR)

    HPP_COL = ROOT.kAzure+10
    H7_COL = ROOT.kBlue
    PY_16_COL = ROOT.kViolet-2
    PY_18_COL = ROOT.kMagenta+3

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, 
    #                  "label": label("hpp", "all", "L1L2L3"),
    #                  "colour": HPP_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_L1, 
    #                  "label": label("h7", "all", "L1L2L3"),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='nPV',
    #             output_filename=os.path.join(OUTPUT_DIR, "nPV_16vs18.pdf"),
    #             plot_kwargs=dict(xtitle="# PV", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, 
    #                  "label": label("hpp", "all", "L1L2L3"),
    #                  "colour": HPP_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_L1, 
    #                  "label": label("h7", "all", "L1L2L3"),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='nPU',
    #             output_filename=os.path.join(OUTPUT_DIR, "nPU_16vs18.pdf"),
    #             plot_kwargs=dict(xtitle="# Interactions", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, 
    #                  "label": label("hpp", "all", "L1L2L3"),
    #                  "colour": HPP_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_L1, 
    #                  "label": label("h7", "all", "L1L2L3"),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='nTrueInt',
    #             output_filename=os.path.join(OUTPUT_DIR, "nTrueInt_16vs18.pdf"),
    #             plot_kwargs=dict(xtitle="# True Interactions", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, 
    #                  "label": label("hpp", "all", "L1L2L3"),
    #                  "colour": HPP_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_L1, 
    #                  "label": label("h7", "all", "L1L2L3"),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='rho',
    #             output_filename=os.path.join(OUTPUT_DIR, "rho_herwig_16vs18.pdf"),
    #             plot_kwargs=dict(xtitle="rho", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": pythia_18_all_L1L2L3, 
    #                  "label": label("py18", "all", "L1L2L3"),
    #                  "colour": PY_18_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_L1, 
    #                  "label": label("h7", "all", "L1L2L3"),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='rho',
    #             output_filename=os.path.join(OUTPUT_DIR, "rho_pyvsherwig_18.pdf"),
    #             plot_kwargs=dict(xtitle="rho", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7/PY8", subplot_limits=None)
    # )

    pt_ylim = [1.E-13, 5.]

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, 
    #                  "label": label("hpp", "all", "L1L2L3"),
    #                  "colour": HPP_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_L1, 
    #                  "label": label("h7", "all", "L1L2L3"),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='pThat', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "pthat_herwig_16vs18_L1L2L3.pdf"),
    #             plot_kwargs=dict(xtitle="#hat{p}_{T} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1, 
    #                  "label": label("hpp", "all", "L1"),
    #                  "colour": HPP_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_L1, 
    #                  "label": label("h7", "all", "L1"),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='pThat', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "pthat_herwig_16vs18_L1.pdf"),
    #             plot_kwargs=dict(xtitle="#hat{p}_{T} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_NoPU, 
    #                  "label": label("hpp", "all", ""),
    #                  "colour": HPP_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": herwig7_all_NoPU, 
    #                  "label": label("h7", "all", ""),
    #                  "colour": H7_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='pThat', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "pthat_herwig_16vs18_NoPU.pdf"),
    #             plot_kwargs=dict(xtitle="#hat{p}_{T} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": pythia_16_all_L1L2L3, 
    #                  "label": label("py16", "all", "L1L2L3"),
    #                  "colour": PY_16_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": pythia_18_all_L1L2L3, 
    #                  "label": label("py18", "all", "L1L2L3"),
    #                  "colour": PY_18_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='pThat', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "pthat_pythia_16vs18_L1L2L3.pdf"),
    #             plot_kwargs=dict(xtitle="#hat{p}_{T} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="CP5 / CUETP8M1", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": pythia_16_all_L1, 
    #                  "label": label("py16", "all", "L1"),
    #                  "colour": PY_16_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  }, 
    #                 {"file": pythia_18_all_L1, 
    #                  "label": label("py18", "all", "L1"),
    #                  "colour": PY_18_COL,
    #                  'line_style': 1,
    #                  'line_width': lw,
    #                  "subplot_ind": 0,
    #                  }
    #             ],
    #             histname='pThat', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "pthat_pythia_16vs18_L1.pdf"),
    #             plot_kwargs=dict(xtitle="#hat{p}_{T} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="CP5 / CUETP8M1", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, "label": label("hpp", "all", "L1L2L3"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_L1L2L3, "label": label("h7", "all", "L1L2L3"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, "subplot_ind": 0, },
    #                 {"file": pythia_16_all_L1L2L3, "label": label("py16", "all", "L1L2L3"), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": pythia_18_all_L1L2L3, "label": label("py18", "all", "L1L2L3"), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw, "subplot_ind": 2, }
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_PYvsHERWIG_16vs18_L1L2L3.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="2018 / 2016", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, "label": label("hpp", "all", "L1L2L3"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_L1L2L3, "label": label("h7", "all", "L1L2L3"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": pythia_16_all_L1L2L3, "label": label("py16", "all", "L1L2L3"), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #                 {"file": pythia_18_all_L1L2L3, "label": label("py18", "all", "L1L2L3"), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw, "subplot_ind": 1, }
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_PYvsHERWIG_L1L2L3.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="Pythia / Herwig", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1, "label": label("hpp", "all", "L1"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_L1, "label": label("h7", "all", "L1"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, "subplot_ind": 0, },
    #                 {"file": pythia_16_all_L1, "label": label("py16", "all", "L1"), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": pythia_18_all_L1, "label": label("py18", "all", "L1"), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw, "subplot_ind": 2, }
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_PYvsHERWIG_16vs18_L1.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="2018 / 2016", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1, "label": label("hpp", "all", "L1"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_L1, "label": label("h7", "all", "L1"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": pythia_16_all_L1, "label": label("py16", "all", "L1"), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #                 {"file": pythia_18_all_L1, "label": label("py18", "all", "L1"), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw, "subplot_ind": 1, }
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_PYvsHERWIG_L1.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="Pythia / Herwig", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_NoPU, "label": label("hpp", "all", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_NoPU, "label": label("h7", "all", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_HERWIG_NoPU.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1L2L3, "label": label("hpp", "all", "L1L2L3"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_L1L2L3, "label": label("h7", "all", "L1L2L3"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_HERWIG_L1L2L3.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1, "label": label("hpp", "all", "L1"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_L1, "label": label("h7", "all", "L1"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_HERWIG_L1.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )


    # do_1d_plot(entries=[
    #                 {"file": herwigpp_ud_NoPU, "label": label("hpp", "ud", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_ud_NoPU, "label": label("h7", "ud", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_HERWIG_NoPU_ud.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_g_NoPU, "label": label("hpp", "g", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_g_NoPU, "label": label("h7", "g", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_HERWIG_NoPU_g.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_ud_NoPU, "label": label("hpp", "ud", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwigpp_g_NoPU, "label": label("hpp", "g", ""), "colour": HPP_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 0},
    #                 {"file": herwig7_ud_NoPU, "label": label("h7", "ud", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw,  },
    #                 {"file": herwig7_g_NoPU, "label": label("h7", "g", ""), "colour": H7_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 2 },
    #             ],
    #             histname='refPtAll', logx=True, logy=True,
    #             output_filename=os.path.join(OUTPUT_DIR, "refPtAll_HERWIG_NoPU_udVSg.pdf"),
    #             plot_kwargs=dict(xtitle="p_{T}^{GenJet} [GeV]", ytitle="p.d.f", title=None, ylim=pt_ylim,
    #                              subplot_type='ratio', subplot_title="g / ud", subplot_limits=None)
    # )
    

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_L1, "label": label("hpp", "all", "L1"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_L1, "label": label("h7", "all", "L1"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refEtaLeading',
    #             output_filename=os.path.join(OUTPUT_DIR, "refEtaLeading_HERWIG_L1.pdf"),
    #             plot_kwargs=dict(xtitle="#eta^{Leading GenJet}", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )


    # do_1d_plot(entries=[
    #                 {"file": herwigpp_ud_NoPU, "label": label("hpp", "ud", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_ud_NoPU, "label": label("h7", "ud", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refEtaLeading',
    #             output_filename=os.path.join(OUTPUT_DIR, "refEtaLeading_HERWIG_NoPU_ud.pdf"),
    #             plot_kwargs=dict(xtitle="#eta^{Leading GenJet}", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_g_NoPU, "label": label("hpp", "g", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_g_NoPU, "label": label("h7", "g", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'subplot_ind': 0 },
    #             ],
    #             histname='refEtaLeading',
    #             output_filename=os.path.join(OUTPUT_DIR, "refEtaLeading_HERWIG_NoPU_g.pdf"),
    #             plot_kwargs=dict(xtitle="#eta^{Leading GenJet}", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_ud_NoPU, "label": label("hpp", "ud", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwigpp_g_NoPU, "label": label("hpp", "g", ""), "colour": HPP_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 0},
    #                 {"file": herwig7_ud_NoPU, "label": label("h7", "ud", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw,  },
    #                 {"file": herwig7_g_NoPU, "label": label("h7", "g", ""), "colour": H7_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 2 },
    #             ],
    #             histname='refEtaLeading',
    #             output_filename=os.path.join(OUTPUT_DIR, "refEtaLeading_HERWIG_NoPU_udVSg.pdf"),
    #             plot_kwargs=dict(xtitle="#eta^{Leading GenJet}", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="g / ud", subplot_limits=None)
    # )

    # do_1d_plot(entries=[
    #                 {"file": herwigpp_all_NoPU, "label": label("hpp", "all", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": herwig7_all_NoPU, "label": label("h7", "all", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw, "subplot_ind": 0, }
    #             ],
    #             histname='nRef',
    #             output_filename=os.path.join(OUTPUT_DIR, "nRef_16vs18_NoPU.pdf"),
    #             plot_kwargs=dict(xtitle="# GenJets (p_{T}^{Gen} > 10 GeV)", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=None)
    # )
    # 
  
    # do_1d_plot(entries=[
    #                 {"file": pythia_16_ud_L1L2L3, "label": label("py16", "ud", ""), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, },
    #                 {"file": pythia_16_g_L1L2L3, "label": label("py16", "g", ""), "colour": PY_16_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 0},
    #                 {"file": pythia_18_ud_L1L2L3, "label": label("py18", "ud", ""), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw,  },
    #                 {"file": pythia_18_g_L1L2L3, "label": label("py18", "g", ""), "colour": PY_18_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 2 },
    #             ],
    #             histname='refEtaLeading',
    #             output_filename=os.path.join(OUTPUT_DIR, "refEtaLeading_Pythia_L1L2L3_udVSg.pdf"),
    #             plot_kwargs=dict(xtitle="#eta^{Leading GenJet} [GeV]", ytitle="p.d.f", title=None,
    #                              subplot_type='ratio', subplot_title="g / ud", subplot_limits=None)
    # )
    
    for start, end in [(30, 35), (120, 150), (400, 550), (1000, 1500)]:
        rebin = 4
        if start > 200:
            rebin=2
        # make_projection_plot(entries=[
        #                 {"file": herwigpp_ud_NoPU, "label": label("hpp", "ud", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
        #                 {"file": herwigpp_g_NoPU, "label": label("hpp", "g", ""), "colour": HPP_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 0, 'rebin': rebin},
        #                 {"file": herwig7_ud_NoPU, "label": label("h7", "ud", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw,  'rebin': rebin},
        #                 {"file": herwig7_g_NoPU, "label": label("h7", "g", ""), "colour": H7_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 2, 'rebin': rebin},
        #             ],
        #             histname='RelRspVsRefPt', start_val=start, end_val=end,
        #             output_filename=os.path.join(OUTPUT_DIR, "RelRsp_Pt%gto%g_Herwig_gVSud_NoPU.pdf" % (start, end)),
        #             plot_kwargs=dict(xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", xlim=(0.5, 1.5), title="%g < p_{T}^{Gen} < %g GeV" % (start, end),
        #                              subplot_type='ratio', subplot_title="g / ud", subplot_limits=[0, 2])
        # )
        # make_projection_plot(entries=[
        #                 {"file": herwigpp_ud_NoPU, "label": label("hpp", "ud", ""), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
        #                 {"file": herwigpp_g_NoPU, "label": label("hpp", "g", ""), "colour": HPP_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin},
        #                 {"file": herwig7_ud_NoPU, "label": label("h7", "ud", ""), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin, 'subplot_ind': 0},
        #                 {"file": herwig7_g_NoPU, "label": label("h7", "g", ""), "colour": H7_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 1, 'rebin': rebin},
        #             ],
        #             histname='RelRspVsRefPt', start_val=start, end_val=end,
        #             output_filename=os.path.join(OUTPUT_DIR, "RelRsp_Pt%gto%g_Herwig_16vs18_NoPU.pdf" % (start, end)),
        #             plot_kwargs=dict(xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", xlim=(0.5, 1.5), title="%g < p_{T}^{Gen} < %g GeV" % (start, end),
        #                              subplot_type='ratio', subplot_title="H7 / H++", subplot_limits=[0, 2])
        # )

        make_projection_plot(entries=[
                        {"file": herwigpp_ud_L1, "label": label("hpp", "ud", "L1"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
                        {"file": herwigpp_g_L1, "label": label("hpp", "g", "L1"), "colour": HPP_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 0, 'rebin': rebin},
                        {"file": herwig7_ud_L1, "label": label("h7", "ud", "L1"), "colour": H7_COL, 'line_style': 1, 'line_width': lw,  'rebin': rebin},
                        {"file": herwig7_g_L1, "label": label("h7", "g", "L1"), "colour": H7_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 2, 'rebin': rebin},
                        {"file": pythia_16_ud_L1, "label": label("py16", "ud", "L1"), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
                        {"file": pythia_16_g_L1, "label": label("py16", "g", "L1"), "colour": PY_16_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 4, 'rebin': rebin},
                        {"file": pythia_18_ud_L1, "label": label("py18", "ud", "L1"), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw,  'rebin': rebin},
                        {"file": pythia_18_g_L1, "label": label("py18", "g", "L1"), "colour": PY_18_COL, 'line_style': 2, 'line_width': lw, 'subplot_ind': 6, 'rebin': rebin},
                    ],
                    histname='RelRspVsRefPt', start_val=start, end_val=end,
                    output_filename=os.path.join(OUTPUT_DIR, "RelRsp_Pt%gto%g_PythiaHerwig_gVSud_L1_eta2p5to2p964.pdf" % (start, end)),
                    plot_kwargs=dict(xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", xlim=(0.5, 1.5), title="%g < p_{T}^{Gen} < %g GeV" % (start, end),
                                     subplot_type='ratio', subplot_title="g / ud", subplot_limits=[0, 2])
        )
        make_projection_plot(entries=[
                        {"file": herwigpp_g_L1, "label": label("hpp", "g", "L1"), "colour": HPP_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin},
                        {"file": herwig7_g_L1, "label": label("h7", "g", "L1"), "colour": H7_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin, 'subplot_ind': 0},
                        {"file": pythia_16_g_L1, "label": label("py16", "g", "L1"), "colour": PY_16_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin},
                        {"file": pythia_18_g_L1, "label": label("py18", "g", "L1"), "colour": PY_18_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin, 'subplot_ind':2},
                    ],
                    histname='RelRspVsRefPt', start_val=start, end_val=end,
                    output_filename=os.path.join(OUTPUT_DIR, "RelRsp_Pt%gto%g_PythiaHerwig_16vs18_g_L1_eta2p5to2p964.pdf" % (start, end)),
                    plot_kwargs=dict(xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", xlim=(0.5, 1.5), title="%g < p_{T}^{Gen} < %g GeV" % (start, end),
                                     subplot_type='ratio', subplot_title="2018 / 2016", subplot_limits=[0, 2])
        )
        make_projection_plot(entries=[
                        {"file": herwigpp_ud_L1, "label": label("hpp", "ud", "L1"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
                        {"file": herwig7_ud_L1, "label": label("h7", "ud", "L1"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin, 'subplot_ind': 0},
                        {"file": pythia_16_ud_L1, "label": label("py16", "ud", "L1"), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
                        {"file": pythia_18_ud_L1, "label": label("py18", "ud", "L1"), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw,  'rebin': rebin, 'subplot_ind': 2},
                    ],
                    histname='RelRspVsRefPt', start_val=start, end_val=end,
                    output_filename=os.path.join(OUTPUT_DIR, "RelRsp_Pt%gto%g_PythiaHerwig_16vs18_ud_L1_eta2p5to2p964.pdf" % (start, end)),
                    plot_kwargs=dict(xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", xlim=(0.5, 1.5), title="%g < p_{T}^{Gen} < %g GeV" % (start, end),
                                     subplot_type='ratio', subplot_title="2018 / 2016", subplot_limits=[0, 2])
        )
        for ecf in ["chf", "nhf", "cef", "nef"]: #, "chmult", "nmult"]:
            xlim = None
            if "mult" in ecf:
                rebin = 2
                xlim = [0, 80]
            make_projection_plot(entries=[
                            {"file": herwigpp_g_L1, "label": label("hpp", "g", "L1"), "colour": HPP_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin},
                            {"file": herwig7_g_L1, "label": label("h7", "g", "L1"), "colour": H7_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin, 'subplot_ind': 0},
                            {"file": pythia_16_g_L1, "label": label("py16", "g", "L1"), "colour": PY_16_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin},
                            {"file": pythia_18_g_L1, "label": label("py18", "g", "L1"), "colour": PY_18_COL, 'line_style': 2, 'line_width': lw, 'rebin': rebin, 'subplot_ind':2},
                        ],
                        histname='Ref%sVsRefPt' % ecf, start_val=start, end_val=end,
                        output_filename=os.path.join(OUTPUT_DIR, "Ref%s_Pt%gto%g_PythiaHerwig_16vs18_g_L1_eta2p5to2p964.pdf" % (ecf, start, end)),
                        plot_kwargs=dict(xtitle="%s" % ecf.upper(), title="%g < p_{T}^{Gen} < %g GeV" % (start, end), xlim=xlim,
                                         subplot_type='ratio', subplot_title="2018 / 2016", subplot_limits=[0, 2])
            )
            make_projection_plot(entries=[
                            {"file": herwigpp_ud_L1, "label": label("hpp", "ud", "L1"), "colour": HPP_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
                            {"file": herwig7_ud_L1, "label": label("h7", "ud", "L1"), "colour": H7_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin, 'subplot_ind': 0},
                            {"file": pythia_16_ud_L1, "label": label("py16", "ud", "L1"), "colour": PY_16_COL, 'line_style': 1, 'line_width': lw, 'rebin': rebin},
                            {"file": pythia_18_ud_L1, "label": label("py18", "ud", "L1"), "colour": PY_18_COL, 'line_style': 1, 'line_width': lw,  'rebin': rebin, 'subplot_ind': 2},
                        ],
                        histname='Ref%sVsRefPt' % ecf, start_val=start, end_val=end,
                        output_filename=os.path.join(OUTPUT_DIR, "Rel%s_Pt%gto%g_PythiaHerwig_16vs18_ud_L1_eta2p5to2p964.pdf" % (ecf, start, end)),
                        plot_kwargs=dict(xtitle="%s" % ecf.upper(), title="%g < p_{T}^{Gen} < %g GeV" % (start, end), xlim=xlim,
                                         subplot_type='ratio', subplot_title="2018 / 2016", subplot_limits=[0, 2])
            )