#!/usr/bin/env python

"""Do H7 vs H++ vs Py8 plots"""


from __future__ import print_function
import os
import sys
import argparse
from copy import deepcopy

import ROOT
sys.path.append("..")
from MyStyle import My_Style
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
My_Style.cd()


FONT_SIZE = 0.032

lw = 2


def do_2d_plot(h2d, output_filename, renorm=None, logx=False, logy=False, logz=False, xtitle=None, ytitle=None, rebinx=1, rebiny=1):
    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetLogx(logx)
    canv.SetLogy(logy)
    canv.SetLogz(logz)
    h2d_new = h2d.Clone(cu.get_unique_str())
    h2d_new.RebinX(rebinx)
    h2d_new.RebinY(rebiny)
    if renorm.lower() in ['x', 'y']:
        h2d_new = cu.make_normalised_TH2(h2d_new, renorm, True)
    if xtitle:
        h2d_new.GetXaxis().SetTitle(xtitle);
    if ytitle:
        h2d_new.GetYaxis().SetTitle(ytitle);
    h2d_new.Draw("COLZ")
    dirname = os.path.dirname(os.path.abspath(output_filename))
    cu.check_dir_exists_create(dirname)
    canv.SaveAs(output_filename)


def do_all_2d_plots(tfile, output_dir):
    rebin = 2
    hist_dicts = [
    {
    "hname": "JtcefVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jtcef",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtchfVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jtchf",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtchmultVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jtchmult",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JthfefVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jthfef",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JthfhfVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jthfhf",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtmufVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jtmuf",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtnefVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jtnef",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtnhfVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jtnhf",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtnmultVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Jtnmult",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtRefcefRatioVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "JtRefcefRatio",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtRefchfRatioVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "JtRefchfRatio",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtRefchmultRatioVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "JtRefchmultRatio",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtRefmufRatioVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "JtRefmufRatio",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtRefnefRatioVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "JtRefnefRatio",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtRefnhfRatioVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "JtRefnhfRatio",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "JtRefnmultRatioVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "JtRefnmultRatio",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "RefcefVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Refcef",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "RefchfVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Refchf",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "RefchmultVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Refchmult",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "RefmufVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Refmuf",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "RefnefVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Refnef",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "RefnhfVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Refnhf",
    "rebinx": rebin,
    "rebiny": rebin,
    },
    {
    "hname": "RefnmultVsRsp_JetEta0to0.783",
    "renorm": "X",
    "xtitle": "Response",
    "ytitle": "Refnmult",
    "rebinx": rebin,
    "rebiny": rebin,
    },



    ]
    for hdict in hist_dicts:
        for renorm in ['', 'X', 'Y']:
            h2d = cu.get_from_tfile(tfile, "ak4pfchs/"+hdict['hname'])
            ofile = output_dir+"/"+hdict.get('ofile', hdict['hname'])
            if renorm == "X":
                ofile += "_renormX"
            elif renorm == "Y":
                ofile += "_renormY"
            ofile +=".pdf"
            do_2d_plot(h2d, output_filename=ofile,
                       renorm=renorm,
                       logz=hdict.get('logz', False),
                       xtitle=hdict.get('xtitle', ''), ytitle=hdict.get('ytitle', ''),
                       rebinx=hdict.get('rebinx', 1), rebiny=hdict.get('rebiny', 1),
                       )



if __name__ == "__main__":
    PY8_CP5_DIR = "QCD_Pt_18_NoJEC_relPtHatCut2p5_jtptmin4_PhysicsParton_nbinsrelrsp_10k_Autumn18_V3_MC_L1FastJetL2L3"
    do_all_2d_plots(
        cu.open_root_file(os.path.join(PY8_CP5_DIR, "Closure_ak4pfchs_QCD_Pt_15toInf_Py8_CP5_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin750_ptGenMax1000_nrefmax2_b.root")),
        os.path.join(PY8_CP5_DIR, "Closure_ak4pfchs_QCD_Pt_15toInf_Py8_CP5_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin750_ptGenMax1000_nrefmax2_b_plots_2d")
    )
    do_all_2d_plots(
        cu.open_root_file(os.path.join(PY8_CP5_DIR, "Closure_ak4pfchs_QCD_Pt_15toInf_Py8_CP5_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin72_ptGenMax90_nrefmax2_b.root")),
        os.path.join(PY8_CP5_DIR, "Closure_ak4pfchs_QCD_Pt_15toInf_Py8_CP5_NoJEC_newFlav__relRspMax2_relPtHatMax2p5_ptGenMin72_ptGenMax90_nrefmax2_b_plots_2d")
    )

    H7_DIR = "QCD_Pt_Herwig7_NoJEC_relPtHatCut2p5_jtptmin4_PhysicsParton_nbinsrelrsp_10k_Autumn18_V3_MC_L1FastJetL2L3"
    do_all_2d_plots(
        cu.open_root_file(os.path.join(H7_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin750_ptGenMax1000_nrefmax2_b_0.root")),
        os.path.join(H7_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin750_ptGenMax1000_nrefmax2_b_plots_2d")
    )
    do_all_2d_plots(
        cu.open_root_file(os.path.join(H7_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin72_ptGenMax90_nrefmax2_b_0.root")),
        os.path.join(H7_DIR, "Closure_ak4pfchsQCD_Pt_15to7000_Herwig7_NoJEC_miniaod__relRspMax2_relPtHatMax2p5_ptGenMin72_ptGenMax90_nrefmax2_b_plots_2d")
    )
