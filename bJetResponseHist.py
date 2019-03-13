#!/usr/bin/env python


"""Quickly dump response hists for specific flavour, pt, eta bin"""


from __future__ import print_function
import os
import sys
import argparse
from copy import deepcopy

import ROOT
from MyStyle import My_Style
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
My_Style.cd()


FONT_SIZE = 0.032


f_hpp = cu.open_root_file("QCD_Pt_Herwigpp_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_L1FastJet_Summer16_07Aug2017_V20_MC/jra_QCD_Pt_15to7000_Herwigpp_NoJEC_newFlav_ak4pfchs_rspRangeLarge_absEta_hadronParton.root")
f_hpp = cu.open_root_file("QCD_Pt_Pythia8_CP5_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_L1FastJet_Autumn18_V3_MC/jra_QCD_Pt_15toInf_NoJEC_newFlav_ak4pfchs_rspRangeLarge_absEta_hadronParton.root")
f_h7 = cu.open_root_file("QCD_Pt_Herwig7_NoJEC_relPtHatCut2p5_jtptmin4_HadronParton_nbinsrelrsp_10k_L1FastJet_Autumn18_V3_MC/jra_QCD_Pt_15to7000_Herwig7_NoJEC_newFlav_ak4pfchs_rspRangeLarge_absEta_hadronParton.root")

for flav in ['b', 'ud', 'g', 'c', 's']:
    for ptmin, ptmax in [(750, 1000), (72, 90)]:
        ptbin = "%dto%d" % (ptmin, ptmax)
        histname = "ak4pfchsl1/%s_RelRsp_JetEta0to0.783_RefPt%s" % (flav, ptbin)
        hpp = cu.get_from_tfile(f_hpp, histname)
        h7 = cu.get_from_tfile(f_h7, histname)

        entries = [
        (hpp, {"label":"Pythia8 (CP5)", "line_color": ROOT.kAzure+10, "line_width": 2}),
        (h7, {"label":"Herwig7", "line_color": ROOT.kBlue, "line_width": 2, "subplot": hpp})
        ]

        plot_kwargs = {
            "title": "0 < |#eta| < 0.783\n%d < p_{T}^{Gen} < %d GeV\n%s jets (+ L1FastJet JEC)" %(ptmin, ptmax, flav),
            "xlim": [0.5, 1.5],
            "xtitle": "Response (p_{T}^{Reco} / p_{T}^{Gen})",
            "subplot_type": "ratio",
            # "subplot_title": "H7 / H++",
            "subplot_title": "H7 / Py8",
        }
        cu.do_comparison_plot(entries, output_filename="h7_vs_py8_jec_%s_%s.pdf" % (ptbin, flav), rebin=25, **plot_kwargs)