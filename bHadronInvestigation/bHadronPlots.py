#!/usr/bin/env python

"""Do b hadron plots"""


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

PDGID_STR = {
    11: "e",
    12: "#nu_{e}",
    13: "#mu",
    14: "#nu_{#mu}",
    15: "#tau",
    16: "#nu_{#tau}",
    21: "g",
    22: "#gamma",
    111: "#pi^{0}",
    113: "#rho^{0}",
    213: "#rho^{+}",
    221: "#eta",
    223: "#omega",
    310: "K^{0}_{S}",
    311: "K^{0}",
    313: "K^{*0}",
    323: "K^{*+}",
    331: "#eta^{'}",
    333: "#phi",
    411: "D^{+}",
    413: "D^{*+}",
    415: "D_{2}^{*+}",
    421: "D^{0}",
    423: "D^{*0}",
    425: "D_{2}^{*0}",
    431: "D_{s}^{+}",
    433: "D_{s}^{*+}",
    435: "D_{s2}^{*+}",
    441: "#eta_{C}",
    443: "J/#psi",
    445: "#chi_{c2}",
    2212: "p",
    2214: "#Delta^{+}",
    3212: "#Sigma^{0}",
    4112: "#Sigma_{c}^{0}",
    4114: "#Sigma_{c}^{*0}",
    4122: "#Lambda_{c}^{+}",
    4124: "4124",
    4132: "#Xi_{c}^{0}",
    4212: "#Sigma_{c}^{+}",
    4214: "#Sigma_{c}^{*+}",
    4222: "#Sigma_{c}^{++}",
    4224: "#Sigma_{c}^{*++}",
    4232: "#Xi_{c}^{+}",
    4312: "#Xi_{c}^{'0}",
    4314: "#Xi_{c}^{*0}",
    4322: "#Xi_{c}^{'+}",
    4324: "#Xi_{c}^{*0}",
    4332: "#Xi_{c}^{*+}",
    4334: "#Omega_{c}^{*0}",
    511: "B^{0}",
    521: "B^{+}",
    531: "B_{s}^{0}",
    541: "B_{c}^{+}",
    551: "#eta_{b}",
    553: "#Upsilon",
    555: "#chi_{b2}",
    5122: "#Lambda_{b}^{0}",
    5132: "#Xi_{b}^{-}",
    5232: "#Xi_{b}^{0}",
    5332: "#Omega_{b}^{-}",
}


def make_plot(entries, output_filename, pt_bin, plot_kwargs, is_pdgid_plot=False, logy=False):
    start_val, end_val = 750, 1000
    pt_bin = pt_bin.lower()
    if pt_bin == "low":
        start_val, end_val = 57, 90
    elif pt_bin != "high":
        raise RuntimeError("pt_bin must be low or high")
    
    conts = []
    if not is_pdgid_plot:
        for ent in entries:
            h2d = cu.get_from_tfile(ent['file'], "ak4pfchsl1/"+ent['histname']+"_JetEta0to0.783")
            hist = cu.get_projection_plot(h2d, start_val, end_val, 'x')
            contrib = Contribution(hist, label=ent['label'],
                                   line_width=lw, line_color=ent['colour'], line_style=ent.get('line_style', 1),
                                   marker_size=0, marker_color=ent['colour'], marker_style=1,
                                   normalise_hist=True, rebin_hist=ent.get('rebin', None))
            conts.append(contrib)
    
    else:
        custom_bins = []
        for ent in entries:
            h2d = cu.get_from_tfile(ent['file'], "ak4pfchsl1/"+ent['histname']+"_JetEta0to0.783")
            hist = cu.get_projection_plot(h2d, start_val, end_val, 'x')
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
                                   marker_size=0, marker_color=ent['colour'], marker_style=1,
                                   normalise_hist=True, rebin_hist=ent.get('rebin', None))
            conts.append(contrib)

    for ind, ent in enumerate(entries):
        if ent.get('subplot_ind', -1) >= 0:
            conts[ind].subplot = conts[ent['subplot_ind']].obj

    plot = Plot(conts, what="hist", ytitle="p.d.f.", has_data=False, **plot_kwargs)
    if is_pdgid_plot and conts[0].obj.GetNbinsX() > 10:
        plot.default_canvas_size = (800, 800)
    plot.legend.SetX1(0.5)
    plot.legend.SetX2(0.97)
    if len(entries) > 4:
        plot.legend.SetY1(0.6)
        plot.legend.SetNColumns(2)
    else:
        plot.legend.SetY1(0.75)
    plot.legend.SetY2(0.9)
    plot.plot("HISTE NOSTACK")
    if logy:
        plot.set_logy()
    plot.save(output_filename)



if __name__ == "__main__":
    pythia_highpt_file = cu.open_root_file("BHadron_ak4pfchsl1_pythia_miniaod_highPt.root")
    pythia_lowpt_file = cu.open_root_file("BHadron_ak4pfchsl1_pythia_miniaod_lowPt.root")
    herwig_highpt_file = cu.open_root_file("BHadron_ak4pfchsl1_herwig_miniaod_all_highPt.root")
    herwig_lowpt_file = cu.open_root_file("BHadron_ak4pfchsl1_herwig_miniaod_all_lowPt.root")
    herwig_highpt_file = cu.open_root_file("BHadron_ak4pfchsl1_herwig_miniaod_ext1234_highPt.root")
    herwig_lowpt_file = cu.open_root_file("BHadron_ak4pfchsl1_herwig_miniaod_ext1234_lowPt.root")

    hist_dicts = [
        {
            "histname": "NHadronsVsRefPt", 
            "title": "# B hadron in jet"
        },
        {
            "histname": "BJetRefVsRefPt_AtLeast2Hadron", 
            "title": "B jet index (2+ B hadrons)"
        },
        {
            "histname": "BJetRefVsRefPt_SingleHadron", 
            "title": "B jet index (1 B hadron)"
        },
        {
            "histname": "JtcefVsRefPt", 
            "title": "Electron EF",
            "rebin": 5,
            "logy": True,
        },
        {
            "histname": "JtchfVsRefPt", 
            "title": "Charged Hadron EF",
            "rebin": 5,
            "logy": True,
        },
        {
            "histname": "JthfefVsRefPt", 
            "title": "HFEF",
            "rebin": 5,
            "logy": True,
        },
        {
            "histname": "JthfhfVsRefPt", 
            "title": "HFHF",
            "rebin": 5,
            "logy": True,
        },
        {
            "histname": "JtmufVsRefPt", 
            "title": "Muon EF",
            "rebin": 5,
            "logy": True,
        },
        {
            "histname": "JtnefVsRefPt", 
            "title": "Photon EF",
            "rebin": 5,
            "logy": True,
        },
        {
            "histname": "JtnhfVsRefPt", 
            "title": "Neutral Hadron EF",
            "rebin": 5,
            "logy": True,
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_FirstHadron", 
            "title": "Hadron decay PDGID [1st B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_FirstHadron_Lepton", 
            "title": "Hadron decay PDGID [1st B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_FirstHadron_Neutrino", 
            "title": "Hadron decay PDGID [1st B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_SecondHadron", 
            "title": "Hadron decay PDGID [2nd B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_SecondHadron_Lepton", 
            "title": "Hadron decay PDGID [2nd B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_SecondHadron_Neutrino", 
            "title": "Hadron decay PDGID [2nd B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_SingleHadron", 
            "title": "Hadron decay PDGID [single B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_SingleHadron_Lepton", 
            "title": "Hadron decay PDGID [single B]"
        },
        {
            "histname": "RefHadronDecayPdgidVsRefPt_SingleHadron_Neutrino", 
            "title": "Hadron decay PDGID [single B]"
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_FirstHadron", 
            "title": "Hadron decay p_{T}/p_{T}^{B} [1st B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_FirstHadron_Lepton", 
            "title": "Hadron decay p_{T}^{l}/p_{T}^{B} [1st B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_FirstHadron_Neutrino", 
            "title": "Hadron decay p_{T}^{#nu}/p_{T}^{B} [1st B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_SecondHadron", 
            "title": "Hadron decay p_{T}/p_{T}^{B} [2nd B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_SecondHadron_Lepton", 
            "title": "Hadron decay p_{T}^{l}/p_{T}^{B} [2nd B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_SecondHadron_Neutrino", 
            "title": "Hadron decay p_{T}^{#nu}/p_{T}^{B} [2nd B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_SingleHadron", 
            "title": "Hadron decay p_{T}/p_{T}^{B} [single B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_SingleHadron_Lepton", 
            "title": "Hadron decay p_{T}^{l}/p_{T}^{B} [single B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioVsRefPt_SingleHadron_Neutrino", 
            "title": "Hadron decay p_{T}^{#nu}/p_{T}^{B} [single B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronNdecayVsRefPt_FirstHadron", 
            "title": "# hadron decay products [1st B]"
        },
        {
            "histname": "RefHadronNdecayVsRefPt_SecondHadron", 
            "title": "# hadron decay products [2nd B]"
        },
        {
            "histname": "RefHadronNdecayVsRefPt_SingleHadron", 
            "title": "# hadron decay products [single B]"
        },
        {
            "histname": "RefHadronPdgidVsRefPt_FirstHadron", 
            "title": "Hadron PDGID [1st B]"
        },
        {
            "histname": "RefHadronPdgidVsRefPt_SecondHadron", 
            "title": "Hadron PDGID [2nd B]"
        },
        {
            "histname": "RefHadronPdgidVsRefPt_SingleHadron", 
            "title": "Hadron PDGID [single B]"
        },
        {
            "histname": "RefHadronPtRatioVsRefPt_FirstHadron", 
            "title": "Hadron p_{T}/p_{T}^{Gen} [1st B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronPtRatioVsRefPt_SecondHadron", 
            "title": "Hadron p_{T}/p_{T}^{Gen} [2nd B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronPtRatioVsRefPt_SingleHadron", 
            "title": "Hadron p_{T}/p_{T}^{Gen} [single B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronSldecayVsRefPt_FirstHadron", 
            "title": "Hadron semileptonic decay [1st B]"
        },
        {
            "histname": "RefHadronSldecayVsRefPt_SecondHadron", 
            "title": "Hadron semileptonic decay [2nd B]"
        },
        {
            "histname": "RefHadronSldecayVsRefPt_SingleHadron", 
            "title": "Hadron semileptonic decay [single B]"
        },
        {
            "histname": "RelRspVsRefPt_AtLeast2Hadron_HadDecay", 
            "title": "p_{T}^{reco}/p_{T}^{Gen} [2+ B hadrons, hadronic decay]",
            "rebin": 5,
        },
        {
            "histname": "RelRspVsRefPt_AtLeast2Hadron", 
            "title": "p_{T}^{reco}/p_{T}^{Gen} [2+ B hadrons]",
            "rebin": 5,
        },
        {
            "histname": "RelRspVsRefPt_AtLeast2Hadron_SLDecay", 
            "title": "p_{T}^{reco}/p_{T}^{Gen} [2+ B hadrons, semileptonic decay]",
            "rebin": 5,
        },
        {
            "histname": "RelRspVsRefPt", 
            "title": "p_{T}^{reco}/p_{T}^{Gen}",
            "rebin": 5,
        },
        {
            "histname": "RelRspVsRefPt_SingleHadron_HadDecay", 
            "title": "p_{T}^{reco}/p_{T}^{Gen} [1 B hadron, hadronic decay]",
            "rebin": 5,
        },
        {
            "histname": "RelRspVsRefPt_SingleHadron", 
            "title": "p_{T}^{reco}/p_{T}^{Gen} [1 B hadron]",
            "rebin": 5,
        },
        {
            "histname": "RelRspVsRefPt_SingleHadron_SLDecay", 
            "title": "p_{T}^{reco}/p_{T}^{Gen} [1 B hadron, semileptonic decay]",
            "rebin": 5,
        },
    ]
    
    hist_dicts = [
        {
            "histname": "RefHadronDecayPtRatioJetVsRefPt_SingleHadron", 
            "title": "Hadron decay p_{T}/p_{T}^{GenJet} [single B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioJetVsRefPt_SingleHadron_Lepton", 
            "title": "Hadron decay p_{T}^{l}/p_{T}^{GenJet} [single B]",
            "rebin": 5,
        },
        {
            "histname": "RefHadronDecayPtRatioJetVsRefPt_SingleHadron_Neutrino", 
            "title": "Hadron decay p_{T}^{#nu}/p_{T}^{GenJet} [single B]",
            "rebin": 5,
        },
    ]

    high_pt_setup = [pythia_highpt_file, herwig_highpt_file, True]
    low_pt_setup = [pythia_lowpt_file, herwig_lowpt_file, False]
    
    for hdict in hist_dicts:
        # if 'pdgid' not in hdict['histname'].lower():
        #     continue
        for pythia_file, herwig_file, is_highPt in [high_pt_setup, low_pt_setup]:

            pt_bin_str = "high" if is_highPt else "low"
            entries_dicts = [
                {
                    "file": pythia_file,
                    "histname": hdict['histname'],
                    "label": "PYTHIA8 [%s p_{T}]" % pt_bin_str,
                    "colour": ROOT.kBlue,
                    "line_style": 1,
                    'rebin': hdict.get('rebin', None),
                },
                {
                    "file": herwig_file,
                    "histname": hdict['histname'],
                    "label": "HERWIG++ [%s p_{T}]" % pt_bin_str,
                    "colour": ROOT.kRed,
                    "line_style": 1,
                    "subplot_ind": 0,
                    'rebin': hdict.get('rebin', None),
                }
            ]
            plot_kwargs = dict(
                title=None, xtitle=hdict['title'],
                xlim=None, ylim=None,
                subplot_type="ratio", 
                subplot_title="Herwig++/Pythia8", subplot_limits=None
            )
            make_plot(entries_dicts, is_pdgid_plot="pdgid" in hdict['histname'].lower(),
                      output_filename=hdict['histname']+"_%sPt.pdf" % pt_bin_str, 
                      pt_bin=pt_bin_str, logy=hdict.get('logy', False), plot_kwargs=plot_kwargs)

