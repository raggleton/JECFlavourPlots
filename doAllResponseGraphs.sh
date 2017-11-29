#!/bin/bash -e
set -u

PYTHIA_TEXT="PYTHIA8 Flat QCD"
./plotFlavourResponseGraphs.py --input jra_QCD_FLAT_NoJEC_f_g.root --outputDir QCD_FLAT_NoJEC --title "Without JEC" --sampleName "$PYTHIA_TEXT"
./plotFlavourResponseGraphs.py --input jra_QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4_f_g.root --outputDir QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4 --title "Summer16_23Sep2016V4" --sampleName "$PYTHIA_TEXT"
./plotFlavourResponseGraphs.py --input jra_QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8_f_g.root --outputDir QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8 --title "Summer16_03Feb2017_V8" --sampleName "$PYTHIA_TEXT"

HERWIG_TEXT="HERWIG++ Flat QCD"
# ./plotFlavourResponseGraphs.py --input jra_QCD_FLAT_HERWIG_NoJEC_ak4pfchs_newFlav_f_g.root --outputDir QCD_FLAT_HERWIG_NoJEC --title "Without JEC" --sampleName "$HERWIG_TEXT"
# ./plotFlavourResponseGraphs.py --input jra_QCD_FLAT_HERWIG_withL1L2L3_Summer16_03Feb2017_V8_ak4pfchs_newFlav_f_g.root --outputDir QCD_FLAT_HERWIG_withL1L2L3_Summer16_03Feb2017_V8 --title "Summer16_03Feb2017_V8" --sampleName "$HERWIG_TEXT"
