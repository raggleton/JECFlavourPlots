#!/bin/bash -e
set -u

PYTHIA8_TEXT="PYTHIA8 Flat QCD"
./plotFlavourResponseHists.py --input jra_QCD_FLAT_NoJEC_f.root --outputDir QCD_FLAT_NoJEC --title "Without JEC" --sampleName "$PYTHIA8_TEXT"
./plotFlavourResponseHists.py --input jra_QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4_f.root --outputDir QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4 --title "Summer16_23Sep2016V4" --sampleName "$PYTHIA8_TEXT"
./plotFlavourResponseHists.py --input jra_QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8_f.root --outputDir QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8 --title "Summer16_03Feb2017_V8" --sampleName "$PYTHIA8_TEXT"

HERWIG_TEXT="HERWIG++ Flat QCD"
./plotFlavourResponseHists.py --input jra_QCD_FLAT_HERWIG_NoJEC_ak4pfchs_newFlav_f.root --outputDir QCD_FLAT_HERWIG_NoJEC --title "Without JEC" --sampleName "$HERWIG_TEXT"
./plotFlavourResponseHists.py --input jra_QCD_FLAT_HERWIG_withL1L2L3_Summer16_03Feb2017_V8_ak4pfchs_newFlav_f.root --outputDir QCD_FLAT_HERWIG_withL1L2L3_Summer16_03Feb2017_V8 --title "Summer16_03Feb2017_V8" --sampleName "$HERWIG_TEXT"
