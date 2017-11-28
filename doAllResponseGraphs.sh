#!/bin/bash -e

./plotFlavourResponseGraphs.py --inputGraphs jra_QCD_FLAT_NoJEC_f_g.root --outputDir QCD_FLAT_NoJEC --title "Without JEC"
./plotFlavourResponseGraphs.py --inputGraphs jra_QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4_f_g.root --outputDir QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4 --title "Summer16_23Sep2016V4"
./plotFlavourResponseGraphs.py --inputGraphs jra_QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8_f_g.root --outputDir QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8 --title "Summer16_03Feb2017_V8"
