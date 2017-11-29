#!/bin/bash -e
set -u

# Do all old/new flavour comparison plots

./compareFlavourResponseGraphs.py \
--input jra_QCD_FLAT_NoJEC_f_g.root --label "New flavour" \
--input jra_QCD_FLAT_NoJEC_oldFlav_f_g.root --label "Old flavour" \
--outputDir QCD_FLAT_NoJEC/comparisonOldNewFlav \
--title "Without JEC"

./compareFlavourResponseGraphs.py \
--input jra_QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4_f_g.root --label "New flavour" \
--input jra_QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4_oldFlav_f_g.root --label "Old flavour" \
--outputDir QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4/comparisonOldNewFlav \
--title "Summer16_23Sep2016V4"

./compareFlavourResponseGraphs.py \
--input jra_QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8_f_g.root --label "New flavour" \
--input jra_QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8_oldFlav_f_g.root --label "Old flavour" \
--outputDir QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8/comparisonOldNewFlav \
--title "Summer16_03Feb2017_V8"