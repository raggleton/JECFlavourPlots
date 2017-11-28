#!/bin/bash -e
# So all old/new corrections comparison plots

./compareFlavourResponseGraphs.py \
--input jra_QCD_FLAT_withL1L2L3_Summer16_23Sep2016V4_f_g.root --label "23Sep2016V4" \
--input jra_QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8_f_g.root --label "03Feb2017_V8" \
--outputDir QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8/comparison_Summer16_23Sep2016V4 \
--title "With Summer16_* JEC"
