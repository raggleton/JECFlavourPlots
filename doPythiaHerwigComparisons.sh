#!/bin/bash -e

# Compare response
./compareFlavourResponseGraphs.py \
--input jra_QCD_FLAT_NoJEC_f_g.root --label "PYTHIA8" \
--input jra_QCD_FLAT_HERWIG_NoJEC_ak4pfchs_newFlav_f_g.root --label "HERWIG++" \
--outputDir QCD_FLAT_HERWIG_NoJEC/comparisonPythiaHerwig \
--title "Without JEC"

./compareFlavourResponseGraphs.py \
--input jra_QCD_FLAT_withL1L2L3_Summer16_03Feb2017_V8_f_g.root --label "PYTHIA8" \
--input jra_QCD_FLAT_HERWIG_withL1L2L3_Summer16_03Feb2017_V8_ak4pfchs_newFlav_f_g.root --label "HERWIG++" \
--outputDir QCD_FLAT_HERWIG_withL1L2L3_Summer16_03Feb2017_V8NoJEC/comparisonPythiaHerwig \
--title "Summer_03Feb2017_V8"
