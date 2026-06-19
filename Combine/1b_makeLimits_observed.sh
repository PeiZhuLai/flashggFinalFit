#!/bin/bash
# Observed (UNBLINDED) AsymptoticLimits for all mA.
# Same combine config as 1_makeLimits.sh (freeze MH=125.38, rMax 2) but WITHOUT
# --run blind, so combine returns observed + expected quantiles using real SR data.
# Output goes to output_combine_results_observed/ so the blind-expected results in
# output_combine_results/ are NOT overwritten.
cmsenv

ALPmassList=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )

outdir=output_combine_results_observed
mkdir -p "$outdir"

for mA in "${ALPmassList[@]}"; do
    echo "Processing OBSERVED for mA = ${mA}"
    path_datacard="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/${mA}_pruned_datacard_leptons.txt"
    combine -M AsymptoticLimits "$path_datacard" --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n "${mA}"
    mv "higgsCombine${mA}.AsymptoticLimits.mH125.38.root" "$outdir"/ 2>/dev/null
done
