#!/bin/bash
# Observed LOCAL significance (discovery q0) for each mA. Uses real data (unblinded).
# Same MH/r treatment as the limits. Output -> output_significance/.
cmsenv

ALPmassList=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
outdir=output_significance
mkdir -p "$outdir"

for mA in "${ALPmassList[@]}"; do
    echo "Processing SIGNIFICANCE for mA = ${mA}"
    path_datacard="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/${mA}_pruned_datacard_leptons.txt"
    combine -M Significance "$path_datacard" --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 -n "${mA}"
    mv "higgsCombine${mA}.Significance.mH125.38.root" "$outdir"/ 2>/dev/null
done
