#!/bin/bash
# Scan the mA=2 BDT working point: for each candidate cut, build MC pseudo-data, refit the
# background envelope, and run the closure (observed-on-MC significance + obs/exp 95% CL limit).
# Shows how the fake-signal significance and sensitivity trade off vs the cut, to pick a WP.
# Usage: bash scan_mA2_wp.sh "0.975 0.978 0.980 0.982 0.985 0.988 0.990 0.992"
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
source $FFIT/setup.sh >/dev/null 2>&1
CUTS="${1:-0.975 0.978 0.980 0.982 0.985 0.988 0.990 0.992}"
realDC=$FFIT/Datacard/output_Datacard_leptons/2_pruned_datacard_leptons.txt
SCAN=$PC/scan_mA2; mkdir -p $SCAN
printf "%-8s %-8s %-9s %-9s %-9s %-7s\n" cut R_mc obs_lim exp_lim signifZ pseudoN | tee $SCAN/scan_summary.txt

for cut in $CUTS; do
    tag=$(echo $cut | tr -d '.')
    OB=$SCAN/cut_${tag}
    # 1) pseudo-data ws MUST be pre-built (conda ROOT, NO cmsenv -- cling clashes otherwise)
    WS=$OB/mA_M2/ws/run3.root
    if [[ ! -f "$WS" ]]; then echo "[SKIP $cut] missing pre-built ws $WS"; continue; fi
    N=$(grep -oE "pseudo-data N=[0-9]+" $OB.build.log 2>/dev/null | grep -oE "[0-9]+" | head -1)
    # 2) refit envelope
    FIT=$OB/fit; mkdir -p $FIT/ft
    ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $WS \
        --saveMultiPdf $FIT/CMS-HGG_mva_13p6TeV_multipdf.root -D $FIT/ft \
        --mass_ALP 2 -c 1 --isFlashgg 0 --isData 0 -f data, \
        --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 30 ) > $OB.ftest.log 2>&1
    # 3) datacard swap -> combine
    psDC=$OB/datacard.txt
    sed "s#/Background/ALP_BkgModel_ReReco/fit_results_run3/2/#/pseudodata_closure/scan_mA2/cut_${tag}/fit/#g; s#CMS-HGG_mva_13p6TeV_multipdf.root multipdf:roohist_data_mass_cat0#CMS-HGG_mva_13p6TeV_multipdf.root multipdf:roohist_data_mass_cat0#g" "$realDC" > "$psDC"
    text2workspace.py "$psDC" -o $OB/ws_t2w.root > $OB.t2w.log 2>&1
    combine -M AsymptoticLimits $OB/ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n Scan${tag} > $OB.lim.log 2>&1
    combine -M Significance $OB/ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n Scan${tag} > $OB.sig.log 2>&1
    rm -f higgsCombineScan${tag}.*.root 2>/dev/null
    obs=$(grep "Observed Limit" $OB.lim.log | sed -E 's/.*r < //')
    exp=$(grep "Expected 50.0%" $OB.lim.log | sed -E 's/.*r < //')
    Z=$(grep -E "Significance:" $OB.sig.log | sed -E 's/.*Significance: //')
    printf "%-8s %-8s %-9s %-9s %-9s %-7s\n" "$cut" "-" "${obs:-NA}" "${exp:-NA}" "${Z:-NA}" "${N:-NA}" | tee -a $SCAN/scan_summary.txt
done
echo "[scan] done -> $SCAN/scan_summary.txt"
