#!/bin/bash
# Background-fit plots (same style as Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults)
# for the DY-MC pseudo-data closure test.
#
# Mirrors step 2 of shellScripts/bkg/fit_bkg.sh (makeBkgPlots_ALP), but reads the
# pseudo-data multipdf refit by scripts/run_ftest_pseudodata.sh:
#   input  = pseudodata_closure/fit_results_run3_pseudodata/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root
#   output = pseudodata_closure/fit_results_run3_pseudodata/{AllFitResults_blind,AllFitResults_unblind}/
#            {bkgplot,allPdfs}_<mA>.pdf   (+ per-mass BkgPlots_{blind,unblind}/ and .root)
# Reuses the already-built ./bin/makeBkgPlots_ALP (no rebuild). Parallel NPROC<=6 (CERN rule).
#   ALPMASS="1 2 3" to restrict.
#
# Note on the fit window: all mass points are built and plotted on [95,180]
# (massLow=95 for all mA, select_lib.py; mA1,2 restored from 105 to 95 on 2026-07-30 now
# that the low-mass sculpting is controlled). This matches the real-data AllFitResults.
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
source $FFIT/setup.sh >/dev/null 2>&1
cd $FFIT/Background
PC=$FFIT/pseudodata_closure
LUMI="${LUMI:-172.13}"
# Band sampling step. Default 2.5 matches the real-data AllFitResults; it can be
# coarsened (e.g. MASSSTEP=10) for a mass point whose band scan is prohibitively
# slow. This affects only how finely the uncertainty band is sampled, not the fit.
MASSSTEP="${MASSSTEP:-2.5}"
NPROC="${NPROC:-6}"
# Legend text for the points: these are DY simulation treated as pseudo-data, not real data.
# Needs the --dataLegend option of makeBkgPlots_ALP (added 2026-07-28; default is "Data",
# so real-data output is unchanged).
DATALEG="${DATALEG:-MC pseudo-data}"
RATIOTITLE="${RATIOTITLE:-MC - Fit}"

if [[ ! -x ./bin/makeBkgPlots_ALP ]]; then
    echo "[ERROR] ./bin/makeBkgPlots_ALP missing -- build Background once first"; exit 1
fi

RES=$PC/fit_results_run3_pseudodata
OUT_BLIND=$RES/AllFitResults_blind
OUT_UNBLIND=$RES/AllFitResults_unblind
mkdir -p $OUT_BLIND $OUT_UNBLIND $PC/logs

massList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && massList=( ${ALPMASS} )

run_one() {
    local mA="$1"
    local DIR=$RES/${mA}
    local MPDF=$DIR/CMS-HGG_mva_13p6TeV_multipdf.root
    if [[ ! -f "$MPDF" ]]; then echo "[SKIP mA=$mA] missing $MPDF"; return 1; fi

    local mhLow=95

    local common=( -b "$MPDF" --sqrts 13p6TeV --isMultiPdf --useBinnedData
                   --massStep "$MASSSTEP" --mhVal 125.0 --maVal "$mA"
                   --mhLow "$mhLow" --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135
                   --intLumi "$LUMI" -c 0 --isFlashgg 0 --doBands
                   --dataLegend "$DATALEG" --ratioYTitle "$RATIOTITLE" )

    ./bin/makeBkgPlots_ALP "${common[@]}" \
        -d $DIR/BkgPlots_blind   --total_OutDir $OUT_BLIND   -o $DIR/BkgPlots_blind.root \
        > $DIR/makeBkgPlots.stdout.log 2>&1
    local st1=$?
    ./bin/makeBkgPlots_ALP "${common[@]}" --unblind \
        -d $DIR/BkgPlots_unblind --total_OutDir $OUT_UNBLIND -o $DIR/BkgPlots_unblind.root \
        >> $DIR/makeBkgPlots.stdout.log 2>&1
    local st2=$?
    echo "[mA=$mA] mhLow=$mhLow blind_exit=$st1 unblind_exit=$st2"
}
export -f run_one
export RES OUT_BLIND OUT_UNBLIND LUMI DATALEG RATIOTITLE

for mA in "${massList[@]}"; do
    run_one "$mA" > $PC/logs/bkgplots_${mA}.log 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 3; done
done
wait

echo "[bkgplots-pseudodata] done"
echo "  blind   -> $OUT_BLIND"
echo "  unblind -> $OUT_UNBLIND"
grep -h "\[mA=" $PC/logs/bkgplots_*.log 2>/dev/null | sort -t= -k2 -n
