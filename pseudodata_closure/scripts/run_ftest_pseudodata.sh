#!/bin/bash
# Refit the background envelope on the DY-MC pseudo-data workspaces (sculpting closure).
# Mirrors shellScripts/bkg/fit_bkg.sh's fTest_ALP_turnOn call, but:
#   input  = pseudodata_closure/pseudodata/mA_M<mA>/ws/run3.root  (MC-as-data)
#   output = pseudodata_closure/fit_results_run3_pseudodata/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root
# Reuses the already-built ./bin/fTest_ALP_turnOn (no rebuild). Parallel NPROC<=6.
#   ALPMASS="1 2 3" to restrict; GOF_TOYS=<N> (default 100, same as real).
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
source $FFIT/setup.sh >/dev/null 2>&1
cd $FFIT/Background
PC=$FFIT/pseudodata_closure
GOF_TOYS="${GOF_TOYS:-100}"
NPROC="${NPROC:-6}"

massList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && massList=( ${ALPMASS} )
mkdir -p $PC/logs

run_one() {
    local mA="$1"
    local WS=$PC/pseudodata/mA_M${mA}/ws/run3.root
    local OUT=$PC/fit_results_run3_pseudodata/${mA}
    if [[ ! -f "$WS" ]]; then echo "[SKIP mA=$mA] missing $WS"; return 1; fi
    mkdir -p $OUT/HZAmassInde_fTest
    ./bin/fTest_ALP_turnOn -i $WS \
        --saveMultiPdf $OUT/CMS-HGG_mva_13p6TeV_multipdf.root -D $OUT/HZAmassInde_fTest \
        --mass_ALP ${mA} -c 1 --isFlashgg 0 --isData 0 -f data, \
        --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys "$GOF_TOYS" \
        > $OUT/ftest.stdout.log 2>&1
    local st=$?
    local bf=$(grep -E "Best Fit Pdf =" $OUT/ftest.stdout.log | tail -1 | sed 's/.*Best Fit Pdf = //')
    echo "[mA=$mA] exit=$st  bestFitPdf=${bf:-NA}"
}
export -f run_one
export PC

for mA in "${massList[@]}"; do
    run_one "$mA" > $PC/logs/ftest_${mA}.log 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 3; done
done
wait
echo "[ftest-pseudodata] done -> $PC/fit_results_run3_pseudodata/*/CMS-HGG_mva_13p6TeV_multipdf.root"
grep -h "\[mA=" $PC/logs/ftest_*.log 2>/dev/null | sort -t= -k2 -n
