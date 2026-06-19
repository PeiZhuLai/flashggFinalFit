#!/bin/bash
# UNBLIND step 1b: observed LOCAL significance (discovery q0) per mA on real data.
#   combine -M Significance, MH frozen at 125.38, r in [-5,5].
#   Output -> Combine/output_significance/higgsCombine{mA}.Significance.mH125.38.root
# mA1 uses the R=1 working point datacard (its spurious 5.78 sigma is removed).
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
DC=$FFIT/Datacard/output_Datacard_leptons
cd $FFIT/Combine
mkdir -p output_significance logs_parallel

ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"

run_one() {
    local mA="$1"
    echo "=== LOCAL significance mA = ${mA} ==="
    combine -M Significance $DC/${mA}_pruned_datacard_leptons.txt \
        --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 -n "${mA}"
    mv "higgsCombine${mA}.Significance.mH125.38.root" output_significance/ 2>/dev/null || true
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > logs_parallel/significance_${mA}.log 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait
echo "[1b] local significance -> $FFIT/Combine/output_significance/ (per-mA logs: Combine/logs_parallel/significance_*.log)"
