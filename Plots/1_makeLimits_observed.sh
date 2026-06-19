#!/bin/bash
# UNBLIND step 1a: observed (UNBLINDED) 95% CL upper limits per mA.
#   text2workspace -> combine -M AsymptoticLimits (observed, real SR data).
#   MH frozen at 125.38, rMax 2 (avoids the spurious rMax/MH interplay; see
#   doc/HZa/limit_mA2_spurious_rMax_bug.md). Output -> Combine/output_combine_results_observed/.
# mA1 uses the R=1 working point datacard; mA2..30 the production cut.
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
DC=$FFIT/Datacard/output_Datacard_leptons
cd $FFIT/Combine

ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"
mkdir -p root_t2w output_combine_results_observed logs_parallel

run_one() {
    local mA="$1"
    echo "=== OBSERVED limit mA = ${mA} ==="
    text2workspace.py $DC/${mA}_pruned_datacard_leptons.txt -o root_t2w/${mA}_Datacard_leptons.root
    combine -M AsymptoticLimits root_t2w/${mA}_Datacard_leptons.root \
        --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n "${mA}"
    mv "higgsCombine${mA}.AsymptoticLimits.mH125.38.root" output_combine_results_observed/ 2>/dev/null || true
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > logs_parallel/makeLimits_${mA}.log 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait
echo "[1a] observed limits -> $FFIT/Combine/output_combine_results_observed/ (per-mA logs: Combine/logs_parallel/makeLimits_*.log)"
