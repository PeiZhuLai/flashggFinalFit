#!/bin/bash
# UNBLIND step 1d: per-mA goodness-of-fit on the UNBLINDED (observed) data.
#   combine -M GoodnessOfFit --algo saturated : observed test statistic + toy distribution.
#   GOF p-value = fraction of toys with t_toy >= t_obs.
# Output combine roots -> Plots/unblind_GOF/roots/ ; plots via plot_gof_unblind.py
# -> Plots/unblind_GOF/.  NTOYS overridable (default 500).
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T2W=$FFIT/Combine/root_t2w
OUT=$FFIT/Plots/unblind_GOF/roots
mkdir -p "$OUT" "$OUT/../logs_parallel"
cd "$OUT"

NTOYS="${NTOYS:-500}"
ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"

run_one() {
    local mA="$1"
    local ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1a first)"; return 0; }
    echo "=== GOF mA = ${mA} (saturated, ${NTOYS} toys) ==="
    combine -M GoodnessOfFit "$ws" --algo saturated -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH -n ${mA}_obs
    combine -M GoodnessOfFit "$ws" --algo saturated -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH -t ${NTOYS} -s 12345 -n ${mA}_toys
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > "$OUT/../logs_parallel/gof_${mA}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait

echo "[1d] GOF roots -> $OUT (per-mA logs: Plots/unblind_GOF/logs_parallel/) ; now plotting..."
cd $FFIT/Plots
source /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh
conda activate higgs-alp-ana
python3 plot_gof_unblind.py
