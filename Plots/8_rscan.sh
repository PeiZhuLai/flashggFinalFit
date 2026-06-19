#!/bin/bash
# UNBLIND step: 1D profile-likelihood scan of the signal strength r (OBSERVED) per mA.
#   combine -M MultiDimFit --algo grid : -2 dLL vs r over [-5,5]; plot1DScan.py -> pdf
#   (gives r-hat and the 68%/95% CL intervals from the scan). MH frozen at 125.38.
#   Needs root_t2w (built by 1_makeLimits_observed.sh). mA1 = R=1 working-point datacard.
# Output -> Combine/output_rscan_observed/{mA}_rscan.{root,pdf}
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
# Each mA runs in its own tmp workdir (avoids concurrent combine_logger.out clobber).
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T2W=$FFIT/Combine/root_t2w
OUT=$FFIT/Combine/output_rscan_observed
mkdir -p "$OUT" "$OUT/logs_parallel"

POINTS="${POINTS:-50}"
ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"

run_one() {
    local mA="$1"
    local ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1_makeLimits_observed.sh first)"; return 0; }
    echo "=== r scan mA = ${mA} (${POINTS} points) ==="
    local wd; wd=$(mktemp -d "$OUT/work_${mA}_XXXX")
    ( cd "$wd"
      combine -M MultiDimFit "$ws" --algo grid --points "$POINTS" -m 125.38 \
          --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 \
          --cminDefaultMinimizerStrategy 0 -n ${mA}_rscan
      cp -f "higgsCombine${mA}_rscan.MultiDimFit.mH125.38.root" "$OUT/" 2>/dev/null || true
      plot1DScan.py "higgsCombine${mA}_rscan.MultiDimFit.mH125.38.root" --POI r -o "$OUT/${mA}_rscan" \
          || echo "[warn] plot1DScan failed for mA=${mA} (scan root still saved)" )
    rm -rf "$wd"
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > "$OUT/logs_parallel/rscan_${mA}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait
echo "[scan] r profile-likelihood scans -> $OUT/{mA}_rscan.{root,pdf} (logs: $OUT/logs_parallel/)"
