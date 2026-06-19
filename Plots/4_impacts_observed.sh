#!/bin/bash
# UNBLIND step 1e: OBSERVED impact plots + best-fit r-hat per mA (real SR data).
#   combineTool.py -M Impacts: doInitialFit -> doFits (one fit per nuisance) -> collect
#   json -> plotImpacts.py. The header of the impact plot shows r = observed (the best-fit
#   signal strength r-hat with its uncertainty), and the ranked nuisance impacts on r.
#   MH frozen at 125.38, r in [-5,5], robustFit. NO -t / --expectSignal (this is OBSERVED
#   data, not an Asimov/expected dataset). mA1 uses the R=1 working-point datacard.
# Output -> Combine/output_impacts_observed/{mA}_impacts.{json,pdf}
# TWO parallel axes: NPROC mA-workers (default 2) x PARALLEL nuisance-fits each (default 2).
# Keep NPROC*PARALLEL modest on shared lxplus. ALPMASS="1 2 3" to restrict.
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T2W=$FFIT/Combine/root_t2w
OUT=$FFIT/Combine/output_impacts_observed
mkdir -p "$OUT" "$OUT/logs_parallel"

PARALLEL="${PARALLEL:-2}"   # combineTool --parallel (nuisance fits per mA)
NPROC="${NPROC:-2}"          # mA workers
ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )

common="-m 125.38 --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 --robustFit 1 --cminDefaultMinimizerStrategy 0"

run_one() {
    local mA="$1"
    local ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1_makeLimits_observed.sh first)"; return 0; }
    echo "=== OBSERVED impacts mA = ${mA} ==="
    local workdir; workdir=$(mktemp -d "$OUT/work_${mA}_XXXX")
    ( cd "$workdir"
      combineTool.py -M Impacts -d "$ws" $common --doInitialFit
      combineTool.py -M Impacts -d "$ws" $common --doFits --parallel "$PARALLEL"
      combineTool.py -M Impacts -d "$ws" $common -o "$OUT/${mA}_impacts.json"
      plotImpacts.py -i "$OUT/${mA}_impacts.json" -o "$OUT/${mA}_impacts" )
    rm -rf "$workdir"
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > "$OUT/logs_parallel/impacts_${mA}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait
echo "[1e] observed impacts (r=observed) -> $OUT/{mA}_impacts.{json,pdf} (logs: $OUT/logs_parallel/)"
