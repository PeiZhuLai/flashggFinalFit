#!/bin/bash
# UNBLIND step 1f: S+B (and B-only) fit to OBSERVED data per mA via FitDiagnostics.
#   Gives the best-fit signal strength r-hat (S+B fit) with its uncertainty, plus the
#   post-fit signal+background shapes and per-bin uncertainties (saveShapes). The B-only
#   fit (shapes_fit_b) and S+B fit (shapes_fit_s) are both stored in fitDiagnostics*.root.
#   MH frozen at 125.38, r in [-5,5]. mA1 uses the R=1 working-point datacard.
# Output -> Combine/output_splusb_fit_observed/fitDiagnostics{mA}_sb.root
#           Combine/output_splusb_fit_observed/{mA}_rhat.txt   (best-fit r summary)
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
# Each mA runs in its own tmp workdir so the concurrent combines don't clobber combine_logger.out.
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T2W=$FFIT/Combine/root_t2w
OUT=$FFIT/Combine/output_splusb_fit_observed
mkdir -p "$OUT" "$OUT/logs_parallel"

ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"

run_one() {
    local mA="$1"
    local ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1_makeLimits_observed.sh first)"; return 0; }
    echo "=== S+B FitDiagnostics mA = ${mA} ==="
    local wd; wd=$(mktemp -d "$OUT/work_${mA}_XXXX")
    ( cd "$wd"
      combine -M FitDiagnostics "$ws" -m 125.38 \
          --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 \
          --cminDefaultMinimizerStrategy 0 \
          --saveShapes --saveWithUncertainties --saveNormalizations -n ${mA}_sb
      mv "fitDiagnostics${mA}_sb.root" "$OUT/" 2>/dev/null || true )
    python3 - "$OUT/fitDiagnostics${mA}_sb.root" "${mA}" > "$OUT/${mA}_rhat.txt" <<'PY'
import sys, ROOT
f=ROOT.TFile(sys.argv[1]); mA=sys.argv[2]
t=f.Get("tree_fit_sb")
if t and t.GetEntries()>0:
    t.GetEntry(0)
    print(f"mA={mA}: r-hat = {t.r:+.4f}  +{abs(t.rHiErr):.4f} / -{abs(t.rLoErr):.4f}  (fit_status={t.fit_status})")
else:
    print(f"mA={mA}: tree_fit_sb missing/empty")
PY
    cat "$OUT/${mA}_rhat.txt"
    rm -rf "$wd"
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > "$OUT/logs_parallel/splusb_${mA}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait
echo "[1f] S+B fits -> $OUT/fitDiagnostics{mA}_sb.root ; r-hat -> $OUT/{mA}_rhat.txt (logs: $OUT/logs_parallel/)"
cat "$OUT"/*_rhat.txt 2>/dev/null | sort -t= -k2 -V
