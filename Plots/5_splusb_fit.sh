#!/bin/bash
# UNBLIND step 1f: S+B (and B-only) fit to OBSERVED data per mA via FitDiagnostics.
#   Gives the best-fit signal strength r-hat (S+B fit) with its uncertainty, plus the
#   post-fit signal+background shapes and per-bin uncertainties (saveShapes). The B-only
#   fit (shapes_fit_b) and S+B fit (shapes_fit_s) are both stored in fitDiagnostics*.root.
#   MH frozen at 125.38, r in [-3,3]. mA1 uses the R=1 working-point datacard.
#   NOTE: r range is [-3,3] not [-5,5]: with the discrete-profiling multipdf bkg a
#   few mA (e.g. mA11/mA12) find a spurious S+B minimum at the |r|=5 boundary
#   (mA11 status=-1, mA12 r=+5.000); narrowing to ±3 recovers the true r-hat~=0
#   (identical at ±2 and ±3). All 30 mA have |r-hat|<0.11, so ±3 is ample.
# Output -> roots: Combine/output_splusb_fit_observed/fitDiagnostics{mA}_sb.root
#           r-hat: Plots/plot_limits/6_splusb_fit/mA{NN}_rhat.txt   (best-fit r summary)
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
# Each mA runs in its own tmp workdir so the concurrent combines don't clobber combine_logger.out.
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T2W=$FFIT/Combine/root_t2w
OUT=$FFIT/Combine/output_splusb_fit_observed           # FitDiagnostics roots stay here
PLOTDIR=$FFIT/Plots/plot_limits/6_splusb_fit           # human-facing r-hat txt (execution-order #6)
mkdir -p "$OUT" "$OUT/logs_parallel" "$PLOTDIR"

ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"
PDFIDX="${PDFIDX:-pdfindex_cat0_13p6TeV}"

# exit 0 iff the FitDiagnostics root has a usable S+B fit (tree_fit_sb + shapes_fit_s)
has_fit_s() {
    python3 - "$1" <<'PY'
import ROOT, sys
ROOT.gErrorIgnoreLevel = ROOT.kError
try:
    f = ROOT.TFile(sys.argv[1]); t = f.Get("tree_fit_sb")
    sys.exit(0 if (t and t.GetEntries() > 0 and f.Get("shapes_fit_s")) else 1)
except Exception:
    sys.exit(1)
PY
}

run_one() {
    local mA="$1"
    local ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1_makeLimits_observed.sh first)"; return 0; }
    echo "=== S+B FitDiagnostics mA = ${mA} ==="
    local mApad; printf -v mApad "%02d" "$mA"
    local wd; wd=$(mktemp -d "$OUT/work_${mA}_XXXX")
    ( cd "$wd"
      combine -M FitDiagnostics "$ws" -m 125.38 \
          --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 \
          --cminDefaultMinimizerStrategy 0 \
          --saveShapes --saveWithUncertainties --saveNormalizations -n ${mA}_sb
      # Fallback: a few mA (e.g. mA4) fail the profiled multipdf S+B fit ("Fit failed",
      # no tree_fit_sb/shapes_fit_s). Freeze the envelope index to its best-fit function
      # and refit. (--robustHesse does NOT help -- it makes the fit fail outright.)
      if ! has_fit_s "fitDiagnostics${mA}_sb.root"; then
          echo "[fallback] mA=${mA}: profiled S+B fit failed; refit with frozen best pdfindex"
          combine -M MultiDimFit "$ws" -m 125.38 --setParameters MH=125.38 --freezeParameters MH \
              --rMin -3 --rMax 3 --cminDefaultMinimizerStrategy 0 --saveSpecifiedIndex "$PDFIDX" \
              -n ${mA}_sbbest >/dev/null 2>&1 || true
          bi=$(python3 - "higgsCombine${mA}_sbbest.MultiDimFit.mH125.38.root" "$PDFIDX" <<'PY'
import ROOT, sys
ROOT.gErrorIgnoreLevel = ROOT.kError
try:
    f = ROOT.TFile(sys.argv[1]); t = f.Get("limit"); t.GetEntry(0)
    print(int(getattr(t, sys.argv[2])))
except Exception:
    print(0)
PY
)
          bi=${bi:-0}
          combine -M FitDiagnostics "$ws" -m 125.38 \
              --setParameters MH=125.38,${PDFIDX}=${bi} --freezeParameters MH,${PDFIDX} --rMin -3 --rMax 3 \
              --cminDefaultMinimizerStrategy 0 \
              --saveShapes --saveWithUncertainties --saveNormalizations -n ${mA}_sb
      fi
      mv "fitDiagnostics${mA}_sb.root" "$OUT/" 2>/dev/null || true )
    python3 - "$OUT/fitDiagnostics${mA}_sb.root" "${mA}" > "$PLOTDIR/mA${mApad}_rhat.txt" <<'PY'
import sys, ROOT
f=ROOT.TFile(sys.argv[1]); mA=sys.argv[2]
t=f.Get("tree_fit_sb")
if t and t.GetEntries()>0:
    t.GetEntry(0)
    print(f"mA={mA}: r-hat = {t.r:+.4f}  +{abs(t.rHiErr):.4f} / -{abs(t.rLoErr):.4f}  (fit_status={t.fit_status})")
else:
    print(f"mA={mA}: tree_fit_sb missing/empty")
PY
    cat "$PLOTDIR/mA${mApad}_rhat.txt"
    rm -rf "$wd"
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > "$OUT/logs_parallel/splusb_${mA}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait
echo "[1f] S+B fits -> $OUT/fitDiagnostics{mA}_sb.root ; r-hat -> $PLOTDIR/mA{NN}_rhat.txt (logs: $OUT/logs_parallel/)"
cat "$PLOTDIR"/*_rhat.txt 2>/dev/null | sort -t= -k2 -V
