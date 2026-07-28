#!/bin/bash
# UNBLIND step: 1D profile-likelihood scan of the signal strength r (OBSERVED) per mA.
#   combine -M MultiDimFit --algo grid : -2 dLL vs r over [-5,5]; plot1DScan.py -> pdf
#   (gives r-hat and the 68%/95% CL intervals from the scan). MH frozen at 125.38.
#   Needs root_t2w (built by 1_makeLimits_observed.sh). mA1 = R=1 working-point datacard.
# Output -> root: Combine/output_rscan_observed/{mA}_rscan.root ; pdf: Plots/plot_limits/8_rscan/mA{NN}_rscan.pdf
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
# Each mA runs in its own tmp workdir (avoids concurrent combine_logger.out clobber).
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T2W=$FFIT/Combine/root_t2w
OUT=$FFIT/Combine/output_rscan_observed                # MultiDimFit scan roots stay here
PLOTDIR=$FFIT/Plots/plot_limits/8_rscan                # human-facing scan pdf (execution-order #8)
mkdir -p "$OUT" "$OUT/logs_parallel" "$PLOTDIR"

POINTS="${POINTS:-50}"
ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"

# RooMultiPdf index (envelope discrete nuisance). Frozen to its best-fit function
# during the scan: discrete profiling at fixed r is unstable and returns combine's
# 9990 fail-sentinel (2*dNLL=19980) at every grid point -> degenerate 1-point scan.
PDFIDX="${PDFIDX:-pdfindex_cat0_13p6TeV}"
OBSDIR=$FFIT/Combine/output_combine_results_observed

run_one() {
    local mA="$1"
    local ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1_makeLimits_observed.sh first)"; return 0; }
    local mApad; printf -v mApad "%02d" "$mA"
    local wd; wd=$(mktemp -d "$OUT/work_${mA}_XXXX")
    # Adaptive r window from the observed 95% CL limit. These signals are tightly
    # constrained (sigma_r ~ limit/2), so a fixed [-5,5] grid never samples the
    # minimum -> the scan is uninformative. rlo=-1.5*limit, rhi=3.5*limit spans a
    # few sigma each side of r_hat~0 with grid spacing << sigma_r.
    local obsroot=$OBSDIR/higgsCombine${mA}.AsymptoticLimits.mH125.38.root
    local rlo rhi
    read rlo rhi < <(python3 - "$obsroot" <<'PY'
import ROOT,sys,os
ROOT.gErrorIgnoreLevel=ROOT.kError
fn=sys.argv[1]; obs=0.05
if os.path.exists(fn):
    f=ROOT.TFile(fn); t=f.Get("limit")
    v=[e.limit for e in t] if t else []
    f.Close()
    if v: obs=v[5] if len(v)>5 else v[-1]
if obs<=0: obs=0.05
print(f"{-1.5*obs:.5f} {3.5*obs:.5f}")
PY
)
    echo "=== r scan mA = ${mA} (${POINTS} pts, r in [${rlo}, ${rhi}]) ==="
    ( cd "$wd"
      # best-fit pdf index (float r over the adaptive window)
      combine -M MultiDimFit "$ws" -m 125.38 --setParameters MH=125.38 --freezeParameters MH \
          --rMin "$rlo" --rMax "$rhi" --cminDefaultMinimizerStrategy 0 \
          --saveSpecifiedIndex "$PDFIDX" -n ${mA}_best >/dev/null 2>&1
      local bestidx
      bestidx=$(python3 - "higgsCombine${mA}_best.MultiDimFit.mH125.38.root" "$PDFIDX" <<'PY'
import ROOT,sys
ROOT.gErrorIgnoreLevel=ROOT.kError
try:
    f=ROOT.TFile(sys.argv[1]); t=f.Get("limit"); t.GetEntry(0)
    print(int(getattr(t, sys.argv[2])))
except Exception:
    print(0)
PY
)
      bestidx=${bestidx:-0}
      combine -M MultiDimFit "$ws" --algo grid --points "$POINTS" -m 125.38 \
          --setParameters MH=125.38,${PDFIDX}=${bestidx} --freezeParameters MH,${PDFIDX} \
          --rMin "$rlo" --rMax "$rhi" --cminDefaultMinimizerStrategy 0 -n ${mA}_rscan
      cp -f "higgsCombine${mA}_rscan.MultiDimFit.mH125.38.root" "$OUT/" 2>/dev/null || true
      # plot1DScan.py mangles an absolute -o into ".//afs/..." (unwritable from the tmp
      # cwd); write with a relative name here, then copy the pdf/png to PLOTDIR.
      if plot1DScan.py "higgsCombine${mA}_rscan.MultiDimFit.mH125.38.root" --POI r -o "mA${mApad}_rscan"; then
          cp -f "mA${mApad}_rscan.pdf" "mA${mApad}_rscan.png" "$PLOTDIR/" 2>/dev/null || true
      else
          echo "[warn] plot1DScan failed for mA=${mA} (scan root still saved)"
      fi )
    rm -rf "$wd"
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > "$OUT/logs_parallel/rscan_${mA}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait
echo "[scan] r scans -> root: $OUT/{mA}_rscan.root ; pdf: $PLOTDIR/mA{NN}_rscan.pdf (logs: $OUT/logs_parallel/)"
