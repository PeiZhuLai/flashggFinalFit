#!/bin/bash
# UNBLIND step: S+B post-fit MODEL PLOT per mA (data + S+B curve + B component) built from
# the FitDiagnostics output of 5_splusb_fit.sh (shapes_fit_s). PyROOT plotter. Two sets:
#   _unblind (data shown everywhere) and _blind (data blanked in 115-135 GeV).
#   NOTE: 5_splusb_fit.sh's Hesse cannot force a pos-def covariance with the discrete-
#   profiling multipdf bkg (fit_status=300), so --saveWithUncertainties is dropped and these
#   curves carry NO post-fit uncertainty band (central S+B / B curves only).
# Output -> Plots/plot_limits/7_splusb_modelplot/mA{NN}_splusb_model_{unblind,blind}.pdf
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SB=$FFIT/Combine/output_splusb_fit_observed              # FitDiagnostics roots (input) stay here
T2W=$FFIT/Combine/root_t2w                               # combine workspaces (post-fit pdf source)
PLOTDIR=$FFIT/Plots/plot_limits/7_splusb_modelplot       # human-facing model pdf (execution-order #7)
mkdir -p "$PLOTDIR"
LUMI="172.13 fb^{-1}"
ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )

cd "$FFIT/Plots"
for mA in "${ALPmassList[@]}"; do
    fd=$SB/fitDiagnostics${mA}_sb.root
    ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$fd" ] || { echo "[skip] missing $fd (run 5_splusb_fit.sh first)"; continue; }
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1_makeLimits_observed.sh first)"; continue; }
    printf -v mApad "%02d" "$mA"
    echo "=== S+B model plot mA = ${mA} (blind + unblind) ==="
    python3 plot_splusb_model.py "$mA" "$ws" "$fd" "$PLOTDIR/mA${mApad}_splusb_model_unblind.pdf" "$LUMI" 0
    python3 plot_splusb_model.py "$mA" "$ws" "$fd" "$PLOTDIR/mA${mApad}_splusb_model_blind.pdf"   "$LUMI" 1 115 135
done
echo "[plot] S+B post-fit model plots -> $PLOTDIR/mA{NN}_splusb_model_{unblind,blind}.pdf"
