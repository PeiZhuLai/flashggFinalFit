#!/bin/bash
# Update mA=1 unblind diagnostics (GOF / impacts / S+B fit + model plot) after the
# mA1 signal-model core-window change. Scoped to mA=1 only via ALPMASS.
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
cd "$FFIT/Plots"
export ALPMASS="1" NPROC=1 PARALLEL=1 NTOYS=500

echo "############ [1/4] GOF (mA1) ############"
bash 3_gof_unblind.sh 2>&1 | tail -8
echo "############ [2/4] Impacts (mA1) ############"
bash 4_impacts_observed.sh 2>&1 | tail -6
echo "############ [3/4] S+B FitDiagnostics (mA1) ############"
bash 5_splusb_fit.sh 2>&1 | tail -6
echo "############ [4/4] S+B model plot (mA1) ############"
bash 9_splusb_modelplot.sh 2>&1 | tail -6
echo "############ ALL DONE ############"
