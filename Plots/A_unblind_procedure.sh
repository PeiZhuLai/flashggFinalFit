#!/bin/bash
# HZa low-mass ALP (mA = 1..30) -- UNBLINDING procedure driver (PARALLEL).
# Each step self-contains cmsenv + paths and operates on the OBSERVED signal-region data.
# mA1 uses the R=1 working-point datacard. Step 1 builds Combine/root_t2w + the observed
# limits that the later steps depend on, so it MUST run (and finish) first.
#
# PARALLELISM: every per-mA step now parallelizes its mA loop with NPROC workers
# (default 4; impacts uses NPROC mA-workers x PARALLEL nuisance-fits). Steps run
# sequentially (one heavy combine step at a time) to stay friendly on shared lxplus --
# the big win is the ~NPROC speedup *within* each step. Tune with NPROC=.. / PARALLEL=..
# To go further on a dedicated node, the independent steps (3,4,5,6,8) can be launched
# concurrently (see the commented block at the bottom).
#
#   1 1_makeLimits_observed.sh  observed 95% CL upper limits (+ BUILDS root_t2w)   [must be first]
#   2 2_runLimitsPlot.sh        limit-vs-mA plots (expected + observed overlay)
#   3 3_gof_unblind.sh          goodness-of-fit on observed data (+ plot)
#   4 4_impacts_observed.sh     observed impact plots + best-fit r-hat (r = observed)
#   5 5_splusb_fit.sh           S+B (and B-only) FitDiagnostics: r-hat + post-fit shapes
#   6 6_significance_local.sh   observed local significance (discovery q0) per mA   [-> 7]
#   7 7_significance_global.sh  global significance (look-elsewhere, Gross-Vitells)
#   8 8_rscan.sh                1D profile-likelihood scan of r (observed)
#   9 9_splusb_modelplot.sh     S+B post-fit model plot (blinded + unblinded)       [needs 5]

baseDir=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Plots
export NPROC="${NPROC:-4}"
export PARALLEL="${PARALLEL:-2}"

cd $baseDir
bash 1_makeLimits_observed.sh    # gate: builds root_t2w + observed limits
bash 2_runLimitsPlot.sh
bash 6_significance_local.sh
bash 7_significance_global.sh    # needs 6
bash 3_gof_unblind.sh
bash 5_splusb_fit.sh
bash 9_splusb_modelplot.sh       # needs 5
bash 8_rscan.sh
bash 4_impacts_observed.sh       # long pole (2 mA-workers x PARALLEL nuisance-fits)

# ---- Optional: concurrent independent steps (dedicated node only; oversubscribes shared lxplus) ----
# bash 1_makeLimits_observed.sh
# bash 2_runLimitsPlot.sh
# NPROC=2 bash 3_gof_unblind.sh        &
# NPROC=2 bash 5_splusb_fit.sh         &
# NPROC=2 bash 8_rscan.sh              &
# NPROC=2 bash 6_significance_local.sh &
# NPROC=2 bash 4_impacts_observed.sh   &
# wait
# bash 7_significance_global.sh
# bash 9_splusb_modelplot.sh
