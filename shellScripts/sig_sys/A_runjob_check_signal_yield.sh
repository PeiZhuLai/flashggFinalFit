#!/bin/bash
pwd
source /cvmfs/cms.cern.ch/cmsset_default.sh
CMSSW_TOP="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4"
cd "${CMSSW_TOP}/src"
cmsenv

dir_sig="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal"
source /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/tools

# ------------- Open for Condor -------------
# python3 $dir_sig/RunPlotter.py --mass_ALP $1 --years $2 --channel $3

# ------------- Open for Dry Run -------------
# python3 $dir_sig/RunPlotter.py --mass_ALP 5 --years 2022preEE --channel ele

# ------------- Open for Local Run -------------
python3 $dir_sig/scripts/audit_signal_yields.py \
  --masses 20 \
  --years 2024 \
  --channels ele \
  --csv $dir_sig/scripts/audit_m20_2024_ele.csv
