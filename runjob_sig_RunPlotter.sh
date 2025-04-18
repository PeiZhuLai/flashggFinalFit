#!/bin/bash
/bin/hostname
gcc -v
pwd
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv

dir_sig="/afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal"
source /afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/tools

python3 $dir_sig/RunPlotter.py --mass_ALP $1 --years $2 --channel $3


