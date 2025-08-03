#!/bin/bash
# /bin/hostname
# gcc -v
pwd
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv

dir_sig="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/scripts"
source /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/tools

python3 $dir_sig/fTest.py --mass_ALP $1 --year $2 --mass $3 --channel $4 --inputWSDir $5

# Comissioning
# BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
# python3 $dir_sig/fTest.py --mass_ALP 5 --year 2022preEE --mass 125 --channel ele --inputWSDir $BaseDir/MVAcut/run3_ReReco/output/sig/ele