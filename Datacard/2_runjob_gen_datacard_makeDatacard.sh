#!/bin/bash
# /bin/hostname
# gcc -v
pwd
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv

ulimit -s unlimited
set -e
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
export SCRAM_ARCH=el9_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/tools

path_makeDatacard="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard"

# python3 $path_makeDatacard/makeDatacard.py --mass_ALP $1 --years $2 --channel $3

# Test
ALPmassList=( 5 15 30 )

for iALPmass in "${!ALPmassList[@]}"; do
    python3 $path_makeDatacard/makeDatacard.py --mass_ALP ${ALPmassList[iALPmass]} --years 2022preEE --channel leptons
done