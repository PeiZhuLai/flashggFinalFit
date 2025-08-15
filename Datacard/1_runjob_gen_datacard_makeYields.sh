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

path_makeYields="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard"

# python3 $path_makeYields/makeYields.py --inputWSDirMap $1 --mass_ALP $2 --year $3 --channel $4

ALPmassList=( 5 15 30 )

# Test
InputWSDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run3_ReReco/output"
for iALPmass in "${!ALPmassList[@]}"; do
    python3 $path_makeYields/makeYields.py --inputWSDirMap 2022preEE=$InputWSDir --mass_ALP ${ALPmassList[iALPmass]} --channel leptons
done