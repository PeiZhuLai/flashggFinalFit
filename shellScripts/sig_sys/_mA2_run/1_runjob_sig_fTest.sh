#!/bin/bash
pwd
source /cvmfs/cms.cern.ch/cmsset_default.sh
CMSSW_TOP="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4"
cd "${CMSSW_TOP}/src"
cmsenv

dir_sig="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/scripts"
EosDir="/eos/home-p/pelai/HZa/root_MVAcut/sig"
source /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/tools

# ------------- Open for Condor -------------
# python3 $dir_sig/fTest.py --mass_ALP $1 --year $2 --mass $3 --channel $4 --inputWSDir $5

# ------------- Open for Dry Run -------------
BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
# python3 $dir_sig/fTest.py --mass_ALP 5 --year 2022preEE --mass 125 --channel ele --inputWSDir $EosDir/ALP_M5/ws_Tree2WS

# ------------- Open for Local Run -------------
mAs=( 2 )
channels=( ele mu )
years=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024)

for mA in "${mAs[@]}"; do
    for channel in "${channels[@]}"; do
        for year in "${years[@]}"; do
            python3 $dir_sig/fTest.py --mass_ALP ${mA} --year ${year} --mass 125 --channel ${channel} --inputWSDir $EosDir/mA_M${mA}/ws_Tree2WS
        done
    done
done
