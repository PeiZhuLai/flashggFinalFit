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

# python3 $dir_sig/RunPlotter.py --mass_ALP 30 --channel mu --years '2022preEE,2022postEE,2023preBPix,2023postBPix,2024'

# ------------- Open for Local Run -------------
mAs=( 2 )
channels=( ele mu )
years=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024)

for mA in "${mAs[@]}"; do
    for channel in "${channels[@]}"; do
        for year in "${years[@]}"; do
            python3 $dir_sig/RunPlotter.py --mass_ALP ${mA} --years ${year} --channel ${channel}
        done
        python3 $dir_sig/RunPlotter.py --mass_ALP ${mA} --channel ${channel} --years '2022preEE,2022postEE,2023preBPix,2023postBPix,2024'
    done
done