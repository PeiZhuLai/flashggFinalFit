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
mAs=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
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

# ------------- Sub-GeV (merged low-mA) points: 2024 only -------------
# 0pX signalFit output exists for 2024 only; other eras/combined would raise OSError.
# Requires RunPlotter.py --mass_ALP as type='string' (0p5 is not an int).
mAs_lowMA=( 0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9 )
for mA in "${mAs_lowMA[@]}"; do
    for channel in "${channels[@]}"; do
        python3 $dir_sig/RunPlotter.py --mass_ALP ${mA} --years 2024 --channel ${channel}
    done
done