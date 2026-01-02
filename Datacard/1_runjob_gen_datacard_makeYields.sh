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

mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
years=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024)

InputWSDir="/eos/home-p/pelai/HZa/root_MVAcut"
for mA in "${mAs[@]}"; do
    python3 $path_makeYields/makeYields.py \
        --inputWSDirMap 2022preEE=$InputWSDir,2022postEE=$InputWSDir,2023preBPix=$InputWSDir,2023postBPix=$InputWSDir,2024=$InputWSDir \
        --mass_ALP "$mA" \
        --channel leptons \
        --doSystematics
done
