#!/bin/bash
/bin/hostname
gcc -v
pwd
# export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/
source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsenv

/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/bin/fTest_ALP_turnOn -i $1 --saveMultiPdf $2 -D $3 --mass_ALP $4 -c $5 --isFlashgg $6 --isData $7 -f $8

