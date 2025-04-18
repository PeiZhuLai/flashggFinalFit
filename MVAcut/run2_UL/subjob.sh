#!/bin/bash

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# massList=( 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 )

nMass=${#massList[@]}

path_code='/afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL'

for ((iBin=0; iBin<$nMass; iBin++))
do
    echo "Submitting job for mass ${massList[$iBin]}"
    nohup $path_code/runjob.sh ${massList[$iBin]} > ./log_file/job${massList[$iBin]}.log 2>&1 &
    # nohup runjob.sh -o ./log_file/job${massList[$iBin]}.log -e ./log_file/job${massList[$iBin]}.err -argu ${massList[$iBin]}
done