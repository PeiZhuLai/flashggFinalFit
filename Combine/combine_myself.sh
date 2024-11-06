#!/bin/bash

cmsenv

ChannelList=( ele mu )
# ChannelList=( ele )
nChannel=${#ChannelList[@]}

ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# ALPmassList=( 1 )
nALPmass=${#ALPmassList[@]}

YearsList=( 16 16APV 17 18 )
# YearsList=( 16 )
nYear=${#YearsList[@]}
# combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m 125 --rAbsAcc 0.00000001 --rMin -200 --rMax 200 --freezeParameters MH


# Expected (--doObserved False)
path_datacard="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine"
combine -M AsymptoticLimits $path_datacard/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root -v 3 --cminDefaultMinimizerStrategy 0 --setParameterRanges MH=120,130:r=-100,100

# path_datacard="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_ele/1_pruned_datacard_16_ele.txt"
# combine -M AsymptoticLimits -v 3 $path_datacard --cminDefaultMinimizerStrategy 0 --setParameterRanges MH=120,130 -m 125
# combine -M AsymptoticLimits -v 3 $path_datacard --cminDefaultMinimizerStrategy 0 -t -1 --setParameterRanges MH=120,130 -m 125

# First step to debug 
# path_datacard="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_ele/1_pruned_datacard_16_ele.txt"
# combine -M MultiDimFit -v 3 $path_datacard

# path_datacard="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine"
# combine -M GenerateOnly $path_datacard/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root -t -1 --saveToys --setParameters r=1
# combineTool.py -M FastScan -w $path_datacard/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root:w -d higgsCombineTest.GenerateOnly.mH120.123456.root:toys/toy_asimov

# combineTool.py -M FastScan -w $path_datacard/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root:w

#other people 
# path_datacard=/publicfs/cms/user/houbaorui/WorkSpace/combine/CMSSW_14_1_0_pre4/src/Limit/datacard_RH2800LowmassSRElectron2022PreEERegionA.txt
# combine -M AsymptoticLimits $path_datacard

# official counting experiment OK
# path_datacard=realistic-counting-experiment.txt
# combine -M AsymptoticLimits $path_datacard

# path_datacard=./Shape_tutorial/simple-shapes-TH1.txt
# combine -M AsymptoticLimits $path_datacard

# Data (--doObserved True)



# for ((iChannel=0; iChannel<$nChannel; iChannel++))
#   do
#   for ((iALPmass=0; iALPmass<$nALPmass; iALPmass++))
#     do
#     for ((iYear=0; iYear<$nYear; iYear++))
#       do
#       python3 RunText2Workspace.py --mode mu_inclusive --mass_ALP ${ALPmassList[$iALPmass]} --year ${YearsList[$iYear]} --channel ${ChannelList[$iChannel]}
#       done
#     done
#   done

