#!/bin/bash

cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine

eval `scramv1 runtime -sh`

text2workspace.py /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_ele/20_pruned_datacard_16APV_ele.txt -o ./output_datacard_rootfile_ele/20_Datacard_16APV_ele_mu_inclusive.root -m 125 higgsMassRange=122,128 