#!/bin/bash

cmsenv

# ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
mAs=( 5 15 30)
# mAs=( 5 )

mkdir -p output_impacts

# combineTool.py -M Impacts -d root_t2w/5_Datacard_leptons.root -n _mA5 -m 125.38 --doInitialFit --robustFit 1 -t -1 --expectSignal 1
# combineTool.py -M Impacts -d root_t2w/5_Datacard_leptons.root -m 125.38 --doFits --robustFit 1  -t -1 --expectSignal 1
# combineTool.py -M Impacts -d root_t2w/5_Datacard_leptons.root -m 125.38 -o output_impacts/5_impacts.json
# plotImpacts.py -i output_impacts/5_impacts.json -o output_impacts/5_Expected_Impacts


for mA in "${mAs[@]}"; do
    # ---------- Expected Impacts ---------- 
    # 1. 生成workspace文件
    # python3 RunText2Workspace.py --batch local --queue workday --mA ${mA}
    
    CARD=root_t2w/${mA}_Datacard_leptons.root
    MH=125.38

    if [ $mA -lt 5 ]; then
        PDFIDX=2
    elif [ $mA -lt 15 ]; then
        PDFIDX=3
    elif [ $mA -lt 30 ]; then
        PDFIDX=2
    fi

    # 2. 运行initial fit
    combineTool.py -M Impacts -d $CARD -m $MH \
    --doInitialFit \
    --redefineSignalPOIs r --rMin -5 --rMax 5 \
    -t -1 --expectSignal 1 \
    --robustFit 1 --cminPreScan \
    --setParameters pdfindex_cat0_13p6TeV=$PDFIDX \
    --freezeParameters=pdfindex_cat0_13p6TeV,MH \
    --setParameterRanges 'Exp_turnon_p1=95,125:Exp_width_p1=0.1,50:Exp_sigma_p1=0.05,20:Pow_turnon_p1=95,125:Pow_width_p1=0.1,50:Pow_sigma_p1=0.05,20:Lau_turnon_p1=95,125:Lau_width_p1=0.1,50:Lau_sigma_p1=0.05,20:Bern_gsigma=0.05,20:Bern_step=90,130:Bern_stepWidth=0.1,50:lumi_13p6TeV_Uncorrelated_2022preEE=-5,5' \
    --cminDefaultMinimizerType Minuit2 \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    --cminFallbackAlgo Minuit2,0:0.1 \
    -v 1
    
    # 3. 运行所有systematic的fit
    combineTool.py -M Impacts -d $CARD -m $MH \
    --doFits --parallel 4 \
    --redefineSignalPOIs r --rMin -5 --rMax 5 \
    -t -1 --expectSignal 1 \
    --robustFit 1 --cminPreScan \
    --setParameters pdfindex_cat0_13p6TeV=$PDFIDX \
    --freezeParameters=pdfindex_cat0_13p6TeV,MH \
    --setParameterRanges 'Exp_turnon_p1=95,125:Exp_width_p1=0.1,50:Exp_sigma_p1=0.05,20:Pow_turnon_p1=95,125:Pow_width_p1=0.1,50:Pow_sigma_p1=0.05,20:Lau_turnon_p1=95,125:Lau_width_p1=0.1,50:Lau_sigma_p1=0.05,20:Bern_gsigma=0.05,20:Bern_step=90,130:Bern_stepWidth=0.1,50:lumi_13p6TeV_Uncorrelated_2022preEE=-5,5' \
    --cminDefaultMinimizerType Minuit2 \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    --cminFallbackAlgo Minuit2,0:0.1 \
    -v 1

    # 4. 生成json文件
    combineTool.py -M Impacts -d $CARD -m 125.38 \
    --redefineSignalPOIs r \
    -o output_impacts/${mA}_impacts.json

    # 5. 画图
    plotImpacts.py -i output_impacts/${mA}_impacts.json -o output_impacts/${mA}_Expected_Impacts --impact-xrange -0.13 0.13

    # ---------- Expected Impacts ----------
    # combineTool.py -M Impacts -d root_t2w/${mA}_Datacard_leptons.root -n mA${mA} -m 125.38 --doInitialFit --robustFit 1 -t -1 --expectSignal 1
    # combineTool.py -M Impacts -d root_t2w/${mA}_Datacard_leptons.root -m 125.38 --doFits --robustFit 1  -t -1 --expectSignal 1
    # combineTool.py -M Impacts -d root_t2w/${mA}_Datacard_leptons.root -m 125.38 -o output_impacts/${mA}_impacts.json
    # plotImpacts.py -i output_impacts/${mA}_impacts.json -o output_impacts/${mA}_Expected_Impacts

done