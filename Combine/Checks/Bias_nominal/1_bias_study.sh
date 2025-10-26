#!/bin/bash

cmsenv

# ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
mAs=( 5 15 30 )
# mAs=( 5 )

root_datacard_path="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/root_t2w"

for mA in "${mAs[@]}"; do

    rm -rf BiasFits BiasToys 2>/dev/null
    
    ./RunBiasStudy.py -d $root_datacard_path/${mA}_Datacard_leptons.root --mA ${mA} -t 
    ./RunBiasStudy.py -d $root_datacard_path/${mA}_Datacard_leptons.root --mA ${mA} -f -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --freezeParameters MH" 
    ./RunBiasStudy.py -d $root_datacard_path/${mA}_Datacard_leptons.root --mA ${mA} -p --gaussianFit

    python3 plot_bias.py --mA ${mA}
done