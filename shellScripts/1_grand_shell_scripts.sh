
baseDir=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit

### ----------- MVA Cut 
# python3 $baseDir/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py
# python3 $baseDir/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py

### ----------- Tree2WS
# cd $baseDir/Trees2WS
# sh $baseDir/Trees2WS/run_tree2ws.sh

### ----------- Background 
# cd $baseDir/shellScripts
# sh $baseDir/shellScripts/bkg/fit_bkg.sh
### ----------- Background Condor
# bash $baseDir/shellScripts/bkg/Condor/subjob_bkg.sh
# bash $baseDir/shellScripts/bkg/Condor/collect_bkg_results.sh
### ----------- Background Condor run locally
# bash bkg/Condor/subjob_bkg.sh
# bash bkg/Condor/collect_bkg_results.sh

### ----------- Signal
# cd $baseDir/shellScripts
# sh $baseDir/shellScripts/sig_sys/1_runjob_sig_fTest.sh
# sh $baseDir/shellScripts/sig_sys/2_runjob_sig_calcPhotonSyst.sh
# sh $baseDir/shellScripts/sig_sys/3_runjob_sig_signalFit.sh
# sh $baseDir/shellScripts/sig_sys/4_runjob_sig_RunPlotter.sh
# sh $baseDir/shellScripts/sig_sys/5_runjob_sig_plotEffSigma.sh

### ----------- Datacard
# cd $baseDir/Datacard
# sh 1_runjob_gen_datacard_makeYields.sh
# sh 2_runjob_gen_datacard_makeDatacard.sh
# sh 3_rysn_datacard.sh

# ### ----------- Combine Limits
# cd $baseDir/Combine
# sh 1_makeLimits.sh

# ### ----------- Plot Limts
# cd $baseDir/Plots
# sh 1_runLimitsPlot.sh

### ----------- Impact Plot
# cd $baseDir/Combine
# sh 2_expectedImpact.sh
### ----------- Impact Plot Condor (Generate ws for bias study)
# bash $baseDir/shellScripts/impact/Condor/subjob_expectedImpact.sh

### ----------- Bias Study
# cd $baseDir/Combine/Checks/Bias_nominal
# sh 1_bias_study.sh
### ----------- Bias Study Condor
bash $baseDir/shellScripts/bias/Condor/subjob_bias_study.sh

cd $baseDir/shellScripts