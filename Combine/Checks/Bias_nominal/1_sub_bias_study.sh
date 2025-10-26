#!/bin/bash

# mAs=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
mAs=( 5 15 30 )
BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal"
Executable="$BaseDir/1_run_bias_study.sh"
LogDir="$BaseDir/logs"
SubmitFile="1_sub_bias_study.submit"

# Create log directory if it doesn't exist
mkdir -p "$LogDir"

# Generate HTCondor submit file
cat > $SubmitFile << EOF
# subjob_sig_fTest
universe              = vanilla
executable            = $Executable
getenv                = True
request_memory        = 1000
+JobFlavour           = "workday"

EOF


# Generate job combinations and append to the submit file
for mA in "${mAs[@]}"; do

    log_file="$LogDir/${mA}_BiasStudy.log"
    output_file="$LogDir/${mA}_BiasStudy.out"
    error_file="$LogDir/${mA}_BiasStudy.err"
    arguments="${mA}"

    echo "log = $log_file" >> $SubmitFile
    echo "output = $output_file" >> $SubmitFile
    echo "error = $error_file" >> $SubmitFile
    echo "arguments = $arguments" >> $SubmitFile
    echo "queue" >> $SubmitFile
    echo "" >> $SubmitFile

done

# Submit the HTCondor job
condor_submit $SubmitFile
