#!/bin/bash

mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
# mAs=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# mAs=( 5 15 30 )
BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal"
FlashggBase="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
CmsswTop="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4"
Executable="$BaseDir/1_bias_study.sh"
LogDir="$BaseDir/logs"
SubmitFile="1_sub_bias_study.submit"

# Clean logs from the previous run so they don't accumulate and fill the work quota.
# Done here at submit time (not after the run): Condor writes these asynchronously
# while the jobs are still running, so they can only be safely wiped before the next launch.
echo ">>> Cleaning logs from previous run (logs/, condor_logs/, bias_jobs/)..."
rm -rf "$LogDir" "$BaseDir/condor_logs" "$BaseDir/bias_jobs"

# Create log directory if it doesn't exist
mkdir -p "$LogDir"

# Generate HTCondor submit file
cat > $SubmitFile << EOF
# subjob_sig_fTest
universe              = vanilla
executable            = $Executable
getenv                = True
environment           = "BASE_DIR=$FlashggBase CMSSW_TOP=$CmsswTop"
request_memory        = 2500
transfer_output_files = ""
+JobFlavour           = "tomorrow"
on_exit_hold          = (ExitBySignal == True) || (ExitCode != 0)
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

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
