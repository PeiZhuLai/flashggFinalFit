#!/bin/bash

# Configuration variables
ChannelList=( ele mu )
# ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
ALPmassList=( 5 15 30 )
YearsList=( 2022preEE )
BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
EosDir="/eos/home-p/pelai/HZa/root_MVAcut/sig"
Executable="$BaseDir/shellScripts/sig_sys/2_runjob_sig_calcPhotonSyst.sh"
LogDir="$BaseDir/Signal/log_files"
SubmitFile="2_subjob_sig_calcPhotonSyst.submit"

# Create log directory if it doesn't exist
mkdir -p "$LogDir"

# Generate HTCondor submit file
cat > $SubmitFile << EOF
# subjob_sig_signalFit
universe              = vanilla
executable            = $Executable
getenv                = True
request_memory        = 1000
+JobFlavour           = "workday"

EOF

# Generate job combinations and append to the submit file
for iChannel in "${!ChannelList[@]}"; do
  for iALPmass in "${!ALPmassList[@]}"; do
    for iYear in "${!YearsList[@]}"; do
      log_file="$LogDir/${ChannelList[iChannel]}/${ALPmassList[iALPmass]}_calcPhotonSys_job_${YearsList[iYear]}_${ChannelList[iChannel]}.log"
      output_file="$LogDir/${ChannelList[iChannel]}/${ALPmassList[iALPmass]}_calcPhotonSys_job_${YearsList[iYear]}_${ChannelList[iChannel]}.out"
      error_file="$LogDir/${ChannelList[iChannel]}/${ALPmassList[iALPmass]}_calcPhotonSys_job_${YearsList[iYear]}_${ChannelList[iChannel]}.err"
      arguments="${ALPmassList[iALPmass]} ${YearsList[iYear]} ${ChannelList[iChannel]} $EosDir/ALP_M${ALPmassList[iALPmass]}/ws_Tree2WS"

      echo "log = $log_file" >> $SubmitFile
      echo "output = $output_file" >> $SubmitFile
      echo "error = $error_file" >> $SubmitFile
      echo "arguments = $arguments" >> $SubmitFile
      echo "queue" >> $SubmitFile
      echo "" >> $SubmitFile
    done
  done
done

# Submit the HTCondor job
condor_submit $SubmitFile
