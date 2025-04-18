#!/bin/bash

# Configuration variables
ChannelList=( ele mu )
ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
YearsList=( 16 16APV 17 18 )
BaseDir="/afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard"
Executable="$BaseDir/runjob_gen_datacard_makeYields.sh"
InputWSDir="/afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output"
LogDir="$BaseDir/logfiles"
SubmitFile="subjob_gen_datacard_makeYields.submit"

# Create log directory if it doesn't exist
mkdir -p "$LogDir"

# Generate HTCondor submit file
cat > $SubmitFile << EOF
# subjob_gen_datacard_makeYields
universe              = vanilla
executable            = $Executable
getenv                = True
request_memory        = 2000
+JobFlavour           = "tomorrow"

EOF

# Generate job combinations and append to the submit file
for iChannel in "${!ChannelList[@]}"; do
  for iALPmass in "${!ALPmassList[@]}"; do
    for iYear in "${!YearsList[@]}"; do
      log_file="$LogDir/${ChannelList[iChannel]}/${ALPmassList[iALPmass]}_makeYields_job_${YearsList[iYear]}_${ChannelList[iChannel]}.log"
      output_file="$LogDir/${ChannelList[iChannel]}/${ALPmassList[iALPmass]}_makeYields_job_${YearsList[iYear]}_${ChannelList[iChannel]}.out"
      error_file="$LogDir/${ChannelList[iChannel]}/${ALPmassList[iALPmass]}_makeYields_job_${YearsList[iYear]}_${ChannelList[iChannel]}.err"
      arguments="${YearsList[iYear]}=$InputWSDir ${ALPmassList[iALPmass]} ${YearsList[iYear]} ${ChannelList[iChannel]}"

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
