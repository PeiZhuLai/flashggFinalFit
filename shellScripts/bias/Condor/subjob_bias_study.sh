#!/bin/bash

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -eo pipefail

BASE_DIR="${BASE_DIR:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit}"
CMSSW_TOP="${CMSSW_TOP:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}"

if [[ $# -gt 0 ]]; then
  mAs=( "$@" )
else
  mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
fi

format_duration() {
  local total_seconds="${1:-0}"
  local hours=$(( total_seconds / 3600 ))
  local minutes=$(( (total_seconds % 3600) / 60 ))
  local seconds=$(( total_seconds % 60 ))
  printf "%02dh:%02dm:%02ds" "$hours" "$minutes" "$seconds"
}

script_start_time=$(date +%s)

echo ">>> Preparing bias study Condor submission"
echo "  BASE_DIR  = ${BASE_DIR}"
echo "  CMSSW_TOP = ${CMSSW_TOP}"
echo "  masses    = ${mAs[*]}"

submit_dir="${BASE_DIR}/shellScripts/bias/Condor"
bias_dir="${BASE_DIR}/Combine/Checks/Bias_nominal"
log_dir="${bias_dir}/condor_logs"
submit_file="${submit_dir}/subjob_bias_study.submit"
executable="${submit_dir}/runjob_bias_study.sh"

mkdir -p "$log_dir" "${bias_dir}/bias_jobs" "${bias_dir}/bias_outputs"

cat > "$submit_file" << EOF
# subjob_bias_study
universe              = vanilla
executable            = ${executable}
getenv                = True
request_memory        = 4000
request_cpus          = 1
+JobFlavour           = "tomorrow"

# Send the job to Held state on failure.
on_exit_hold          = (ExitBySignal == True) || (ExitCode != 0)

# Periodically retry failed jobs every 10 minutes, up to a maximum of 3 starts.
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

EOF

for mA in "${mAs[@]}"; do
  cat >> "$submit_file" << EOF
log                   = ${log_dir}/bias_mA${mA}.\$(ClusterId).log
output                = ${log_dir}/bias_mA${mA}.\$(ClusterId).\$(ProcId).out
error                 = ${log_dir}/bias_mA${mA}.\$(ClusterId).\$(ProcId).err
arguments             = ${mA} ${BASE_DIR} ${CMSSW_TOP}
queue

EOF
done

echo ">>> Generated ${submit_file}"
condor_submit "$submit_file"
echo "[Timer] Total subjob_bias_study.sh runtime: $(format_duration $(( $(date +%s) - script_start_time )))"
