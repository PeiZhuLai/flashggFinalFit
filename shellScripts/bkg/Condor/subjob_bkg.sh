#!/bin/bash

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -eo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_base_dir="$(cd "${script_dir}/../../.." && pwd)"

BASE_DIR="${BASE_DIR:-${default_base_dir}}"
CMSSW_TOP="${CMSSW_TOP:-${CMSSW_BASE:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}}"
lable="${LABLE:-run3}"
version="${VERSION:-ReReco}"
Lumi_run3="${LUMI_RUN3:-170.84}"
dir_input="${DIR_INPUT:-/eos/home-p/pelai/HZa/root_MVAcut/data}"

if [[ $# -gt 0 ]]; then
  massList=( "$@" )
else
  massList=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
fi

format_duration() {
  local total_seconds="${1:-0}"
  local hours=$(( total_seconds / 3600 ))
  local minutes=$(( (total_seconds % 3600) / 60 ))
  local seconds=$(( total_seconds % 60 ))
  printf "%02dh:%02dm:%02ds" "$hours" "$minutes" "$seconds"
}

script_start_time=$(date +%s)

echo ">>> Preparing background Condor submission"
echo "  BASE_DIR  = ${BASE_DIR}"
echo "  CMSSW_TOP = ${CMSSW_TOP}"
echo "  version   = ${version}"
echo "  lable     = ${lable}"
echo "  dir_input = ${dir_input}"
echo "  masses    = ${massList[*]}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_TOP}/src"
cmsenv

source "${BASE_DIR}/setup.sh"

cd "${BASE_DIR}/Background"
build_start_time=$(date +%s)
make -f makefile clean
make -f makefile -B
build_end_time=$(date +%s)
echo "[Timer] Background build finished in $(format_duration $((build_end_time - build_start_time)))"

if [[ ! -x ./bin/fTest_ALP_turnOn || ! -x ./bin/makeBkgPlots_ALP ]]; then
  echo "[ERROR] build failed or binaries missing in ${BASE_DIR}/Background/bin"
  ls -l ./bin || true
  exit 1
fi

dir_output="${BASE_DIR}/Background/ALP_BkgModel_${version}"
path_out_bkg="${dir_output}/fit_results_${lable}"
total_OutDir="${path_out_bkg}/AllFitResults"
log_dir="${path_out_bkg}/condor_logs"
submit_file="${BASE_DIR}/shellScripts/bkg/Condor/subjob_bkg.submit"
executable="${BASE_DIR}/shellScripts/bkg/Condor/runjob_bkg.sh"

mkdir -p "$path_out_bkg" "$total_OutDir" "$log_dir"
for mass in "${massList[@]}"; do
  mkdir -p "${path_out_bkg}/${mass}"
done

cat > "$submit_file" << EOF
# subjob_bkg
universe              = vanilla
executable            = ${executable}
getenv                = True
request_memory        = 2000
request_cpus          = 1
+JobFlavour           = "workday"
batch_name            = Fit_bkg

# Send the job to Held state on failure.
on_exit_hold          = (ExitBySignal == True) || (ExitCode != 0)

# Periodically retry failed jobs every 10 minutes, up to a maximum of 3 starts.
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

EOF

for mass in "${massList[@]}"; do
  cat >> "$submit_file" << EOF
log                   = ${log_dir}/bkg_mA${mass}.\$(ClusterId).log
output                = ${log_dir}/bkg_mA${mass}.\$(ClusterId).\$(ProcId).out
error                 = ${log_dir}/bkg_mA${mass}.\$(ClusterId).\$(ProcId).err
arguments             = ${mass} ${lable} ${version} ${Lumi_run3} ${dir_input} ${path_out_bkg} ${total_OutDir} ${BASE_DIR} ${CMSSW_TOP}
queue

EOF
done

echo ">>> Generated ${submit_file}"
condor_submit "$submit_file"

echo ">>> After jobs finish, merge per-mass summaries with:"
echo "    bash ${BASE_DIR}/shellScripts/bkg/Condor/collect_bkg_results.sh"
echo "[Timer] Total subjob_bkg.sh runtime: $(format_duration $(( $(date +%s) - script_start_time )))"
