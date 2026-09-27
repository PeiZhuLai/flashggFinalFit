#!/bin/bash
#
# Submit the high-statistics bias study: NCHUNKS parallel chunks per mass point.
#
# Usage: subjob_bias_highstat.sh MASS [MASS ...]
#   NCHUNKS=10            number of parallel chunks per mass point
#   NTOYS_PER_CHUNK=1000  toys per chunk (total = NCHUNKS x this)

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -eo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_base_dir="$(cd "${script_dir}/../../.." && pwd)"
BASE_DIR="${BASE_DIR:-${default_base_dir}}"
CMSSW_TOP="${CMSSW_TOP:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}"
NCHUNKS="${NCHUNKS:-10}"
NTOYS_PER_CHUNK="${NTOYS_PER_CHUNK:-1000}"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 MASS [MASS ...]"
  exit 2
fi
mAs=( "$@" )

bias_dir="${BASE_DIR}/Combine/Checks/Bias_nominal"
log_dir="${bias_dir}/condor_logs_highstat"
submit_file="${script_dir}/subjob_bias_highstat.submit"
executable="${script_dir}/runjob_bias_highstat.sh"

mkdir -p "$log_dir" "${bias_dir}/bias_jobs_highstat" "${bias_dir}/bias_outputs_highstat"

total=$(( ${#mAs[@]} * NCHUNKS ))
echo ">>> High-statistics bias submission"
echo "  masses          = ${mAs[*]}"
echo "  chunks per mass = ${NCHUNKS}  (${NTOYS_PER_CHUNK} toys each)"
echo "  total toys/mass = $(( NCHUNKS * NTOYS_PER_CHUNK ))"
echo "  total jobs      = ${total}"
echo "  output          = ${bias_dir}/bias_outputs_highstat/mA_<m>/chunk<i>/BiasFits"
echo "  logs            = ${log_dir}"

cat > "$submit_file" << EOF
# subjob_bias_highstat
universe              = vanilla
executable            = ${executable}
getenv                = True
request_memory        = 4000
request_cpus          = 1
+JobFlavour           = "tomorrow"
batch_name            = Bias_HighStat
environment           = "NTOYS_PER_CHUNK=${NTOYS_PER_CHUNK}"

on_exit_hold          = (ExitBySignal == True) || (ExitCode != 0)
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

EOF

for mA in "${mAs[@]}"; do
  for (( c=0; c<NCHUNKS; c++ )); do
    cat >> "$submit_file" << EOF
log                   = ${log_dir}/hs_mA${mA}_c${c}.\$(ClusterId).log
output                = /dev/null
error                 = ${log_dir}/hs_mA${mA}_c${c}.\$(ClusterId).\$(ProcId).err
arguments             = ${mA} ${c} ${BASE_DIR} ${CMSSW_TOP}
queue

EOF
  done
done

echo ">>> Generated ${submit_file}"
condor_submit "$submit_file"
