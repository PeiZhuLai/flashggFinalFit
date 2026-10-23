#!/bin/bash
#
# Submit the high-statistics bias study on the rebuilt envelopes: 10 chunks x
# 1000 toys per mass point, same job script and seeds as the production run
# (seed = 100000 + 100*mA + chunk, split seeds distinct), so truth functions
# shared with the production envelope are fitted on the same toys and the
# effect of the restored member is a paired comparison.
#
# Outputs : outputs/mA_<m>/chunk<i>/BiasFits   (BIAS_HS_ROOT)
# Cards   : t2w/<m>_Datacard_leptons.root       (ROOT_DATACARD_PATH)
# Logs    : logs/condor/
#
# Usage: bash 3_submit.sh [--run]   (default: write the submit file only)

set -eo pipefail
EC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CMSSW_TOP=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4
exe="${FF}/shellScripts/bias/Condor/runjob_bias_highstat.sh"
masses=( ${MASSES:-2 3 4 20 24 25} )
NCHUNKS="${NCHUNKS:-10}"
log_dir="${EC}/logs/condor"
sub="${EC}/envcheck_bias.submit"
mkdir -p "$log_dir" "${EC}/outputs"

cat > "$sub" << EOF
universe              = vanilla
executable            = ${exe}
getenv                = True
request_memory        = 4000
request_cpus          = 1
+JobFlavour           = "tomorrow"
batch_name            = Bias_EnvCheck
environment           = "NTOYS_PER_CHUNK=1000 BIAS_HS_ROOT=${EC}/outputs ROOT_DATACARD_PATH=${EC}/t2w"

on_exit_hold          = (ExitBySignal == True) || (ExitCode != 0)
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

EOF
for m in "${masses[@]}"; do
  [[ -f "${EC}/t2w/${m}_Datacard_leptons.root" ]] || { echo "[ERROR] missing t2w for mA=${m}"; exit 1; }
  for (( c=0; c<NCHUNKS; c++ )); do
    cat >> "$sub" << EOF
log                   = ${log_dir}/ec_mA${m}_c${c}.\$(ClusterId).log
output                = ${log_dir}/ec_mA${m}_c${c}.\$(ClusterId).out
error                 = ${log_dir}/ec_mA${m}_c${c}.\$(ClusterId).err
arguments             = ${m} ${c} ${FF} ${CMSSW_TOP}
queue

EOF
  done
done
echo "jobs: $(( ${#masses[@]} * NCHUNKS ))  submit file: ${sub}"
[[ "${1:-}" == "--run" ]] && condor_submit "$sub" || echo "(not submitted; pass --run)"
