#!/bin/bash

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -eo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 MASS [BASE_DIR] [CMSSW_TOP]"
  exit 2
fi

mA="$1"
base_dir="${2:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit}"
cmssw_top="${3:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}"
combine_dir="${base_dir}/Combine"
mh="125.38"

format_duration() {
  local total_seconds="${1:-0}"
  local hours=$(( total_seconds / 3600 ))
  local minutes=$(( (total_seconds % 3600) / 60 ))
  local seconds=$(( total_seconds % 60 ))
  printf "%02dh:%02dm:%02ds" "$hours" "$minutes" "$seconds"
}

script_start_time=$(date +%s)

echo ">>> Running expected impacts Condor job"
echo "  mA        = ${mA}"
echo "  BASE_DIR  = ${base_dir}"
echo "  CMSSW_TOP = ${cmssw_top}"

ulimit -s unlimited

if [[ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]]; then
  source /cvmfs/cms.cern.ch/cmsset_default.sh
else
  echo "[ERROR] Missing /cvmfs/cms.cern.ch/cmsset_default.sh"
  exit 1
fi

cd "${cmssw_top}/src" || exit 1
cmsenv

source_crab_env() {
  set --
  source /cvmfs/cms.cern.ch/crab3/crab.sh
}

if [[ -r /cvmfs/cms.cern.ch/crab3/crab.sh ]]; then
  source_crab_env
fi

cd "$combine_dir" || exit 1
mkdir -p root_t2w output_impacts impact_jobs

echo ">>> Processing Text2Workspace for mA = ${mA}"
python3 RunText2Workspace.py --batch local --queue workday --mA "$mA"

card="${combine_dir}/root_t2w/${mA}_Datacard_leptons.root"
if [[ ! -f "$card" ]]; then
  echo "[ERROR] Missing workspace: ${card}"
  exit 1
fi

work_dir="${combine_dir}/impact_jobs/mA_${mA}_$(date +%Y%m%d_%H%M%S)_${RANDOM}"
mkdir -p "$work_dir"
cd "$work_dir" || exit 1

echo ">>> Running expected impacts for mA = ${mA}"
echo "  card = ${card}"

common_opts=(
  -M Impacts
  -d "$card"
  -m "$mh"
  --redefineSignalPOIs r
  --rMin -5
  --rMax 5
  -t -1
  --expectSignal 1
  --robustFit 1
  --cminPreScan
  --freezeParameters=MH
  --setParameterRanges 'Exp_turnon_p1=95,125:Exp_width_p1=0.1,50:Exp_sigma_p1=0.05,20:Pow_turnon_p1=95,125:Pow_width_p1=0.1,50:Pow_sigma_p1=0.05,20:Lau_turnon_p1=95,125:Lau_width_p1=0.1,50:Lau_sigma_p1=0.05,20:Bern_gsigma=0.05,20:Bern_step=90,130:Bern_stepWidth=0.1,50:lumi_13p6TeV_Uncorrelated_2022preEE=-5,5'
  --cminDefaultMinimizerType Minuit2
  --cminDefaultMinimizerStrategy 0
  --cminDefaultMinimizerTolerance 0.01
  --cminFallbackAlgo Minuit2,0:0.1
  -v 1
)

initial_start_time=$(date +%s)
combineTool.py "${common_opts[@]}" --doInitialFit
initial_end_time=$(date +%s)
echo "[Timer] mA=${mA} initial fit finished in $(format_duration $((initial_end_time - initial_start_time)))"

fits_start_time=$(date +%s)
combineTool.py "${common_opts[@]}" --doFits --parallel 4
fits_end_time=$(date +%s)
echo "[Timer] mA=${mA} nuisance fits finished in $(format_duration $((fits_end_time - fits_start_time)))"

json_file="${combine_dir}/output_impacts/${mA}_impacts.json"
plot_file="${combine_dir}/output_impacts/${mA}_Expected_Impacts"

json_start_time=$(date +%s)
combineTool.py -M Impacts -d "$card" -m "$mh" --redefineSignalPOIs r -o "$json_file"
json_end_time=$(date +%s)
echo "[Timer] mA=${mA} impacts json finished in $(format_duration $((json_end_time - json_start_time)))"

plot_start_time=$(date +%s)
plotImpacts.py -i "$json_file" -o "$plot_file"
plot_end_time=$(date +%s)
echo "[Timer] mA=${mA} impact plot finished in $(format_duration $((plot_end_time - plot_start_time)))"

echo "[OK] mA=${mA}"
echo "[Timer] Total runjob_expectedImpact.sh runtime: $(format_duration $(( $(date +%s) - script_start_time )))"
