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
bias_dir="${base_dir}/Combine/Checks/Bias_nominal"
root_datacard_path="${ROOT_DATACARD_PATH:-${base_dir}/Combine/root_t2w}"
fit_config="${bias_dir}/bias_fit_config.sh"

format_duration() {
  local total_seconds="${1:-0}"
  local hours=$(( total_seconds / 3600 ))
  local minutes=$(( (total_seconds % 3600) / 60 ))
  local seconds=$(( total_seconds % 60 ))
  printf "%02dh:%02dm:%02ds" "$hours" "$minutes" "$seconds"
}

copy_outputs() {
  local src_dir="$1"
  local dst_dir="$2"

  mkdir -p "$dst_dir"
  for subdir in BiasToys BiasFits BiasPlots BiasJson plots_bias; do
    if [[ -d "${src_dir}/${subdir}" ]]; then
      mkdir -p "${dst_dir}/${subdir}"
      cp -p "${src_dir}/${subdir}"/* "${dst_dir}/${subdir}/" 2>/dev/null || true
    fi
  done
}

require_glob() {
  local pattern="$1"
  local label="$2"

  if ! compgen -G "$pattern" >/dev/null; then
    echo "[ERROR] Missing ${label}: ${pattern}"
    exit 1
  fi
}

script_start_time=$(date +%s)

echo ">>> Running bias study Condor job"
echo "  mA        = ${mA}"
echo "  BASE_DIR  = ${base_dir}"
echo "  CMSSW_TOP = ${cmssw_top}"
echo "  BIAS_DIR  = ${bias_dir}"
echo "  CARD_DIR  = ${root_datacard_path}"

ulimit -s unlimited

if [[ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]]; then
  source /cvmfs/cms.cern.ch/cmsset_default.sh
else
  echo "[ERROR] Missing /cvmfs/cms.cern.ch/cmsset_default.sh"
  exit 1
fi

cd "${cmssw_top}/src" || exit 1
cmsenv

if [[ -r "${base_dir}/setup.sh" ]]; then
  source "${base_dir}/setup.sh"
fi

if [[ ! -r "${fit_config}" ]]; then
  echo "[ERROR] Missing ${fit_config}"
  exit 1
fi

source "${fit_config}"
bias_fit_opts="$(resolve_bias_fit_combine_options)"
echo "  fit opts  = ${bias_fit_opts}"

source_crab_env() {
  set --
  source /cvmfs/cms.cern.ch/crab3/crab.sh
}

if [[ -r /cvmfs/cms.cern.ch/crab3/crab.sh ]]; then
  source_crab_env
fi

card="${root_datacard_path}/${mA}_Datacard_leptons.root"
if [[ ! -f "$card" ]]; then
  echo "[ERROR] Missing datacard workspace: ${card}"
  exit 1
fi

mkdir -p "${bias_dir}/bias_jobs" "${bias_dir}/bias_outputs"
work_dir="${bias_dir}/bias_jobs/mA_${mA}_$(date +%Y%m%d_%H%M%S)_${RANDOM}"
output_dir="${bias_dir}/bias_outputs/mA_${mA}"
mkdir -p "$work_dir" "$output_dir"
cd "$work_dir" || exit 1

echo ">>> Work dir: ${work_dir}"
echo ">>> Output snapshot dir: ${output_dir}"

toys_start_time=$(date +%s)
python3 "${bias_dir}/RunBiasStudy.py" -d "$card" --mA "$mA" -t
require_glob "BiasToys/biasStudy_*_toys.root" "bias toys"
toys_end_time=$(date +%s)
echo "[Timer] mA=${mA} toys finished in $(format_duration $((toys_end_time - toys_start_time)))"

fits_start_time=$(date +%s)
python3 "${bias_dir}/RunBiasStudy.py" \
  -d "$card" \
  --mA "$mA" \
  -f \
  -c "${bias_fit_opts}"
require_glob "BiasFits/biasStudy_*_fits.root" "bias fits"
fits_end_time=$(date +%s)
echo "[Timer] mA=${mA} fits finished in $(format_duration $((fits_end_time - fits_start_time)))"

plots_start_time=$(date +%s)
python3 "${bias_dir}/RunBiasStudy.py" -d "$card" --mA "$mA" -p --gaussianFit
require_glob "BiasPlots/biasStudy_*_pulls.pdf" "per-pdf bias plots"
require_glob "BiasJson/${mA}_gaussfit.json" "bias gaussian-fit json"
BIAS_BASE_OVERRIDE="$work_dir" python3 "${bias_dir}/plot_bias.py" --mA "$mA"
require_glob "plots_bias/${mA}_pull_overlay_distribution.pdf" "bias overlay plot"
plots_end_time=$(date +%s)
echo "[Timer] mA=${mA} plots finished in $(format_duration $((plots_end_time - plots_start_time)))"

copy_outputs "$work_dir" "$output_dir"

echo "[OK] mA=${mA}"
echo "[Timer] Total runjob_bias_study.sh runtime: $(format_duration $(( $(date +%s) - script_start_time )))"
