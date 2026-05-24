#!/bin/bash

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -eo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_dir="${BASE_DIR:-$(cd "${script_dir}/../../.." && pwd)}"
cmssw_top="${CMSSW_TOP:-${CMSSW_BASE:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}}"
bias_dir="${base_dir}/Combine/Checks/Bias_nominal"
root_datacard_path="${ROOT_DATACARD_PATH:-${base_dir}/Combine/root_t2w}"
fit_config="${bias_dir}/bias_fit_config.sh"

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

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 MASS"
  exit 2
fi

if [[ ! -r "${fit_config}" ]]; then
  echo "[ERROR] Missing ${fit_config}"
  exit 1
fi

MA="$1"
card="${root_datacard_path}/${MA}_Datacard_leptons.root"

if [[ ! -f "${card}" ]]; then
  echo "[ERROR] Missing datacard workspace: ${card}"
  exit 1
fi

source "${fit_config}"
bias_fit_opts="$(resolve_bias_fit_combine_options)"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${cmssw_top}/src"
cmsenv

source "${base_dir}/setup.sh"
export PYTHONPATH="${PYTHONPATH}:${base_dir}/tools:${base_dir}/Signal/tools"

mkdir -p "${bias_dir}/bias_jobs" "${bias_dir}/bias_outputs"
work_dir="${bias_dir}/bias_jobs/mA_${MA}_local_$(date +%Y%m%d_%H%M%S)_${RANDOM}"
output_dir="${bias_dir}/bias_outputs/mA_${MA}"
mkdir -p "${work_dir}" "${output_dir}"
cd "${work_dir}"

echo ">>> Running bias study for mA = ${MA}"
echo "  BASE_DIR           = ${base_dir}"
echo "  CMSSW_TOP          = ${cmssw_top}"
echo "  ROOT_DATACARD_PATH = ${root_datacard_path}"
echo "  fit opts           = ${bias_fit_opts}"
echo "  work dir           = ${work_dir}"
echo "  output dir         = ${output_dir}"

python3 "${bias_dir}/RunBiasStudy.py" -d "${card}" --mA "${MA}" -t
python3 "${bias_dir}/RunBiasStudy.py" -d "${card}" --mA "${MA}" -f -c "${bias_fit_opts}"
python3 "${bias_dir}/RunBiasStudy.py" -d "${card}" --mA "${MA}" -p --gaussianFit
BIAS_BASE_OVERRIDE="${work_dir}" python3 "${bias_dir}/plot_bias.py" --mA "${MA}"

copy_outputs "${work_dir}" "${output_dir}"
