#!/bin/bash
#
# One chunk of a high-statistics bias study.
#
# The nominal bias study (runjob_bias_study.sh) runs 1000 toys per mass point in
# a single job, which takes ~3 h. At m_a = 8 and 11 the measured bias is around
# the 0.2 threshold and repeat runs of the *same* configuration moved it by
# ~0.13 -- far more than the ~0.035 uncertainty quoted from a single run -- so
# the number is dominated by toy statistics rather than by a stable bias.
#
# This script runs one chunk of toys with its own seed and writes only the
# per-split fit files. Several chunks are submitted in parallel and merged
# afterwards by merge_bias_highstat.sh, which hadds the splits and runs the
# Gaussian fit once over the full sample. Ten chunks give 10000 toys with the
# wall time of a single 1000-toy job.
#
# Plot generation is skipped here on purpose: the pulls are only meaningful
# once all chunks are merged.
#
# Usage: runjob_bias_highstat.sh MASS CHUNK [BASE_DIR] [CMSSW_TOP]

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -eo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 MASS CHUNK [BASE_DIR] [CMSSW_TOP]"
  exit 2
fi

mA="$1"
chunk="$2"
base_dir="${3:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit}"
cmssw_top="${4:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}"

bias_dir="${base_dir}/Combine/Checks/Bias_nominal"
root_datacard_path="${ROOT_DATACARD_PATH:-${base_dir}/Combine/root_t2w}"
fit_config="${bias_dir}/bias_fit_config.sh"

# Toys per chunk. Keep at or below the nominal 1000 so one chunk costs about the
# same as one nominal job; the total is NCHUNKS x NTOYS_PER_CHUNK.
ntoys="${NTOYS_PER_CHUNK:-1000}"
split="${TOY_SPLIT:-500}"
# Distinct, reproducible seed per (mass, chunk).
seed=$(( 100000 + mA * 100 + chunk ))

format_duration() {
  local total_seconds="${1:-0}"
  printf "%02dh:%02dm:%02ds" $(( total_seconds / 3600 )) $(( (total_seconds % 3600) / 60 )) $(( total_seconds % 60 ))
}

script_start_time=$(date +%s)

echo ">>> High-statistics bias chunk"
echo "  mA        = ${mA}"
echo "  chunk     = ${chunk}"
echo "  ntoys     = ${ntoys} (split ${split})"
echo "  seed      = ${seed}"
echo "  BASE_DIR  = ${base_dir}"

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

if [[ -r /cvmfs/cms.cern.ch/crab3/crab.sh ]]; then
  ( set --; source /cvmfs/cms.cern.ch/crab3/crab.sh )
fi

card="${root_datacard_path}/${mA}_Datacard_leptons.root"
if [[ ! -f "$card" ]]; then
  echo "[ERROR] Missing datacard workspace: ${card}"
  exit 1
fi

out_root="${BIAS_HS_ROOT:-${bias_dir}/bias_outputs_highstat}/mA_${mA}/chunk${chunk}"
# Toys and combine logs (~8 MB per job) stay on the worker's local scratch; only
# the split fit files are copied back to AFS. Left in AFS, 300 jobs per rerun
# accumulated to 12 GB and filled the work quota. Local runs fall back to AFS.
work_parent="${_CONDOR_SCRATCH_DIR:-${bias_dir}/bias_jobs_highstat}"
work_dir="${work_parent}/mA_${mA}_chunk${chunk}_$(date +%Y%m%d_%H%M%S)_${RANDOM}"
mkdir -p "$work_dir" "$out_root"
cd "$work_dir" || exit 1

echo ">>> Work dir : ${work_dir}"
echo ">>> Output   : ${out_root}"

t0=$(date +%s)
python3 "${bias_dir}/RunBiasStudy.py" -d "$card" --mA "$mA" -t -n "$ntoys" --split "$split" -s "$seed"
if ! compgen -G "BiasToys/biasStudy_*_toys.root" >/dev/null; then
  echo "[ERROR] no toys produced"
  exit 1
fi
echo "[Timer] toys $(format_duration $(( $(date +%s) - t0 )))"

t0=$(date +%s)
python3 "${bias_dir}/RunBiasStudy.py" -d "$card" --mA "$mA" -f -n "$ntoys" --split "$split" -s "$seed" -c "${bias_fit_opts}"
if ! compgen -G "BiasFits/biasStudy_*_fits.root" >/dev/null; then
  echo "[ERROR] no fits produced"
  exit 1
fi
echo "[Timer] fits $(format_duration $(( $(date +%s) - t0 )))"

# Keep only the per-split files: the merge step re-hadds across all chunks, and
# shipping the chunk-local hadd would double-count those toys.
mkdir -p "${out_root}/BiasFits"
n=0
for f in BiasFits/biasStudy_*_split*_fits.root; do
  [[ -f "$f" ]] || continue
  cp -p "$f" "${out_root}/BiasFits/$(basename "$f")"
  n=$(( n + 1 ))
done
echo "[OK] mA=${mA} chunk=${chunk}: copied ${n} split fit files"
echo "[Timer] total $(format_duration $(( $(date +%s) - script_start_time )))"
