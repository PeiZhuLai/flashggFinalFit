#!/usr/bin/env bash
#
# Collect per-mass bias-study plots into the flat directory that
# sync_figures.sh (BIAS_SRC) reads.
#
# The bias jobs (runjob_bias_study.sh) write their plots into per-mass subdirs
#   Combine/Checks/Bias_nominal/bias_outputs/mA_<m>/plots_bias/<m>_*.pdf
# but the AN sync source is the FLAT directory
#   Combine/Checks/Bias_nominal/plots_bias/
# so without this collect step the flat dir stays empty and the AN keeps the
# stale figures from a previous run. Run this after the bias jobs finish.
#
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "${script_dir}/../../.." && pwd)}"   # .../flashggFinalFit
BN="${BASE_DIR}/Combine/Checks/Bias_nominal"
OUTPUTS="${BN}/bias_outputs"
OUTPUTS_HS="${BN}/bias_outputs_highstat"
FLAT="${BN}/plots_bias"

mkdir -p "${FLAT}"

shopt -s nullglob
n=0
hs=0

# Nominal (1000-toy) plots first, so that any mass point that only exists there
# is still represented.
for d in "${OUTPUTS}"/mA_*/plots_bias; do
  [[ -d "${d}" ]] || continue
  for f in "${d}"/*.pdf "${d}"/*.png; do
    cp -p "${f}" "${FLAT}/" && n=$((n+1))
  done
done

# High-statistics (10000-toy) plots overwrite them where available. These are
# the numbers quoted in the AN, so they must win; the nominal run is only a
# fallback for mass points that have not been reprocessed.
for d in "${OUTPUTS_HS}"/mA_*/merged/plots_bias; do
  [[ -d "${d}" ]] || continue
  for f in "${d}"/*.pdf "${d}"/*.png; do
    cp -p "${f}" "${FLAT}/" && hs=$((hs+1))
  done
done

echo "[collect_bias] copied ${n} nominal + ${hs} high-statistics bias plot files into ${FLAT}"
