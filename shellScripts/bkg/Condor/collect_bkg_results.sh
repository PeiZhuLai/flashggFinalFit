#!/bin/bash

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -euo pipefail

BASE_DIR="${BASE_DIR:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit}"
lable="${LABLE:-run3}"
version="${VERSION:-ReReco}"

path_out_bkg="${BASE_DIR}/Background/ALP_BkgModel_${version}/fit_results_${lable}"
summary_log="${path_out_bkg}/background_fit_summary.csv"
failed_log="${path_out_bkg}/failed_mass_points.log"

mkdir -p "$path_out_bkg"

echo "mass,ftest_status,best_fit_pdf,n_pdfs,ftest_results,envelope_results,ftest_stdout,bkgplots_status,bkgplots_stdout" > "$summary_log"
: > "$failed_log"

while IFS= read -r per_mass_summary; do
  awk 'NR > 1' "$per_mass_summary" >> "$summary_log"
done < <(find "$path_out_bkg" -mindepth 2 -maxdepth 2 -name 'background_fit_summary_*.csv' | sort -V)

while IFS= read -r per_mass_failed_log; do
  if [[ -s "$per_mass_failed_log" ]]; then
    cat "$per_mass_failed_log" >> "$failed_log"
  fi
done < <(find "$path_out_bkg" -mindepth 2 -maxdepth 2 -name 'failed_mass_*.log' | sort -V)

echo "Merged summary: ${summary_log}"
echo "Merged failed log: ${failed_log}"
