#!/bin/bash

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -o pipefail

if [[ $# -lt 7 ]]; then
  echo "Usage: $0 MASS LABEL VERSION LUMI DIR_INPUT PATH_OUT_BKG TOTAL_OUTDIR [BASE_DIR] [CMSSW_TOP]"
  exit 2
fi

mass="$1"
lable="$2"
version="$3"
int_lumi="$4"
dir_input="$5"
path_out_bkg="$6"
total_outdir="$7"
base_dir="${8:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit}"
cmssw_top="${9:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}"
gof_toys="${GOF_TOYS:-100}"

format_duration() {
  local total_seconds="${1:-0}"
  local hours=$(( total_seconds / 3600 ))
  local minutes=$(( (total_seconds % 3600) / 60 ))
  local seconds=$(( total_seconds % 60 ))
  printf "%02dh:%02dm:%02ds" "$hours" "$minutes" "$seconds"
}

write_summary() {
  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "$mass" "$ftest_status" "$best_fit_pdf" "$n_pdfs" "$ftest_results" \
    "$envelope_results" "$ftest_stdout" "$bkgplots_status" "$bkgplots_stdout" >> "$summary_log"
}

finish_skipped() {
  write_summary
  echo "[Timer] mA=${mass} finished in $(format_duration $(( $(date +%s) - mass_start_time )))" | tee -a "$failed_log"
  exit 0
}

finish_failed() {
  local step="$1"
  write_summary
  echo "[FAIL][${step}] mA=${mass}" | tee -a "$failed_log"
  echo "[Timer] mA=${mass} finished in $(format_duration $(( $(date +%s) - mass_start_time )))" | tee -a "$failed_log"
  exit 1
}

script_start_time=$(date +%s)
mass_start_time="$script_start_time"

echo ">>> Running background fit Condor job"
echo "  mass      = ${mass}"
echo "  lable     = ${lable}"
echo "  version   = ${version}"
echo "  int_lumi  = ${int_lumi}"
echo "  dir_input = ${dir_input}"
echo "  output    = ${path_out_bkg}"
echo "  GOF toys  = ${gof_toys}"

ulimit -s unlimited

if [[ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]]; then
  source /cvmfs/cms.cern.ch/cmsset_default.sh
else
  echo "[ERROR] Missing /cvmfs/cms.cern.ch/cmsset_default.sh"
  exit 1
fi

cd "${cmssw_top}/src" || exit 1
cmsenv

source "${base_dir}/setup.sh"
cd "${base_dir}/Background" || exit 1

if [[ ! -x ./bin/fTest_ALP_turnOn || ! -x ./bin/makeBkgPlots_ALP ]]; then
  echo "[ERROR] Background binaries missing. Run subjob_bkg.sh once to build before submitting."
  ls -l ./bin || true
  exit 1
fi

mkdir -p "${path_out_bkg}/${mass}" "$total_outdir"
path_bkg="${path_out_bkg}/${mass}"
ftest_outdir="${path_bkg}/HZAmassInde_fTest"
ftest_stdout="${path_bkg}/ftest.stdout.log"
bkgplots_stdout="${path_bkg}/makeBkgPlots.stdout.log"
ftest_results="${ftest_outdir}/fTestResults.txt"
envelope_results="${ftest_outdir}/EnvelopeResults.txt"
failed_log="${path_bkg}/failed_mass_${mass}.log"
summary_log="${path_bkg}/background_fit_summary_${mass}.csv"

: > "$failed_log"
echo "mass,ftest_status,best_fit_pdf,n_pdfs,ftest_results,envelope_results,ftest_stdout,bkgplots_status,bkgplots_stdout" > "$summary_log"

best_fit_pdf="NA"
n_pdfs="NA"
ftest_status="OK"
bkgplots_status="OK"

tree_input="${dir_input}/mA_M${mass}/run3.root"
ws_input="${dir_input}/mA_M${mass}/ws/run3.root"

if [[ ! -f "$tree_input" ]]; then
  ftest_status="MISSING_TREE"
  bkgplots_status="SKIP"
  echo "[FAIL][InputTreeMissing] mA=${mass} missing ${tree_input}" | tee -a "$failed_log"
  finish_skipped
fi

if [[ ! -f "$ws_input" ]]; then
  ftest_status="MISSING_WS"
  bkgplots_status="SKIP"
  echo "[FAIL][WorkspaceMissing] mA=${mass} missing ${ws_input}" | tee -a "$failed_log"
  echo "  Regenerate workspaces first, e.g. run ${base_dir}/Trees2WS/run_tree2ws.sh" | tee -a "$failed_log"
  finish_skipped
fi

if [[ "$ws_input" -ot "$tree_input" ]]; then
  ftest_status="STALE_WS"
  bkgplots_status="SKIP"
  echo "[FAIL][WorkspaceStale] mA=${mass} workspace is older than input tree" | tee -a "$failed_log"
  echo "  tree: ${tree_input}" | tee -a "$failed_log"
  echo "  ws  : ${ws_input}" | tee -a "$failed_log"
  echo "  Regenerate workspaces first, e.g. run ${base_dir}/Trees2WS/run_tree2ws.sh" | tee -a "$failed_log"
  finish_skipped
fi

ftest_start_time=$(date +%s)
./bin/fTest_ALP_turnOn \
  -i "$ws_input" \
  --saveMultiPdf "${path_bkg}/CMS-HGG_mva_13p6TeV_multipdf.root" \
  -D "$ftest_outdir" \
  --mass_ALP "$mass" \
  -c 1 \
  --isFlashgg 0 \
  --isData 0 \
  -f data, \
  --mhLow 95 \
  --mhHigh 180 \
  --mhLowBlind 115 \
  --mhHighBlind 135 \
  --gtoys "$gof_toys" \
  > "$ftest_stdout" 2>&1
ftest_cmd_status=$?
ftest_end_time=$(date +%s)
echo "[Timer] mA=${mass} fTest finished in $(format_duration $((ftest_end_time - ftest_start_time)))"

if [[ $ftest_cmd_status -ne 0 ]]; then
  ftest_status="FAIL"
  finish_failed "fTest"
fi

if [[ -f "$ftest_stdout" ]]; then
  best_fit_pdf=$(grep -E "Best Fit Pdf =" "$ftest_stdout" | tail -1 | sed 's/.*Best Fit Pdf = //' | sed 's/, /:/' || true)
  n_pdfs=$(grep -E "with a total of [0-9]+ pdfs" "$ftest_stdout" | tail -1 | sed -E 's/.*with a total of ([0-9]+) pdfs.*/\1/' || true)
  [[ -n "$best_fit_pdf" ]] || best_fit_pdf="NA"
  [[ -n "$n_pdfs" ]] || n_pdfs="NA"
fi

bkgplots_start_time=$(date +%s)
./bin/makeBkgPlots_ALP \
  -b "${path_bkg}/CMS-HGG_mva_13p6TeV_multipdf.root" \
  -d "${path_bkg}/BkgPlots" \
  --total_OutDir "$total_outdir" \
  -o "${path_bkg}/BkgPlots.root" \
  --sqrts 13p6TeV \
  --isMultiPdf \
  --useBinnedData \
  --massStep 2.5 \
  --mhVal 125.0 \
  --maVal "$mass" \
  --mhLow 95 \
  --mhHigh 180 \
  --mhLowBlind 115 \
  --mhHighBlind 135 \
  --intLumi "$int_lumi" \
  -c 0 \
  --isFlashgg 0 \
  --doBands \
  > "$bkgplots_stdout" 2>&1
bkgplots_cmd_status=$?
bkgplots_end_time=$(date +%s)
echo "[Timer] mA=${mass} makeBkgPlots finished in $(format_duration $((bkgplots_end_time - bkgplots_start_time)))"

if [[ $bkgplots_cmd_status -ne 0 ]]; then
  bkgplots_status="FAIL"
  finish_failed "BkgPlots"
fi

write_summary

echo "[OK] mA=${mass}"
echo "[Timer] mA=${mass} finished in $(format_duration $(( $(date +%s) - mass_start_time )))"
echo "[Timer] Total runjob_bkg.sh runtime: $(format_duration $(( $(date +%s) - script_start_time )))"
