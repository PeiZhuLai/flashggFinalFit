#!/usr/bin/env bash
set -euo pipefail

baseDir="${baseDir:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit}"
script="$baseDir/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py"
logDir="$baseDir/MVAcut/run3_ReReco_Sys/logs/apply_bdt_sig"

mkdir -p "$logDir"

sample_chunks=(
  "mA_M1,mA_M2,mA_M3"
  "mA_M4,mA_M5,mA_M6"
  "mA_M7,mA_M8"
  "mA_M9,mA_M10"
  "mA_M15,mA_M20"
  "mA_M25,mA_M30"
)

pids=()
for i in "${!sample_chunks[@]}"; do
  job_id=$((i + 1))
  samples="${sample_chunks[$i]}"
  log_file="$logDir/job_${job_id}.log"
  echo "[INFO] Start apply_bdt_sig job ${job_id}/6: ${samples}"
  python3 "$script" --samples "$samples" "$@" > "$log_file" 2>&1 &
  pids+=("$!")
done

status=0
for i in "${!pids[@]}"; do
  job_id=$((i + 1))
  if wait "${pids[$i]}"; then
    echo "[INFO] apply_bdt_sig job ${job_id}/6 finished"
  else
    echo "[ERROR] apply_bdt_sig job ${job_id}/6 failed; see $logDir/job_${job_id}.log" >&2
    status=1
  fi
done

exit "$status"
