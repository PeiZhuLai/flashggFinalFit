#!/usr/bin/env bash
set -euo pipefail

baseDir="${baseDir:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit}"
script="$baseDir/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py"
logDir="$baseDir/MVAcut/run3_ReReco_Sys/logs/apply_bdt_sig"

mkdir -p "$logDir"

samples=(
  "mA_M1" "mA_M2" "mA_M3" "mA_M4" "mA_M5" "mA_M6" "mA_M7"
  "mA_M8" "mA_M9" "mA_M10" "mA_M15" "mA_M20" "mA_M25" "mA_M30"
)

years=(
  "2022preEE" "2022postEE" "2023preBPix" "2023postBPix" "2024"
)

max_parallel=6
batch_id=1
job_in_batch=0
total_jobs=0
pids=()
labels=()

status=0
wait_batch() {
  if [ "${#pids[@]}" -eq 0 ]; then
    return
  fi

  echo "[INFO] Waiting for batch ${batch_id} (${#pids[@]} jobs)"
  for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
      echo "[INFO] Finished ${labels[$i]}"
    else
      echo "[ERROR] Failed ${labels[$i]}; see $logDir/${labels[$i]}.log" >&2
      status=1
    fi
  done

  pids=()
  labels=()
  job_in_batch=0
  batch_id=$((batch_id + 1))
}

for sample in "${samples[@]}"; do
  for year in "${years[@]}"; do
    total_jobs=$((total_jobs + 1))
    job_in_batch=$((job_in_batch + 1))
    label="${sample}_${year}"
    log_file="$logDir/${label}.log"

    echo "[INFO] Start batch ${batch_id} job ${job_in_batch}/${max_parallel}: ${label}"
    python3 "$script" --samples "$sample" --years "$year" "$@" > "$log_file" 2>&1 &
    pids+=("$!")
    labels+=("$label")

    if [ "$job_in_batch" -eq "$max_parallel" ]; then
      wait_batch
    fi
  done
done

wait_batch

echo "[INFO] Submitted ${total_jobs} apply_bdt_sig jobs in ordered batches of ${max_parallel}."

exit "$status"
