#!/usr/bin/env bash
#
# Watchdog for the envelope recheck, condor cluster 13052752
# (mA 2,3,4,20,24,25 x 10 chunks x 1000 toys = 60 jobs, restored pre-change envelopes).
#
# Polls every 10 min and writes a heartbeat. A chunk counts as done when its
# BiasFits dir holds split0 and split1 of exactly the envelope functions listed
# for that mA in expected_functions.txt, each newer than the submission. The
# job's [OK] line cannot be used: the submit file sends stdout to /dev/null.
#
# When the cluster leaves the queue it checks every chunk and every exit code
# (condor_history); if all are good it runs merge_bias_highstat.sh one mass at
# a time (so one failing mass does not stop the rest) and appends DONE to the
# registry. Otherwise it appends INCOMPLETE with the bad (mA, chunk) list.
#
# Only numbers are written here; reading the bias values is left for later.

set -uo pipefail

CLUSTER=13052752
WANT=60
MASSES="2 3 4 20 24 25"
SUBMITTED="2026-10-22 14:06"
FFF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
BN=${FFF}/Combine/Checks/Bias_nominal
OUT=${BN}/envcheck_20261022/outputs
HERE=${BN}/envcheck_20261022
EXPECT=${HERE}/expected_functions.txt
STATUS=${HERE}/status.txt
MERGE_LOG=${HERE}/merge.log
REG=${BN}/PENDING_RESULTS.txt
NAME="HZa bias envcheck ${CLUSTER}"
MERGE=${FFF}/shellScripts/bias/Condor/merge_bias_highstat.sh

# Prints nothing for a complete chunk, otherwise a short reason.
chunk_problem() {
  local m=$1 c=$2 d="${OUT}/mA_$1/chunk$2/BiasFits" want got fn s
  want=$(awk -v m="$m" '$1==m{for(i=2;i<=NF;i++)print $i}' "$EXPECT" | sort | tr '\n' ' ')
  [ -d "$d" ] || { echo "nodir"; return; }
  got=$(find "$d" -name 'biasStudy_*_split*_fits.root' -newermt "$SUBMITTED" -size +1k -printf '%f\n' \
        | sed 's/^biasStudy_//; s/_split.*//' | sort -u | tr '\n' ' ')
  [ "$got" = "$want" ] || { echo "funcs[${got}]"; return; }
  for fn in $want; do
    for s in 0 1; do
      [ -n "$(find "$d" -name "biasStudy_${fn}_split${s}_fits.root" -newermt "$SUBMITTED" -size +1k)" ] \
        || { echo "missing_${fn}_split${s}"; return; }
    done
  done
}

count_ok() {
  local n=0 m c
  for m in $MASSES; do for c in $(seq 0 9); do
    [ -z "$(chunk_problem "$m" "$c")" ] && n=$((n+1))
  done; done
  echo "$n"
}

held_reported=0
for _ in $(seq 1 300); do
  q=$(condor_q "$CLUSTER" -af JobStatus NumJobStarts 2>/dev/null) ; qrc=$?
  n_ok=$(count_ok)
  if [ "$qrc" -ne 0 ]; then
    echo "poll $(date '+%F %H:%M:%S') condor_q FAILED rc=$qrc chunks_ok=${n_ok}/${WANT}" > "$STATUS"
    sleep 600; continue
  fi
  n_q=$(printf '%s\n' "$q" | grep -c . || true)
  n_run=$(printf '%s\n' "$q" | awk '$1==2' | wc -l)
  n_idle=$(printf '%s\n' "$q" | awk '$1==1' | wc -l)
  n_held=$(printf '%s\n' "$q" | awk '$1==5' | wc -l)
  n_dead=$(printf '%s\n' "$q" | awk '$1==5 && $2>=3' | wc -l)
  echo "poll $(date '+%F %H:%M:%S') queue=${n_q} R=${n_run} I=${n_idle} H=${n_held} (perm ${n_dead}) chunks_ok=${n_ok}/${WANT}" > "$STATUS"
  if [ "$n_dead" -gt 0 ] && [ "$held_reported" -eq 0 ]; then
    echo "$(date '+%F %H:%M') HELD  ${NAME}  ${n_dead} jobs held after 3 starts" >> "$REG"
    held_reported=1
  fi
  [ "$n_q" -eq 0 ] && break
  sleep 600
done

if [ "$n_q" -ne 0 ]; then
  echo "$(date '+%F %H:%M') TIMEOUT  ${NAME}  queue still ${n_q}, chunks_ok=${n_ok}/${WANT}" >> "$REG"
  exit 1
fi

# Per-chunk reconciliation and exit codes.
bad=()
for m in $MASSES; do
  for c in $(seq 0 9); do
    p=$(chunk_problem "$m" "$c")
    [ -n "$p" ] && bad+=("mA${m}c${c}:${p}")
  done
done
nonzero=$(condor_history "$CLUSTER" -limit 60 -af ProcId ExitCode 2>/dev/null | awk '$2!=0' | tr '\n' ' ')
nhist=$(condor_history "$CLUSTER" -limit 60 -af ProcId 2>/dev/null | sort -u | wc -l)
if [ "${#bad[@]}" -gt 0 ] || [ -n "$nonzero" ]; then
  echo "$(date '+%F %H:%M') INCOMPLETE  ${NAME}  bad chunks: ${bad[*]:-none}; nonzero exit: ${nonzero:-none}; history ${nhist}/${WANT}" >> "$REG"
  exit 1
fi

: > "$MERGE_LOG"
fail=()
for m in $MASSES; do
  echo "##### merge mA=${m} $(date '+%F %H:%M:%S')" >> "$MERGE_LOG"
  BIAS_HS_ROOT="$OUT" ROOT_DATACARD_PATH="${HERE}/t2w" bash "$MERGE" "$m" >> "$MERGE_LOG" 2>&1 || fail+=("$m")
  [ -n "$(find "${OUT}/mA_${m}/merged/BiasJson" -name "${m}_gaussfit.json" -newermt "$SUBMITTED" -size +0 2>/dev/null)" ] \
    || fail+=("${m}nojson")
done
echo "##### MERGE_END $(date '+%F %H:%M:%S') fail=${fail[*]:-none}" >> "$MERGE_LOG"

if [ "${#fail[@]}" -gt 0 ]; then
  echo "$(date '+%F %H:%M') MERGEFAIL  ${NAME}  60/60 chunks ok, merge failed for: ${fail[*]} -> ${MERGE_LOG}" >> "$REG"
  exit 1
fi
echo "$(date '+%F %H:%M') DONE  ${NAME}  60/60 chunks ok (history ${nhist}), merged 6/6 -> ${MERGE_LOG}" >> "$REG"
