#!/bin/bash
# Re-run signalFit for the combinations whose fTest nGauss (nRV) changed during the
# 2026-08-09 full fTest re-run. signalFit reads
#   outdir_<ch>/nGaussScan/json/<mA>_nGauss_<era>_<ch>_Hm125.json
# at signalFit.py:561 whenever useDCB is False (the default), so a changed nRV makes
# the previously-produced signal workspace stale.
#
# calcPhotonSyst is NOT re-run: calcPhotonSyst.py has no nGauss dependency and the
# pkl files it produces already exist for all of these combinations.
#
# All 9 belong to the HZa project (mA = 1 and 8), so the input workspaces all come
# from /eos/home-p/pelai/HZa/root_MVAcut/sig.

set -uo pipefail

BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
CMSSW_TOP="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4"
SigScripts="${BaseDir}/Signal/scripts"
EosDir="/eos/home-p/pelai/HZa/root_MVAcut/sig"
LogDir="${BaseDir}/Signal/logs_sigfit_rerun_20260809"

NPAR="${NPAR:-6}"   # lxplus foreground parallelism cap

mkdir -p "$LogDir"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_TOP}/src"
eval "$(scramv1 runtime -sh)"
source "${BaseDir}/setup.sh"
export PYTHONPATH="${PYTHONPATH}:${BaseDir}/tools:${BaseDir}/Signal/tools"

COMBOS=(
  "1|2022postEE|ele"
  "1|2022preEE|ele"
  "1|2023postBPix|ele"
  "1|2023preBPix|ele"
  "1|2024|ele"
  "8|2022preEE|ele"
  "1|2023postBPix|mu"
  "1|2023preBPix|mu"
  "8|2022preEE|mu"
)

TASKS="${LogDir}/tasklist.txt"
printf "%s\n" "${COMBOS[@]}" > "$TASKS"

echo "[INFO] ${#COMBOS[@]} signalFit jobs, ${NPAR} at a time"
echo "[INFO] logs -> ${LogDir}"

run_one() {
  IFS='|' read -r mA yr ch <<< "$1"
  local lg="${LogDir}/${mA}_${yr}_${ch}.log"
  cd "${SigScripts}" || return 1
  if python3 "${SigScripts}/signalFit.py" --mass_ALP "${mA}" --year "${yr}" \
       --channel "${ch}" --inputWSDir "${EosDir}/mA_M${mA}/ws_Tree2WS" \
       --doSystematics > "${lg}" 2>&1; then
    echo "[OK]   ${mA} ${yr} ${ch}"
  else
    echo "[FAIL] ${mA} ${yr} ${ch}  (see ${lg})"
  fi
}
export -f run_one
export LogDir SigScripts EosDir

xargs -a "$TASKS" -d '\n' -I{} -P "${NPAR}" bash -c 'run_one "$@"' _ {}

echo "[DONE] $(date)"
