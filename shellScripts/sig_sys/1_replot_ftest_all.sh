#!/bin/bash
# Re-run fTest over every (mA, era, channel) combination so the chi2-vs-nGauss
# plots are regenerated with the adaptive y-axis range (plottingTools.py).
#
# The plot can only be produced by re-fitting: plotFTestResults() takes live
# SimultaneousFit objects, nothing is cached to disk.
#
# Two projects share this Plots directory and have DIFFERENT input workspaces:
#   HZa        : mA = 1..30            x 5 eras  -> /eos/home-p/pelai/HZa/root_MVAcut
#   HZa_merged : mA = 0p1..0p9 (sub-GeV, merged photon), 2024 only
#                                                -> /eos/cms/store/group/phys_susy/pelai/HZa_merged/root_MVAcut
#
# NOTE: fTest.py also rewrites outdir_<ch>/nGaussScan/json/<mA>_nGauss_<era>_<ch>_Hm125.json,
# which signalFit consumes. Back that directory up before running and diff it after.

set -uo pipefail

BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
CMSSW_TOP="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4"
SigScripts="${BaseDir}/Signal/scripts"
LogDir="${BaseDir}/Signal/logs_replot_20260809"

EOS_HZA="/eos/home-p/pelai/HZa/root_MVAcut/sig"
EOS_MERGED="/eos/cms/store/group/phys_susy/pelai/HZa_merged/root_MVAcut/sig"

NPAR="${NPAR:-6}"   # lxplus foreground parallelism cap

mkdir -p "$LogDir"

# ---------------------------------------------------------------- environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_TOP}/src"
eval "$(scramv1 runtime -sh)"
source "${BaseDir}/setup.sh"
export PYTHONPATH="${PYTHONPATH}:${BaseDir}/tools:${BaseDir}/Signal/tools"

# ---------------------------------------------------------------- task list
TASKS="${LogDir}/tasklist.txt"
: > "$TASKS"

mAs_hza=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
years_hza=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024 )
mAs_merged=( 0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9 )
channels=( ele mu )

for mA in "${mAs_hza[@]}"; do
  for yr in "${years_hza[@]}"; do
    for ch in "${channels[@]}"; do
      echo "${mA}|${yr}|${ch}|${EOS_HZA}/mA_M${mA}/ws_Tree2WS" >> "$TASKS"
    done
  done
done
for mA in "${mAs_merged[@]}"; do
  for ch in "${channels[@]}"; do
    echo "${mA}|2024|${ch}|${EOS_MERGED}/mA_M${mA}/ws_Tree2WS" >> "$TASKS"
  done
done

echo "[INFO] $(wc -l < "$TASKS") fTest jobs, ${NPAR} at a time"
echo "[INFO] logs -> ${LogDir}"

# ---------------------------------------------------------------- run
run_one() {
  IFS='|' read -r mA yr ch ws <<< "$1"
  local lg="${LogDir}/${mA}_${yr}_${ch}.log"
  cd "${SigScripts}" || return 1
  if python3 "${SigScripts}/fTest.py" --mass_ALP "${mA}" --year "${yr}" --mass 125 \
       --channel "${ch}" --inputWSDir "${ws}" > "${lg}" 2>&1; then
    echo "[OK]   ${mA} ${yr} ${ch}"
  else
    echo "[FAIL] ${mA} ${yr} ${ch}  (see ${lg})"
  fi
}
export -f run_one
export LogDir SigScripts

xargs -a "$TASKS" -d '\n' -I{} -P "${NPAR}" bash -c 'run_one "$@"' _ {}

echo "[DONE] $(date)"
