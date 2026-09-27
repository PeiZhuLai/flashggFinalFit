#!/bin/bash
# Stage 4+5 for the combinations re-fitted by 3_rerun_sigfit_changed_ngauss.sh.
#
#   RunPlotter.py       -> outdir_<ch>/signalFit/Plots/effSigma_<mA>_<era>_<ch>.json (+ signal model pdf)
#   plot_effisigma.py   -> outdir_<ch>/signalFit/Plots/effSigmaVmA_<ch>.pdf
#
# Per-era runs cover only the 9 changed combinations; the run3 combined JSON is
# rebuilt for every (mA, channel) that had any era change, i.e. mA = 1 and 8 in
# both channels, because the combined value sums all five era workspaces.
#
# plot_effisigma.py reads EVERY effSigma json and rebuilds both channel plots, so
# it runs once at the end, after all RunPlotter jobs have landed.

set -uo pipefail

BaseDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
CMSSW_TOP="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4"
SigDir="${BaseDir}/Signal"
LogDir="${SigDir}/logs_plotter_rerun_20260809"
ALLYEARS='2022preEE,2022postEE,2023preBPix,2023postBPix,2024'

NPAR="${NPAR:-6}"

mkdir -p "$LogDir"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_TOP}/src"
eval "$(scramv1 runtime -sh)"
source "${BaseDir}/setup.sh"
export PYTHONPATH="${PYTHONPATH}:${BaseDir}/tools:${SigDir}/tools"

# mA|years|channel|label
TASKS="${LogDir}/tasklist.txt"
cat > "$TASKS" <<EOF
1|2022postEE|ele
1|2022preEE|ele
1|2023postBPix|ele
1|2023preBPix|ele
1|2024|ele
8|2022preEE|ele
1|2023postBPix|mu
1|2023preBPix|mu
8|2022preEE|mu
1|${ALLYEARS}|ele
1|${ALLYEARS}|mu
8|${ALLYEARS}|ele
8|${ALLYEARS}|mu
EOF

echo "[INFO] $(wc -l < "$TASKS") RunPlotter jobs, ${NPAR} at a time"
echo "[INFO] logs -> ${LogDir}"

run_one() {
  IFS='|' read -r mA yr ch <<< "$1"
  local tag="${yr//,/_}"
  [[ "$yr" == *","* ]] && tag="run3"
  local lg="${LogDir}/${mA}_${tag}_${ch}.log"
  cd "${SigDir}" || return 1
  if python3 "${SigDir}/RunPlotter.py" --mass_ALP "${mA}" --years "${yr}" \
       --channel "${ch}" > "${lg}" 2>&1; then
    echo "[OK]   RunPlotter ${mA} ${tag} ${ch}"
  else
    echo "[FAIL] RunPlotter ${mA} ${tag} ${ch}  (see ${lg})"
  fi
}
export -f run_one
export LogDir SigDir

xargs -a "$TASKS" -d '\n' -I{} -P "${NPAR}" bash -c 'run_one "$@"' _ {}

echo "[INFO] RunPlotter done, running plot_effisigma.py"
cd "${SigDir}"
if python3 "${SigDir}/scripts/plot_effisigma.py" > "${LogDir}/plot_effisigma.log" 2>&1; then
  echo "[OK]   plot_effisigma"
else
  echo "[FAIL] plot_effisigma  (see ${LogDir}/plot_effisigma.log)"
fi

echo "[DONE] $(date)"
