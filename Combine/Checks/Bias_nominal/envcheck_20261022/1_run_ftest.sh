#!/bin/bash
#
# Rebuild the background envelopes as they were before the bias-driven changes
# of AN Sec 7.4, with the current data, so that those decisions can be rechecked
# with independent toys (the earlier high-statistics runs reused the same seed
# for both splits of a chunk, i.e. only 5000 distinct toys).
#
#   mA2  : Bernstein back in, up to Bern2   (removed: Bern2 -0.222)
#   mA3  : Bernstein up to Bern2            (capped at 1: Bern2 -0.211)
#   mA4  : Bernstein back in, up to Bern4   (removed: Bern4 -0.221)
#   mA20 : Bernstein up to Bern6            (capped at 5: Bern6 +0.208)
#   mA24 : Laurent back in                  (removed: Lau1 -0.232)
#   mA25 : Bernstein up to Bern6 + Laurent  (capped/removed: Bern6 +0.212, Lau1 -0.215)
#
# Uses the patched copy of fTest_ALP_turnOn in build/ (env overrides only); the
# production binary and ALP_BkgModel_ReReco are not touched. Same options as
# shellScripts/bkg/fit_bkg.sh. One background process per mass point (6).
#
# Usage: bash 1_run_ftest.sh [--run]   (default: dry-run, print the commands)

set -eo pipefail
EC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CMSSW_TOP=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4
dir_input=/eos/home-p/pelai/HZa/root_MVAcut/data
gof_toys=100

declare -A BERN=( [2]=2 [3]=2 [4]=4 [20]=6 [25]=6 )
declare -A LAU=( [24]=1 [25]=1 )
masses=( 2 3 4 20 24 25 )

run=0; [[ "${1:-}" == "--run" ]] && run=1

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_TOP}/src" && eval "$(scramv1 runtime -sh)"

for m in "${masses[@]}"; do
  out="${EC}/bkg/${m}"
  bern="${BERN[$m]:-}"; lau="${LAU[$m]:-}"
  envs="BERN_ORDER_OVERRIDE=${bern:+${m}:${bern}} KEEP_LAURENT_MASSES=${lau:+${m}}"
  cmd="${EC}/build/bin/fTest_ALP_turnOn -i ${dir_input}/mA_M${m}/ws/run3.root --saveMultiPdf ${out}/CMS-HGG_mva_13p6TeV_multipdf.root -D ${out}/HZAmassInde_fTest --mass_ALP ${m} -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys ${gof_toys} -d ${out}/dat/fTest.dat"
  echo "[mA=${m}] env ${envs}"
  echo "  ${cmd}"
  [[ $run -eq 1 ]] || continue
  mkdir -p "${out}/dat"
  ( cd "$out" && env ${envs} setsid nohup ${cmd} > "${out}/ftest.stdout.log" 2>&1 & )
done
[[ $run -eq 1 ]] || echo "(dry-run; pass --run to start)"
