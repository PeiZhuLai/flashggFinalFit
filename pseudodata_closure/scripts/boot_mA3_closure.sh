#!/usr/bin/env bash
# MC-statistics bootstrap of the mA3 pseudo-data closure (2026-09-27).
#
# Why: in the mA3 working-point re-optimization (scan_mA3_reopt.sh) the closure fake signal
# jumps non-monotonically with the cut (Z = 2.3, 2.2, 0.8, 1.4, 2.1 for 0.988..0.996) while the
# background MC in the 120-130 window has only n_eff ~ 36 (0.988) / 20 (0.992) events. Each
# replica rebuilds the pseudo-data with Poisson(1) weights per MC event and reruns the same
# envelope + combine, so the replica spread of r_hat is the MC-statistics uncertainty of the
# closure. If it is comparable to the observed r_hat, the closure cannot tell sculpting from
# MC fluctuation at these cuts.
#
# Usage: boot_mA3_closure.sh "<cuts>" <n_replicas>     (default "0.988 0.992" 40)
# Output: pseudodata_closure/reopt_mA3/boot/cut_<tag>/rep_<seed>/ and boot_summary.txt
# Parallel: 6 replicas at a time (lxplus foreground limit).
set -u
CUTS="${1:-0.988 0.992}"; NREP="${2:-40}"
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
BOOT=$PC/reopt_mA3/boot; mkdir -p $BOOT
realDC=$FFIT/Datacard/output_Datacard_leptons/3_pruned_datacard_leptons.txt

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $FFIT/.. && eval `scramv1 runtime -sh`
source $FFIT/setup.sh >/dev/null 2>&1
set -u

one_rep() {
  local cut=$1 seed=$2 tag=$(echo $1 | tr -d '.')
  local R=$BOOT/cut_$tag/rep_$seed; mkdir -p $R/fit/ft
  ( cd $PC/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 \
      $CPY build_pseudodata_ws.py --mA 3 --cut $cut --bootstrap $seed --outbase $R ) > $R/build.log 2>&1
  ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $R/mA_M3/ws/run3.root \
      --saveMultiPdf $R/fit/CMS-HGG_mva_13p6TeV_multipdf.root -D $R/fit/ft \
      --mass_ALP 3 -c 1 --isFlashgg 0 --isData 0 -f data, \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 30 ) > $R/ftest.log 2>&1
  sed "s#/Background/ALP_BkgModel_ReReco/fit_results_run3/3/#/pseudodata_closure/reopt_mA3/boot/cut_${tag}/rep_${seed}/fit/#g" "$realDC" > $R/datacard.txt
  grep -q "boot/cut_${tag}/rep_${seed}/fit/CMS-HGG" $R/datacard.txt || { echo "$cut $seed FAIL_datacard" > $R/result.txt; return; }
  ( cd $R && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
    combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n B > sig.log 2>&1
    combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n B > rhat.log 2>&1
    rm -f ws_t2w.root )
  local Z=$(grep -E 'Significance:' $R/sig.log | sed -E 's/.*Significance: //')
  local rh=$(grep -E '^\s+r :' $R/rhat.log | awk '{print $3}')
  local N=$(grep -oE 'pseudo-data N=[0-9]+' $R/build.log | grep -oE '[0-9]+')
  echo "$cut $seed N=${N:-NA} Z=${Z:-NA} rhat=${rh:-NA}" > $R/result.txt
}
export -f one_rep; export BOOT PC CPY FFIT realDC

echo "=== bootstrap: cuts [$CUTS], $NREP replicas each, 6 parallel ($(date))"
for cut in $CUTS; do for s in $(seq 1 $NREP); do echo "$cut $s"; done; done \
  | xargs -P 6 -n 2 bash -c 'one_rep "$0" "$1"'

echo "=== summary ($(date))"
{
  for cut in $CUTS; do
    tag=$(echo $cut | tr -d '.')
    cat $BOOT/cut_$tag/rep_*/result.txt 2>/dev/null | $CPY -c "
import sys, numpy as np
rows=[l.split() for l in sys.stdin if 'rhat=' in l]
r=np.array([float(x[4][5:]) for x in rows if x[4][5:] not in ('NA','')])
z=np.array([float(x[3][2:]) for x in rows if x[3][2:] not in ('NA','')])
print('cut $cut: replicas ok %d/%d  rhat mean %+.4f  sd %.4f  median %+.4f  |  Z mean %.2f sd %.2f  frac(Z>2) %.2f'
      % (r.size, len(rows), r.mean(), r.std(ddof=1), np.median(r), z.mean(), z.std(ddof=1), (z>2).mean()))"
  done
} | tee $BOOT/boot_summary.txt
