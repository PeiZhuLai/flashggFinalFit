#!/usr/bin/env bash
# mA3 working-point closure scan (tab:closure_mA3_scan) with the current signal model: the signal is
# reweighted like the background MC and the datacard is the adopted one at 0.988. The background-only
# pseudo-data envelopes of reopt_mA3/cut_<tag>/ps/fit depend only on the background and are reused;
# only text2workspace, Significance and MultiDimFit (r_hat) are redone, with the options of
# scan_mA3_reopt.sh. Step [0] reruns the old datacard at 0.988 and must reproduce the old result.
# Output: pseudodata_closure/reopt_mA3_sigrw/{cut_<tag>/, summary.txt}; reopt_mA3/ is not modified.
set -u
CUTS="${1:-0.980 0.984 0.986 0.988 0.990 0.992 0.994 0.996}"
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure
OLD=$PC/reopt_mA3
OUT=$PC/reopt_mA3_sigrw; mkdir -p $OUT
realDC=$FFIT/Datacard/output_Datacard_leptons/3_pruned_datacard_leptons.txt
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $FFIT/.. && eval `scramv1 runtime -sh`
source $FFIT/setup.sh >/dev/null 2>&1
set -u

fit() {   # $1 dir with datacard.txt, $2 tag
  ( cd $1 && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
    combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n R$2 > sig.log 2>&1
    combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n R$2 > rhat.log 2>&1 )
  echo "closureZ=$(grep -E 'Significance:' $1/sig.log | sed -E 's/.*Significance: //') rhat=$(grep -E '^\s+r :' $1/rhat.log | sed -E 's/^\s+r :\s+//; s/ \(68%\)//')"
}

echo "=== [0] reproduce the old result at 0.988 with the old datacard"
V=$OUT/check_old_0988; mkdir -p $V && cp -p $OLD/cut_0988/ps/datacard.txt $V/
r=$(fit $V check); echo "  $r   (old: $(grep -oE 'closureZ=[^ ]+ rhat=[^|]+' $OLD/cut_0988/result.txt))" | tee $OUT/check_old_0988.txt

echo "=== [1] current signal model"
one() {
  local cut=$1 tag=$(echo $1 | tr -d '.') OB=$OUT/cut_$(echo $1 | tr -d '.')
  mkdir -p $OB
  [ -s $OLD/cut_$tag/ps/fit/CMS-HGG_mva_13p6TeV_multipdf.root ] || { echo "$cut FAIL no pseudo-data envelope" > $OB/result.txt; return; }
  sed "s#/Background/ALP_BkgModel_ReReco/fit_results_run3/3/#/pseudodata_closure/reopt_mA3/cut_${tag}/ps/fit/#g" "$realDC" > $OB/datacard.txt
  grep -q "reopt_mA3/cut_${tag}/ps/fit/CMS-HGG_mva_13p6TeV_multipdf.root" $OB/datacard.txt || { echo "$cut FAIL datacard" > $OB/result.txt; return; }
  echo "$cut $(fit $OB $tag)" > $OB/result.txt
}
n=0
for c in $CUTS; do one $c & n=$((n+1)); [ $((n % 3)) -eq 0 ] && wait; done; wait
for c in $CUTS; do cat $OUT/cut_$(echo $c | tr -d '.')/result.txt; done | tee $OUT/summary.txt
