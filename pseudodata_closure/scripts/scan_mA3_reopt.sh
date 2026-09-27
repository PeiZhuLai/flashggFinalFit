#!/usr/bin/env bash
# mA3 working-point re-optimization on the FSR-fix re-production (2026-09-27, user request).
#
# Two constraints failed so far: at 0.940 (R 1.10) no background family describes the data
# (GOF ~ 0); at the adopted 0.988 (R 1.21) the background-only pseudo-data closure shows a fake
# signal r_hat = +0.063 (Z = 2.3 sigma, ~90% of the expected limit). Looser cuts are worse on
# both. This scan measures, per candidate cut:
#   R, Asimov Z          scan_score_R_significance.py machinery (stored scores, raw weights)
#   data N, data GOF     envelope fitted to real data at that cut (build_data_ws_cut.py + fTest,
#                        100 GOF toys as in the official fit_bkg.sh)
#   closure Z, r_hat     background-only MC pseudo-data (build_pseudodata_ws.py) with the real
#                        mA3 signal model; same combine settings as the closure
# Proposed criterion (to be confirmed by the user): best Asimov Z with closure Z <~ 1.5 sigma
# (the level accepted for mA2 in 2026-07) and all data-envelope GOF > 0.01.
#
# Output: pseudodata_closure/reopt_mA3/{cut_<tag>/, summary.txt}; nothing in the analysis chain
# is modified. At most 3 cuts in parallel, one fTest each (<= 6 foreground processes).
set -u
CUTS="${1:-0.980 0.984 0.986 0.988 0.990 0.992 0.994 0.996}"
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
OUT=$PC/reopt_mA3; mkdir -p $OUT
realDC=$FFIT/Datacard/output_Datacard_leptons/3_pruned_datacard_leptons.txt

echo "=== [0] R and Asimov Z per cut"
( cd $PLOT && env -i HOME=$HOME PATH=/usr/bin:/bin PYTHONPATH=$PLOT/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts \
  $CPY - $CUTS <<'PY'
import sys, numpy as np
src = open("scripts/scan_score_R_significance.py").read()
g = {"__name__": "scan_mod"}; exec(compile(src[:src.index("def main")], "scan", "exec"), g)
pts, bk = g["scan"](3)
for c in map(float, sys.argv[1:]):
    print("RZ %.3f %.3f %.2f" % (c, g["_R_direct"](bk, c), float(np.interp(c, pts[:, 0], pts[:, 2]))))
PY
) 2>&1 | grep '^RZ' | tee $OUT/rz.txt

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $FFIT/.. && eval `scramv1 runtime -sh`
source $FFIT/setup.sh >/dev/null 2>&1
set -u

one_cut() {
  local cut=$1 tag=$(echo $1 | tr -d '.') OB=$OUT/cut_$(echo $1 | tr -d '.')
  mkdir -p $OB
  # --- real data envelope at this cut
  ( cd $PC/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 \
      $CPY build_data_ws_cut.py --mA 3 --cut $cut --outbase $OB/data ) > $OB/data.build.log 2>&1
  mkdir -p $OB/data/fit/ft
  ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $OB/data/mA_M3/ws/run3.root \
      --saveMultiPdf $OB/data/fit/CMS-HGG_mva_13p6TeV_multipdf.root -D $OB/data/fit/ft \
      --mass_ALP 3 -c 1 --isFlashgg 0 --isData 0 -f data, \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $OB/data.ftest.log 2>&1
  # --- background-only MC pseudo-data closure at this cut
  ( cd $PC/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 \
      $CPY build_pseudodata_ws.py --mA 3 --cut $cut --outbase $OB/ps ) > $OB/ps.build.log 2>&1
  mkdir -p $OB/ps/fit/ft
  ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $OB/ps/mA_M3/ws/run3.root \
      --saveMultiPdf $OB/ps/fit/CMS-HGG_mva_13p6TeV_multipdf.root -D $OB/ps/fit/ft \
      --mass_ALP 3 -c 1 --isFlashgg 0 --isData 0 -f data, \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 30 ) > $OB/ps.ftest.log 2>&1
  sed "s#/Background/ALP_BkgModel_ReReco/fit_results_run3/3/#/pseudodata_closure/reopt_mA3/cut_${tag}/ps/fit/#g" "$realDC" > $OB/ps/datacard.txt
  grep -q "reopt_mA3/cut_${tag}/ps/fit/CMS-HGG_mva_13p6TeV_multipdf.root" $OB/ps/datacard.txt || { echo "FAIL datacard $cut" > $OB/result.txt; return; }
  ( cd $OB/ps && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
    combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n R$tag > sig.log 2>&1
    combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n R$tag > rhat.log 2>&1 )
  local nd=$(grep -oE 'data N=[0-9]+' $OB/data.build.log | grep -oE '[0-9]+')
  local np=$(grep -oE 'pseudo-data N=[0-9]+' $OB/ps.build.log | grep -oE '[0-9]+')
  local Z=$(grep -E 'Significance:' $OB/ps/sig.log | sed -E 's/.*Significance: //')
  local rh=$(grep -E '^\s+r :' $OB/ps/rhat.log | sed -E 's/^\s+r :\s+//; s/ \(68%\)//')
  local dg=$(grep -oE 'pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+' $OB/data/fit/ft/EnvelopeResults.txt 2>/dev/null | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')
  local fb=$(grep -c 'fallback' $OB/data/fit/ft/EnvelopeResults.txt 2>/dev/null)
  echo "$cut dataN=${nd:-NA} psN=${np:-NA} closureZ=${Z:-NA} rhat=${rh:-NA} data_fallback=${fb:-NA} | $dg" > $OB/result.txt
}

echo "=== [1] per-cut data envelope + closure (3 in parallel)"
running=0
for cut in $CUTS; do
  one_cut $cut & running=$((running+1))
  if [ $running -ge 3 ]; then wait -n; running=$((running-1)); fi
done
wait

echo "=== summary"
{
  printf "%-6s %-6s %-6s | %s\n" cut R Z_asimov "data N, closure Z, r_hat, data-envelope fallbacks | data GOF per member"
  for cut in $CUTS; do
    rz=$(awk -v c=$cut '$2==c{printf "%-6s %-6s", $3, $4}' $OUT/rz.txt)
    printf "%-6s %s | %s\n" "$cut" "$rz" "$(cut -d' ' -f2- $OUT/cut_$(echo $cut | tr -d '.')/result.txt)"
  done
} | tee $OUT/summary.txt
