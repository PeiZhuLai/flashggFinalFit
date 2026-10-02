#!/usr/bin/env bash
# Background-only pseudo-data closure for all 30 mass points on the FSR-fix re-production
# (user request 2026-09-27; the figure in AN Sec 6.4.2 was still the pre-FSR-fix one).
#
# Same procedure as the July closure (build_pseudodata_ws -> fTest envelope -> real datacard with
# the pseudo-data multipdf -> AsymptoticLimits / Significance / FitDiagnostics), with
#   CLOSURE_FSRFIX=1  : run3_bdt_scored_fsrfix inputs, 2022 DY+gamma = PTG-10to100, and the
#                       working points of the current MVAcut_points_run3.json (nearest anchor for
#                       interpolated points; validated against the real data roohists)
#   the fTest binary built 2026-09-27 (toy GOF everywhere, open step-Gaussian ranges at mA24-29)
# It waits until the mA24-30 real-data refit (unit hza-open2429) has finished its fits, so the two
# never exceed the lxplus limit of 6 parallel processes.
#
# Output: pseudodata_closure/{pseudodata,fit_results_run3_pseudodata,datacards,root_t2w,
#         output_limits,output_significance,output_fitdiag,closure_results.json,plots}
#         + figure_ALP/run3/pseudodata_bkg_alias/{bkgplot,allPdfs}_<mA>.pdf in the AN.
# Old (pre-FSR-fix) results are moved to pseudodata_closure/preFSRfix_20260927/.
set -u
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure
S=$PC/scripts
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
L=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix
AN=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172
export CLOSURE_FSRFIX=1
die() { echo "STOPPING: $*"; exit 1; }

START="${START:-0}"
if [ "$START" -le 2 ]; then
echo "=== [0] wait for the mA24-30 refit chain (fits + limits) to finish ($(date '+%T'))"
while systemctl --user is-active -q hza-open2429; do sleep 60; done
grep -q STOPPING $L/open2429.log && die "open2429 chain stopped -- fix it first"
echo "    go ($(date '+%T'))"

echo "=== [1] move the pre-FSR-fix closure results aside"
B=$PC/preFSRfix_20260927; mkdir -p $B
for d in pseudodata fit_results_run3_pseudodata datacards root_t2w output_limits output_significance output_fitdiag plots closure_results.json; do
  [ -e $PC/$d ] && mv $PC/$d $B/
done
mv $PC/higgsCombinePS*.root $B/ 2>/dev/null
mkdir -p $PC/plots

echo "=== [2] pseudo-data workspaces (30 points)"
( cd $S && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 $CPY build_pseudodata_ws.py --mA 1-30 ) > $L/closure30_build.log 2>&1 || die "build failed"
n=$(ls $PC/pseudodata/mA_M*/ws/run3.root 2>/dev/null | wc -l); [ $n -eq 30 ] || die "only $n pseudo-data workspaces"
grep -E "^mA" $L/closure30_build.log | sed 's/ -> .*//; s/^/    /'

fi
if [ "$START" -le 6 ]; then
echo "=== [3] envelope fits (fTest, 100 GOF toys, 6 parallel)"
NPROC=6 GOF_TOYS=100 bash $S/run_ftest_pseudodata.sh > $L/closure30_ftest.log 2>&1
n=$(ls $PC/fit_results_run3_pseudodata/*/CMS-HGG_mva_13p6TeV_multipdf.root 2>/dev/null | wc -l); [ $n -eq 30 ] || die "only $n multipdf files"

echo "=== [4] limits / significance / fitdiag"
NPROC=6 bash $S/run_limits_pseudodata.sh > $L/closure30_limits.log 2>&1
n=$(ls $PC/output_limits/*.root 2>/dev/null | wc -l); [ $n -eq 30 ] || die "only $n limit files"

echo "=== [5] results table + summary figure"
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd $FFIT/.. && eval `scramv1 runtime -sh`; set -u
( cd $S && python3 extract_results.py ) > $L/closure30_extract.log 2>&1 || die "extract_results failed"
( cd $S && env -i HOME=$HOME PATH=/usr/bin:/bin $CPY plot_closure.py ) > $L/closure30_plot.log 2>&1 || die "plot_closure failed"
tail -40 $L/closure30_extract.log

echo "=== [6] AN figure copies (summary plot now; bkg-fit panels after step 7)"
mkdir -p $AN/figure_ALP/run3/HiggsZaAna/flashggFinalFit/pseudodata_closure/plots
cp -p $PC/plots/sculpting_closure_pseudodata.{pdf,png} $AN/figure_ALP/run3/HiggsZaAna/flashggFinalFit/pseudodata_closure/plots/

fi
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd $FFIT/.. && eval `scramv1 runtime -sh`; set -u
echo "=== [7] background-fit panels (unblind only, 4 parallel, 45 min cap per mass then no bands)"
# The AN uses only the unblind panels. --doBands profiles the envelope in every m_llgg bin and
# diverged for mA25 on the rounded low-statistics pseudo-data in August (> 8 h), hence the cap.
set +u; source $FFIT/setup.sh >/dev/null 2>&1; set -u
U=$PC/fit_results_run3_pseudodata/AllFitResults_unblind; mkdir -p $U
panel() {
  local m=$1 D=$PC/fit_results_run3_pseudodata/$1
  local c=( -b $D/CMS-HGG_mva_13p6TeV_multipdf.root --sqrts 13p6TeV --isMultiPdf --useBinnedData
            --massStep 2.5 --mhVal 125.0 --maVal $m --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135
            --intLumi 172.13 -c 0 --isFlashgg 0 --dataLegend "MC pseudo-data" --ratioYTitle "MC - Fit" --unblind
            -d $D/BkgPlots_unblind --total_OutDir $U -o $D/BkgPlots_unblind.root )
  ( cd $FFIT/Background && timeout 45m ./bin/makeBkgPlots_ALP "${c[@]}" --doBands > $D/makeBkgPlots.stdout.log 2>&1 ) \
    || { echo "[mA=$m] bands timed out/failed -> no bands"; ( cd $FFIT/Background && ./bin/makeBkgPlots_ALP "${c[@]}" >> $D/makeBkgPlots.stdout.log 2>&1 ); }
  echo "[mA=$m] panels done"
}
export -f panel; export PC FFIT U
seq 1 30 | xargs -P 4 -I{} bash -c 'panel {}' > $L/closure30_bkgplots.log 2>&1
n=0; for m in $(seq 1 30); do
  for k in bkgplot allPdfs; do f=$(ls $U/${k}_${m}.pdf 2>/dev/null); [ -n "$f" ] && { cp -p $f $AN/figure_ALP/run3/pseudodata_bkg_alias/; n=$((n+1)); }; done
done
echo "    copied $n / 60 panels into the AN"
echo "DONE $(date '+%F %T')"
