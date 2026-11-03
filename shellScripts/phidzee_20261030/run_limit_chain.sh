#!/usr/bin/env bash
# =============================================================================
# run_limit_chain.sh — photon ID（自產 Z->ee）SF 更新後：datacard → text2workspace → expected limit → expected limit 圖
# =============================================================================
# 為什麼：signal 模型與效率（run_sig_chain.sh）重做後，datacard 的 signal yield 與 SF nuisance 要更新，
# expected limit 與圖跟著重做。背景模型、bias、closure 不動。
# 刻意不做的：
#   * Datacard/3_rysn_datacard.sh：會同步到 ~/public/hza_datacard（給別人看的位置），等使用者決定；
#   * observed limit 與 observed 疊圖：沒有使用者明確要求不 unblind；舊的 output_combine_results_observed* 不動；
#   * impacts：另外用 condor 跑。
# 舊產物（cp -a 成 *_prephidzee_20261030）：Datacard/output_Datacard_leptons、Datacard/yields_leptons、
#   Combine/root_t2w_{xsec,wilson}、Combine/output_combine_results_{xsec,wilson}、/eos/home-p/pelai/HZa/output_plots/plot_limits/2_runLimitsPlot
# log ：shellScripts/phidzee_20261030/logs/limit_*.log；每步完成寫 logs/limit_<step>.done
# 用法：bash run_limit_chain.sh [--go] [--from STEP]   STEP ∈ datacard t2w limits plots
set -o pipefail
GO=0; FROM=datacard
while [ $# -gt 0 ]; do case "$1" in --go) GO=1;; --from) FROM="$2"; shift;; esac; shift; done
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
D=$FF/shellScripts/phidzee_20261030; L=$D/logs; mkdir -p "$L"
run=0; want() { [ "$1" = "$FROM" ] && run=1; [ "$run" -eq 1 ]; }
cmsenv_() { source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1; cd "$FF/.." && eval "$(scramv1 runtime -sh)"; }
step() { echo "[$(date '+%F %T')] === $1"; }
die() { echo "[$(date '+%F %T')] !! $*"; exit 1; }
bk() { [ -e "$1" ] || return 0; [ -e "${1}_prephidzee_20261030" ] || cp -a "$1" "${1}_prephidzee_20261030" || die "backup $1"; }

[ -e "$L/sigfit.done" ] || die "signal chain not finished (no logs/sigfit.done)"
if [ "$GO" -eq 0 ]; then
  echo "DRY-RUN from=$FROM: backups; Datacard 1+2; Combine 2_text2ws + 1_makeLimits; expected limit plots"; exit 0
fi

if want datacard; then
  step datacard
  bk "$FF/Datacard/output_Datacard_leptons"; bk "$FF/Datacard/yields_leptons"
  ( cd "$FF/Datacard" && bash 1_runjob_gen_datacard_makeYields.sh ) > "$L/limit_makeYields.log" 2>&1 || die "makeYields failed"
  ( cd "$FF/Datacard" && bash 2_runjob_gen_datacard_makeDatacard.sh ) > "$L/limit_makeDatacard.log" 2>&1 || die "makeDatacard failed"
  for m in $(seq 1 30); do
    for s in "" _xsec; do
      f="$FF/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons${s}.txt"
      [ "$f" -nt "$L/sigfit.done" ] || die "datacard not rewritten: $f"
    done
  done
  touch "$L/limit_datacard.done"
fi

if want t2w; then
  step t2w
  bk "$FF/Combine/root_t2w_xsec"; bk "$FF/Combine/root_t2w_wilson"
  ( cmsenv_; cd "$FF/Combine" && bash 2_text2ws.sh ) > "$L/limit_t2w.log" 2>&1 || die "text2workspace failed"
  touch "$L/limit_t2w.done"
fi

if want limits; then
  step limits
  bk "$FF/Combine/output_combine_results_xsec"; bk "$FF/Combine/output_combine_results_wilson"
  ( cmsenv_; cd "$FF/Combine" && bash 1_makeLimits.sh ) > "$L/limit_combine.log" 2>&1 || die "combine failed"
  for c in xsec wilson; do
    n=$(ls "$FF/Combine/output_combine_results_$c"/higgsCombine*.AsymptoticLimits.mH125.38.root 2>/dev/null | wc -l)
    [ "$n" -eq 30 ] || die "limits $c: expected 30 files, found $n"
  done
  touch "$L/limit_limits.done"
fi

if want plots; then
  step plots
  bk "/eos/home-p/pelai/HZa/output_plots/plot_limits/2_runLimitsPlot"
  ( cmsenv_; cd "$FF/Plots"
    python3 makeLimitsPlot.py --outdir /eos/home-p/pelai/HZa/output_plots/plot_limits/2_runLimitsPlot --assume-xs 100 --higgs-xs 59703.914 --lumi 172.13 \
      --formats pdf --masses 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 --only xs,br,wilson,compare || exit 1
    python3 makeLimitsPlot_full.py --outdir /eos/home-p/pelai/HZa/output_plots/plot_limits/2_runLimitsPlot --assume-xs 100 --lumi 172.13 --formats pdf,png --tag full_0p1_30 || exit 2
  ) > "$L/limit_plots.log" 2>&1 || die "limit plots failed"
  touch "$L/limit_plots.done"
fi
echo "[$(date '+%F %T')] LIMIT CHAIN DONE"
