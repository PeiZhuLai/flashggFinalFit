#!/usr/bin/env bash
# 基本使用 (與舊版相容)
# 質量點現在固定於 makeLimitsPlot.py 內的 massPoints 列表

# makeLimitsPlot.py uses PyROOT -> need cmsenv (else: ModuleNotFoundError: No module named 'ROOT').
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Plots

# 假设乘以某个 xs（若不需要，设为 0.1 pb）
# 1.0 pb = 1000.0 fb
# 0.1 pb = 100.0 fb

# ggF xs for 14 TeV from https://arxiv.org/pdf/2402.09955
# VBF xs for 14 TeV from https://arxiv.org/pdf/2402.09955
# ggF = 51960 fb for 125.38 GeV
# VBF = 4067 fb for 125.38 GeV
# WH = 1442 fb for 125.38 GeV
# ZH = 936.1 fb for 125.38 GeV
# ttH = 563.4 fb t-chan for 125.38 GeV
# ttH = 3.044 fb s-chan for 125.38 GeV
# ttH = 17.20 fb tWH-chan for 125.38 GeV
# tH = 83.17 fb for 125.38 GeV
# bbH = 632 fb for 125.40 GeV
# 13.6 Tot = 59703.914 fb for 125.38 GeV

# Run2
# pp to H 56000 Zebing

  # --masses 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 \

MASSES_RES="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30"
OBSDIR="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results_observed"

# resolved-only (m_a = 1..30) EXPECTED: XS / BR / Wilson / compare
python3 makeLimitsPlot.py \
  --outdir plot_limits \
  --assume-xs 100 \
  --ggf-xs 59703.914 \
  --lumi 172.13 \
  --formats pdf \
  --masses $MASSES_RES \
  --only xs,br,wilson,compare

# resolved-only (m_a = 1..30) with OBSERVED line (reads output_combine_results_observed)
python3 makeLimitsPlot.py \
  --outdir plot_limits \
  --assume-xs 100 \
  --ggf-xs 59703.914 \
  --lumi 172.13 \
  --formats pdf,png \
  --masses $MASSES_RES \
  --only xs,br \
  --results-dir "$OBSDIR" \
  --tag-suffix observed

# Full-range (m_a = 0.1..30 GeV) XS limit: merged(0.1-0.9) + resolved(1-30) on one log-x canvas.
# Expected (blind) reads Combine/output_combine_results; observed reads output_combine_results_observed.
OBSDIR="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results_observed"

# expected-only
python3 makeLimitsPlot_full.py \
  --outdir plot_limits \
  --assume-xs 100 \
  --lumi 172.13 \
  --formats pdf,png \
  --tag full_0p1_30

# with observed overlay
python3 makeLimitsPlot_full.py \
  --outdir plot_limits \
  --assume-xs 100 \
  --lumi 172.13 \
  --formats pdf,png \
  --resolved-dir "$OBSDIR" \
  --draw-observed \
  --tag full_0p1_30_observed

# python3 makeLimitsPlot.py \
#   --outdir limitPlots \
#   --assume-xs 100 \
#   --ggf-xs 56000 \
#   --lumi 61.89 \
#   --formats png,pdf \
#   --masses 5,15,30 \
#   --only xs,br,wilson

# 只畫 Wilson:
# python3 makeWilsonLimitsPlot.py --outdir limitPlots --assume-xs 100 --ggf-xs 52170 --lumi 62.5 --formats png,pdf --only wilson

# 只畫 BR 與指定質量:
# python3 makeWilsonLimitsPlot.py --outdir limitPlots --masses 10,20,40 --only br

# 只做 BR 線性 y 範例 (仍使用檔內 massPoints)
# python3 makeLimitsPlot.py --br --linear-y