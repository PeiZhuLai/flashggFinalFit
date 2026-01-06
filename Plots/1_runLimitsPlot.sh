#!/usr/bin/env bash
# 基本使用 (與舊版相容)
# 質量點現在固定於 makeLimitsPlot.py 內的 massPoints 列表

# 假设乘以某个 xs（若不需要，设为 0.1 pb）
# 1.0 pb = 1000.0 fb
# 0.1 pb = 100.0 fb

# ggF xs for 14 TeV from https://arxiv.org/pdf/2402.09955
# VBF xs for 14 TeV from https://arxiv.org/pdf/2402.09955
# ggF = 51960 fb for 125.38 GeV
# VBF = 4067 fb for 125.38 GeV
# ggF + VBF = 56027 fb for 125.38 GeV

  # --masses 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 \

python3 makeLimitsPlot.py \
  --outdir limitPlots \
  --assume-xs 100 \
  --ggf-xs 56027 \
  --lumi 170.84 \
  --formats png,pdf \
  --masses 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 \
  --only xs,br,wilson

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