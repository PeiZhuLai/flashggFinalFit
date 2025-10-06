#!/usr/bin/env bash
# 基本使用 (與舊版相容)
# 質量點現在固定於 makeLimitsPlot.py 內的 massPoints 列表

# 假设乘以某个 xs（若不需要，设为 0.1 pb）
# 1.0 pb = 1000.0 fb
# 0.1 pb = 100.0 fb

python3 makeLimitsPlot.py \
  --outdir limitPlots \
  --assume-xs 100 \
  --ggf-xs 51960 \
  --lumi 62.5 \
  --formats png,pdf \
  --masses 5,15,30 \
  --only xs,br,wilson

# 只畫 Wilson:
# python3 makeWilsonLimitsPlot.py --outdir limitPlots --assume-xs 100 --ggf-xs 52170 --lumi 62.5 --formats png,pdf --only wilson

# 只畫 BR 與指定質量:
# python3 makeWilsonLimitsPlot.py --outdir limitPlots --masses 10,20,40 --only br

# 只做 BR 線性 y 範例 (仍使用檔內 massPoints)
# python3 makeLimitsPlot.py --br --linear-y