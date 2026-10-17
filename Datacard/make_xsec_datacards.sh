#!/bin/bash
# Two datacards per m_a (user decision 2026-10-15):
#   <m>_pruned_datacard_leptons.txt       Wilson-coefficient card (makeDatacard.py output, with the theory uncertainties
#                                          on the SM Higgs production cross section: QCDscale_ggH, pdf_Higgs_ggH, alphaS_ggH)
#   <m>_pruned_datacard_leptons_xsec.txt  model-independent cross-section card (signal normalized to 0.1 pb, no SM cross
#                                          section involved): the same card without those three lines
# The xsec card is derived here and checked to differ from the Wilson card by exactly those three lines.
# Usage: make_xsec_datacards.sh [mA ...]   (default 1-30)
set -euo pipefail
DC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons
THEORY="QCDscale_ggH pdf_Higgs_ggH alphaS_ggH"
mAs=( "$@" ); [ ${#mAs[@]} -eq 0 ] && mAs=( $(seq 1 30) )
for mA in "${mAs[@]}"; do
  W=$DC/${mA}_pruned_datacard_leptons.txt; X=$DC/${mA}_pruned_datacard_leptons_xsec.txt
  [ -s "$W" ] || { echo "[ERROR] missing $W"; exit 1; }
  for th in $THEORY; do
    [ "$(awk -v t=$th '$1==t' "$W" | wc -l)" -eq 1 ] || { echo "[ERROR] mA$mA: $th not exactly once in $W"; exit 1; }
  done
  awk '$1!="QCDscale_ggH" && $1!="pdf_Higgs_ggH" && $1!="alphaS_ggH"' "$W" > "$X.tmp"
  [ $(( $(wc -l < "$W") - $(wc -l < "$X.tmp") )) -eq 3 ] || { echo "[ERROR] mA$mA: xsec card does not differ by exactly 3 lines"; exit 1; }
  # the line count in the card header ("kmax") is "*" in these cards, so nothing else changes
  grep -q '^kmax[[:space:]]*\*' "$W" || { echo "[ERROR] mA$mA: kmax is not '*', update it for the xsec card"; exit 1; }
  mv "$X.tmp" "$X"
  echo "[OK] mA$mA: $(basename $X) (Wilson card minus $THEORY)"
done
