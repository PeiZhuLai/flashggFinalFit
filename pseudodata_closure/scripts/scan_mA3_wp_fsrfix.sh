#!/usr/bin/env bash
# mA3 pseudo-data sculpting closure on the FSR-fix re-production (2026-09-26).
#
# Why: after the FSR fix the mA3 working point stays at 0.988 (R = 1.21) instead of the R ~ 1.1
# cut 0.940, because at 0.940 no background family describes the 1885-event data spectrum
# (GOF ~ 0). A cut with R > 1 can sculpt a fake peak, so the background-only closure is redone
# for mA3 at the adopted 0.988, with 0.940/0.960/0.975 as comparison points -- the same test
# that set mA2 (scan_mA2_wp.sh, 2026-07-26).
#
# Inputs : run3_bdt_scored_fsrfix (select_lib with CLOSURE_FSRFIX=1: 2022 DY+gamma = 10to100),
#          the real mA3 datacard (signal model + systematics at 0.988).
# Outputs: pseudodata_closure/scan_mA3_fsrfix/cut_<tag>/ (pseudo-data ws, envelope, limits, Z)
#          and scan_summary.txt. The 30-point closure in pseudodata/ is NOT touched.
set -u
CUTS="${1:-0.988 0.975 0.960 0.940}"
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
SCAN=$PC/scan_mA3_fsrfix; mkdir -p $SCAN
realDC=$FFIT/Datacard/output_Datacard_leptons/3_pruned_datacard_leptons.txt
export CLOSURE_FSRFIX=1

echo "=== [0] selection check at 0.988: fsrfix scored Data vs the real mA3 MVA-cut data tree"
( cd $PC/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 $CPY - <<'PY'
import numpy as np, uproot, select_lib as sl
from build_pseudodata_ws import _load_data_mass
m, s = _load_data_mass(3)
n_sel = int(((s > 0.988) & (m > 95) & (m < 180)).sum())
t = uproot.open("/eos/home-p/pelai/HZa/root_MVAcut/data/mA_M3/run3.root")["DiphotonTree/Data_13p6TeV"]
br = [b for b in t.keys() if b in ("CMS_hza_mass", "H_mass", "mass")][0]
mm = t[br].array(library="np"); n_real = int(((mm > 95) & (mm < 180)).sum())
print("  base %s\n  scored Data with score>0.988, 95<m<180: %d   real MVA-cut tree (%s in 95-180): %d   %s"
      % (sl.SCORED_BASE, n_sel, br, n_real, "MATCH" if n_sel == n_real else "MISMATCH"))
PY
) 2>&1 | grep -vE 'cling|TList|LC_ALL|Invoking|Results was|With exit'

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $FFIT/.. && eval `scramv1 runtime -sh`
source $FFIT/setup.sh >/dev/null 2>&1
set -u
printf "%-8s %-9s %-9s %-9s %-8s\n" cut obs_lim exp_lim signifZ pseudoN | tee $SCAN/scan_summary.txt
for cut in $CUTS; do
    tag=$(echo $cut | tr -d '.')
    OB=$SCAN/cut_${tag}
    ( cd $PC/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 \
        $CPY build_pseudodata_ws.py --mA 3 --cut $cut --outbase $OB ) > $OB.build.log 2>&1
    WS=$OB/mA_M3/ws/run3.root
    [ -s "$WS" ] || { echo "[FAIL $cut] no pseudo-data ws (see $OB.build.log)"; continue; }
    N=$(grep -oE "pseudo-data N=[0-9]+" $OB.build.log | grep -oE "[0-9]+" | head -1)
    FIT=$OB/fit; mkdir -p $FIT/ft
    ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $WS \
        --saveMultiPdf $FIT/CMS-HGG_mva_13p6TeV_multipdf.root -D $FIT/ft \
        --mass_ALP 3 -c 1 --isFlashgg 0 --isData 0 -f data, \
        --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 30 ) > $OB.ftest.log 2>&1
    [ -s $FIT/CMS-HGG_mva_13p6TeV_multipdf.root ] || { echo "[FAIL $cut] no envelope"; continue; }
    psDC=$OB/datacard.txt
    sed "s#/Background/ALP_BkgModel_ReReco/fit_results_run3/3/#/pseudodata_closure/scan_mA3_fsrfix/cut_${tag}/fit/#g" "$realDC" > "$psDC"
    grep -q "scan_mA3_fsrfix/cut_${tag}/fit/CMS-HGG_mva_13p6TeV_multipdf.root" "$psDC" || { echo "[FAIL $cut] datacard path not replaced"; continue; }
    text2workspace.py "$psDC" -o $OB/ws_t2w.root > $OB.t2w.log 2>&1
    ( cd $OB && combine -M AsymptoticLimits $OB/ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n Scan${tag} > $OB.lim.log 2>&1
      combine -M Significance $OB/ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n Scan${tag} > $OB.sig.log 2>&1 )
    obs=$(grep "Observed Limit" $OB.lim.log | sed -E 's/.*r < //')
    exp=$(grep "Expected 50.0%" $OB.lim.log | sed -E 's/.*r < //')
    Z=$(grep -E "Significance:" $OB.sig.log | sed -E 's/.*Significance: //')
    printf "%-8s %-9s %-9s %-9s %-8s\n" "$cut" "${obs:-NA}" "${exp:-NA}" "${Z:-NA}" "${N:-NA}" | tee -a $SCAN/scan_summary.txt
    echo "   envelope: $(tr '\n' ' ' < $FIT/ft/EnvelopeResults.txt 2>/dev/null | grep -oE 'pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+' | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')"
done
echo "[scan] done -> $SCAN/scan_summary.txt"
