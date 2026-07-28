#!/bin/bash
# UNBLIND step 1d: per-mA goodness-of-fit on the UNBLINDED (observed) data.
#   combine -M GoodnessOfFit --algo saturated : observed test statistic + toy distribution.
#   GOF p-value = fraction of toys with t_toy >= t_obs.
# Output -> combine roots: Plots/unblind_GOF/roots/ ; plots via plot_gof_unblind.py
# -> Plots/plot_limits/5_gof/mA{NN}_gof.{pdf,png}.  NTOYS overridable (default 500).
# PARALLEL over mA: NPROC workers (default 4; NPROC=1 -> serial). ALPMASS="1 2 3" to restrict.
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T2W=$FFIT/Combine/root_t2w
OUT=$FFIT/Plots/unblind_GOF/roots
mkdir -p "$OUT" "$OUT/../logs_parallel"
cd "$OUT"

NTOYS="${NTOYS:-500}"
ALPmassList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && ALPmassList=( ${ALPMASS} )
NPROC="${NPROC:-4}"

# Freeze the RooMultiPdf envelope index to its best-fit function. Letting the
# discrete index profile across the GOF toys is not a well-defined GOF (each toy
# would pick its own bkg function). NOTE: freezing does NOT rescue the mA1/mA22
# low p-values -- those are the saturated GOF being unreliable at ~3 events/bin
# over the fine 1 GeV binning (75-85 bins), NOT a bkg-mismodel. The authoritative
# bkg-quality metric is the flashgg per-pdf GOF (mA1 Lau1 p=0.3 etc.), which passes.
PDFIDX="${PDFIDX:-pdfindex_cat0_13p6TeV}"

run_one() {
    local mA="$1"
    local ws=$T2W/${mA}_Datacard_leptons.root
    [ -f "$ws" ] || { echo "[skip] missing $ws (run 1a first)"; return 0; }
    echo "=== GOF mA = ${mA} (saturated, ${NTOYS} toys, pdfindex frozen) ==="
    combine -M MultiDimFit "$ws" -m 125.38 --setParameters MH=125.38 --freezeParameters MH \
        --rMin -3 --rMax 3 --cminDefaultMinimizerStrategy 0 --saveSpecifiedIndex "$PDFIDX" \
        -n ${mA}_gofbest >/dev/null 2>&1
    local bi
    bi=$(python3 - "higgsCombine${mA}_gofbest.MultiDimFit.mH125.38.root" "$PDFIDX" <<'PY'
import ROOT, sys
ROOT.gErrorIgnoreLevel = ROOT.kError
try:
    f = ROOT.TFile(sys.argv[1]); t = f.Get("limit"); t.GetEntry(0)
    print(int(getattr(t, sys.argv[2])))
except Exception:
    print(0)
PY
)
    bi=${bi:-0}
    combine -M GoodnessOfFit "$ws" --algo saturated -m 125.38 \
        --setParameters MH=125.38,${PDFIDX}=${bi} --freezeParameters MH,${PDFIDX} -n ${mA}_obs
    combine -M GoodnessOfFit "$ws" --algo saturated -m 125.38 \
        --setParameters MH=125.38,${PDFIDX}=${bi} --freezeParameters MH,${PDFIDX} -t ${NTOYS} -s 12345 -n ${mA}_toys
}

for mA in "${ALPmassList[@]}"; do
    run_one "$mA" > "$OUT/../logs_parallel/gof_${mA}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
done
wait

echo "[1d] GOF roots -> $OUT (per-mA logs: Plots/unblind_GOF/logs_parallel/) ; now plotting..."
cd $FFIT/Plots
source /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh
conda activate higgs-alp-ana
python3 plot_gof_unblind.py
