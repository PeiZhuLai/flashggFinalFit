#!/bin/bash
# Sculpting closure: build pseudo-data datacards (swap bkg+data_obs to the MC-refit
# multipdf, keep the real signal model + systematics) and run combine per mA:
#   - AsymptoticLimits : observed (on MC pseudo-data) + expected(Asimov) quantiles
#   - Significance     : observed local significance (fake-signal probe)
#   - FitDiagnostics   : best-fit r-hat (rMin -3 rMax 3, avoids multipdf boundary)
# MH frozen at 125.38, rMax 2 (same conventions as Plots/1_makeLimits_observed.sh).
# Parallel NPROC<=6. ALPMASS="1 2 3" to restrict.
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure
DCREAL=$FFIT/Datacard/output_Datacard_leptons
NPROC="${NPROC:-6}"

massList=( $(seq 1 30) )
[ -n "${ALPMASS:-}" ] && massList=( ${ALPMASS} )
mkdir -p $PC/datacards $PC/root_t2w $PC/output_limits $PC/output_significance \
         $PC/output_fitdiag $PC/logs

run_one() {
    local mA="$1"
    local realDC=$DCREAL/${mA}_pruned_datacard_leptons.txt
    local psDC=$PC/datacards/${mA}_pseudodata.txt
    local psMPDF=$PC/fit_results_run3_pseudodata/${mA}/CMS-HGG_mva_13p6TeV_multipdf.root
    [[ -f "$realDC" ]] || { echo "[mA=$mA] MISSING real datacard"; return 1; }
    [[ -f "$psMPDF" ]] || { echo "[mA=$mA] MISSING pseudodata multipdf"; return 1; }
    # Swap ONLY the background-model file path (bkg_mass + data_obs lines both point to it);
    # signal shapes and all systematics stay identical to the real analysis.
    sed "s#/Background/ALP_BkgModel_ReReco/fit_results_run3/${mA}/#/pseudodata_closure/fit_results_run3_pseudodata/${mA}/#g" \
        "$realDC" > "$psDC"

    text2workspace.py "$psDC" -o $PC/root_t2w/${mA}_pseudodata.root >/dev/null 2>&1
    local WS=$PC/root_t2w/${mA}_pseudodata.root

    # Limits: observed (MC pseudo-data) + expected bands
    combine -M AsymptoticLimits "$WS" --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n "PS${mA}" \
        > $PC/logs/limit_${mA}.log 2>&1
    mv higgsCombinePS${mA}.AsymptoticLimits.mH125.38.root $PC/output_limits/ 2>/dev/null

    # Observed local significance (fake-signal probe)
    combine -M Significance "$WS" --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n "PS${mA}" \
        > $PC/logs/signif_${mA}.log 2>&1
    mv higgsCombinePS${mA}.Significance.mH125.38.root $PC/output_significance/ 2>/dev/null

    # Best-fit r-hat
    combine -M FitDiagnostics "$WS" --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n "PS${mA}" \
        > $PC/logs/fitdiag_${mA}.log 2>&1
    mv fitDiagnosticsPS${mA}.root $PC/output_fitdiag/ 2>/dev/null

    echo "[mA=$mA] limits+signif+fitdiag done"
}
export -f run_one
export PC DCREAL

cd $PC
for mA in "${massList[@]}"; do
    run_one "$mA" &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 3; done
done
wait
echo "[limits-pseudodata] done."
