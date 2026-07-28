#!/bin/bash
# Scoped Trees2WS for mA=2 ONLY (cut moved 0.98 -> 0.982, 2026-06-19).
# Mirrors run_tree2ws.sh but restricts to mA_M2 so the other 29 (frozen) masses are untouched.
set -uo pipefail

baseDir=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
# CMSSW environment (scram runtime) — self-contained so this can run outside fit_bkg.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
source ${baseDir}/setup.sh
export PYTHONPATH=$PYTHONPATH:${baseDir}/tools:${baseDir}/Signal/tools

cd ${baseDir}/Trees2WS

years=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)

### ---- Signal (ele+mu, with systematics, massLow 105 for mA2) ----
sig_path="/eos/home-p/pelai/HZa/root_MVAcut/sig/mA_M2"
for year in "${years[@]}"; do
    for lep in ele mu; do
        echo "[mA2-tree2ws][sig] year=${year} lep=${lep}"
        python3 trees2ws.py --inputConfig config.py \
            --inputTreeFile "${sig_path}/output_${year}.root" \
            --inputMass 125 --productionMode ggh --year "${year}" --lepton "${lep}" \
            --doSystematics --massLow 105
    done
done

### ---- Data (massLow 105 + applyMassCut for mA2) ----
data_path="/eos/home-p/pelai/HZa/root_MVAcut/data/mA_M2"
echo "[mA2-tree2ws][data]"
python3 trees2ws_data.py --inputConfig config.py \
    --inputTreeFile "${data_path}/run3.root" \
    --massLow 105 --massCutRange 105,180 --applyMassCut

echo "[mA2-tree2ws] DONE"
