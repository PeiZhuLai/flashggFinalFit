#!/bin/bash
#
# Datacards for the envelope recheck: the production cross-section card of each
# mass point with only the background multipdf (bkg_mass and data_obs lines)
# pointed at the envelope rebuilt by 1_run_ftest.sh, then text2workspace with
# the production options (Combine/t2w_jobs/<m>_xsec_t2w.sh).
#
# Usage: bash 2_make_cards.sh [--run]   (default: dry-run)

set -eo pipefail
EC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CMSSW_TOP=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4
masses=( ${MASSES:-2 3 4 20 24 25} )

run=0; [[ "${1:-}" == "--run" ]] && run=1
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_TOP}/src" && eval "$(scramv1 runtime -sh)"

for m in "${masses[@]}"; do
  prod_card="${FF}/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons_xsec.txt"
  prod_bkg="${FF}/Background/ALP_BkgModel_ReReco/fit_results_run3/${m}/CMS-HGG_mva_13p6TeV_multipdf.root"
  new_bkg="${EC}/bkg/${m}/CMS-HGG_mva_13p6TeV_multipdf.root"
  card="${EC}/datacards/${m}_pruned_datacard_leptons_xsec.txt"
  echo "[mA=${m}] ${prod_card} -> ${card}; t2w -> ${EC}/t2w/${m}_Datacard_leptons.root"
  [[ $run -eq 1 ]] || continue
  [[ -f "$new_bkg" ]] || { echo "[ERROR] missing ${new_bkg}"; exit 1; }
  sed "s#${prod_bkg}#${new_bkg}#g" "$prod_card" > "$card"
  n=$(grep -c "${new_bkg}" "$card" || true)
  [[ "$n" -eq 2 ]] || { echo "[ERROR] mA=${m}: expected 2 replaced shape lines, got ${n}"; exit 1; }
  # only the two shape lines may differ from the production card
  nd=$(diff "$prod_card" "$card" | grep -c '^>' || true)
  [[ "$nd" -eq 2 ]] || { echo "[ERROR] mA=${m}: ${nd} lines differ from the production card"; exit 1; }
  text2workspace.py "$card" -o "${EC}/t2w/${m}_Datacard_leptons.root" -m 125.38 higgsMassRange=115,135 \
      > "${EC}/logs/t2w_${m}.log" 2>&1
  echo "  ok"
done
[[ $run -eq 1 ]] || echo "(dry-run; pass --run to start)"
