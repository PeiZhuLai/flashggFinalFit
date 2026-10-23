#!/bin/bash
#
# Merge the chunks of a high-statistics bias study and run the Gaussian fit once
# over the full toy sample.
#
# Each chunk wrote only its per-split fit files. Here they are renamed to a
# global split index (so no chunk overwrites another), hadded per function, and
# handed to RunBiasStudy.py -p --gaussianFit, which produces the usual
# BiasJson/<mA>_gaussfit.json and the pull plots -- but now over
# NCHUNKS x NTOYS_PER_CHUNK toys instead of 1000.
#
# Usage: merge_bias_highstat.sh MASS [MASS ...]

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi

set -eo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "${script_dir}/../../.." && pwd)}"
CMSSW_TOP="${CMSSW_TOP:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4}"
bias_dir="${BASE_DIR}/Combine/Checks/Bias_nominal"
hs_root="${BIAS_HS_ROOT:-${bias_dir}/bias_outputs_highstat}"   # override for side studies (e.g. envcheck_20261022)
root_datacard_path="${ROOT_DATACARD_PATH:-${BASE_DIR}/Combine/root_t2w}"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 MASS [MASS ...]"
  exit 2
fi

if [[ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]]; then
  source /cvmfs/cms.cern.ch/cmsset_default.sh
fi
cd "${CMSSW_TOP}/src" || exit 1
eval "$(scramv1 runtime -sh)"
if [[ -r "${BASE_DIR}/setup.sh" ]]; then
  source "${BASE_DIR}/setup.sh"
fi

for mA in "$@"; do
  src="${hs_root}/mA_${mA}"
  if [[ ! -d "$src" ]]; then
    echo "[skip] no chunks for mA=${mA} under ${src}"
    continue
  fi

  merge_dir="${hs_root}/mA_${mA}/merged"
  rm -rf "$merge_dir"
  mkdir -p "$merge_dir/BiasFits"
  cd "$merge_dir" || exit 1

  echo "=== mA=${mA}"
  gi=0
  declare -A seen=()
  for f in "${src}"/chunk*/BiasFits/biasStudy_*_split*_fits.root; do
    [[ -f "$f" ]] || continue
    b=$(basename "$f")
    fn="${b#biasStudy_}"; fn="${fn%%_split*}"
    cp -p "$f" "BiasFits/biasStudy_${fn}_split${gi}_fits.root"
    seen[$fn]=1
    gi=$(( gi + 1 ))
  done
  if [[ "$gi" -eq 0 ]]; then
    echo "[skip] mA=${mA}: no split fit files found"
    continue
  fi
  echo "  collected ${gi} split files over functions: ${!seen[*]}"

  min_toys=""
  for fn in "${!seen[@]}"; do
    rm -f "BiasFits/biasStudy_${fn}_fits.root"
    hadd -f "BiasFits/biasStudy_${fn}_fits.root" BiasFits/biasStudy_${fn}_split*_fits.root > /dev/null
    n=$(python3 - "$fn" <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True)
f = ROOT.TFile.Open("BiasFits/biasStudy_%s_fits.root" % sys.argv[1])
t = f.Get("limit") if f and not f.IsZombie() else None
print(t.GetEntries() if t else -1)
PY
)
    # MultiDimFit --algo singles writes three entries per toy (best fit, -1sigma,
    # +1sigma), so the number of usable toys is entries/3.
    toys=$(( n / 3 ))
    echo "    ${fn}: merged tree entries = ${n}  (${toys} toys)"
    if [[ -z "$min_toys" || "$toys" -lt "$min_toys" ]]; then min_toys=$toys; fi
  done

  if [[ -z "$min_toys" || "$min_toys" -le 0 ]]; then
    echo "[ERROR] mA=${mA}: could not determine the merged toy count"
    continue
  fi

  # RunBiasStudy's plotting loop runs over range(opts.nToys) and defaults to
  # 1000. Without passing the merged total it would silently analyse only the
  # first 1000 toys of the merged tree and the uncertainty would not shrink.
  # --split is only used when generating/fitting; the plotting loop ignores it.
  # It still has to satisfy RunBiasStudy's "nToys divisible by split" guard, so
  # set it equal to nToys.
  echo "  running Gaussian fit over ${min_toys} toys"
  card="${root_datacard_path}/${mA}_Datacard_leptons.root"
  sel=()   # truth functions restricted as in the jobs (runjob_bias_highstat.sh SELECT_FUNCTION)
  [[ -n "${SELECT_FUNCTION:-}" ]] && sel=( --selectFunction "${SELECT_FUNCTION}" )
  python3 "${bias_dir}/RunBiasStudy.py" -d "$card" --mA "$mA" -p --gaussianFit \
      -n "$min_toys" --split "$min_toys" "${sel[@]}"

  # Also produce the overlay/r-vs-bias plots that the AN shows in Figure 46.
  # plot_bias.py sizes its loop from the tree itself, so it picks up the full
  # merged sample without needing the toy count passed in.
  BIAS_BASE_OVERRIDE="$merge_dir" python3 "${bias_dir}/plot_bias.py" --mA "$mA"
  if compgen -G "plots_bias/${mA}_r_vs_bias.pdf" >/dev/null; then
    echo "  --> ${merge_dir}/plots_bias/${mA}_r_vs_bias.pdf"
  else
    echo "  [WARN] mA=${mA}: plot_bias.py produced no r_vs_bias plot"
  fi
  echo "  --> ${merge_dir}/BiasJson/${mA}_gaussfit.json"
  cat "BiasJson/${mA}_gaussfit.json" 2>/dev/null || echo "  [WARN] no json produced"
  unset seen
done
