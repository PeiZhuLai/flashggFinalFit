#!/usr/bin/env bash
set -euo pipefail

# -------- config (from your snippet) --------
datacardDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons"
publicDir="/afs/cern.ch/user/p/pelai/public/hza_datacard"

mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )

# rsync options
RSYNC_OPTS=( -avh --progress )
# If you want to delete files in destination that no longer exist in source, add:
# RSYNC_OPTS=( -avh --progress --delete )

# -------- helpers --------
die(){ echo "[ERROR] $*" >&2; exit 1; }

# Extract unique .root paths from datacard "shapes" lines.
# Assumes the .root path is the 4th column in the shapes line (typical combine datacard format):
# shapes <proc> <cat> <file.root> <obj>
extract_root_paths() {
  local dc="$1"
  awk '
    $1=="shapes" {
      for (i=1; i<=NF; i++) {
        if ($i ~ /\.root$/) print $i
      }
    }
  ' "$dc" | sort -u
}

# -------- main --------
mkdir -p "$publicDir"

for mA in "${mAs[@]}"; do
  echo "============================================================"
  echo "[INFO] mA = M${mA}"

  # datacard file name pattern (adjust if your naming differs)
  dc_name="${mA}_pruned_datacard_leptons.txt"
  dc_src="${datacardDir}/${dc_name}"

  [[ -f "$dc_src" ]] || die "Datacard not found: $dc_src"

  # destination layout
  dst_base="${publicDir}/mA_M${mA}"
  dst_sig="${dst_base}/root_sig"
  dst_bkg="${dst_base}/root_bkg"
  mkdir -p "$dst_base" "$dst_sig" "$dst_bkg"

  echo "[INFO] Copy datacard -> $dst_base/"
  rsync "${RSYNC_OPTS[@]}" "$dc_src" "${dst_base}/"

  echo "[INFO] Parse root files from datacard..."
  mapfile -t root_files < <(extract_root_paths "$dc_src")

  if [[ "${#root_files[@]}" -eq 0 ]]; then
    die "No .root paths parsed from shapes lines in: $dc_src"
  fi

  echo "[INFO] Found ${#root_files[@]} root files."
  # Classify signal vs bkg by path heuristic; background multipdf usually sits under Background/
  for rf in "${root_files[@]}"; do
    [[ -f "$rf" ]] || die "Referenced ROOT not found: $rf"

    if [[ "$rf" == *"/Background/"* ]] || [[ "$rf" == *"multipdf"* ]]; then
      echo "  [BKG] $rf"
      rsync "${RSYNC_OPTS[@]}" "$rf" "${dst_bkg}/"
    else
      echo "  [SIG] $rf"
      rsync "${RSYNC_OPTS[@]}" "$rf" "${dst_sig}/"
    fi
  done

  echo "[INFO] Done mA=M${mA} -> ${dst_base}"
done

echo "============================================================"
echo "[INFO] All done."
