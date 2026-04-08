#!/usr/bin/env python3
"""
Compare the data yield in one plot bin between:
  1) the plot-side source ntuples (run3_mergedBDT/Data/*.root), with the same
     mA-specific MVA cut used by plot_bkgmcScupltingCheck.py
  2) the apply_bdt_data.py outputs (root_MVAcut/data/mA_Mx/output_*.root + run3.root)

This is meant to debug cases where a bin around 114 GeV looks different between
the sculpting-check plot and the final data ROOT files.
"""

import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import awkward as ak
import uproot


H_M_NBINS = 85
H_M_XMIN = 95.0
H_M_XMAX = 180.0
BLIND_LOW = 115.0
BLIND_HIGH = 135.0
YEARS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

DEFAULT_SOURCE_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_mergedBDT/Data"
DEFAULT_OUTPUT_BASE = "/eos/home-p/pelai/HZa/root_MVAcut/data"
DEFAULT_MVA_CUT_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
LOCAL_MVA_CUT_CANDIDATES = [
    REPO_ROOT / "AN-25-172/figure_ALP/run3/HiggsZaAna/Plot/output/MVAcut_points_run3.json",
    REPO_ROOT / "AN-25-172/figure_ALP/run3/optimize_alias/MVAcut_points_run3.json",
    REPO_ROOT / "AN-25-172/figure_ALP/run3/selection_alias/MVAcut_points_run3.json",
]


def _to_int(value) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return int(round(value))
    if isinstance(value, str):
        match = re.search(r"-?\d+", value)
        return int(match.group(0)) if match else None
    return None


def _to_float(value) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        match = re.search(r"-?\d+(?:\.\d+)?", value)
        return float(match.group(0)) if match else None
    return None


def _parse_ma(name: str) -> Optional[int]:
    match = re.search(r"[Mm](\d+)", name)
    return int(match.group(1)) if match else None


def _entry_to_pair(entry: dict) -> Tuple[Optional[int], Optional[float]]:
    mass_keys = ("mA", "ma", "mass", "massA")
    cut_keys = ("MVAcut", "mvaCut", "mva_cut", "cut", "bdtCut", "bdt_cut", "best_MVAcut")

    ma = None
    for key in mass_keys:
        if key in entry:
            ma = _to_int(entry[key])
            break

    if ma is None:
        for key in ("sample", "name", "label", "title"):
            if key in entry and isinstance(entry[key], str):
                ma = _parse_ma(entry[key])
                if ma is not None:
                    break

    cut = None
    for key in cut_keys:
        if key in entry:
            cut = _to_float(entry[key])
            break

    if cut is None:
        for key in ("best", "opt", "result", "payload"):
            if key in entry and isinstance(entry[key], dict):
                for cut_key in cut_keys:
                    if cut_key in entry[key]:
                        cut = _to_float(entry[key][cut_key])
                        break
            if cut is not None:
                break

    return ma, cut


def parse_mva_cuts(path: Path) -> Dict[int, float]:
    with open(path, "r") as handle:
        data = json.load(handle)

    out: Dict[int, float] = {}

    if isinstance(data, list):
        for entry in data:
            if not isinstance(entry, dict):
                continue
            ma, cut = _entry_to_pair(entry)
            if ma is not None and cut is not None:
                out[ma] = cut

    if not out and isinstance(data, dict):
        for key in ("results", "points", "entries"):
            entries = data.get(key)
            if not isinstance(entries, list):
                continue
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                ma, cut = _entry_to_pair(entry)
                if ma is not None and cut is not None:
                    out[ma] = cut
            if out:
                break

    if not out and isinstance(data, dict):
        for key, value in data.items():
            ma_key = _parse_ma(str(key)) or _to_int(key)
            if isinstance(value, dict):
                ma, cut = _entry_to_pair({"mA": ma_key, **value})
            else:
                ma, cut = ma_key, _to_float(value)
            if ma is not None and cut is not None:
                out[ma] = cut

    return out


def complete_mva_cuts(cuts: Dict[int, float], target_masses: Iterable[int]) -> Dict[int, float]:
    if not cuts:
        raise ValueError("No valid MVA cuts were parsed.")

    points = sorted((int(mass), float(cut)) for mass, cut in cuts.items())
    completed = dict(cuts)

    for mass in target_masses:
        if mass in completed:
            continue
        if mass <= points[0][0]:
            completed[mass] = points[0][1]
        elif mass >= points[-1][0]:
            completed[mass] = points[-1][1]
        else:
            for (m1, c1), (m2, c2) in zip(points[:-1], points[1:]):
                if m1 <= mass <= m2:
                    frac = (mass - m1) / float(m2 - m1) if m2 > m1 else 0.0
                    completed[mass] = c1 + frac * (c2 - c1)
                    break
    return completed


def resolve_mva_cut_json(user_path: str) -> Path:
    candidates = [Path(user_path), Path(DEFAULT_MVA_CUT_JSON), *LOCAL_MVA_CUT_CANDIDATES]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    tried = "\n  - ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"Cannot find MVA cut JSON. Tried:\n  - {tried}")


def strip_cycle(name: str) -> str:
    return "/".join(part.split(";")[0] for part in name.split("/"))


def list_trees(root_file) -> List[str]:
    trees = []
    for name, cls in root_file.classnames(recursive=True).items():
        if cls == "TTree":
            trees.append(strip_cycle(name))
    return trees


def pick_output_tree(root_file) -> str:
    candidates = [
        "DiphotonTree/Data_13p6TeV",
        "Data_13p6TeV",
        "inclusive",
    ]
    trees = list_trees(root_file)
    for candidate in candidates:
        if candidate in trees:
            return candidate
    for tree in trees:
        if tree.endswith("Data_13p6TeV"):
            return tree
    if len(trees) == 1:
        return trees[0]
    raise KeyError(f"Could not determine output tree. Available trees: {trees}")


def resolve_source_mva_branch(tree, ma: int) -> str:
    branches = {str(name) for name in tree.keys()}
    mass_tag = f"M{ma}"
    candidates = [
        f"MVA_Score_mA_{mass_tag}",
        f"MVA_Score_{mass_tag}",
        f"MVA_score_{mass_tag}",
        "MVA_Score",
    ]
    for candidate in candidates:
        if candidate in branches:
            return candidate
    raise KeyError(f"No suitable MVA branch found for mA={ma}. Available branches do not contain {candidates}.")


def resolve_output_mass_branch(tree) -> str:
    branches = {str(name) for name in tree.keys()}
    for candidate in ("CMS_hza_mass", "H_m", "H_mass"):
        if candidate in branches:
            return candidate
    raise KeyError(f"No mass branch found in output tree. Available branches: {sorted(branches)}")


def get_plot_bin_edges(mass_value: float) -> Tuple[int, float, float]:
    width = (H_M_XMAX - H_M_XMIN) / H_M_NBINS
    clipped = min(max(mass_value, H_M_XMIN), math.nextafter(H_M_XMAX, H_M_XMIN))
    index = int((clipped - H_M_XMIN) / width)
    low = H_M_XMIN + index * width
    high = H_M_XMAX if index == H_M_NBINS - 1 else low + width
    return index + 1, low, high


def in_root_bin(values, low: float, high: float) -> ak.Array:
    if high >= H_M_XMAX:
        return (values >= low) & (values <= high)
    return (values >= low) & (values < high)


def count_source_bin(
    path: Path,
    ma: int,
    cut: float,
    bin_low: float,
    bin_high: float,
    only_ele: bool,
    only_mu: bool,
    blind: bool,
) -> Tuple[int, int, str]:
    with uproot.open(path) as root_file:
        tree = root_file["inclusive"]
        mva_branch = resolve_source_mva_branch(tree, ma)
        needed = ["H_m", mva_branch]
        if only_ele:
            needed.append("z_mumu")
        if only_mu:
            needed.append("z_ee")
        arrays = tree.arrays(needed, library="ak")

    h_mass = arrays["H_m"]
    score = arrays[mva_branch]
    mask = (h_mass >= H_M_XMIN) & (h_mass <= H_M_XMAX) & (score >= cut)

    if only_ele:
        mask = mask & (abs(arrays["z_mumu"]) != 1)
    if only_mu:
        mask = mask & (abs(arrays["z_ee"]) != 1)
    if blind:
        mask = mask & ~((h_mass >= BLIND_LOW) & (h_mass <= BLIND_HIGH))

    bin_mask = mask & in_root_bin(h_mass, bin_low, bin_high)
    return int(ak.sum(bin_mask)), int(ak.sum(mask)), mva_branch


def count_output_bin(path: Path, bin_low: float, bin_high: float) -> Tuple[int, int, str, str]:
    with uproot.open(path) as root_file:
        tree_name = pick_output_tree(root_file)
        tree = root_file[tree_name]
        mass_branch = resolve_output_mass_branch(tree)
        arrays = tree.arrays([mass_branch], library="ak")

    mass_values = arrays[mass_branch]
    bin_mask = in_root_bin(mass_values, bin_low, bin_high)
    return int(ak.sum(bin_mask)), int(len(mass_values)), tree_name, mass_branch


def fmt_path(path: Path) -> str:
    return str(path)


def fmt_delta(delta: int) -> str:
    return f"{delta:+d}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Check data bin consistency between plot input and apply_bdt_data output.")
    parser.add_argument("--ma", type=int, default=1, help="Target ALP mass, e.g. 1 for mA_M1.")
    parser.add_argument("--mass", type=float, default=114.0, help="Mass value used to locate the plot bin, e.g. 114.")
    parser.add_argument("--source-base", default=DEFAULT_SOURCE_BASE, help="Directory containing source data ROOT files: {year}.root.")
    parser.add_argument("--output-base", default=DEFAULT_OUTPUT_BASE, help="Directory containing apply_bdt_data outputs.")
    parser.add_argument("--mva-cut-json", default=DEFAULT_MVA_CUT_JSON, help="Path to the MVA cut JSON.")
    parser.add_argument("--year", default="run3", choices=["run3", *YEARS], help="Check one year or all years merged as run3.")
    parser.add_argument("--unblind-source", action="store_true", default=False, help="Do not apply the plot-side 115-135 GeV data blind.")
    parser.add_argument("--ele", action="store_true", default=False, help="Apply the same electron-only filter as plot_bkgmcScupltingCheck.py.")
    parser.add_argument("--mu", action="store_true", default=False, help="Apply the same muon-only filter as plot_bkgmcScupltingCheck.py.")
    args = parser.parse_args()

    if args.ele and args.mu:
        raise ValueError("Choose at most one of --ele or --mu.")

    years = YEARS if args.year == "run3" else [args.year]
    mva_json = resolve_mva_cut_json(args.mva_cut_json)
    mva_cuts = complete_mva_cuts(parse_mva_cuts(mva_json), range(1, 31))
    cut = float(mva_cuts[args.ma])

    bin_idx, bin_low, bin_high = get_plot_bin_edges(args.mass)
    sample_name = f"mA_M{args.ma}"
    source_base = Path(args.source_base)
    output_dir = Path(args.output_base) / sample_name

    print(f"[Config] mA={args.ma}, cut={cut:.6f}")
    print(f"[Config] plot bin #{bin_idx}: [{bin_low:.3f}, {bin_high:.3f}{']' if bin_high >= H_M_XMAX else ')'}")
    print(f"[Config] source base: {source_base}")
    print(f"[Config] output dir : {output_dir}")
    print(f"[Config] source blind applied: {not args.unblind_source}")
    if args.ele:
        print("[Config] channel filter: electron-only")
    elif args.mu:
        print("[Config] channel filter: muon-only")
    else:
        print("[Config] channel filter: inclusive")
    print("")

    source_counts: Dict[str, int] = {}
    output_counts: Dict[str, int] = {}
    source_totals: Dict[str, int] = {}
    output_totals: Dict[str, int] = {}
    source_branch = None
    output_tree = None
    output_mass_branch = None

    header = f"{'Year':<14} {'SourceBin':>10} {'OutputBin':>10} {'Delta':>8} {'SourceSel':>10} {'OutputRows':>10}"
    print(header)
    print("-" * len(header))

    for year in years:
        source_file = source_base / f"{year}.root"
        output_file = output_dir / f"output_{year}.root"

        if not source_file.is_file():
            print(f"{year:<14} {'MISSING':>10} {'-':>10} {'-':>8} {'-':>10} {'-':>10}  source={fmt_path(source_file)}")
            continue
        if not output_file.is_file():
            print(f"{year:<14} {'-':>10} {'MISSING':>10} {'-':>8} {'-':>10} {'-':>10}  output={fmt_path(output_file)}")
            continue

        src_bin, src_total, src_branch = count_source_bin(
            source_file,
            args.ma,
            cut,
            bin_low,
            bin_high,
            only_ele=args.ele,
            only_mu=args.mu,
            blind=not args.unblind_source,
        )
        out_bin, out_total, tree_name, mass_branch = count_output_bin(output_file, bin_low, bin_high)

        source_counts[year] = src_bin
        output_counts[year] = out_bin
        source_totals[year] = src_total
        output_totals[year] = out_total
        source_branch = src_branch
        output_tree = tree_name
        output_mass_branch = mass_branch

        print(f"{year:<14} {src_bin:>10d} {out_bin:>10d} {fmt_delta(out_bin - src_bin):>8} {src_total:>10d} {out_total:>10d}")

    if source_counts:
        src_sum = sum(source_counts.values())
        out_sum = sum(output_counts.values())
        src_total_sum = sum(source_totals.values())
        out_total_sum = sum(output_totals.values())
        print("-" * len(header))
        print(f"{'run3(sum)':<14} {src_sum:>10d} {out_sum:>10d} {fmt_delta(out_sum - src_sum):>8} {src_total_sum:>10d} {out_total_sum:>10d}")

    merged_file = output_dir / "run3.root"
    if args.year == "run3" and merged_file.is_file():
        merged_bin, merged_total, merged_tree, merged_mass_branch = count_output_bin(merged_file, bin_low, bin_high)
        src_sum = sum(source_counts.values()) if source_counts else 0
        print(f"{'run3(root)':<14} {src_sum:>10d} {merged_bin:>10d} {fmt_delta(merged_bin - src_sum):>8} {'-':>10} {merged_total:>10d}")
        output_tree = merged_tree
        output_mass_branch = merged_mass_branch
    elif args.year == "run3":
        print(f"{'run3(root)':<14} {'-':>10} {'MISSING':>10} {'-':>8} {'-':>10} {'-':>10}  output={fmt_path(merged_file)}")

    print("")
    if source_branch:
        print(f"[Source] MVA branch used: {source_branch}")
    if output_tree and output_mass_branch:
        print(f"[Output] tree used: {output_tree}")
        print(f"[Output] mass branch used: {output_mass_branch}")

    if args.year == "run3" and source_counts:
        src_sum = sum(source_counts.values())
        out_sum = sum(output_counts.values())
        if merged_file.is_file():
            merged_bin, _, _, _ = count_output_bin(merged_file, bin_low, bin_high)
            if out_sum == src_sum and merged_bin != src_sum:
                print("[Hint] Year-by-year output matches the source, but run3.root does not. Rebuild the hadd output.")
            elif out_sum != src_sum:
                print("[Hint] Year-by-year output already differs from the source. Check which input folder apply_bdt_data.py used and whether the output files are stale.")


if __name__ == "__main__":
    main()
