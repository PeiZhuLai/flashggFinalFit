#!/usr/bin/env python3

import argparse
import csv
import glob
import math
import os
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import awkward as ak
import pyarrow.parquet as pq
import uproot


DEFAULT_MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
DEFAULT_YEARS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
DEFAULT_CHANNELS = ["ele", "mu"]


def parse_csv_arg(value: str, cast=str) -> List:
    items = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        items.append(cast(part))
    return items


def safe_ratio(num: Optional[float], den: Optional[float]) -> Optional[float]:
    if num is None or den is None:
        return None
    if den == 0:
        return None
    return num / den


def fmt(value: Optional[float]) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    return f"{value:.10g}"


def find_first_existing(patterns: Sequence[str]) -> Optional[str]:
    for pattern in patterns:
        matches = sorted(glob.glob(pattern, recursive=True))
        for match in matches:
            if os.path.isfile(match):
                return match
    return None


def resolve_parquet_file(base: str, mass: int, year: str) -> Optional[str]:
    mass_tag = f"mA_M{mass}_{year}"
    patterns = [
        os.path.join(base, mass_tag, "merged_nominal.parquet"),
        os.path.join(base, "Sig_MC", mass_tag, "merged_nominal.parquet"),
        os.path.join(base, "parquet_DNA", "Sig_MC", mass_tag, "merged_nominal.parquet"),
        os.path.join(base, "parquet_cutflow_DNA", "Sig_MC", mass_tag, "merged_nominal.parquet"),
        os.path.join(base, "**", mass_tag, "merged_nominal.parquet"),
    ]
    return find_first_existing(patterns)


def resolve_root_file(base: str, mass: int, year: str, filename: str) -> Optional[str]:
    sample_tags = [f"mA_M{mass}", f"ALP_M{mass}"]
    patterns: List[str] = []
    for tag in sample_tags:
        patterns.extend(
            [
                os.path.join(base, tag, filename),
                os.path.join(base, "**", tag, filename),
            ]
        )
    return find_first_existing(patterns)


def resolve_workspace_file(base: str, mass: int, year: str, channel: str, ws_dir: str) -> Optional[str]:
    sample_tags = [f"mA_M{mass}", f"ALP_M{mass}"]
    filename = f"ws_{channel}_{year}.root"
    patterns: List[str] = []
    for tag in sample_tags:
        patterns.extend(
            [
                os.path.join(base, tag, ws_dir, filename),
                os.path.join(base, "**", tag, ws_dir, filename),
            ]
        )
    return find_first_existing(patterns)


def sum_parquet_column(path: str, column: str) -> float:
    table = pq.read_table(path, columns=[column])
    arr = table.column(column).combine_chunks().to_numpy(zero_copy_only=False)
    total = 0.0
    for value in arr:
        total += float(value)
    return total


def sum_root_tree_weight(path: str, tree_name: str, weight_branch: str) -> float:
    total = 0.0
    with uproot.open(path) as root_file:
        if tree_name not in root_file:
            raise KeyError(f"Tree '{tree_name}' not found in {path}")
        tree = root_file[tree_name]
        for chunk in tree.iterate([weight_branch], library="ak", step_size="200 MB"):
            total += float(ak.sum(chunk[weight_branch]))
    return total


def sum_workspace_dataset(path: str, dataset_name: str) -> float:
    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError(
            "PyROOT is required for workspace checks. Run this script inside your CMSSW/ROOT environment."
        ) from exc

    def find_workspace(node, hint: Optional[str] = None):
        if not node:
            return None

        obj = node.Get(hint) if hint and hasattr(node, "Get") else None

        def _find_ws(obj_inner):
            if not obj_inner:
                return None
            if obj_inner.InheritsFrom("RooWorkspace"):
                return obj_inner
            if obj_inner.InheritsFrom("TDirectory"):
                keys = obj_inner.GetListOfKeys()
                if keys:
                    for key in keys:
                        child = key.ReadObj()
                        ws = _find_ws(child)
                        if ws:
                            return ws
            return None

        ws = _find_ws(obj)
        if ws:
            return ws

        for candidate in [hint, "CMS_hza_workspace", "cms_hgg_13TeV", "workspace", "w", "CMS_hgg_workspace"]:
            if not candidate:
                continue
            child = node.Get(candidate) if hasattr(node, "Get") else None
            ws = _find_ws(child)
            if ws:
                return ws

        return _find_ws(node)

    root_file = ROOT.TFile.Open(path, "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Failed to open workspace file: {path}")

    try:
        ws = find_workspace(root_file, "CMS_hza_workspace")
        if ws is None:
            raise KeyError(f"RooWorkspace not found in {path}")

        dataset = ws.data(dataset_name)
        if dataset is None:
            raise KeyError(f"Dataset '{dataset_name}' not found in workspace '{ws.GetName()}' from {path}")
        return float(dataset.sumEntries())
    finally:
        root_file.Close()


def collect_row(
    mass: int,
    year: str,
    channel: str,
    parquet_base: Optional[str],
    p2root_base: Optional[str],
    mvacut_base: Optional[str],
    ws_dirname: str,
    production_mode: str,
    skip_workspace: bool,
) -> Dict[str, Optional[float]]:
    row: Dict[str, Optional[float]] = {
        "mass": mass,
        "year": year,
        "channel": channel,
        "parquet_sumw": None,
        "inclusive_sumw": None,
        "test_sumw": None,
        "test_x2_sumw": None,
        "mvacut_sumw": None,
        "mvacut_x2_sumw": None,
        "ws_sumEntries": None,
        "ws_x2_sumEntries": None,
        "inclusive_over_parquet": None,
        "testx2_over_inclusive": None,
        "mvacut_over_test": None,
        "mvacut_over_testx2": None,
        "mvacutx2_over_inclusive": None,
        "ws_over_mvacut": None,
        "wsx2_over_inclusive": None,
    }

    parquet_file = resolve_parquet_file(parquet_base, mass, year) if parquet_base else None
    p2root_file = resolve_root_file(p2root_base, mass, year, f"{year}.root") if p2root_base else None
    mvacut_file = resolve_root_file(mvacut_base, mass, year, f"output_{year}.root") if mvacut_base else None
    ws_file = resolve_workspace_file(mvacut_base, mass, year, channel, ws_dirname) if mvacut_base else None

    row["parquet_file"] = parquet_file
    row["p2root_file"] = p2root_file
    row["mvacut_file"] = mvacut_file
    row["ws_file"] = ws_file

    if parquet_file:
        row["parquet_sumw"] = sum_parquet_column(parquet_file, "weight_central")

    if p2root_file:
        row["inclusive_sumw"] = sum_root_tree_weight(p2root_file, "inclusive", "weight")
        row["test_sumw"] = sum_root_tree_weight(p2root_file, "test", "weight")
        if row["test_sumw"] is not None:
            row["test_x2_sumw"] = 2.0 * row["test_sumw"]

    nominal_tree = f"DiphotonTree/{production_mode}_125_Za_{channel}_13p6TeV_cat0"
    if mvacut_file:
        row["mvacut_sumw"] = sum_root_tree_weight(mvacut_file, nominal_tree, "weight")
        if row["mvacut_sumw"] is not None:
            row["mvacut_x2_sumw"] = 2.0 * row["mvacut_sumw"]

    dataset_name = f"{production_mode}_125_Za_{channel}_13p6TeV_cat0"
    if ws_file and not skip_workspace:
        row["ws_sumEntries"] = sum_workspace_dataset(ws_file, dataset_name)
        if row["ws_sumEntries"] is not None:
            row["ws_x2_sumEntries"] = 2.0 * row["ws_sumEntries"]

    row["inclusive_over_parquet"] = safe_ratio(row["inclusive_sumw"], row["parquet_sumw"])
    row["testx2_over_inclusive"] = safe_ratio(row["test_x2_sumw"], row["inclusive_sumw"])
    row["mvacut_over_test"] = safe_ratio(row["mvacut_sumw"], row["test_sumw"])
    row["mvacut_over_testx2"] = safe_ratio(row["mvacut_sumw"], row["test_x2_sumw"])
    row["mvacutx2_over_inclusive"] = safe_ratio(row["mvacut_x2_sumw"], row["inclusive_sumw"])
    row["ws_over_mvacut"] = safe_ratio(row["ws_sumEntries"], row["mvacut_sumw"])
    row["wsx2_over_inclusive"] = safe_ratio(row["ws_x2_sumEntries"], row["inclusive_sumw"])
    return row


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit signal yield propagation from parquet to P2Root, MVAcut ROOT, and ws_Tree2WS."
    )
    parser.add_argument(
        "--parquet-base",
        default="/eos/home-p/pelai/HZa",
        help="Base directory containing HiggsDNA merged signal parquet outputs.",
    )
    parser.add_argument(
        "--p2root-base",
        default="/eos/home-p/pelai/HZa/root_P2Root/run3_BDT",
        help="Base directory containing split signal ROOT files with inclusive/train/test trees.",
    )
    parser.add_argument(
        "--mvacut-base",
        default="/eos/home-p/pelai/HZa/root_MVAcut/sig",
        help="Base directory containing MVA-cut ROOT files and ws_Tree2WS subdirectories.",
    )
    parser.add_argument(
        "--workspace-dirname",
        default="ws_Tree2WS",
        help="Workspace subdirectory name under each signal sample directory.",
    )
    parser.add_argument(
        "--masses",
        default=",".join(str(x) for x in DEFAULT_MASSES),
        help="Comma-separated ALP masses to check.",
    )
    parser.add_argument(
        "--years",
        default=",".join(DEFAULT_YEARS),
        help="Comma-separated years to check.",
    )
    parser.add_argument(
        "--channels",
        default=",".join(DEFAULT_CHANNELS),
        help="Comma-separated channels to check.",
    )
    parser.add_argument(
        "--production-mode",
        default="ggh",
        help="Production-mode prefix used in ROOT trees and RooDataSet names.",
    )
    parser.add_argument(
        "--csv",
        default="",
        help="Optional CSV output path. If omitted, prints CSV to stdout.",
    )
    parser.add_argument(
        "--skip-workspace",
        action="store_true",
        help="Skip workspace checks when ROOT/PyROOT is unavailable.",
    )
    return parser


def write_rows(rows: Sequence[Dict[str, Optional[float]]], output_path: str = "") -> None:
    fieldnames = [
        "mass",
        "year",
        "channel",
        "parquet_sumw",
        "inclusive_sumw",
        "test_sumw",
        "test_x2_sumw",
        "mvacut_sumw",
        "mvacut_x2_sumw",
        "ws_sumEntries",
        "ws_x2_sumEntries",
        "inclusive_over_parquet",
        "testx2_over_inclusive",
        "mvacut_over_test",
        "mvacut_over_testx2",
        "mvacutx2_over_inclusive",
        "ws_over_mvacut",
        "wsx2_over_inclusive",
        "parquet_file",
        "p2root_file",
        "mvacut_file",
        "ws_file",
    ]

    stream = open(output_path, "w", newline="") if output_path else sys.stdout
    try:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = {}
            for key in fieldnames:
                value = row.get(key)
                if key.endswith("_file"):
                    out[key] = value or ""
                elif isinstance(value, (int, str)):
                    out[key] = value
                else:
                    out[key] = fmt(value)
            writer.writerow(out)
    finally:
        if output_path:
            stream.close()


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    masses = parse_csv_arg(args.masses, int)
    years = parse_csv_arg(args.years, str)
    channels = parse_csv_arg(args.channels, str)

    rows = []
    errors: List[Tuple[int, str, str, str]] = []

    for mass in masses:
        for year in years:
            for channel in channels:
                try:
                    row = collect_row(
                        mass=mass,
                        year=year,
                        channel=channel,
                        parquet_base=args.parquet_base,
                        p2root_base=args.p2root_base,
                        mvacut_base=args.mvacut_base,
                        ws_dirname=args.workspace_dirname,
                        production_mode=args.production_mode,
                        skip_workspace=args.skip_workspace,
                    )
                    rows.append(row)
                except Exception as exc:
                    errors.append((mass, year, channel, str(exc)))

    write_rows(rows, args.csv)

    if errors:
        sys.stderr.write("\n[WARN] Some rows failed:\n")
        for mass, year, channel, err in errors:
            sys.stderr.write(f"  mA={mass}, year={year}, channel={channel}: {err}\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
