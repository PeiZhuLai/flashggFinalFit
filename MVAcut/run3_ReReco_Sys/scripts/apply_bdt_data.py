import os
import json
import re, ast
import logging
from argparse import ArgumentParser
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import uproot
import awkward as ak
import subprocess

optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/Data"
mAs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
procductions = ['Data']
years = ['2022preEE', '2022postEE', '2023preBPix', '2023postBPix', '2024']

# 讀 INPUT_BASE 的樹名與 iterate step
INPUT_BASE_TREE_NAME = "inclusive"
UPROOT_STEP = "200 MB"

def setup_logging(level_str: str = "INFO") -> None:
    level = getattr(logging, level_str.upper(), logging.INFO)
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=level)
    for name in ['uproot', 'fsspec', 'XRootD', 'asyncio', 'urllib3']:
        logging.getLogger(name).setLevel(logging.WARNING if level > logging.DEBUG else level)

def get_args():
    parser = ArgumentParser(description='Apply BDT data MVA cut and write ROOT files')
    parser.add_argument('-i', '--inputFolder', default='/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/Data', help='Path to the input folder')
    parser.add_argument('-o', '--outputFolder', default='/eos/home-p/pelai/HZa/root_MVAcut/data', help='Path to the output folder')
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Logging level')
    return parser.parse_args()

def parse_mva_cuts(txt_path: str) -> Dict[int, float]:
    cuts: Dict[int, float] = {}
    try:
        with open(txt_path, "r") as f:
            data = json.load(f)
        if isinstance(data, dict) and "results" in data and isinstance(data["results"], list):
            for item in data["results"]:
                try:
                    cuts[int(item["mA"])] = float(item["MVAcut"])
                except (KeyError, TypeError, ValueError):
                    continue
        elif isinstance(data, dict):
            for k, v in data.items():
                try:
                    cuts[int(k)] = float(v if not isinstance(v, dict) else v.get("MVAcut"))
                except (TypeError, ValueError):
                    continue
        if cuts:
            logging.info(f"Loaded {len(cuts)} MVA cuts from JSON.")
            return cuts
        logging.warning("JSON parsed but no valid (mA, MVAcut) pairs found, try legacy txt.")
    except Exception as e:
        logging.warning(f"Failed to parse JSON '{txt_path}': {e}, try legacy txt.")
    try:
        with open(txt_path, "r") as f:
            content = f.read()
        m = re.search(r"\{.*\}", content, flags=re.S)
        if m:
            legacy = {int(k): float(v) for k, v in ast.literal_eval(m.group(0)).items()}
            logging.info(f"Loaded {len(legacy)} MVA cuts from legacy txt.")
            return legacy
    except Exception as e:
        logging.error(f"Failed to parse legacy txt format: {e}")
    logging.error("No MVA cuts could be loaded; returning empty dict.")
    return cuts

def parse_ma_from_name(name: str) -> Optional[int]:
    m = re.search(r"[Mm](\d+)", name)
    return int(m.group(1)) if m else None

# 新增：找距離最近的 mA 與其 cut
def find_nearest_cut(ma: Optional[int], cuts: Dict[int, float]) -> Tuple[Optional[int], Optional[float]]:
    if ma is None or not cuts:
        return None, None
    closest_ma = min(cuts.keys(), key=lambda k: abs(k - ma))
    return closest_ma, cuts.get(closest_ma)

def choose_id_columns(keys: List[str]) -> Tuple[str, ...]:
    kset = set(keys)
    if {"run", "luminosityBlock", "event"}.issubset(kset):
        return ("run", "luminosityBlock", "event")
    if {"run", "lumi", "event"}.issubset(kset):
        return ("run", "lumi", "event")
    return ("event",)

def list_root_files(base: str, sample: str, years: List[str]) -> List[str]:
    files: List[str] = []
    base = base.rstrip("/")
    sample_dir = os.path.join(base, sample)
    if not os.path.exists(sample_dir):
        return files
    for y in years:
        single = os.path.join(sample_dir, f"{y}.root")
        if os.path.isfile(single):
            files.append(single)
    if not files:
        for root, _, fnames in os.walk(sample_dir):
            for fn in fnames:
                if not fn.endswith(".root"):
                    continue
                full = os.path.join(root, fn)
                if years and not any((y in fn) or full.endswith(f"{y}.root") or (os.path.sep + y + os.path.sep in full) for y in years):
                    continue
                files.append(full)
    return sorted(list(dict.fromkeys(files)))

def build_pass_event_map(samples: List[str], years_list: List[str], input_base: str, mva_cuts: Dict[int, float]) -> Tuple[Dict[tuple, set], Tuple[str, ...]]:
    pass_map: Dict[tuple, set] = {}
    id_cols: Tuple[str, ...] = ()
    for s in samples:
        ma = parse_ma_from_name(s)
        # 由最近 mA 取得 cut（若完全沒有，略過）
        nearest_ma, cut = find_nearest_cut(ma, mva_cuts)
        if cut is None:
            logging.warning(f"No available MVA cut for sample {s} (ma={ma}); skip building pass map.")
            continue
        if nearest_ma != ma:
            logging.info(f"Use nearest mA for {s}: target mA={ma}, nearest mA={nearest_ma}, cut={cut}")
        files = list_root_files(input_base, s, years_list)
        year_sets = {y: set() for y in years_list}
        for fpath in files:
            try:
                with uproot.open(fpath) as f:
                    if INPUT_BASE_TREE_NAME not in f:
                        continue
                    t = f[INPUT_BASE_TREE_NAME]
                    keys = list(map(str, t.keys()))
                    if "MVA_Score" not in keys or "event" not in keys:
                        continue
                    if not id_cols:
                        id_cols = choose_id_columns(keys)
                        logging.info(f"Using ID columns: {id_cols}")
                    if any(col not in keys for col in id_cols):
                        logging.warning(f"Skip {fpath}: missing ID columns {id_cols}")
                        continue
                    y_match = None
                    for y in years_list:
                        if fpath.endswith(f"/{y}.root") or y in fpath:
                            y_match = y
                            break
                    if y_match is None:
                        if len(years_list) == 1:
                            y_match = years_list[0]
                        else:
                            continue
                    cols = ["MVA_Score"] + list(id_cols)
                    for arrs in t.iterate(cols, library="ak", step_size=UPROOT_STEP):
                        mask = arrs["MVA_Score"] > cut
                        if len(id_cols) == 1:
                            evs = ak.to_numpy(arrs[id_cols[0]][mask]).astype(np.int64, copy=False)
                            year_sets[y_match].update(int(x) for x in evs.tolist())
                        else:
                            comp = [ak.to_numpy(arrs[c][mask]).astype(np.int64, copy=False) for c in id_cols]
                            for tpl in zip(*comp):
                                year_sets[y_match].add(tuple(int(x) for x in tpl))
            except Exception as e:
                logging.warning(f"Failed to scan {fpath}: {e}")
                continue
        for y in years_list:
            pass_map[(s, y)] = year_sets[y]
            logging.info(f"MVA pass events for {s}, {y}: {len(pass_map[(s,y)])}")
    if not id_cols:
        id_cols = ("event",)
    return pass_map, id_cols

def get_mva_col(columns: List[str], syst: str) -> Optional[str]:
    cols = set(columns)
    if syst == 'nominal':
        return 'MVA_Score' if 'MVA_Score' in cols else None
    return 'MVA_Score' if 'MVA_Score' in cols else None

def filter_columns(df: pd.DataFrame) -> pd.DataFrame:
    # 僅保留必要與系統誤差分支（資料僅 nominal 但保留 Up/Down/central 以防萬一）
    columns_to_keep = []
    for col in df.columns:
        if col in ["CMS_hgg_mass", "CMS_hza_mass", "weight", "dZ"] or col.endswith("Up") or col.endswith("Down") or col.endswith("central"):
            columns_to_keep.append(col)
    for essential_col in ["CMS_hgg_mass", "CMS_hza_mass", "weight", "dZ"]:
        if essential_col in df.columns and essential_col not in columns_to_keep:
            columns_to_keep.append(essential_col)
    return df[columns_to_keep]

def process_files(output_folder: str, input_folder: str, pass_map: Dict[tuple, set], id_cols: Tuple[str, ...], mva_cuts: Dict[int, float]) -> None:
    os.makedirs(output_folder, exist_ok=True)
    syst_variations = ["nominal"]

    samples = [f"mA_M{m}" for m in mAs]
    for year in years:
        for s in samples:
            out_dir = os.path.join(output_folder, s)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"output_{year}.root")

            for syst in syst_variations:
                if syst == "nominal":
                    input_path = f"{input_folder}/{s}/{year}.root"
                else:
                    # 資料不處理 syst；保留接口
                    input_path = f"{input_folder}/{s}/{year}.root"
                if not os.path.exists(input_path):
                    logging.warning(f"Missing input: {input_path}")
                    continue

                with uproot.open(input_path) as infile:
                    if INPUT_BASE_TREE_NAME not in infile:
                        logging.error(f"Tree '{INPUT_BASE_TREE_NAME}' not found in: {input_path}")
                        continue
                    df = infile[INPUT_BASE_TREE_NAME].arrays(library='pd')
                    if 'dZ' not in df.columns:
                        df['dZ'] = np.zeros(len(df), dtype=np.float32)

                    before = len(df)
                    # 1) 先用事件鍵與 pass_map 過濾
                    pass_set = pass_map.get((s, year), set())
                    if pass_set:
                        if len(id_cols) == 1 and id_cols[0] == 'event':
                            df['event'] = df['event'].astype(np.int64, copy=False)
                            df = df[df['event'].isin(pass_set)]
                        else:
                            missing = [c for c in id_cols if c not in df.columns]
                            if missing:
                                logging.warning(f"Missing ID columns in {s} {year}: {missing}, skip ID merge.")
                            else:
                                for c in id_cols:
                                    df[c] = df[c].astype(np.int64, copy=False)
                                pass_df = pd.DataFrame(list(pass_set), columns=list(id_cols))
                                df = df.merge(pass_df, on=list(id_cols), how='inner')

                    # 2) 用最近 mA 的 MVA cut 再次保險過濾
                    ma_val = parse_ma_from_name(s)
                    nearest_ma, thr = find_nearest_cut(ma_val, mva_cuts)
                    mva_col = get_mva_col(df.columns.tolist(), 'nominal')
                    if thr is not None and mva_col is not None:
                        if nearest_ma != ma_val:
                            logging.info(f"Use nearest mA for filtering {s} {year}: target mA={ma_val}, nearest mA={nearest_ma}, cut={thr}")
                        df = df[df[mva_col] > thr]
                    logging.info(f"Filtered {s} {year}: {before} -> {len(df)} rows")

                    # 3) 輸出前清理欄位與命名
                    if 'H_mass' in df.columns:
                        df = df.rename(columns={"H_mass": "CMS_hza_mass"})
                    df = filter_columns(df)

                # 4) 寫檔
                with uproot.recreate(out_path) as outfile:
                    outfile['DiphotonTree/Data_13p6TeV'] = df
                logging.info(f"Wrote: {out_path}")

# 新增：將每個 mA 目錄下的各年度輸出合併為 run3.root
def hadd_outputs(output_folder: str, samples: List[str], years_list: List[str]) -> None:
    for s in samples:
        out_dir = os.path.join(output_folder, s)
        inputs = [os.path.join(out_dir, f"output_{y}.root") for y in years_list
                  if os.path.exists(os.path.join(out_dir, f"output_{y}.root"))]
        if not inputs:
            logging.warning(f"No yearly files found to hadd for {s} under {out_dir}")
            continue
        merged = os.path.join(out_dir, "run3.root")
        if os.path.exists(merged):
            try:
                os.remove(merged)
            except Exception:
                pass
        cmd = ["hadd", "-f", merged] + inputs
        logging.info(f"HADD {s}: {' '.join(inputs)} -> {merged}")
        try:
            subprocess.run(cmd, check=True)
        except Exception as e:
            logging.error(f"hadd failed for {s}: {e}")

if __name__ == "__main__":
    args = get_args()
    setup_logging(args.log_level)

    mva_cuts = parse_mva_cuts(optimized_BDT_Cut)
    samples = [f"mA_M{m}" for m in mAs]
    pass_map, id_cols = build_pass_event_map(samples, years, INPUT_BASE, mva_cuts)
    process_files(args.outputFolder, INPUT_BASE, pass_map, id_cols, mva_cuts)
    # 新增：合併每個 mA 的年度輸出
    hadd_outputs(args.outputFolder, samples, years)
