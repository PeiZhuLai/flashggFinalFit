import os
import json
import uproot
import pandas as pd
import numpy as np
import xgboost as xgb
import pickle
from argparse import ArgumentParser
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import logging
import re, ast
from typing import Dict, List, Optional, Tuple  # 新增 Tuple
import awkward as ak

# 新增：統一設定日誌等級與抑制吵雜第三方 logger
def setup_logging(level_str: str = "INFO") -> None:
    level = getattr(logging, level_str.upper(), logging.INFO)
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=level)
    # uproot/fsspec/XRootD 等在 DEBUG 時會非常吵，非 DEBUG 時壓到 WARNING
    noisy_loggers = ['uproot', 'fsspec', 'XRootD', 'asyncio', 'urllib3']
    for name in noisy_loggers:
        logging.getLogger(name).setLevel(logging.WARNING if level > logging.DEBUG else level)

optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

procductions = ["ggh"]
INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal/"
sig_samples = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
years_sig  = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]  # 信号

# 新增：讀取 INPUT_BASE/*.root 用到的設定
INPUT_BASE_TREE_NAME = "test"
UPROOT_STEP = "200 MB"

def get_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description='Apply BDT signal correction to data')
    parser.add_argument('-c', '--config', default='data/training_config_BDT.json', help='Path to the training config file')
    parser.add_argument('-i', '--inputFolder', default='/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal/', help='Path to the input folder')
    parser.add_argument('-o', '--outputFolder', default='/eos/home-p/pelai/HZa/root_MVAcut/sig', help='Path to the output folder')
    # 新增：控制日誌等級，預設 INFO；若需要完整追蹤，指定 --log-level DEBUG
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                        help='Logging level (default: INFO)')
    return parser.parse_args()

def parse_mva_cuts(txt_path: str) -> Dict[int, float]:
    """
    讀取優化後的 MVA cut。
    支援：
      1) JSON 格式：
         {
           "results": [
             {"mA": 5, "MVAcut": 0.975, ...},
             {"mA": 15, "MVAcut": 0.97, ...},
             ...
           ]
         }
         或者 {"5": 0.975, "15": 0.97, ...}
      2) 舊 txt 格式（回退）：以大括號的字典字串表示。
    回傳：{ mA(int): MVAcut(float) }
    """
    cuts: Dict[int, float] = {}

    # 先嘗試 JSON
    try:
        with open(txt_path, "r") as f:
            data = json.load(f)

        if isinstance(data, dict) and "results" in data and isinstance(data["results"], list):
            for item in data["results"]:
                try:
                    ma = int(item["mA"])
                    cut = float(item["MVAcut"])
                    cuts[ma] = cut
                except (KeyError, TypeError, ValueError):
                    continue
        elif isinstance(data, dict):
            # 也支援 {"5": 0.975, "15": 0.97, ...} 或 {"5": {"MVAcut": 0.975}, ...}
            for k, v in data.items():
                try:
                    ma = int(k)
                    cut = float(v if not isinstance(v, dict) else v.get("MVAcut"))
                    cuts[ma] = cut
                except (TypeError, ValueError):
                    continue

        if cuts:
            logging.info(f"Loaded {len(cuts)} MVA cuts from JSON: {sorted(cuts.keys())}")
            return cuts
        else:
            logging.warning("JSON parsed but no valid (mA, MVAcut) pairs found, will try legacy txt parsing.")
    except Exception as e:
        logging.warning(f"Failed to parse JSON '{txt_path}': {e}, will try legacy txt parsing.")

    # 回退：舊的 txt 解析（維持相容）
    try:
        with open(txt_path, "r") as f:
            content = f.read()
        m = re.search(r"\{.*\}", content, flags=re.S)
        if m:
            legacy = {int(k): float(v) for k, v in ast.literal_eval(m.group(0)).items()}
            logging.info(f"Loaded {len(legacy)} MVA cuts from legacy txt format.")
            return legacy
    except Exception as e:
        logging.error(f"Failed to parse legacy txt format: {e}")

    logging.error("No MVA cuts could be loaded; returning empty dict.")
    return cuts

# 新增：從樣本名抓取質量點
def parse_ma_from_name(name: str) -> Optional[int]:
    m = re.search(r"[Mm](\d+)", name)
    return int(m.group(1)) if m else None

# 新增：選擇唯一事件鍵（優先 run, luminosityBlock, event）
def choose_id_columns(keys: List[str]) -> Tuple[str, ...]:
    kset = set(keys)
    if {"run", "luminosityBlock", "event"}.issubset(kset):
        return ("run", "luminosityBlock", "event")
    if {"run", "lumi", "event"}.issubset(kset):
        return ("run", "lumi", "event")
    return ("event",)

# 新增：列出 INPUT_BASE/<sample>/<year>.root（若不存在則掃描目錄）
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

# 新增：由 sample 名稱取得對應的 MVA 分支（era 檔內同時有多個 mA 分支時使用）
def get_mva_branch_for_sample(keys: List[str], sample: str, fallback_syst: str = "nominal") -> Optional[str]:
    """
    e.g. sample='mA_M4' -> branch='MVA_Score_mA_M4'（優先）
    若不存在，回退到既有 get_mva_col(keys, syst)（相容舊格式）
    """
    target = f"MVA_Score_{sample}"
    kset = set(map(str, keys))
    if target in kset:
        return target
    # fallback: 舊檔可能只有 MVA_Score，或 syst 命名不同
    return get_mva_col(list(kset), fallback_syst)

# 新增：建立通過 MVA cut 的 event 鍵集合（用複合鍵）
def build_pass_event_map(samples: List[str], years: List[str], input_base: str, mva_cuts: Dict[int, float]) -> Tuple[Dict[tuple, set], Tuple[str, ...]]:
    pass_map: Dict[tuple, set] = {}
    id_cols: Tuple[str, ...] = ()
    for s in samples:
        ma = parse_ma_from_name(s)
        if ma is None or ma not in mva_cuts:
            continue
        cut = mva_cuts[ma]
        files = list_root_files(input_base, s, years)
        year_sets = {y: set() for y in years}
        for fpath in files:
            try:
                with uproot.open(fpath) as f:
                    if INPUT_BASE_TREE_NAME not in f:
                        continue
                    t = f[INPUT_BASE_TREE_NAME]
                    keys = list(map(str, t.keys()))

                    # 改：用 sample 對應的 MVA branch（例如 MVA_Score_mA_M4）
                    mva_branch = get_mva_branch_for_sample(keys, s, fallback_syst="nominal")
                    if mva_branch is None or "event" not in keys:
                        continue

                    # 在第一個有效檔案決定 id_cols
                    if not id_cols:
                        id_cols = choose_id_columns(keys)
                        logging.info(f"Using ID columns: {id_cols}")
                    if any(col not in keys for col in id_cols):
                        logging.warning(f"Skip {fpath}: missing ID columns {id_cols}")
                        continue

                    # 判斷年份
                    y_match = None
                    for y in years:
                        if fpath.endswith(f"/{y}.root") or y in fpath:
                            y_match = y
                            break
                    if y_match is None:
                        if len(years) == 1:
                            y_match = years[0]
                        else:
                            continue

                    # 迭代讀取：改成讀 mva_branch
                    cols = [mva_branch] + list(id_cols)
                    for arrs in t.iterate(cols, library="ak", step_size=UPROOT_STEP):
                        mask = arrs[mva_branch] > cut
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
        for y in years:
            pass_map[(s, y)] = year_sets[y]
            logging.info(f"MVA pass events for {s}, {y}: {len(pass_map[(s,y)])}")
    if not id_cols:
        id_cols = ("event",)
    return pass_map, id_cols

# 新增：根據 syst 嘗試找對應的 MVA 分支名稱（含多種命名可能）
def get_mva_col(columns: List[str], syst: str) -> Optional[str]:
    cols = set(columns)
    if syst == 'nominal':
        return 'MVA_Score' if 'MVA_Score' in cols else None
    parts = syst.split('_')  # e.g. ['Photon','smear','up']
    if len(parts) >= 2:
        base = parts[0] + ''.join(p.capitalize() for p in parts[1:-1])  # PhotonSmear
        uod = parts[-1].capitalize()  # Up/Down
        candidates = [
            f"MVA_Score_{parts[0]}_{parts[1]}_{uod}",   # MVA_Score_Photon_smear_Up
            f"MVA_Score_{parts[0]}_{parts[1]}{uod}",     # MVA_Score_Photon_smearUp
            f"MVA_Score_{base}{uod}",                    # MVA_Score_PhotonSmearUp
            f"MVA_Score{base}{uod}",                     # MVA_ScorePhotonSmearUp
            f"MVA_Score_{base}_{uod}",                   # MVA_Score_PhotonSmear_Up
        ]
        for c in candidates:
            if c in cols:
                return c
        # 寬鬆比對（忽略大小寫）
        low_parts = [p.lower() for p in parts]
        for c in columns:
            low = c.lower()
            if 'mva_score' in low and all(p in low for p in [low_parts[0], low_parts[-1]]):
                return c
    # 若找不到，回退為無（交由 pass_map 過濾）
    return 'MVA_Score' if 'MVA_Score' in cols else None

def process_files(output_folder, input_folder, pass_map: Dict[tuple, set], id_cols: Tuple[str, ...], mva_cuts: Dict[int, float]):
    """Process input files and write the results to output ROOT files, after MVA_Score filtering."""

    syst_variations = [
        "nominal", 
        "FNUF_up", "FNUF_down", "Material_up", "Material_down",
        "Electron_scale_up", "Electron_scale_down", "Electron_smear_up", "Electron_smear_down", 
        "Muon_scale_up", "Muon_scale_down", "Muon_smear_up", "Muon_smear_down",
        "Photon_scale_up", "Photon_scale_down", "Photon_smear_up", "Photon_smear_down"
    ]

    os.makedirs(output_folder, exist_ok=True)

    # 以 sig_samples 與 years_sig 取代未定義的 procductions/years
    for proc, year in ((p, y) for p in procductions for y in years_sig):
        for mA in sig_samples:
            if not os.path.exists(f"{output_folder}/{mA}/output_{year}.root"):
                os.makedirs(f"{output_folder}/{mA}", exist_ok=True)
            with uproot.recreate(f"{output_folder}/{mA}/output_{year}.root") as outfile:
                logging.info(f"Output file: {output_folder}/{mA}/output_{year}.root")
                for syst in syst_variations:
                    logging.info(f"Processing {year} {mA} {syst}")
                    if syst == 'nominal':
                        input_path = f'{input_folder}/{mA}/{year}.root'
                        syst_suffix = ''
                    else:
                        input_path = f'{input_folder}/{mA}_{syst}/{year}.root'
                        syst_parts = syst.split('_')
                        syst_type = syst_parts[0] + ''.join(part.capitalize() for part in syst_parts[1:-1])
                        syst_uod = syst_parts[-1].capitalize()
                        syst_suffix = f'_{syst_type}{syst_uod}01sigma'
                    if not os.path.exists(input_path):
                        logging.warning(f"Missing input: {input_path}")
                        continue
                    with uproot.open(input_path) as infile:
                        # 新增：檢查與使用指定樹名（預設 'test'）
                        if INPUT_BASE_TREE_NAME not in infile:
                            logging.error(f"Tree '{INPUT_BASE_TREE_NAME}' not found in: {input_path}")
                            continue
                        all_data = infile[INPUT_BASE_TREE_NAME].arrays(library='pd')
                        all_data.columns = [col.replace('up', 'Up').replace('down', 'Down') for col in all_data.columns]
                        if 'dZ' not in all_data.columns:
                            all_data['dZ'] = np.zeros(len(all_data))

                        before = len(all_data)
                        # 先以唯一鍵過濾
                        pass_set = pass_map.get((mA, year), set())
                        if pass_set:
                            if len(id_cols) == 1 and id_cols[0] == 'event':
                                all_data['event'] = all_data['event'].astype(np.int64, copy=False)
                                all_data = all_data[all_data['event'].isin(pass_set)]
                            else:
                                missing = [c for c in id_cols if c not in all_data.columns]
                                if missing:
                                    logging.warning(f"Missing ID columns in {mA} {year} {syst}: {missing}, skip ID merge for this syst.")
                                else:
                                    for c in id_cols:
                                        all_data[c] = all_data[c].astype(np.int64, copy=False)
                                    pass_df = pd.DataFrame(list(pass_set), columns=list(id_cols))
                                    all_data = all_data.merge(pass_df, on=list(id_cols), how='inner')

                        # 改：優先使用對應 sample 的 MVA 分支（MVA_Score_mA_Mx）
                        ma_val = parse_ma_from_name(mA)
                        sample_mva_col = f"MVA_Score_{mA}"
                        mva_col = sample_mva_col if sample_mva_col in all_data.columns else get_mva_col(all_data.columns.tolist(), syst)

                        if ma_val is not None and ma_val in mva_cuts and mva_col is not None and mva_col in all_data.columns:
                            all_data = all_data[all_data[mva_col] > mva_cuts[ma_val]]
                        elif syst != 'nominal':
                            logging.warning(
                                f"No suitable MVA column for {mA} {year} {syst} "
                                f"(tried '{sample_mva_col}'), skip second-stage MVA filter."
                            )

                        logging.info(f"Filtered {mA} {year} {syst}: {before} -> {len(all_data)} rows")

                        for lep in ['ele', 'mu']:
                            data = all_data.query('n_electrons==2' if lep == 'ele' else 'n_muons==2')
                            logging.info(f"Number of events in {syst if syst != 'nominal' else 'nominal'} ({lep}): {len(data)}")
                            data = data.rename(columns={"H_mass": "CMS_hza_mass"})
                            tree_name = f'{proc}_125_Za_{lep}_13p6TeV_cat0{syst_suffix}'
                            outfile[f'DiphotonTree/{tree_name}'] = data

if __name__ == "__main__":
    args = get_args()
    # 新增：在進入主流程前設定日誌等級與抑制第三方刷屏
    setup_logging(args.log_level)
    # 1) 讀取每個 ma 的最佳化 cut
    mva_cuts = parse_mva_cuts(optimized_BDT_Cut)
    # 2) 根據 INPUT_BASE 的 MVA_Score 建立通過事件對照表（回傳唯一鍵欄位）
    pass_map, id_cols = build_pass_event_map(sig_samples, years_sig, INPUT_BASE, mva_cuts)
    # 3) 在主流程中套用唯一鍵與 MVA_Score 雙重篩選
    process_files(args.outputFolder, args.inputFolder, pass_map, id_cols, mva_cuts)
