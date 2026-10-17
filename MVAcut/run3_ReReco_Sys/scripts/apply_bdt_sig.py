import os
import sys
import json
import glob
import hashlib
import math
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

# shared with the scorer (HiggsZaAna/Parquet2Rootfile): local write + xrdcp, and which models produced the scores
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile")
import eos_safe_write as esw
import score_provenance as sp

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
INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/"
sig_samples = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
years_sig  = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]  # 信号

# 新增：讀取 INPUT_BASE/*.root 用到的設定
INPUT_BASE_TREE_NAME = "test"
UPROOT_STEP = "200 MB"
DEFAULT_REWEIGHT_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json"
MVA_REWEIGHT_SYST_NAME = "mva_reweight"
MVA_REWEIGHT_LNN = None  # set from --mva-reweight-lnn-json (S1 scheme 2)
# S1 scheme 2 is the analysis default (user decision 2026-09-27). Until 2026-10-01 it was opt-in through
# HZA_MVA_REWEIGHT_LNN_JSON, and a chain that re-ran one mass point without the variable silently wrote
# the old per-event variation into the datacard (mA3, 2026-10-01). Now the JSON is the default, a missing
# file is an error, and every entry must carry the working point it was derived at ("mvacut"), which has
# to match the working point of this run. --legacy-mva-reweight restores the old per-event variation.
DEFAULT_MVA_REWEIGHT_LNN_JSON = ("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/"
                                 "s1_altsideband_20260927/results/s1_lnn_scheme2_sigrw_dyveto.json")
# [PZ 2026-10-02] Signal reweighting (user decision): the signal is reweighted with the same nominal
# sideband (data/MC) reweight as the background simulation, so that every simulated sample is corrected
# to the data shape. The training implementation (HZaMVA/scripts/sideband_reweight.py) is used, with the
# mass-hypothesis variable at the TRUE ALP mass, param = (ALP_m - m_a)/H_m, as in the BDT training and in
# eval_alt_reweight.py. The reweight is normalized per (m_a, era, channel) so that the preselected yield of
# the nominal test sample is unchanged: it changes the BDT efficiency and the shapes, not the cross section
# normalization. The same per-channel normalization is used for the systematic-variation trees. Absolute
# weight columns are multiplied; scale-factor columns (weight_*_sf_*) are ratios and are left alone.
# The CMS_hza_mva_reweight lnN must then be derived relative to the nominal reweight
# (make_s1_lnn_json.py --reference nominal); a JSON with another reference is refused.
HZAMVA_SCRIPTS = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts"
SIGNAL_RW = None          # SidebandReweighter, None with --no-signal-reweight
SIGNAL_RW_NORM = {}       # (mA, year, lep) -> sum(w) / sum(w R) over the preselected nominal test events
ABS_WEIGHT_COLS = ("weight", "weight_central", "factor", "weight_central_no_lumi")
REWEIGHT_UNCERTAINTY_FRACTION = 1.0

REWEIGHT_VAR_ALIASES = {
    "H_m": ("H_m", "H_mass", "CMS_hza_mass"),
    "ALP_m": ("ALP_m", "ALP_mass"),
    "pho1ECALIso": ("pho1ECALIso", "pho1PIso_noCorr", "ALP_lead_photon_ecalPFClusterIso"),
    "pho2ECALIso": ("pho2ECALIso", "pho2PIso_noCorr", "ALP_sublead_photon_ecalPFClusterIso"),
}

def get_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description='Apply BDT signal correction to data')
    parser.add_argument('-c', '--config', default='data/training_config_BDT.json', help='Path to the training config file')
    parser.add_argument('-i', '--inputFolder', default='/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/', help='Path to the input folder')
    parser.add_argument('-o', '--outputFolder', default='/eos/home-p/pelai/HZa/root_MVAcut/sig', help='Path to the output folder')
    parser.add_argument('--samples', default=",".join(sig_samples),
                        help='Comma-separated signal samples to process, e.g. mA_M1,mA_M2')
    parser.add_argument('--years', default=",".join(years_sig),
                        help='Comma-separated years to process')
    parser.add_argument('--reweight-json', default=DEFAULT_REWEIGHT_JSON,
                        help='Path to sideband reweight JSON, or a directory containing sideband_run3_iterative.json')
    parser.add_argument('--disable-reweight-uncertainty', action='store_true',
                        help='Do not add MVA reweight uncertainty branches to signal output trees')
    # S1 scheme 2 (user decision 2026-09-27): the reweight nuisance is the envelope of the signal
    # yield after the WP under {nominal, z_m&hsb, z_m_narrow&hsb} reweights relative to no
    # reweight, one asymmetric lnN per mA, correlated across eras and channels. When this JSON
    # is given, weight_mva_reweight_Up/Down become those per-mA constants (central 1), so the
    # datacard a_w factory writes kappa_up/kappa_down directly.
    parser.add_argument('--mva-reweight-lnn-json',
                        default=os.environ.get('HZA_MVA_REWEIGHT_LNN_JSON', DEFAULT_MVA_REWEIGHT_LNN_JSON),
                        help='JSON {"values": {"<mA>": {"up": k_up, "down": k_down, "mvacut": cut}}} replacing the per-event reweight variation (default: S1 scheme 2)')
    parser.add_argument('--no-signal-reweight', action='store_true',
                        help='Do NOT apply the nominal sideband reweight to the signal (pre-2026-10-02 behavior; '
                             'needs an lnN JSON with reference "none")')
    parser.add_argument('--legacy-mva-reweight', action='store_true',
                        help='use the old per-event reweight variation instead of the S1 scheme-2 lnN')
    # 新增：控制日誌等級，預設 INFO；若需要完整追蹤，指定 --log-level DEBUG
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                        help='Logging level (default: INFO)')
    return parser.parse_args()

def parse_csv_arg(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]

def resolve_reweight_json(path: str) -> Optional[str]:
    if not path:
        return None
    if os.path.isdir(path):
        preferred = os.path.join(path, "sideband_run3_iterative.json")
        if os.path.isfile(preferred):
            return preferred
        matches = sorted(glob.glob(os.path.join(path, "*.json")))
        return matches[0] if matches else None
    return path if os.path.isfile(path) else None

def _first_present_column(frame: pd.DataFrame, logical_name: str) -> Optional[str]:
    for candidate in REWEIGHT_VAR_ALIASES.get(logical_name, (logical_name,)):
        if candidate in frame.columns:
            return candidate
    return None

def _finite_float(value, default=np.nan):
    try:
        value = float(value)
    except Exception:
        return default
    return value if math.isfinite(value) else default

def _lookup_reweight(values, edges, factors):
    values = np.asarray(values, dtype=float)
    edges = np.asarray(edges, dtype=float)
    factors = np.asarray(factors, dtype=float)
    if len(edges) != len(factors) + 1:
        raise ValueError("Invalid reweight binning: len(edges) must equal len(factors) + 1")
    idx = np.searchsorted(edges, values, side="right") - 1
    idx = np.clip(idx, 0, len(factors) - 1)
    out = factors[idx]
    return np.where(np.isfinite(values) & np.isfinite(out), out, 1.0)

def _mass_from_event(events, seed, masses):
    events = np.asarray(events, dtype=np.int64)
    masses = np.asarray(masses, dtype=float)
    selected = np.empty(events.shape[0], dtype=float)
    for idx, event in enumerate(events):
        payload = f"{seed}:{int(event)}".encode("utf-8", errors="ignore")
        digest = hashlib.blake2b(payload, digest_size=8).digest()
        selected[idx] = masses[int.from_bytes(digest, "little") % len(masses)]
    return selected

def _mass_from_row_order(n_rows, seed, masses):
    rng = np.random.default_rng(seed)
    return rng.choice(np.asarray(masses, dtype=float), size=n_rows, replace=True)

class SidebandReweightVariation:
    """Evaluate configurable variations of the training reweight corrections."""

    def __init__(self, payload: Dict, source_path: str):
        self.payload = payload
        self.source_path = source_path
        settings = payload.get("settings", {})
        param = payload.get("param", {})
        self.reweight_vars = list(settings.get("reweight_vars", []))
        self.seed = int(settings.get("seed", 12345))
        self.mass_hypotheses = list(param.get("mass_hypotheses", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0]))
        self.param_method = str(param.get("method", "event_hash"))
        self.iterations = list(payload.get("iterations", []))

    @classmethod
    def from_json(cls, path: str):
        with open(path, "r") as handle:
            return cls(json.load(handle), source_path=path)

    def _param_values(self, frame: pd.DataFrame):
        h_col = _first_present_column(frame, "H_m")
        alp_col = _first_present_column(frame, "ALP_m")
        if h_col is None or alp_col is None:
            raise KeyError("Need H_m/H_mass/CMS_hza_mass and ALP_m/ALP_mass to build sideband param.")

        h_m = frame[h_col].to_numpy(dtype=float)
        alp_m = frame[alp_col].to_numpy(dtype=float)
        if self.param_method == "event_hash" and "event" in frame.columns:
            mass_hyp = _mass_from_event(frame["event"].fillna(-1).to_numpy(dtype=np.int64), self.seed, self.mass_hypotheses)
        else:
            mass_hyp = _mass_from_row_order(len(frame), self.seed, self.mass_hypotheses)

        with np.errstate(divide="ignore", invalid="ignore"):
            return (alp_m - mass_hyp) / h_m

    def _values(self, frame: pd.DataFrame, var: str):
        if var == "param":
            if "param" in frame.columns:
                return frame["param"].to_numpy(dtype=float)
            return self._param_values(frame)
        column = _first_present_column(frame, var)
        if column is None:
            raise KeyError(f"Missing sideband reweight variable '{var}'.")
        return frame[column].to_numpy(dtype=float)

    def _step_weight(self, frame: pd.DataFrame, step: Dict):
        factors = _lookup_reweight(self._values(frame, step["var"]), step["edges"], step["clipped_factors"])
        norm = _finite_float(step.get("normalization_scale", 1.0), default=1.0)
        return factors * norm

    def weights_for_dataframe(self, frame: pd.DataFrame, varied_var: Optional[str] = None, variation_scale: float = 1.0, skipped_vars: Optional[set] = None):
        weights = np.ones(len(frame), dtype=float)
        for iteration in self.iterations:
            for step in iteration.get("steps", []):
                try:
                    step_weight = self._step_weight(frame, step)
                except KeyError:
                    if skipped_vars is not None:
                        skipped_vars.add(step.get("var", "?"))
                        continue
                    raise
                if varied_var is not None and step.get("var") == varied_var:
                    # Vary only the correction relative to unity: 1 -> 1 has no uncertainty.
                    step_weight = 1.0 + variation_scale * (step_weight - 1.0)
                weights *= step_weight
        return np.where(np.isfinite(weights), weights, 1.0)

def load_reweight_variation(path: str, disabled: bool = False) -> Optional[SidebandReweightVariation]:
    if disabled:
        logging.info("MVA reweight uncertainty disabled by command-line option.")
        return None
    resolved = resolve_reweight_json(path)
    if resolved is None:
        logging.warning(f"No reweight JSON found from '{path}'. MVA reweight uncertainty branches will be unity.")
        return None
    try:
        evaluator = SidebandReweightVariation.from_json(resolved)
        logging.info(f"Loaded MVA reweight JSON for uncertainty: {resolved}")
        return evaluator
    except Exception as e:
        logging.warning(f"Failed to load reweight JSON '{resolved}': {e}. MVA reweight uncertainty branches will be unity.")
        return None

def choose_yield_weight_column(frame: pd.DataFrame) -> Optional[str]:
    for col in ("weight", "weight_central", "factor"):
        if col in frame.columns:
            return col
    return None

def _weighted_ratio(ratio, base_weight):
    ratio = np.asarray(ratio, dtype=float)
    base_weight = np.asarray(base_weight, dtype=float)
    valid = np.isfinite(ratio) & np.isfinite(base_weight)
    if not np.any(valid):
        return 1.0
    denom = np.sum(base_weight[valid])
    if denom == 0.0:
        return 1.0
    out = float(np.sum(base_weight[valid] * ratio[valid]) / denom)
    return out if math.isfinite(out) and out > 0.0 else 1.0

def add_mva_reweight_uncertainty_columns(frame: pd.DataFrame, evaluator: Optional[SidebandReweightVariation], label: str) -> pd.DataFrame:
    central = f"weight_{MVA_REWEIGHT_SYST_NAME}_central"
    up = f"weight_{MVA_REWEIGHT_SYST_NAME}_Up"
    down = f"weight_{MVA_REWEIGHT_SYST_NAME}_Down"

    if len(frame) == 0 or evaluator is None:
        frame[central] = np.ones(len(frame), dtype=np.float32)
        frame[up] = np.ones(len(frame), dtype=np.float32)
        frame[down] = np.ones(len(frame), dtype=np.float32)
        return frame

    weight_col = choose_yield_weight_column(frame)
    base_weight = frame[weight_col].to_numpy(dtype=float) if weight_col else np.ones(len(frame), dtype=float)
    skipped = set()
    try:
        nominal = evaluator.weights_for_dataframe(frame, skipped_vars=skipped)
    except Exception as e:
        logging.warning(f"{label}: cannot evaluate nominal reweight factors ({e}); set MVA reweight uncertainty to unity.")
        frame[central] = np.ones(len(frame), dtype=np.float32)
        frame[up] = np.ones(len(frame), dtype=np.float32)
        frame[down] = np.ones(len(frame), dtype=np.float32)
        return frame

    # Sample-mean of the nominal reweight factor on THIS signal sample.
    # Rescaling each varied weight to preserve sig_mean_nominal removes the
    # global yield-shift component, so the envelope captures only the
    # shape-only uncertainty of the sideband-data/MC reweight.
    sig_mean_nominal = _weighted_ratio(nominal, base_weight)

    ratios = [1.0]
    with np.errstate(divide="ignore", invalid="ignore"):
        for var in evaluator.reweight_vars:
            if var in skipped:
                continue
            for scale in (1.0 + REWEIGHT_UNCERTAINTY_FRACTION, 1.0 - REWEIGHT_UNCERTAINTY_FRACTION):
                try:
                    varied = evaluator.weights_for_dataframe(frame, varied_var=var, variation_scale=scale, skipped_vars=skipped)
                except KeyError:
                    skipped.add(var)
                    continue
                # Rescale so weighted-mean(varied) == weighted-mean(nominal): isolate shape only.
                sig_mean_varied = _weighted_ratio(varied, base_weight)
                if sig_mean_varied > 0.0:
                    varied = varied * (sig_mean_nominal / sig_mean_varied)
                ratio = np.divide(varied, nominal, out=np.ones_like(varied), where=(nominal != 0.0))
                ratios.append(_weighted_ratio(ratio, base_weight))

    if skipped:
        logging.debug(f"{label}: skipped reweight variables absent from tree: {sorted(skipped)}")

    up_ratio = max(ratios)
    down_ratio = min(ratios)
    frame[central] = np.ones(len(frame), dtype=np.float32)
    frame[up] = np.full(len(frame), up_ratio, dtype=np.float32)
    frame[down] = np.full(len(frame), down_ratio, dtype=np.float32)
    logging.info(f"{label}: {MVA_REWEIGHT_SYST_NAME} yield envelope down/up = {down_ratio:.6f}/{up_ratio:.6f}")
    return frame

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

def build_test_key_map(samples: List[str], years: List[str], input_base: str, id_cols: Tuple[str, ...]) -> Dict[tuple, pd.DataFrame]:
    """(sample, year) -> DataFrame of the event keys of the NOMINAL test tree (all events, before
    any BDT cut). Systematic samples are restricted to these keys (see process_files)."""
    out: Dict[tuple, pd.DataFrame] = {}
    for s in samples:
        for y in years:
            p = os.path.join(input_base, s, f"{y}.root")
            if not os.path.exists(p):
                continue
            with uproot.open(p) as f:
                if INPUT_BASE_TREE_NAME not in f:
                    continue
                df = f[INPUT_BASE_TREE_NAME].arrays(list(id_cols), library='pd')
            for c in id_cols:
                df[c] = df[c].astype(np.int64, copy=False)
            out[(s, y)] = df.drop_duplicates()
            logging.info(f"Nominal test keys for {s}, {y}: {len(out[(s, y)])}")
    return out

def _signal_rw_frame(frame: pd.DataFrame, ma: float) -> pd.DataFrame:
    fr = frame.copy()
    for logical in ("pho1ECALIso", "pho2ECALIso", "H_m", "ALP_m"):
        col = _first_present_column(fr, logical)
        if col is not None and col != logical:
            fr[logical] = fr[col]
    fr["param"] = (fr["ALP_m"].to_numpy(dtype=float) - float(ma)) / fr["H_m"].to_numpy(dtype=float)
    return fr

def signal_rw_factors(frame: pd.DataFrame, ma: float) -> np.ndarray:
    if len(frame) == 0:
        return np.ones(0, dtype=float)
    return np.asarray(SIGNAL_RW.weights_for_dataframe(_signal_rw_frame(frame, ma)), dtype=float)

def compute_signal_rw_norm(input_folder: str, mA: str, year: str, ma: float):
    """Per-channel factor that preserves the preselected nominal test-sample yield."""
    path = f"{input_folder}/{mA}/{year}.root"
    with uproot.open(path) as f:
        frame = f[INPUT_BASE_TREE_NAME].arrays(library="pd")
    w = frame["weight"].to_numpy(dtype=float)
    r = signal_rw_factors(frame, ma)
    for lep, col in (("ele", "n_electrons"), ("mu", "n_muons")):
        sel = frame[col].to_numpy() == 2
        den = float(np.sum(w[sel] * r[sel]))
        norm = float(np.sum(w[sel])) / den if den > 0 else 1.0
        SIGNAL_RW_NORM[(mA, year, lep)] = norm
        logging.info(f"signal reweight {mA} {year} {lep}: <R> = {den / max(float(np.sum(w[sel])), 1e-30):.5f}, "
                     f"normalization {norm:.5f} ({int(sel.sum())} preselected test events)")

SYST_VARIATIONS = [
    "nominal",
    "FNUF_up", "FNUF_down", "Material_up", "Material_down",
    "Electron_scale_up", "Electron_scale_down", "Electron_smear_up", "Electron_smear_down",
    "Muon_scale_up", "Muon_scale_down", "Muon_smear_up", "Muon_smear_down",
    "Photon_scale_up", "Photon_scale_down", "Photon_smear_up", "Photon_smear_down"
]

def process_files(output_folder, input_folder, samples: List[str], years: List[str], pass_map: Dict[tuple, set], id_cols: Tuple[str, ...], mva_cuts: Dict[int, float], reweight_variation: Optional[SidebandReweightVariation]):
    """Process input files and write the results to output ROOT files, after MVA_Score filtering."""

    syst_variations = SYST_VARIATIONS

    os.makedirs(output_folder, exist_ok=True)

    # Process only the selected sample/year slice so multiple shell jobs can run safely.
    for proc, year in ((p, y) for p in procductions for y in years):
        for mA in samples:
            if not os.path.exists(f"{output_folder}/{mA}/output_{year}.root"):
                os.makedirs(f"{output_folder}/{mA}", exist_ok=True)
            if SIGNAL_RW is not None:
                compute_signal_rw_norm(input_folder, mA, year, parse_ma_from_name(mA))
            # written to a local file, verified (every branch, every entry), then xrdcp'd to EOS through .partial and
            # renamed: writing through the EOS FUSE mount damaged files, and a job killed mid-way must not leave a
            # readable output that silently lacks the later systematic trees
            out_final = f"{output_folder}/{mA}/output_{year}.root"
            out_tmp = esw.stage_path(out_final)
            expect = {}
            with uproot.recreate(out_tmp) as outfile:
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
                    # [PZ 2026-09-27] Systematic samples: the train/validation/test split in p2root
                    # (Parque2Root_BDT.py) is row-index % 10 per file, so once a variation moves one
                    # event across the preselection every later row lands in a different bucket and
                    # the syst 'test' tree is a different random 30% than the nominal one. Joining
                    # it with the nominal pass set kept only the overlap (median 48% of the nominal
                    # yield, down to 3%), which faked 5%-capped rate terms in calcPhotonSyst. Now:
                    # read the syst 'inclusive' tree, keep the events whose key is in the NOMINAL
                    # test tree (same split as the nominal, no training events), and apply the
                    # syst's own BDT score cut below (no join with the nominal pass set, so events
                    # can migrate both ways). The nominal branch is unchanged.
                    tree_name_in = INPUT_BASE_TREE_NAME if syst == 'nominal' else "inclusive"
                    with uproot.open(input_path) as infile:
                        if tree_name_in not in infile:
                            logging.error(f"Tree '{tree_name_in}' not found in: {input_path}")
                            continue
                        all_data = infile[tree_name_in].arrays(library='pd')
                        all_data.columns = [col.replace('up', 'Up').replace('down', 'Down') for col in all_data.columns]
                        if 'dZ' not in all_data.columns:
                            all_data['dZ'] = np.zeros(len(all_data))

                        before = len(all_data)
                        # 先以唯一鍵過濾
                        pass_set = pass_map.get((mA, year), set())
                        if syst != 'nominal':
                            test_keys = test_key_map.get((mA, year))
                            missing = [c for c in id_cols if c not in all_data.columns]
                            if test_keys is None or missing:
                                logging.error(f"No nominal test keys / ID columns for {mA} {year} {syst}; skip")
                                continue
                            for c in id_cols:
                                all_data[c] = all_data[c].astype(np.int64, copy=False)
                            all_data = all_data.merge(test_keys, on=list(id_cols), how='inner')
                            pass_set = set()   # no join with the nominal pass set for systematics
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
                        if SIGNAL_RW is not None:
                            all_data["_sigrw_raw"] = signal_rw_factors(all_data, ma_val)

                        for lep in ['ele', 'mu']:
                            data = all_data.query('n_electrons==2' if lep == 'ele' else 'n_muons==2').copy()
                            logging.info(f"Number of events in {syst if syst != 'nominal' else 'nominal'} ({lep}): {len(data)}")
                            data = add_mva_reweight_uncertainty_columns(
                                data,
                                reweight_variation,
                                label=f"{mA} {year} {syst} {lep}",
                            )
                            if MVA_REWEIGHT_LNN is not None:
                                k = MVA_REWEIGHT_LNN[ma_val]
                                data[f"weight_{MVA_REWEIGHT_SYST_NAME}_central"] = np.ones(len(data), dtype=np.float32)
                                data[f"weight_{MVA_REWEIGHT_SYST_NAME}_Up"] = np.full(len(data), k["up"], dtype=np.float32)
                                data[f"weight_{MVA_REWEIGHT_SYST_NAME}_Down"] = np.full(len(data), k["down"], dtype=np.float32)
                            if SIGNAL_RW is not None:
                                rw = data.pop("_sigrw_raw").to_numpy(dtype=float) * SIGNAL_RW_NORM[(mA, year, lep)]
                                for col in ABS_WEIGHT_COLS:
                                    if col in data.columns:
                                        data[col] = data[col].to_numpy(dtype=float) * rw
                                data["weight_sideband_rwgt"] = rw.astype(np.float32)
                            data = data.rename(columns={"H_mass": "CMS_hza_mass"})
                            tree_name = f'{proc}_125_Za_{lep}_13p6TeV_cat0{syst_suffix}'
                            outfile[f'DiphotonTree/{tree_name}'] = data
                            esw.record(expect, f'DiphotonTree/{tree_name}', data)
            esw.verify_local(out_tmp, expect)
            esw.publish(out_tmp, out_final)

if __name__ == "__main__":
    args = get_args()
    # 新增：在進入主流程前設定日誌等級與抑制第三方刷屏
    setup_logging(args.log_level)
    samples = parse_csv_arg(args.samples)
    years = parse_csv_arg(args.years)
    logging.info(f"Signal samples selected: {samples}")
    logging.info(f"Years selected: {years}")
    reweight_variation = load_reweight_variation(args.reweight_json, disabled=args.disable_reweight_uncertainty)
    # 1) 讀取每個 ma 的最佳化 cut
    mva_cuts = parse_mva_cuts(optimized_BDT_Cut)
    if not args.no_signal_reweight:
        sys.path.insert(0, HZAMVA_SCRIPTS)
        from sideband_reweight import SidebandReweighter
        _rw_path = resolve_reweight_json(args.reweight_json)
        if _rw_path is None:
            raise SystemExit(f"Signal reweight JSON not found from '{args.reweight_json}'")
        SIGNAL_RW = SidebandReweighter.from_json(_rw_path)
        logging.info(f"Signal reweighted with {_rw_path} (true-mass param, preselected yield preserved per channel)")
    else:
        logging.warning("--no-signal-reweight: the signal is NOT reweighted")
    MVA_REWEIGHT_LNN = None
    if args.legacy_mva_reweight:
        logging.warning("--legacy-mva-reweight: per-event reweight variation (NOT the S1 scheme-2 lnN)")
    else:
        lnn_path = args.mva_reweight_lnn_json
        if not lnn_path or not os.path.exists(lnn_path):
            raise SystemExit(f"S1 lnN JSON not found: {lnn_path} (use --legacy-mva-reweight for the old variation)")
        _lnn = json.load(open(lnn_path))
        _ref = _lnn.get("reference", "none")
        _want = "none" if args.no_signal_reweight else "nominal"
        if _ref != _want:
            raise SystemExit(f"S1 lnN JSON {lnn_path} is relative to the '{_ref}' reweight, but this run "
                             f"{'does not reweight' if _want == 'none' else 'reweights'} the signal "
                             f"(needs reference '{_want}'; make_s1_lnn_json.py --reference {_want})")
        _vals = _lnn["values"]
        MVA_REWEIGHT_LNN = {int(m): {"up": float(v["up"]), "down": float(v["down"]),
                                     "mvacut": (float(v["mvacut"]) if "mvacut" in v else None)} for m, v in _vals.items()}
        problems = []
        for x in samples:
            m = parse_ma_from_name(x)
            if m not in MVA_REWEIGHT_LNN:
                problems.append(f"mA{m}: no entry"); continue
            c = MVA_REWEIGHT_LNN[m]["mvacut"]
            if c is None:
                problems.append(f"mA{m}: entry has no 'mvacut' (cannot verify it matches the working point)")
            elif m in mva_cuts and abs(c - mva_cuts[m]) > 1e-6:
                problems.append(f"mA{m}: lnN derived at cut {c} but the working point is {mva_cuts[m]}")
        if problems:
            raise SystemExit("S1 lnN JSON %s inconsistent with this run -- rerun eval_alt_reweight.py + "
                             "make_s1_lnn_json.py at the current working points:\n  " % lnn_path + "\n  ".join(problems))
        logging.info(f"MVA reweight nuisance = per-mA lnN from {lnn_path}: {MVA_REWEIGHT_LNN}")
    # the scores must come from the installed models (fingerprint, or recomputed for files written before it):
    # every nominal and systematic input of this slice, and the INPUT_BASE files of the pass map
    to_check = []
    for base in dict.fromkeys([args.inputFolder.rstrip("/"), INPUT_BASE.rstrip("/")]):
        for x in samples:
            for syst in SYST_VARIATIONS:
                for y in years:
                    p = f"{base}/{x}/{y}.root" if syst == "nominal" else f"{base}/{x}_{syst}/{y}.root"
                    if os.path.isfile(p):
                        to_check.append(p)
    for p, how in sp.check_scored_files(to_check).items():
        logging.info(f"Scores of {p}: installed models ({how})")
    # 2) 根據 INPUT_BASE 的 MVA_Score 建立通過事件對照表（回傳唯一鍵欄位）
    pass_map, id_cols = build_pass_event_map(samples, years, INPUT_BASE, mva_cuts)
    # 3) 在主流程中套用唯一鍵與 MVA_Score 雙重篩選
    test_key_map = build_test_key_map(samples, years, args.inputFolder, id_cols)
    try:
        process_files(args.outputFolder, args.inputFolder, samples, years, pass_map, id_cols, mva_cuts, reweight_variation)
    finally:
        esw.cleanup()
