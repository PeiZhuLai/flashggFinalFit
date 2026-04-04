# Hold defs of functions for calculating systematics and adding to dataframe
import os, sys, re, json
import ROOT
from commonTools import *
from commonObjects import *
from .wsUtils import fetchWorkspace

# sd = "systematics dataframe"

# For constant systematics:
def addConstantSyst(sd,_syst,options):
  # -------------- 新增：偵錯旗標 --------------
  debug = _syst.get('debug_constant', False)

  def _constant_value_for_year(value, year, syst_name):
    if isinstance(value, dict):
      if year not in value:
        return '-'
      return value[year]
    return value

  fromJson = False
  if "json" in str(_syst['value']):
    fromJson = True
    try:
      with open(_syst['value'], "r" ) as jsonfile:
        uval = json.load(jsonfile)
      if debug: print(f"[addConstantSyst] Loaded JSON for {_syst['name']}")
    except Exception as e:
      print(f"[addConstantSyst][ERROR] Cannot read json '{_syst['value']}' for {_syst['name']}: {e}")
      return sd

  years = options.years.split(",")
  # -------------- 新增：若需要 dict 但給單值，自動展開 --------------
  if _syst['correlateAcrossYears'] in (0,-1) and not fromJson:
    if not isinstance(_syst['value'], dict):
      _syst['value'] = {y: _syst['value'] for y in years}
      if debug:
        print(f"[addConstantSyst] Auto-expand single value to dict for {_syst['name']}: {_syst['value']}")

  # -------------- correlateAcrossYears == 1 完全相關 --------------
  if _syst['correlateAcrossYears'] == 1:
    sd[_syst['name']] = '-'
    mask_sig = (sd['type']=='sig') & (~sd['cat'].str.contains("NOTAG"))
    if fromJson:
      # 新增：收集未匹配樣本
      unmatched_samples = []
      def _json_apply(row):
        val = getValueFromJson(row,uval,_syst['name'],_syst)
        if val == '-':
          if len(unmatched_samples)<10: unmatched_samples.append(row.get('proc','?'))
        return val
      sd.loc[mask_sig,_syst['name']] = sd[mask_sig].apply(_json_apply, axis=1)
      if debug:
        filled = sd.loc[mask_sig, _syst['name']].ne('-').sum()
        total = mask_sig.sum()
        print(f"[addConstantSyst] {_syst['name']} correlate=1 filled rows={filled}/{total}")
        if unmatched_samples:
          print(f"[addConstantSyst][DEBUG] {_syst['name']} first unmatched proc examples: {unmatched_samples}")
          print(f"[addConstantSyst][DEBUG] Available JSON keys: {list(uval.keys())[:8]}{' ...' if len(uval)>8 else ''}")
    else:
      if isinstance(_syst['value'], dict):
        missing_years = sorted(set(sd.loc[mask_sig, 'year'].unique()) - set(_syst['value'].keys()))
        if missing_years:
          print(f"[addConstantSyst][WARNING] {_syst['name']} missing values for years {missing_years}")
        sd.loc[mask_sig, _syst['name']] = sd.loc[mask_sig].apply(
          lambda row: _constant_value_for_year(_syst['value'], row['year'], _syst['name']),
          axis=1
        )
      else:
        sd.loc[mask_sig, _syst['name']] = _syst['value']
    if debug:
      filled = sd.loc[mask_sig, _syst['name']].ne('-').sum()
      print(f"[addConstantSyst] {_syst['name']} correlate=1 filled rows={filled}")

  # -------------- correlateAcrossYears == -1 部分相關 (單欄位依年份填不同值) --------------
  elif _syst['correlateAcrossYears'] == -1:
    sd[_syst['name']] = '-'
    for year in years:
      if fromJson:
        # JSON 走 getValueFromJson 流程（但這裡語意少見，保留最小行為）
        mask_sig_year = (sd['type']=='sig') & (~sd['cat'].str.contains("NOTAG")) & (sd['year']==year)
        sd.loc[mask_sig_year,_syst['name']] = sd[mask_sig_year].apply(lambda x: getValueFromJson(x,uval,_syst['name'],_syst), axis=1)
      else:
        if year not in _syst['value']:
          print(f"[addConstantSyst][WARNING] {_syst['name']} missing value for year {year}")
          continue
        mask_sig_year = (sd['type']=='sig') & (~sd['cat'].str.contains("NOTAG")) & (sd['year']==year)
        sd.loc[mask_sig_year,_syst['name']] = _syst['value'][year]
    if debug:
      filled = sd[_syst['name']].ne('-').sum()
      print(f"[addConstantSyst] {_syst['name']} correlate=-1 filled rows={filled}")

  # -------------- correlateAcrossYears == 0 不相關（每年獨立欄位） --------------
  else:
    for year in years:
      col = f"{_syst['name']}_{year}"
      sd[col] = '-'
      if fromJson:
        mask_sig_year = (sd['type']=='sig') & (~sd['cat'].str.contains("NOTAG")) & (sd['year']==year)
        # 這裡 JSON 仍假設 getValueFromJson 回傳 list 或 '-'，若需要 per-year JSON key 可自行擴充
        sd.loc[mask_sig_year, col] = sd[mask_sig_year].apply(lambda x: getValueFromJson(x,uval,_syst['name'],_syst), axis=1)
      else:
        if year not in _syst['value']:
          print(f"[addConstantSyst][WARNING] {_syst['name']} missing value for year {year}")
          continue
        mask_sig_year = (sd['type']=='sig') & (~sd['cat'].str.contains("NOTAG")) & (sd['year']==year)
        sd.loc[mask_sig_year, col] = _syst['value'][year]
      if debug:
        filled = sd[col].ne('-').sum()
        print(f"[addConstantSyst] {_syst['name']} year={year} filled rows={filled}")

  return sd

# ---------------- 新增：通用 JSON process key 匹配 ----------------
def _match_json_process_key(proc_name, uncertainties, syst_dict=None):
  """
  回傳 uncertainties 中最佳匹配的 key：
    1. 若 syst_dict 指定 json_key_field，直接用該欄位值（若存在於 uncertainties）
    2. 否則做 longest match：
         exact > startswith > substring
    3. 都找不到則回 None
  """
  if proc_name is None:
    return None
  # 允許使用其他欄位（例如 procOriginal）：
  if syst_dict:
    alt_field = syst_dict.get('json_key_field')
    if alt_field and alt_field in proc_name:
      # (此處 alt_field 是欄位名，不是內容；真正欄位內容需在呼叫前提供，這裡保持簡單)
      pass
  keys = list(uncertainties.keys())
  # exact
  if proc_name in uncertainties:
    return proc_name
  candidates = []
  for k in keys:
    if proc_name.startswith(k):
      candidates.append((len(k), 'startswith', k))
    elif k in proc_name:
      candidates.append((len(k), 'substr', k))
  if candidates:
    candidates.sort(reverse=True)
    return candidates[0][2]
  # 嘗試移除常見尾碼（年份 / decay / channel）後再一次
  stripped = re.sub(r'_(20(16|17|18|19|22|23)\w*)$', '', proc_name)
  if stripped in uncertainties:
    return stripped
  for k in keys:
    if stripped.startswith(k):
      candidates.append((len(k), 'startswith2', k))
    elif k in stripped:
      candidates.append((len(k), 'substr2', k))
  if candidates:
    candidates.sort(reverse=True)
    return candidates[0][2]
  return None

def getValueFromJson(row,uncertainties,sname,syst_dict=None):
  """
  改寫：改用彈性 process key 匹配；若找到：
    * list -> 原樣回傳
    * 單值 -> 包成 list 以符合後續 lnN 多/單參數格式
  """
  # 嘗試多個候選欄位
  proc_fields = ['proc','procOriginal','proc_s0']
  p0 = None
  for f in proc_fields:
    if f in row and isinstance(row[f], str):
      p0 = row[f]
      break
  if p0 is None:
    return '-'
  match_key = _match_json_process_key(p0, uncertainties, syst_dict)
  if match_key is None:
    return '-'
  # 取值
  if match_key not in uncertainties:  # safety
    return '-'
  entry = uncertainties[match_key].get(sname)
  if entry is None:
    return '-'
  if isinstance(entry, list):
    return entry
  else:
    return [entry]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to return type of systematic: to be used by factory functions
# a) Anti-symmetric weight in nominal RooDataSet: "a_w"
# b) Symmetric weight in nominal RooDataSet: "s_w"
# c) Anti-symmetric shifts in RooDataHist: "a_h"
def _find_updown_objects(ws, base):
    """嘗試尋找與 base 對應之 Up/Down 形狀或資料集名稱"""
    candidates = []
    patterns = [
        f"{base}Up", f"{base}Down",
        f"{base}_Up", f"{base}_Down",
        f"{base}up", f"{base}down",
        f"{base}_up", f"{base}_down"
    ]
    for p in patterns:
        obj = ws.obj(p)
        if obj:
            candidates.append(p)
    return candidates

def _extract_dataset_varnames(rdata):
    """從 RooDataSet 第一筆條目解析所有變數名稱 (contentsString)"""
    names = set()
    if rdata and rdata.numEntries()>0:
        try:
            first = rdata.get(0)
            contents = first.contentsString()
            for kv in contents.split(':'):
                if '=' in kv:
                    names.add(kv.split('=')[0])
        except Exception:
            pass
    return names

def _detect_hist_pair(ws, nominalDataName, sname):
    """
    檢查是否存在對應 Up/Down RooDataHist:
      <nominal>_<syst>Up01sigma / <nominal>_<syst>Down01sigma
    允許名稱中大小寫精確匹配；若存在即回 True
    """
    up = f"{nominalDataName}_{sname}Up01sigma"
    down = f"{nominalDataName}_{sname}Down01sigma"
    if ws.data(up) is not None and ws.data(down) is not None:
        return True
    return False

def _weight_name_candidates(sname):
    """
    回傳該 systematic 可能的權重變數名稱 (central, up, down)
    支援兩種常見前綴：
      weight_<syst>_(central|Up|Down)
      <syst>(Up|Down)  / <syst> (central)
    """
    return {
        'central': [
            f"weight_{sname}_central",
            f"weight_{sname}",
            f"{sname}",
        ],
        'up': [
            f"weight_{sname}_Up",
            f"weight_{sname}Up",
            f"weight_{sname}_up",
            f"{sname}Up",
            f"{sname}_Up",
            f"{sname}_up",
        ],
        'down': [
            f"weight_{sname}_Down",
            f"weight_{sname}Down",
            f"weight_{sname}_down",
            f"{sname}Down",
            f"{sname}_Down",
            f"{sname}_down",
        ]
    }

def _detect_weight_triplet(varnames, sname):
    """
    判斷是否擁有 (central + up + down) 權重組；
    只要 central/Up/Down 清單中任一名稱命中各一個即可。
    """
    cands = _weight_name_candidates(sname)
    has_c = any(v in varnames for v in cands['central'])
    has_u = any(v in varnames for v in cands['up'])
    has_d = any(v in varnames for v in cands['down'])
    return has_c and has_u and has_d, has_c, has_u, has_d

def _pick_first_existing_name(cands, contents):
    """
    挑選第一個存在於 contentsString 的變數名稱
    """
    for n in cands:
        if n in contents:
            return n
    return None

def factoryType(d,s):

  # 強制指定型態映射（新增：所有以 weight_<name>_{central,Up,Down} 形式存在的實驗系統誤差）
  forcedTypeMap = {
    'hlt_sf':'a_w',
    'pu_reweight_sf':'a_w',
    'electron_iso_sf_SelectedElectron':'a_w',
    'electron_reco_sf_SelectedElectron':'a_w',
    'electron_wplid_sf_SelectedElectron':'a_w',
    'electron_wplid_sf_nomatch_SelectedGenNoRecoElectron':'a_w',
    'muon_reco_sf_SelectedMuon':'a_w',
    'muon_looseid_sf_SelectedMuon':'a_w',
    'muon_looseid_sf_nomatch_SelectedGenNoRecoMuon':'a_w',
    'muon_iso_sf_SelectedMuon':'a_w',
  }
  if s['name'] in forcedTypeMap:
    print(f" [factoryType] 使用強制映射: {s['name']} -> {forcedTypeMap[s['name']]}")
    return forcedTypeMap[s['name']]

  if('weight_LHEPd' in s['name']): return "s_w"

  missing_ws_logs = []
  scanned_any_workspace = False
  found_candidates = False

  # ---------------- 新增：工作區 / dataset 解析快取 ----------------
  if not hasattr(factoryType, "_cache"):
    factoryType._cache = {}  # key=(inputWSFile, nominalDataName) -> {'vars':set(), 'hasHist:<sname>':bool}

  # 為提升效率：先收集所有 signal rows，逐一打開直到能判斷
  for ir, r in d[d['type']=='sig'].iterrows():
    key = (r.inputWSFile, r.nominalDataName)
    # 若已有在 cache，直接使用
    if key in factoryType._cache and f"hasHist:{s['name']}" in factoryType._cache[key]:
        cache_entry = factoryType._cache[key]
        if cache_entry[f"hasHist:{s['name']}"]:
            print(f" [factoryType] {s['name']} => a_h (cache, hist pair)")
            return "a_h"
        else:
            varnames = cache_entry['vars']
            has_triplet, has_c, has_u, has_d = _detect_weight_triplet(varnames, s['name'])
            if has_triplet:
                print(f" [factoryType] {s['name']} => a_w (cache, weight triplet)")
                return "a_w"
            # 單一權重 (s_w) ：若只有 central 或 up/down 其一
            if any([has_c, has_u, has_d]):
                print(f" [factoryType] {s['name']} => s_w (cache, single weight)")
                return "s_w"
            continue  # cache 中依然無法判斷，試下一個檔案

    # 開檔
    f = ROOT.TFile.Open(r.inputWSFile)
    if (not f) or f.IsZombie():
      missing_ws_logs.append(f"[FILE-ERROR] 無法開啟: {r.inputWSFile}")
      if f: f.Close()
      continue
    scanned_any_workspace = True
    ws, logs = fetchWorkspace(f, inputWSName__)
    missing_ws_logs.extend([f"[WS-SCAN] {r.inputWSFile} {lg}" for lg in logs])
    if ws is None:
      f.Close()
      continue

    # 取得 / 建立 cache entry
    if key not in factoryType._cache:
      rdata = ws.data(r.nominalDataName)
      varnames = _extract_dataset_varnames(rdata)
      factoryType._cache[key] = {'vars': varnames}
    else:
      varnames = factoryType._cache[key]['vars']

    # 1) 判斷是否有 hist pair
    has_hist = _detect_hist_pair(ws, r.nominalDataName, s['name'])
    factoryType._cache[key][f"hasHist:{s['name']}"] = has_hist

    if has_hist:
      # 若同時又有權重 triplet，可列印提醒（預設仍以 a_h 為準）
      has_triplet, *_ = _detect_weight_triplet(varnames, s['name'])
      if has_triplet:
        print(f" [factoryType][INFO] {s['name']} 同時發現 RooDataHist 與權重欄位，優先採用 a_h (shape)")
      f.Close()
      print(f" [factoryType] {s['name']} => a_h (hist pair)")
      return "a_h"

    # 2) 權重 triplet (a_w)
    has_triplet, has_c, has_u, has_d = _detect_weight_triplet(varnames, s['name'])
    if has_triplet:
      f.Close()
      print(f" [factoryType] {s['name']} => a_w (weight triplet)")
      return "a_w"

    # 3) 單權重 / 部分權重 (視為 symmetric weight)
    if any([has_c, has_u, has_d]):
      f.Close()
      print(f" [factoryType] {s['name']} => s_w (partial weight match)")
      return "s_w"

    # 4) 傳統舊邏輯 (以名稱匹配 workspace 變數 / 函式)
    #    若有需要仍可擴充；這裡僅在上述判斷失敗時執行
    varsMatch = ws.allVars().selectByName("%s*"%(s['name']))
    if varsMatch.getSize():
      nWeights = varsMatch.getSize()
      f.Close()
      if nWeights == 2:
        print(f" [factoryType] {s['name']} => a_w (legacy varsMatch==2)")
        return "a_w"
      elif nWeights == 1:
        print(f" [factoryType] {s['name']} => s_w (legacy varsMatch==1)")
        return "s_w"
      else:
        print(f" --> [ERROR] systematic {s['name']}: >2 matches ({nWeights}). Leaving...")
        sys.exit(1)

    # 再檢查 functions
    funcs = ws.allFunctions()
    if funcs and funcs.getSize():
      func_obj = funcs.find(s['name'])
      if func_obj:
        f.Close()
        print(f" [factoryType] {s['name']} => s_w (found in allFunctions)")
        return "s_w"

    f.Close()
    # 嘗試下一個 signal workspace

  if not scanned_any_workspace:
    print(f" --> [WARNING] 未能掃描任何 signal workspace，預設 {s['name']} = s_w")
    return "s_w"

  # 仍無法判斷：預設 s_w
  print(f" --> [WARNING] 無法判斷 systematic 類型，預設使用 s_w : {s['name']}")
  if missing_ws_logs:
    print(" ---- 掃描紀錄 ----")
    for m in missing_ws_logs:
      print("    ", m)
    print(" ------------------")
  return "s_w"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to extract yield variations for signal row in dataFrame
def calcSystYields(_nominalDataName,_nominalDataContents,_inputWS,_systFactoryTypes,skipCOWCorr=True,proc="ggH",year='2016',ignoreWarnings=False):

  errMessage = "WARNING" if ignoreWarnings else "ERROR"
  errString = "Using nominal yield" if ignoreWarnings else ""

  # Define dictionary to store systematic yield counters
  systYields = {}
  # Loop over systematics and create counter in dict
  for s, f in _systFactoryTypes.items():
    if f in ["a_h","a_w"]:
      for direction in ['up','down']: 
        systYields["%s_%s"%(s,direction)] = 0
        if not skipCOWCorr: systYields["%s_%s_COWCorr"%(s,direction)] = 0
    else: 
      systYields[s] = 0
      if not skipCOWCorr: systYields["%s_COWCorr"%s] = 0

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # For systematics stored as weights (a_w,s_w) in nominal RooDataSets
  # Extract nominal dataset
  data_nominal = _inputWS.data(_nominalDataName)
  # CHECK: is weight in contents: if not then add syst to systToSkip container + print warning
  systToSkip = []
  for s,f in _systFactoryTypes.items():
    if f == "a_h":
      continue
    elif f == "a_w":
      # 需同時有 up / down 至少一個命中
      cands = _weight_name_candidates(s)
      has_up = any(u in _nominalDataContents for u in cands['up'])
      has_down = any(dn in _nominalDataContents for dn in cands['down'])
      if not (has_up and has_down):
        systToSkip.append(s)
        print(f" --> [{errMessage}] Weight triplet (up/down) for systematic ({s}) 不存在於 ({proc},{year}). {errString}")
        if not ignoreWarnings: sys.exit(1)
    else:
      # s_w：任一 central/單一權重存在即可，否則跳過
      cands = _weight_name_candidates(s)
      has_any = any(v in _nominalDataContents for v in (cands['central']+cands['up']+cands['down']))
      if not has_any:
        systToSkip.append(s)
        print(f" --> [{errMessage}] Symmetric weight for systematic ({s}) 不存在於 ({proc},{year}). {errString}")
        if not ignoreWarnings: sys.exit(1)

  # Loop over events and extract reweighted yields
  for i in range(0,data_nominal.numEntries()):
    p = data_nominal.get(i)
    w = data_nominal.weight()
    f_COWCorr = p.getRealValue("centralObjectWeight") if "centralObjectWeight" in _nominalDataContents else 1.
    f_NNLOPS = abs(p.getRealValue("NNLOPSweight")) if "NNLOPSweight" in _nominalDataContents else 1.
    # Loop over systematics:
    for s, f in _systFactoryTypes.items():

      if f == "a_h": continue

      # If asymmetric weights:
      elif f == "a_w":

        if s in systToSkip: 
          systYields["%s_up"%s] += w
          systYields["%s_down"%s] += w
          if not skipCOWCorr:
            if f_COWCorr != 0:
              systYields["%s_up_COWCorr"%s] += w*(f_NNLOPS/f_COWCorr)
              systYields["%s_down_COWCorr"%s] += w*(f_NNLOPS/f_COWCorr)

        else:
          # 解析實際變數名稱
          cands = _weight_name_candidates(s)
          centralVar = _pick_first_existing_name(cands['central'], _nominalDataContents)
          upVar = _pick_first_existing_name(cands['up'], _nominalDataContents)
          downVar = _pick_first_existing_name(cands['down'], _nominalDataContents)

          # 若 centralVar 缺失，嘗試使用 'weight_central'；再不行則當作 1
          if centralVar is None:
              centralVar = "weight_central" if "weight_central" in _nominalDataContents else None

          f_central = p.getRealValue(centralVar) if centralVar else 1.
          f_up = p.getRealValue(upVar) if upVar else f_central
          f_down = p.getRealValue(downVar) if downVar else f_central

          if f_central == 0: continue
          if (f_up == f_down): w_up, w_down = w, w
          else:
            w_up = w * (f_up / f_central)
            w_down = w * (f_down / f_central)

          systYields[f"{s}_up"] += w_up
          systYields[f"{s}_down"] += w_down
          if not skipCOWCorr and f_COWCorr != 0:
            systYields[f"{s}_up_COWCorr"] += w_up*(f_NNLOPS/f_COWCorr)
            systYields[f"{s}_down_COWCorr"] += w_down*(f_NNLOPS/f_COWCorr)

      # If symmetric weights
      else:

        if s in systToSkip:
          systYields[s] += w
          if not skipCOWCorr:
            if f_COWCorr != 0:
              systYields["%s_COWCorr"%s] += w*(f_NNLOPS/f_COWCorr)

        else:
          cands = _weight_name_candidates(s)
          centralVar = _pick_first_existing_name(cands['central'], _nominalDataContents)
          # 若只有一個 up/down 其中之一存在，把存在的當作 shift
          shiftVar = (_pick_first_existing_name(cands['up'], _nominalDataContents) or
                      _pick_first_existing_name(cands['down'], _nominalDataContents) or
                      centralVar)

          if centralVar is None:
            centralVar = "weight_central" if "weight_central" in _nominalDataContents else None

          f_central = p.getRealValue(centralVar) if centralVar else 1.
          f_shift = p.getRealValue(shiftVar) if shiftVar else f_central

          if (f_central == 0) and (f_shift == 0):
            systYields[s] += w
            if not skipCOWCorr and f_COWCorr != 0:
              systYields[f"{s}_COWCorr"] += w*(f_NNLOPS/f_COWCorr)
          elif f_central == 0: continue
          else:
            w_shift = w * (f_shift / f_central)
            systYields[s] += w_shift
            if not skipCOWCorr and f_COWCorr != 0:
              systYields[f"{s}_COWCorr"] += w_shift*(f_NNLOPS/f_COWCorr)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # For systematics stored as separate RooDataHists
  for s, f in _systFactoryTypes.items():
    if f == "a_h":
      data_hist_up, data_hist_down = _inputWS.data("%s_%sUp01sigma"%(_nominalDataName,s)), _inputWS.data("%s_%sDown01sigma"%(_nominalDataName,s))
      # Check if datasets exist: if not print warning message and set to nominal weight
      if( data_hist_up == None )|( data_hist_down == None ):
        print(" --> [%s] RooDataHist for systematic (%s) does not exist for (%s,%s). %s"%(errMessage,s,proc,year,errString))
        if not ignoreWarnings: sys.exit(1)
        systYields["%s_up"%s] = data_nominal.sumEntries()
        systYields["%s_down"%s] = data_nominal.sumEntries()
      else:
        systYields["%s_up"%s] = data_hist_up.sumEntries()
        systYields["%s_down"%s] = data_hist_down.sumEntries()

        
  # Add variations to dataFrame
  return systYields
  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPERIMENTAL SYSTEMATICS FACTORY:
# d - dataFrame, systs - dict of systematics, ftype - dict of factoryTypes
def experimentalSystFactory(d,systs,ftype,options,_removal=False):

  # Loop over systematics and add new column in dataFrame
  for s in systs:
    if s['type'] == 'constant': continue
    if s['correlateAcrossYears']: d[s['name']] = '-'
    else:
      for year in options.years.split(","): d['%s_%s'%(s['name'],year)] = '-'

  # Loop over systematics and fill entries for rows which satisfy mask
  for s in systs:
    if s['type'] == 'constant': continue
    # Extract factory type
    f = ftype[s['name']]
    if s['correlateAcrossYears']:
      mask = (d['type']=='sig')&(~d['cat'].str.contains("NOTAG"))
      d.loc[mask,s['name']] = d[mask].apply(lambda x: compareYield(x,f,s['name']), axis=1)
    else:
      for year in options.years.split(","):
        mask = (d['type']=='sig')&(~d['cat'].str.contains("NOTAG"))&(d['year']==year)
        d.loc[mask,'%s_%s'%(s['name'],year)] = d[mask].apply(lambda x: compareYield(x,f,s['name']), axis=1)

    # Remove yield columns from dataFrame
    if _removal:
      if f in ['a_h','a_w']: 
        for direction in ['up','down']: d.drop(['%s_%s_yield'%(s['name'],direction)], axis=1, inplace=True)
      else: d.drop(['%s_yield'%s['name']], axis=1, inplace=True)

  return d

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THEORY SYSTEMATICS FACTORY:
def theorySystFactory(d,systs,ftype,options,stxsMergeScheme=None,_removal=False):

  # For process yields: sum central object weight corrected (remove experimental effects)
  corrExt = "_COWCorr" if not options.skipCOWCorr else ''
   
  # Calculate the per-production mode (per-year) yield variation: add as column in dataFrame
  def _loop_sig_proc_s0(d):
    for proc_s0 in d[d['type']=='sig'].proc_s0.unique():
      for year in options.years.split(","):
        mask = (d['proc_s0']==proc_s0)&(d['year']==year)
        d.loc[mask,'proc_s0_nominal_yield'] = d[mask]['nominal_yield%s'%corrExt].sum()
        for s in systs:
          if s['type'] == 'constant': continue
          f = ftype[s['name']]
          if f in ['a_w','a_h']: 
            for direction in ['up','down']: 
              d.loc[mask,'proc_s0_%s_%s_yield'%(s['name'],direction)] = d[mask]['%s_%s_yield%s'%(s['name'],direction,corrExt)].sum()
          else:
            d.loc[mask,'proc_s0_%s_yield'%s['name']] = d[mask]['%s_yield%s'%(s['name'],corrExt)].sum()

  # Calculate the per-STXS bin (per-year already in proc name) yield variations: add as column in dataFrame
  for proc in d[d['type']=='sig'].proc.unique():
    mask = (d['proc']==proc)
    d.loc[mask,'proc_nominal_yield'] = d[mask]['nominal_yield%s'%corrExt].sum() 
    for s in systs:
      if s['type'] == 'constant': continue
      mask = (d['proc']==proc)
      f = ftype[s['name']]
      if f in ['a_w','a_h']: 
        for direction in ['up','down']: 
          d.loc[mask,'proc_%s_%s_yield'%(s['name'],direction)] = d[mask]['%s_%s_yield%s'%(s['name'],direction,corrExt)].sum()
      else: 
        d.loc[mask,'proc_%s_yield'%s['name']] = d[mask]['%s_yield%s'%(s['name'],corrExt)].sum()

  # For merging STXS bins in parameter scheme:
  if options.doSTXSMerging:
    for mergeName, mergeBins in stxsMergeScheme.items():
      for year in options.years.split(","):
        mBins = [] # add full name (inc year and and decay)
        for mb in mergeBins: mBins.append("%s_%s_hgg"%(mb,year)) 
        mask = (d['type']=='sig')&(d.apply(lambda x: x['proc'] in mBins, axis=1))
        d.loc[mask,'merge_%s_nominal_yield'%mergeName] = d[mask]['nominal_yield%s'%corrExt].sum()
        # Loop over systematics
        for s in systs:
          if s['type'] == 'constant': continue
          f = ftype[s['name']]
          if f in ['a_w','a_h']:
            for direction in ['up','down']:
              d.loc[mask,'merge_%s_%s_%s_yield'%(mergeName,s['name'],direction)] = d[mask]['%s_%s_yield%s'%(s['name'],direction,corrExt)].sum()
          else:
            d.loc[mask,'merge_%s_%s_yield'%(mergeName,s['name'])] = d[mask]['%s_yield%s'%(s['name'],corrExt)].sum()

  # Loop over systematics and add new column in dataFrame for each tier
  for s in systs:
    if s['type'] == 'constant': continue
    for tier in s['tiers']: 
      if tier == 'mnorm': 
        if options.doSTXSMerging:
          for mergeName in stxsMergeScheme: d["%s_%s_mnorm"%(s['name'],mergeName)] = '-'
      else: d["%s_%s"%(s['name'],tier)] = '-'

  # Loop over systematics and fill entries for rows which satisfy mask
  for s in systs:
    if s['type'] == 'constant': continue
    # Extract factory type
    f = ftype[s['name']]
    # For ggH theory uncertainties: require proc contains "ggH"
    if "THU_ggH" in s['name']: mask = (d['type']=='sig')&(d['nominal_yield']!=0)&(d['proc'].str.contains('ggH'))
    else: mask = (d['type']=='sig')&(d['nominal_yield']!=0)
    # Loop over tiers and use appropriate mode for compareYield function: skip mnorm as treated separately below
    for tier in s['tiers']: 
      if tier == 'mnorm': continue
      d.loc[mask,"%s_%s"%(s['name'],tier)] = d[mask].apply(lambda x: compareYield(x,f,s['name'],mode=tier), axis=1)

  # For merging STXS bins in parameter scheme: calculate mnorm systematics (merged-STXS-normalisation)
  # One nuisance per merge
  if options.doSTXSMerging:
    for mergeName in stxsMergeScheme:
      for s in systs:
        if s['type'] == 'constant': continue
        elif 'mnorm' not in s['tiers']: continue
        for year in options.years.split(","):
          # Remove NaN entries and require specific year
          mask = (d['merge_%s_nominal_yield'%mergeName]==d['merge_%s_nominal_yield'%mergeName])&(d['year']==year)&(d['nominal_yield']!=0)
          d.loc[mask,"%s_%s_mnorm"%(s['name'],mergeName)] = d[mask].apply(lambda x: compareYield(x,f,s['name'],mode='mnorm',mname=mergeName), axis=1)

  # Removal: remove yield columns from dataFrame
  if _removal:
    ids_ = ['','proc_','proc_s0_']
    if options.doSTXSMerging:
      for mergeName in stxsMergeScheme: ids_.append("merge_%s_"%mergeName)
    # Loop over systematics
    for s in systs:
      if s['type'] == 'constant': continue
      # Extract factory type
      f = ftype[s['name']]
      if f in ['a_h','a_w']: 
        for direction in ['up','down']: 
          for id_ in ids_:
            d.drop(['%s%s_%s_yield'%(id_,s['name'],direction)], axis=1, inplace=True)
      else: 
        for id_ in ids_: 
          d.drop(['%s%s_yield'%(id_,s['name'])], axis=1, inplace=True)

    # Remove also nominal yields for all combinations
    ids_.remove('')
    for id_ in ids_: d.drop(['%snominal_yield'%id_], axis=1, inplace=True)

  return d
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function for extracting systematic factors:
#   * mode == treatment of theory systematic  
def compareYield(row,factoryType,sname,mode='default',mname=None):
  # 新增：安全取得欄位工具
  def _get(row, key):
    return row[key] if key in row else None

  # 新增：只對缺失欄位警告一次
  if not hasattr(compareYield, "_missing_warned"):
    compareYield._missing_warned = set()
  def _warn_once(tag, msg):
    if tag not in compareYield._missing_warned:
      print(msg)
      compareYield._missing_warned.add(tag)

  # Catch: if any yields in denominators are zero: return 1
  if row['nominal_yield']==0:
    if factoryType in ["a_w","a_h"]: return [1.,1.]
    else: return [1.]
  if mode != 'default': #only for theory uncertainties...
    if row['proc_nominal_yield']==0:
      if factoryType in ["a_w","a_h"]: return [1.,1.]
      else: return [1.]
    if factoryType in ["a_w","a_h"]:
      for direction in ['up','down']:
        if _get(row,f"proc_{sname}_{direction}_yield") in [None,0]: return [1.,1.]
    else:
      if _get(row,f"proc_{sname}_yield") in [None,0]: return [1.]

  # 需要的欄位名稱 (依型態)
  up_y = f"{sname}_up_yield"
  down_y = f"{sname}_down_yield"
  single_y = f"{sname}_yield"

  # 檢查缺失並回傳 1
  if factoryType == "a_h":
    if (_get(row, up_y) is None) or (_get(row, down_y) is None):
      _warn_once(f"miss_{sname}", f" --> [WARNING] Missing hist yield columns for {sname}, 使用單位效應。")
      return [1.,1.]
  elif factoryType == "a_w":
    if (_get(row, up_y) is None) or (_get(row, down_y) is None):
      _warn_once(f"miss_{sname}", f" --> [WARNING] Missing weight up/down yield columns for {sname}, 使用單位效應。")
      return [1.,1.]
  else:
    if _get(row, single_y) is None:
      _warn_once(f"miss_{sname}", f" --> [WARNING] Missing symmetric yield column for {sname}, 使用單位效應。")
      return [1.]

  if( mode == 'default' )|( mode == 'ishape' ):
    if factoryType == "a_h":
      midpoint_yield = 0.5*(_get(row,down_y)+_get(row,up_y))
      if midpoint_yield == 0: return [1.,1.]
      else: return [(_get(row,down_y)/midpoint_yield),(_get(row,up_y)/midpoint_yield)]
    elif factoryType == "a_w":
      return [(_get(row,down_y)/row['nominal_yield']),(_get(row,up_y)/row['nominal_yield'])]
    else:
      return [(_get(row,single_y)/row['nominal_yield'])]

  elif mode=='shape':
    if factoryType in ["a_w","a_h"]:
      shape_up = (_get(row,up_y)/row['nominal_yield'])/(_get(row,f"proc_{sname}_up_yield")/row["proc_nominal_yield"])
      shape_down = (_get(row,down_y)/row['nominal_yield'])/(_get(row,f"proc_{sname}_down_yield")/row["proc_nominal_yield"])
      return [shape_down,shape_up]
    else:
      shape = (_get(row,single_y)/row['nominal_yield'])/(_get(row,f"proc_{sname}_yield")/row["proc_nominal_yield"])
      return [shape]

  elif mode == 'norm':
    if factoryType in ["a_w","a_h"]:
      norm_up = (_get(row,f"proc_{sname}_up_yield")/row["proc_nominal_yield"])/(_get(row,f"proc_s0_{sname}_up_yield")/row["proc_s0_nominal_yield"])
      norm_down = (_get(row,f"proc_{sname}_down_yield")/row["proc_nominal_yield"])/(_get(row,f"proc_s0_{sname}_down_yield")/row["proc_s0_nominal_yield"])
      return [norm_down,norm_up]
    else:
      norm = (_get(row,f"proc_{sname}_yield")/row["proc_nominal_yield"])/(_get(row,f"proc_s0_{sname}_yield")/row["proc_s0_nominal_yield"])
      return [norm]

  elif mode == 'mnorm':
    if factoryType in ["a_w","a_h"]:
      mnorm_up = (_get(row,f"proc_{sname}_up_yield")/row["proc_nominal_yield"])/(_get(row,f"merge_{mname}_{sname}_up_yield")/row[f"merge_{mname}_nominal_yield"])
      mnorm_down = (_get(row,f"proc_{sname}_down_yield")/row["proc_nominal_yield"])/(_get(row,f"merge_{mname}_{sname}_down_yield")/row[f"merge_{mname}_nominal_yield"])
      return [mnorm_down,mnorm_up]
    else:
      mnorm = (_get(row,f"proc_{sname}_yield")/row["proc_nominal_yield"])/(_get(row,f"merge_{mname}_{sname}_yield")/row[f"merge_{mname}_nominal_yield"])
      return [mnorm]

  elif mode == 'inorm':
    # 修正未定義變數：使用 proc_nominal_yield
    if factoryType in ["a_w","a_h"]:
      inorm_up = (_get(row,f"proc_{sname}_up_yield")/row["proc_nominal_yield"])
      inorm_down = (_get(row,f"proc_{sname}_down_yield")/row["proc_nominal_yield"])
      return [inorm_down,inorm_up]
    else:
      inorm = (_get(row,f"proc_{sname}_yield")/row["proc_nominal_yield"])
      return [inorm]

  elif mode == 'inc':
    if factoryType in ["a_w","a_h"]:
      inc_up = (_get(row,f"proc_s0_{sname}_up_yield")/row["proc_s0_nominal_yield"])
      inc_down = (_get(row,f"proc_s0_{sname}_down_yield")/row["proc_s0_nominal_yield"])
      return [inc_down,inc_up]
    else:
      inc = (_get(row,f"proc_s0_{sname}_yield")/row["proc_s0_nominal_yield"])
      return [inc]

  else:
    print(" --> [ERROR] theory systematic tier %s is not supported. Leaving"%mode)
    sys.exit(1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to group systematics: e.g. for scaleWeight where up/down = [1,2],[3,6] etc
def groupSystematics(d,systs,options,prefix="scaleWeight",groupings=[],stxsMergeScheme=None,_removal=False):
  
  # Loop over groupings
  for group_idx in range(len(groupings)): 
    gr = groupings[group_idx]

    # Extract systematic from systs
    s0, s1 = None, None
    for s in systs:
      if s['name'] == "%s_%g"%(prefix,gr[0]): s0 = s
      elif s['name'] == "%s_%g"%(prefix,gr[1]): s1 = s

    skipGroup = False
    if( s0 == None )|( s1 == None ):
      print(" --> [WARNING] No systematic exists for prefix %s and group %s. Skipping"%(prefix,gr))
      skipGroup = True
    if skipGroup: continue

    # Loop over systematic tiers
    for tier in s0['tiers']:
      if tier == 'mnorm':
        if options.doSTXSMerging:
          # Loop over merging schemes
          for mergeName in stxsMergeScheme:
            s0_name = "%s_%s_%s"%(s0['name'],mergeName,tier)
            s1_name = "%s_%s_%s"%(s1['name'],mergeName,tier)
            gr_name = "%s_gr%g_%s_%s"%(prefix,group_idx,mergeName,tier)
            d[gr_name] = '-'
            # Define mask as all entries where d[s0_name]!='-'
            mask = (d[s0_name]!='-')&(d[s1_name]!='-')
            d.loc[mask,gr_name] = d[mask].apply(lambda x: [x[s0_name][0],x[s1_name][0]],axis=1)
            # Remove original columns from dataFrame
            if _removal:
              for i in gr: d.drop( ['%s_%g_%s_%s'%(prefix,i,mergeName,tier)], axis=1, inplace=True )
        else: continue
      else:
        s0_name = "%s_%s"%(s0['name'],tier)
        s1_name = "%s_%s"%(s1['name'],tier)
        gr_name = "%s_gr%g_%s"%(prefix,group_idx,tier)
        d[gr_name] = '-'
        # Define mask as all entries where d[s0_name]!='-' and d[s1_name]!='-'
        mask = (d[s0_name]!='-')&(d[s1_name]!='-')
        d.loc[mask,gr_name] = d[mask].apply(lambda x: [x[s0_name][0],x[s1_name][0]],axis=1)
        # Remove columns from dataFrame
        if _removal:
          for i in gr: d.drop( ['%s_%g_%s'%(prefix,i,tier)], axis=1, inplace=True )

    # Replace individual systs in dict with grouped syst
    for s in systs:
      # Change w0_name and remove w1_name
      if s['name'] == "%s_%g"%(prefix,gr[0]):
        s['name'] = re.sub("%s_%g"%(prefix,gr[0]),"%s_gr%g"%(prefix,group_idx),s['name'])
        s['title'] = re.sub("%s_%g"%(prefix,gr[0]),"%s_gr%g"%(prefix,group_idx),s['title'])
      elif s['name'] == "%s_%g"%(prefix,gr[1]): systs.remove(s)

  return d,systs
      
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate the envelope of systematics with a regexp
def envelopeSystematics(d,systs,options,regexp=None,stxsMergeScheme=None,_removal=False):

  if regexp is None:
    print(" --> [WARNING] No systematics with regexp (None). Cannot form envelope")
    return d, systs

  # Extract systematics with regexp
  s_regexp = []
  for s in systs:
    if regexp in s['name']: s_regexp.append(s)
  if len(s_regexp) == 0:
    print(" --> [WARNING] No systematics with regexp (%s). Cannot form envelope"%regexp)
    return d, systs

  # Determine properties of envelope from first entry: remove "group" tag if in name
  env = {}
  env['name'] = "%s_env"%("_".join(s_regexp[0]['name'].split("_")[:-1]))
  env['title'] = "%s_env"%("_".join(s_regexp[0]['title'].split("_")[:-1]))
  env['tiers'] = s_regexp[0]['tiers'] 
  env['prior'] = s_regexp[0]['prior'] 
  env['correlateAcrossYears'] = s_regexp[0]['correlateAcrossYears'] 
  env['type'] = s_regexp[0]['type'] 

  # Loop over systematic tiers: all enveloped systematics must have the same tiers
  for s in s_regexp:
    if s['tiers'] != env['tiers']:
      print(" --> [WARNING] Systematics in envelope have different tiers. Cannot form envelope")
  for tier in env['tiers']:
    if tier == 'mnorm':
      if options.doSTXSMerging:
        # Loop over merging schemes
        for mergeName in stxsMergeScheme:
          env_name = "%s_%s_%s"%(env['name'],mergeName,tier)
          d[env_name] = '-'      
          # Define mask as all entries when first entry in envelope is set
          mask = (d['%s_%s_%s'%(s_regexp[0]['name'],mergeName,tier)]!='-')
          d.loc[mask,env_name] = d[mask].apply(lambda x: compareSystForEnvelope(x,s_regexp,tier,mname=mergeName) ,axis=1)
          # Remove original columns from dataFrame
          if _removal:
            for s in s_regexp: d.drop( ["%s_%s_%s"%(s['name'],mergeName,tier)], axis=1, inplace=True ) 
      else: continue
    else:
     env_name = "%s_%s"%(env['name'],tier)
     d[env_name] = '-'      
     # Define mask as all entries when first entry in envelope is set
     mask = (d['%s_%s'%(s_regexp[0]['name'],tier)]!='-') 
     d.loc[mask,env_name] = d[mask].apply(lambda x: compareSystForEnvelope(x,s_regexp,tier) ,axis=1)
     # Remove original columns from dataFrame
     if _removal:
       for s in s_regexp: d.drop( ["%s_%s"%(s['name'],tier)], axis=1, inplace=True ) 
 
  # Add envelope to syst dictionary
  systs.append(env)
  if _removal:
    for s in s_regexp: systs.remove(s)
 
  return d, systs

# Function to compare systematic variation for envelope
def compareSystForEnvelope(row,systs,stier,mname=None):
  e_symm_max = 0.
  env_idx = 0
  for sidx in range(len(systs)):
    s = systs[sidx]
    if mname is not None: sname = "%s_%s_%s"%(s['name'],mname,stier)
    else: sname = "%s_%s"%(s['name'],stier)
    if len(row[sname]) == 2: e_symm = 0.5*(abs(row[sname][0]-1)+abs(row[sname][1]-1))
    else: e_symm = abs(row[sname][0]-1)
    if e_symm > e_symm_max: 
      e_symm_max = e_symm
      env_idx = sidx
  env_s = systs[env_idx]
  if mname is not None: env_sname = "%s_%s_%s"%(env_s['name'],mname,stier)
  else: env_sname = "%s_%s"%(env_s['name'],stier)
  return row[env_sname]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to change syst title
def renameSyst(t,oldexp,newexp):
  # 修正: 原本遺漏第三個參數 (string)，導致 TypeError
  return re.sub(oldexp, newexp, t)
