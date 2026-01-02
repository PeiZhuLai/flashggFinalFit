# Script to calculate yields from input flashgg workspaces
#  * Uses Pandas dataframe to store all proc x cat yields
#  * Including systematic variations
#  * Output to be used for creating datacard

import os, sys
import re
from optparse import OptionParser
import ROOT
from ROOT import RooArgList, RooFormulaVar
import pandas as pd
import glob
import pickle
import math
from collections import OrderedDict
from systematics_HToZa import theory_systematics, experimental_systematics, signal_shape_systematics

from commonObjects import *
from commonTools import *

ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
interploate_ma_list = [11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29]

def _nearest_anchor_mass(mass, anchors):
  """給定 mass(int/float)，回傳 anchors 中距離最近者。"""
  m = int(round(float(mass)))
  return min(anchors, key=lambda a: abs(a - m))

def resolve_mass_for_io(mass_alp, anchors, interpolate_list):
  """
  若 mass_alp 在 interpolate_list，回傳最近鄰 anchor (用於讀檔/取模型)；
  否則回傳 mass_alp 本身。
  """
  m = int(round(float(mass_alp)))
  if m in interpolate_list:
    return _nearest_anchor_mass(m, anchors)
  return m

# 讓 od() 可用（有些環境沒在 commonTools 定義）
od = OrderedDict

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZa DATACARD MAKER RUN III ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
def leave():
  print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZa DATACARD MAKER RUN III (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit(0)

def get_options():
  parser = OptionParser()

  parser.add_option('--mass_ALP', dest='mass_ALP', default=5, type='int', help="ALP mass") # PZ
  parser.add_option('--year', dest='year', default='16', help="year") # PZ
  parser.add_option("--channel", dest='channel', default='', help="ele, mu, or leptons") # PZ

  parser.add_option('--inputWSDirMap', dest='inputWSDirMap', default='2022preEE:/eos/home-p/pelai/HZa/root_MVAcut', help="Map. Format: year=inputWSDir (separate years by comma)")
  parser.add_option('--cat', dest='cat', default='cat0', help='Analysis category')
  parser.add_option('--procs', dest='procs', default='GG2H', help='Comma separated list of signal processes. auto = automatically inferred from input workspaces')
  parser.add_option('--ext', dest='ext', default='', help='Extension for saving') 
  parser.add_option('--mass', dest='mass', default='125', help='Input workspace mass')
  parser.add_option('--mergeYears', dest='mergeYears', default=True, action="store_true", help="Merge category across years")
  parser.add_option('--skipBkg', dest='skipBkg', default=False, action="store_true", help="Only add signal processes to datacard")
  parser.add_option('--bkgScaler', dest='bkgScaler', default=1., type="float", help="Add overall scale factor for background")
  parser.add_option('--sigModelWSDir', dest='sigModelWSDir', default='/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal', help='Input signal model WS directory') 
  parser.add_option('--sigModelExt', dest='sigModelExt', default='packaged', help='Extension used when saving signal model') 
  parser.add_option('--bkgModelWSDir', dest='bkgModelWSDir', default='/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3', help='Input background model WS directory') 
  parser.add_option('--bkgModelExt', dest='bkgModelExt', default='multipdf', help='Extension used when saving background model') 
  # For yields calculations:
  parser.add_option('--skipZeroes', dest='skipZeroes', default=False, action="store_true", help="Skip signal processes with 0 sum of weights")
  parser.add_option('--skipCOWCorr', dest='skipCOWCorr', default=False, action="store_true", help="Skip centralObjectWeight correction for events in acceptance. Use if no centralObjectWeight in workspace")
  # For systematics:
  parser.add_option('--doSystematics', dest='doSystematics', default=False, action="store_true", help="Include systematics calculations and add to datacard")
  parser.add_option('--ignore-warnings', dest='ignore_warnings', default=False, action="store_true", help="Skip errors for missing systematics. Instead output warning message")
  parser.add_option('--debugNames', dest='debugNames', default=False, action="store_true", help="列印所有組合字串以除錯名稱")
  parser.add_option('--missingDatasetAction', dest='missingDatasetAction', default='skip',
                    help="當找不到 nominal RooDataSet 時的動作: error|skip|guess (guess=嘗試模糊匹配)")
  parser.add_option('--disableAutoWeightFix', dest='disableAutoWeightFix', default=False, action="store_true",
                    help="停用：自動把缺 RooDataHist 但具 weight_* 變數的系統誤差由 a_h 轉為 a_w")
  return parser.parse_args()

(opt,args) = get_options()

# 決定「實際用來讀檔/取模型」的 mass
mass_for_io = resolve_mass_for_io(opt.mass_ALP, ma_list, interploate_ma_list)
if mass_for_io != int(opt.mass_ALP):
  print(f" --> [INFO] mass_ALP={opt.mass_ALP} 不在 anchor 清單，改用最近鄰 mA={mass_for_io} 的檔案進行計算")

# Extract years and inputWSDir
inputWSDirMap = od()
for i in opt.inputWSDirMap.split(","): 
  print(" --> Taking %s input workspaces from: %s"%(i.split("=")[0],i.split("=")[1]) )
  if not os.path.isdir( i.split("=")[1] ):
    print(" --> [ERROR] Directory %s does not exist. Leaving..."%i.split("=")[1])
    leave()
  inputWSDirMap[i.split("=")[0]] = i.split("=")[1]
years = list(inputWSDirMap.keys())
# whether change years Pei-Zhu

procsMap = od()
if opt.procs == 'auto':
  for y,iWSDir in inputWSDirMap.items():
    WSFileNames = extractWSFileNames(iWSDir)
    procsMap[y] = extractListOfProcs(WSFileNames)
  # Require common procs for each year
  for i,iy in enumerate(years):
    for j,jy in enumerate(years):
      if j > i:
        if set(procsMap[iy].split(",")) != set(procsMap[jy].split(",")):
          print(" --> [ERROR] Mis-match in list of process for %s and %s. Intersection = %s"%(iy,jy,(set(procsMap[jy]).symmetric_difference(set(procsMap[iy])))))
          leave()
  # Define list of procs (alphabetically ordered)
  procs = procsMap[years[0]].split(",")
else:
  procs = opt.procs.split(",")
procs.sort()

# Initiate pandas dataframe
columns_data = ['year','type','procOriginal','proc','proc_s0','cat','inputWSFile','nominalDataName','modelWSFile','model','rate']
data = pd.DataFrame( columns=columns_data )

# 工具：列出 workspace / 檔案中的 dataset 名稱（無 ws.dir() 依賴）
def list_workspace_datasets(ws, tfile=None):
  """列出 workspace 中的 RooDataSet 名稱；若提供 tfile，則遞迴掃描整個 ROOT 檔找出 RooData(Set/Hist) 名稱一起回傳。"""
  names = set()
  # 1) 從 RooWorkspace 取得
  try:
    it = ws.allData()  # RooAbsCollection
    for obj in it:
      names.add(obj.GetName())
  except Exception:
    pass

  # 2) 也可從整個 ROOT 檔遞迴掃描（可選）
  if tfile is not None:
    def _recurse(dir_obj):
      lst = getattr(dir_obj, "GetListOfKeys", lambda: None)()
      if not lst:
        return
      for key in lst:
        try:
          obj = key.ReadObj()
        except Exception:
          continue
        try:
          cname = obj.ClassName()
        except Exception:
          cname = ""
        if "RooDataSet" in cname or "RooDataHist" in cname:
          names.add(key.GetName())
        # 深入子目錄（TDirectory*）
        try:
          if hasattr(obj, "GetListOfKeys") or obj.InheritsFrom("TDirectory"):
            _recurse(obj)
        except Exception:
          pass
    try:
      _recurse(tfile)
    except Exception:
      pass

  return sorted(names)

# 工具：模糊猜測 dataset 名稱
def guess_closest_name(target, candidates):
  if not candidates: return None
  best, best_score = None, -1.0
  for c in candidates:
    common = sum(1 for ch in target if ch in c)
    denom = float(len(target) if len(target)>0 else 1)
    score = common / denom
    if target in c: score += 0.3
    if score > best_score:
      best, best_score = c, score
  return best

# 新增：偵測某 systematic 是否以 weight 形式存在 (central/Up/Down)
def detect_weight_based_systematic(syst_name, contents_string):
  """回傳 True 若 dataset 內容字串含有 weight_<syst_name>_{central,Up,Down}。"""
  base = f"weight_{syst_name}_"
  required = [base + "central", base + "Up", base + "Down"]
  return all(tok in contents_string for tok in required)

# 新增：為 weight_*_central / _Up / _Down 建立別名 (若缺乏無 _central 版本)
def ensure_weight_aliases(rdata, syst_name, debug=False):
  """
  確保 dataset 內存在常見多種權重命名:
    weight_<syst_name>
    weight_<syst_name>Up,   weight_<syst_name>Down
    weight_<syst_name>_up,  weight_<syst_name>_down
  (原本的 *_central / *_Up / *_Down 會被當作來源)
  回傳 True 若有新增任一 alias。
  """
  if rdata is None or rdata.numEntries() == 0:
    return False
  first = rdata.get(0)
  def _has(var): 
    return first.find(var) is not None

  made_any = False

  central_src = f"weight_{syst_name}_central"
  up_src      = f"weight_{syst_name}_Up"
  down_src    = f"weight_{syst_name}_Down"

  alias_targets = [
    (central_src, f"weight_{syst_name}"),
    (up_src,      f"weight_{syst_name}Up"),
    (down_src,    f"weight_{syst_name}Down"),
    (up_src,      f"weight_{syst_name}_up"),
    (down_src,    f"weight_{syst_name}_down"),
  ]

  for src, dest in alias_targets:
    if _has(src) and not _has(dest):
      var_src = first.find(src)
      aliasF = RooFormulaVar(dest, "@0", RooArgList(var_src))
      rdata.addColumn(aliasF)
      made_any = True

  if made_any and debug:
    print(f" [Alias] 新增 {syst_name} 權重別名 (central/Up/Down 多種格式)")

  return made_any

# (可選) 小工具：排版輸出
def _print_block(title, dct):
  print("[DEBUG] ================== %s ==================" % title)
  for k,v in dct.items():
    print("  %s = %s" % (k,v))
  print("[DEBUG] ----------------------------------------")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILL DATAFRAME: all processes
print(" ..........................................................................................")

# Signal processes
for year in years:
  # PZ: define lep channel
  if opt.channel == 'leptons':
    leps = ['ele', 'mu']
  else:
    leps = [opt.channel]

  for proc in procs:
    for lep_channel in leps:

      # Identifier
      _id = "%s_%s_%s_%s_%s"%(proc,year,lep_channel,opt.cat,sqrts__)
      origin_id = "%s_%s_%s_%s"%(proc,year,opt.cat,sqrts__)
      
      # Mapping to STXS definition here
      _procOriginal = proc
      _proc = "%s_%s_%s"%(procToDatacardName(proc),year,lep_channel)
      _proc_s0 = procToData(proc.split("_")[0])

      # Define category: add year tag if not merging
      if opt.mergeYears: _cat = opt.cat
      else: _cat = "%s_%s"%(opt.cat,year)

      # Input Signal flashgg ws 
      _inputWS_pattern = f"{inputWSDirMap[year]}/sig/mA_M{mass_for_io}/ws_Tree2WS/ws_{lep_channel}_{year}.root"
      _inputWSFile_list = glob.glob(_inputWS_pattern)
      if len(_inputWSFile_list) == 0:
        print(f" --> [WARNING] 找不到工作區檔案: pattern={_inputWS_pattern} (skip)")
        continue
      if len(_inputWSFile_list) > 1:
        print(f" --> [WARNING] 匹配到多個檔案(取第一個): {_inputWSFile_list}")
      _inputWSFile = _inputWSFile_list[0]

      if opt.debugNames:
        if not os.path.isfile(_inputWSFile):
          print(f"[DEBUG][ERROR] 檔案不存在: {_inputWSFile}")
        else:
          rchk = ROOT.TFile.Open(_inputWSFile)
          obj = rchk.Get(inputWSName__)
          print(f"[DEBUG] 檢查 workspace  '{inputWSName__}' in {_inputWSFile} -> type={type(obj)}")
          rchk.Close()

      _nominalDataName = "%s_%s_Za_%s_%s_%s"%(_proc_s0,opt.mass,lep_channel,sqrts__,opt.cat)

      # If opt.skipZeroes check nominal yield if 0 then do not add
      skipProc = False
      if opt.skipZeroes:
        f = ROOT.TFile(_inputWSFile)
        w = f.Get(inputWSName__)
        sumw = w.data(_nominalDataName).sumEntries()
        if sumw == 0.: skipProc = True
        w.Delete()
        f.Close()
      if skipProc: continue

      # Input Signal model ws 
      if opt.cat == "NOTAG": 
        _modelWSFile, _model = '-', '-'
      else:
        _modelWSFile = f"{opt.sigModelWSDir}/outdir_{lep_channel}/signalFit/output/{mass_for_io}_CMS-HGG_sigfit_{year}_{lep_channel}_Hm125.root"
        _model = "%s_%s:%s_%s"%(outputWSName__,sqrts__,outputWSObjectTitle__,origin_id)

      # Extract rate from lumi
      _rate = float(lumiMap[year])*1000

      if opt.debugNames:
        _print_block("Signal 名稱組合",
          {
            "_id (唯一識別: proc_year_lep_cat_sqrts)" : _id,
            "origin_id (原始信號組合 key)" : origin_id,
            "_procOriginal (輸入工作區原始流程名)" : _procOriginal,
            "_proc (轉成 datacard 使用之流程+年份+lepton)" : _proc,
            "_proc_s0 (基礎歸一化流程別)" : _proc_s0,
            "_cat (類別, 可能含年分或合併)" : _cat,
            "_inputWSFile_pattern (glob 模式)" : _inputWS_pattern,
            "_inputWSFile_list (實際匹配列表)" : _inputWSFile_list,
            "_inputWSFile (存入 DataFrame 的字串)" : _inputWSFile,
            "_nominalDataName (RooDataSet 名稱)" : _nominalDataName,
            "_modelWSFile (信號模型工作區路徑)" : _modelWSFile if opt.cat != "NOTAG" else "(NOTAG 無模型)",
            "_model (模型對象 full spec)" : _model if opt.cat != "NOTAG" else "(NOTAG 無模型)",
            "_rate (由 lumiMap[year]*1000)" : _rate
          }
        )

      # Add signal process to dataFrame:
      print(" --> Adding to dataFrame: (proc,cat) = (%s,%s)"%(_proc,_cat))
      data.loc[len(data)] = [year,'sig',_procOriginal,_proc,_proc_s0,_cat,_inputWSFile,_nominalDataName,_modelWSFile,_model,_rate]

# Background and data processes
if (not opt.skipBkg) & (opt.cat != "NOTAG"):
  _proc_bkg = "bkg_mass"
  _proc_data = "data_obs"
  if opt.mergeYears:
    _cat = opt.cat
    _modelWSFile = "%s/%s/CMS-HGG_mva_13p6TeV_multipdf.root"%(opt.bkgModelWSDir, mass_for_io)
    _model_bkg = "%s:CMS_%s_%s_%s_bkgshape"%(bkgWSName__,decayMode,_cat,sqrts__)
    _model_data = "%s:roohist_data_mass_%s"%(bkgWSName__,_cat)
    _proc_s0 = 'ggH' # not needed for data/bkg
    # 原碼 year 未定義；用第一個 year 避免 NameError（僅選來源，不改邏輯）
    year_for_data = years[0]
    _inputWSFile = "%s/data/mA_M%s/ws/run3.root"%(inputWSDirMap[year_for_data], mass_for_io)
    _nominalDataName = "Data_13p6TeV" # Pei-Zhu
    if opt.debugNames:
      _print_block("Background/Data (合併年分)",
        {
          "_proc_bkg" : _proc_bkg,
          "_proc_data" : _proc_data,
          "_cat (合併後類別)" : _cat,
          "_modelWSFile (背景 multipdf)" : _modelWSFile,
          "_model_bkg (bkg pdf 名稱)" : _model_bkg,
          "_model_data (data hist 名稱)" : _model_data,
          "_inputWSFile (原始 data/bkg 載入, 供參考)" : _inputWSFile,
          "_nominalDataName (此處不使用, 佔位)" : _nominalDataName
        }
      )
    print(" --> Adding to dataFrame: (proc,cat) = (%s,%s)"%(_proc_bkg,_cat))
    print(" --> Adding to dataFrame: (proc,cat) = (%s,%s)"%(_proc_data,_cat))
    data.loc[len(data)] = ["merged",'bkg',_proc_bkg,_proc_bkg,'-',_cat,_inputWSFile,_nominalDataName,_modelWSFile,_model_bkg,opt.bkgScaler]
    data.loc[len(data)] = ["merged",'data',_proc_data,_proc_data,'-',_cat,_inputWSFile,_nominalDataName,_modelWSFile,_model_data,-1]

  else:
    for year in years:
      _cat = "%s_%s"%(opt.cat,year)
      _modelWSFile = "%s/%s/CMS-HGG_mva_13TeV_multipdf.root"%(opt.bkgModelWSDir, mass_for_io)
      _model_bkg = "%s:CMS_%s_%s_%s_bkgshape"%(bkgWSName__,decayMode,_cat,sqrts__)
      _model_data = "%s:roohist_data_mass_%s"%(bkgWSName__,_cat) # Pei-Zhu 
      _proc_s0 = 'ggH' # not needed for data/bkg
      _inputWSFile = "%s/data/ALP_data_bkg_Am%s_workspace.root"%(inputWSDirMap[year], mass_for_io) # Pei-Zhu 
      _nominalDataName = 'ggh_125_13TeV_cat0' # Pei-Zhu
      if opt.debugNames:
        _print_block("Background/Data (逐年)",
          {
            "year" : year,
            "_proc_bkg" : _proc_bkg,
            "_proc_data" : _proc_data,
            "_cat (含年份)" : _cat,
            "_modelWSFile" : _modelWSFile,
            "_model_bkg" : _model_bkg,
            "_model_data" : _model_data,
            "_inputWSFile" : _inputWSFile,
            "_nominalDataName (佔位)" : _nominalDataName
          }
        )
      print(" --> Adding to dataFrame: (proc,cat) = (%s,%s)"%(_proc_bkg,_cat))
      print(" --> Adding to dataFrame: (proc,cat) = (%s,%s)"%(_proc_data,_cat))
      data.loc[len(data)] = ["year",'bkg',_proc_bkg,_proc_bkg,'-',_cat,_inputWSFile,_nominalDataName,_modelWSFile,_model_bkg,opt.bkgScaler]
      data.loc[len(data)] = ["year",'data',_proc_data,_proc_data,'-',_cat,_inputWSFile,_nominalDataName,_modelWSFile,_model_data,-1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Yields: for each signal row in dataFrame extract the yield
print(" ..........................................................................................")
#   * if systematics=True: also extract reweighted yields for each uncertainty source
from tools.calcSystematics import factoryType, calcSystYields
from tools.wsUtils import fetchWorkspace

# Create columns in dataFrame to store yields
data['nominal_yield'] = '-'
data['sumw2'] = '-'
if not opt.skipCOWCorr:
  data['nominal_yield_COWCorr'] = '-'

# Add columns in dataFrame for systematic yield variations
if opt.doSystematics:
  # Extract type of systematic using factoryType function (defined in tools.calcSystematics)
  #  * a_h: anti-symmetric RooDataHist (2 columns in dataframe)
  #  * a_w: anti-symmetric weight in nominal RooDataSet (2 columns in dataframe)
  #  * s_w: symmetric (single) weight in nominal RooDataSet (1 column in dataframe)
  experimentalFactoryType = {}
  theoryFactoryType = {}
  # No experimental systematics for NOTAG
  if opt.cat != "NOTAG":
    for s in experimental_systematics: 
      if s['type'] == 'factory': 
        # Fix for HEM as only in 2018 workspaces
        if s['name'] == 'JetHEM':
          experimentalFactoryType[s['name']] = "a_h"
        else:
          experimentalFactoryType[s['name']] = factoryType(data,s)
        if experimentalFactoryType[s['name']] in ["a_w","a_h"]:
          data['%s_up_yield'%s['name']] = '-'
          data['%s_down_yield'%s['name']] = '-'
        else:
          data['%s_yield'%s['name']] = '-'
  for s in theory_systematics: 
    if s['type'] == 'factory': 
      theoryFactoryType[s['name']] = factoryType(data,s)
      if theoryFactoryType[s['name']] in ["a_w","a_h"]:
        data['%s_up_yield'%s['name']] = '-'
        data['%s_down_yield'%s['name']] = '-'
        if not opt.skipCOWCorr:
          data['%s_up_yield_COWCorr'%s['name']] = '-'
          data['%s_down_yield_COWCorr'%s['name']] = '-'
      else: 
        data['%s_yield'%s['name']] = '-'
        if not opt.skipCOWCorr:
          data['%s_yield_COWCorr'%s['name']] = '-'

# Loop over signal rows in dataFrame: extract yields (nominal & systematic variations)
totalSignalRows = float(data[data['type']=='sig'].shape[0])
for ir,r in data[data['type']=='sig'].iterrows():

  print(" --> Extracting yields: (%s,%s) [%.1f%%]"%(r['proc'],r['cat'],100*(float(ir)/totalSignalRows)))

  # Open input WS file and extract workspace
  if opt.debugNames:
    print("[DEBUG] 開始讀取信號工作區: file=%s, workspace=%s, nominalData=%s" %
          (r.inputWSFile, inputWSName__, r.nominalDataName))
  f_in = ROOT.TFile(r.inputWSFile)
  if (not f_in) or f_in.IsZombie():
    print(f"[ERROR] 無法開啟工作區檔案: {r.inputWSFile}")
    sys.exit(1)
  ws, logs = fetchWorkspace(f_in, inputWSName__)
  for lg in logs:
    print(f"[DEBUG][WS] {r.inputWSFile} {lg}")
  if ws is None:
    print(f"[ERROR] 找不到 RooWorkspace '{inputWSName__}' 於檔案: {r.inputWSFile}")
    sys.exit(1)
  inputWS = ws

  # ---------------- Robustly fetch nominal RooDataSet ----------------
  ds_name = r['nominalDataName']
  rdata_nominal = inputWS.data(ds_name)

  def _safe_num_entries(ds):
    """Return number of entries or -1 if ds is an invalid/null PyROOT proxy."""
    try:
      return int(ds.numEntries())
    except ReferenceError:
      return -1
    except Exception:
      return -1

  invalid_ds = (rdata_nominal is None) or (_safe_num_entries(rdata_nominal) < 0)
  if invalid_ds:
    # 列出可用名稱（方便你比對），並依據 --missingDatasetAction 處理
    all_ds = list_workspace_datasets(inputWS, f_in)
    print(f" [WARNING] 找不到或無效的 RooDataSet: {ds_name}")
    if opt.debugNames:
      print(" [DEBUG] 可用資料集清單 (%d): %s" % (len(all_ds), ", ".join(all_ds)))
    action = opt.missingDatasetAction.lower()

    if action == 'guess':
      guess = guess_closest_name(ds_name, all_ds)
      if guess:
        print(f" [INFO] 使用猜測資料集名稱: {guess} (原: {ds_name})")
        rdata_nominal = inputWS.data(guess)
        if _safe_num_entries(rdata_nominal) < 0:
          print(" [WARNING] 猜測到的資料集仍是無效指標 -> 視同找不到")
          rdata_nominal = None
      else:
        print(" [WARNING] 無法猜測相近資料集名稱 -> 視同找不到")

    if (rdata_nominal is None) or (_safe_num_entries(rdata_nominal) < 0):
      if action == 'error':
        print(" [ERROR] 依使用者設定 --missingDatasetAction=error -> 停止")
        sys.exit(2)
      # skip / guess失敗：填入 0 並繼續
      data.at[ir,'nominal_yield'] = 0.
      data.at[ir,'sumw2'] = 0.
      if not opt.skipCOWCorr: data.at[ir,'nominal_yield_COWCorr'] = 0.
      if opt.doSystematics:
        for c in data.columns:
          if c.endswith("_yield") or c.endswith("_yield_COWCorr"):
            if data.at[ir,c] == '-': data.at[ir,c] = 0.
      # 清理並跳下一列
      inputWS.Delete()
      f_in.Close()
      continue

  # ---------------- Calculate nominal yield (safe) ----------------
  contents = ""
  y, y_COWCorr = 0.0, 0.0
  sumw2 = 0.0

  n_entries = _safe_num_entries(rdata_nominal)
  if n_entries < 0:
    # 理論上不會到這（上面已處理），但再保護一次
    data.at[ir,'nominal_yield'] = 0.
    data.at[ir,'sumw2'] = 0.
    if not opt.skipCOWCorr:
      data.at[ir,'nominal_yield_COWCorr'] = 0.
    inputWS.Delete()
    f_in.Close()
    continue

  for i in range(n_entries):
    p = rdata_nominal.get(i)
    w = rdata_nominal.weight()
    y += w
    sumw2 += w*w
    if i == 0:
      contents = p.contentsString()
    if not opt.skipCOWCorr:
      f_COWCorr = p.getRealValue("centralObjectWeight") if "centralObjectWeight" in contents else 1.
      f_NNLOPS = abs(p.getRealValue("NNLOPSweight")) if "NNLOPSweight" in contents else 1.
      if f_COWCorr != 0:
        y_COWCorr += w*(f_NNLOPS/f_COWCorr)

  # Need to times 2 on sumw b/c only use half of the dataset
  # Times 4 on sumw2
  data.at[ir,'nominal_yield'] = 2.0*y
  data.at[ir,'sumw2'] = 4.0*sumw2
  if not opt.skipCOWCorr:
    # 也需乘 2.0，與 nominal_yield 一致（只用半個 dataset）
    data.at[ir,'nominal_yield_COWCorr'] = 2.0*y_COWCorr

  if opt.debugNames:
    print("[DEBUG] Nominal yield 結果: proc=%s cat=%s yield=%.6f sumw2=%.6f%s" %
          (r['proc'], r['cat'], y, sumw2,
           ("" if opt.skipCOWCorr else " yield_COWCorr=%.6f"%y_COWCorr)))

  # 在進入系統誤差計算前，動態修正 factoryType：若原判 a_h 但其實只有權重
  if opt.doSystematics and (not opt.disableAutoWeightFix):
    # 只需檢查一次 dataset 內容字串 contents 已取得
    # (1) 處理 experimental
    if "NOTAG" not in r['cat']:
      for sname, ftype in list(experimentalFactoryType.items()):
        if ftype == 'a_h':
          if detect_weight_based_systematic(sname, contents):
            experimentalFactoryType[sname] = 'a_w'
            print(f" [AutoFix] 將系統誤差 {sname} 由 a_h 轉為 a_w (使用權重欄位) for ({r['proc']},{r['year']})")
    # (2) 處理 theory
    for sname, ftype in list(theoryFactoryType.items()):
      if ftype == 'a_h':
        if detect_weight_based_systematic(sname, contents):
          theoryFactoryType[sname] = 'a_w'
          print(f" [AutoFix] 將理論系統誤差 {sname} 由 a_h 轉為 a_w (使用權重欄位) for ({r['proc']},{r['year']})")

  # 在呼叫 calcSystYields 前插入 alias 建立 (只針對 a_w)
  if opt.doSystematics:
    alias_added = False
    # experimental
    if "NOTAG" not in r['cat']:
      for sname, ftype in experimentalFactoryType.items():
        if ftype == 'a_w':
          if ensure_weight_aliases(rdata_nominal, sname, debug=opt.debugNames):
            alias_added = True
    # theory
    for sname, ftype in theoryFactoryType.items():
      if ftype == 'a_w':
        if ensure_weight_aliases(rdata_nominal, sname, debug=opt.debugNames):
          alias_added = True
    # 若有新增 alias，更新 contents 以供後續 calcSystYields 解析
    if alias_added:
      contents = rdata_nominal.get(0).contentsString()

  # Systematics: loop over systematics and use function to extract yield variations
  if opt.doSystematics:

    # For experimental systematics: skip NOTAG events
    if "NOTAG" not in r['cat']:
      # Skip centralObjectWeight correction as concerns events in acceptance
      experimentalSystYields = calcSystYields(
        r['nominalDataName'], contents, inputWS, experimentalFactoryType,
        skipCOWCorr=True, proc=r['proc'], year=r['year'], ignoreWarnings=opt.ignore_warnings
      )
      for s,f in experimentalFactoryType.items():
        if f in ['a_w','a_h']: 
          for direction in ['up','down']: 
            data.at[ir,"%s_%s_yield"%(s,direction)] = experimentalSystYields["%s_%s"%(s,direction)]
        else:
          data.at[ir,"%s_yield"%s] = experimentalSystYields[s]

    # For theoretical systematics:
    theorySystYields = calcSystYields(
      r['nominalDataName'], contents, inputWS, theoryFactoryType,
      skipCOWCorr=opt.skipCOWCorr, proc=r['proc'], year=r['year'], ignoreWarnings=opt.ignore_warnings
    )
    for s,f in theoryFactoryType.items():
      if f in ['a_w','a_h']: 
        for direction in ['up','down']: 
          data.at[ir,"%s_%s_yield"%(s,direction)] = theorySystYields["%s_%s"%(s,direction)]
          if not opt.skipCOWCorr:
            data.at[ir,"%s_%s_yield_COWCorr"%(s,direction)] = theorySystYields["%s_%s_COWCorr"%(s,direction)]
      else:
        data.at[ir,"%s_yield"%s] = theorySystYields[s]
        if not opt.skipCOWCorr:
          data.at[ir,"%s_yield_COWCorr"%s] = theorySystYields["%s_COWCorr"%s]

    # === 新增：系統誤差 yield 全部乘 2.0（因使用半 dataset） ===
    for col in data.columns:
      if (col.endswith("_yield") or col.endswith("_yield_COWCorr")) and not col.startswith("nominal_"):
        val = data.at[ir, col]
        if val != '-' and val is not None:
          try:
            data.at[ir, col] = 2.0 * float(val)
          except Exception:
            pass
    if opt.debugNames:
      print(f"[DEBUG] 系統誤差 yield 已統一乘 2.0: proc={r['proc']} cat={r['cat']}")

  if opt.doSystematics and opt.debugNames:
    print("[DEBUG] 已填入系統誤差變動: proc=%s cat=%s" % (r['proc'], r['cat']))

  # Remove the workspace and file from heap
  inputWS.Delete()
  f_in.Close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE YIELDS DATAFRAME
print(" ..........................................................................................")
extStr = "_%s"%opt.channel if opt.channel != '' else ''
print(" --> Saving yields dataframe: ./yields%s/%s_datacard_%s.pkl"%(extStr,opt.mass_ALP,opt.channel))
if not os.path.isdir("./yields%s"%extStr):
  os.system("mkdir ./yields%s"%extStr)
if opt.debugNames:
  print("[DEBUG] 最終 DataFrame 欄位: %s" % list(data.columns))
  print("[DEBUG] 儲存路徑: ./yields%s/%s_datacard_%s.pkl" %
        (extStr,opt.mass_ALP,opt.channel))
with open("./yields%s/%s_datacard_%s.pkl"%(extStr,opt.mass_ALP,opt.channel),"wb") as fD:
  pickle.dump(data,fD)
