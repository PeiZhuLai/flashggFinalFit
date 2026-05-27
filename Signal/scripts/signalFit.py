import ROOT
# 在載入 pandas/numpy 之前抑制 numpy.core.getlimits 相關的 UserWarning
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module=r"numpy\.core\.getlimits")
import pandas as pd
import pickle
import math
import os, sys
import json
from optparse import OptionParser
import glob
import re
from collections import OrderedDict as od

from commonTools import *
from commonObjects import *
from signalTools import *
from replacementMap import globalReplacementMap
from XSBRMap import *
from simultaneousFit import *
from finalModel import *
from plottingTools import *

# Constant
# MHLow, MHHigh = '120', '130'
MHLow, MHHigh = '100', '180' # In this way, the result will be as same as the fTest
MHNominal = '125'

def _resolve_fit_range(mass_alp, default_low, default_high):
    try:
        mA = int(mass_alp)
    except Exception:
        return default_low, default_high
    if mA <= 1:
        return '118', '135'
    if mA == 2:
        return '115', '135'
    if mA <= 4:
        return '110', '140'
    return default_low, default_high

def _force_dcb_for_low_ma(mass_alp, current_useDCB):
    """Low ma has a real bremsstrahlung/merged-photon low-side tail; DCB describes it
    with a power-law tail instead of inflating sigma like nGaussians do."""
    try:
        mA = int(mass_alp)
    except Exception:
        return current_useDCB
    if mA <= 4:
        return True
    return current_useDCB

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZallgg SIGNAL FITTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZallgg SIGNAL FITTER (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit()

# 新增：穩健取得 RooWorkspace 的工具函式
def get_workspace(tfile, ws_hint):
  # 優先用 hint 嘗試抓取
  obj = tfile.Get(ws_hint) if ws_hint else None
  def _find_ws(node):
    if not node: 
      return None
    # 直接是 RooWorkspace
    if hasattr(ROOT, "RooWorkspace") and node.InheritsFrom("RooWorkspace"):
      return node
    # 在目錄裡遞迴找
    if node.InheritsFrom("TDirectory"):
      keys = node.GetListOfKeys()
      if keys:
        for k in keys:
          o = k.ReadObj()
          ws = _find_ws(o)
          if ws:
            return ws
    return None

  ws = _find_ws(obj)
  if ws:
    return ws

  # 常見名稱的回退列表
  candidates = [
    ws_hint,
    "CMS_hza_workspace",
    "cms_hgg_13TeV",
    "workspace",
    "w",
    "CMS_hgg_workspace"
  ]
  for name in candidates:
    if not name:
      continue
    o = tfile.Get(name)
    ws = _find_ws(o)
    if ws:
      return ws

  # 全檔遞迴搜尋
  ws = _find_ws(tfile)
  if ws:
    return ws

  print(" --> [ERROR] Could not retrieve RooWorkspace from file: %s" % tfile.GetName())
  print("     Please check inputWSName__ and file structure. tfile content:")
  tfile.ls()
  leave()

# 新增：從 workspace 以名稱穩健取得 RooAbsData 的工具函式
def get_dataset(ws, name):
  ds = ws.data(name)
  if ds:
    return ds
  print(" --> [ERROR] Dataset '%s' not found in workspace '%s'." % (name, ws.GetName()))
  try:
    print("     Workspace content:")
    ws.Print()
  except Exception:
    pass
  leave()

# 新增：systematics 合併工具
def _to_list(val):
  if not val: return []
  if isinstance(val, (list, tuple, set)): return [str(x) for x in val if str(x)]
  return [x for x in re.split(r'[,\s]+', str(val).strip()) if x]

def _uniq(seq):
  seen = set(); out = []
  for x in seq:
    if x not in seen:
      out.append(x); seen.add(x)
  return out

def build_channel_systematics(channel, skip, scales, scalesCorr, scalesGlobal, smears):
  if skip:
    return {
      'scales'      : scales,
      'scalesCorr'  : scalesCorr,
      'scalesGlobal': scalesGlobal,
      'smears'      : smears,
    }

  base_scales      = ['PhotonScale']
  base_scalesCorr  = ['FNUF','Material']
  base_scalesGlobal= _to_list(scalesGlobal)  # 保持使用者自定義，不做預設
  base_smears      = ['PhotonSmear']

  user_scales     = _to_list(scales)
  user_scalesCorr = _to_list(scalesCorr)
  user_smears     = _to_list(smears)

  ele_scales  = ['ElectronScale']
  ele_smears  = ['ElectronSmear']
  mu_scales   = ['MuonScale']
  mu_smears   = ['MuonSmear']

  ch = (channel or '').lower()
  if ch in ('ele', 'electron', 'e'):
    keep_scales = base_scales + ele_scales + user_scales
    keep_smears = base_smears + ele_smears + user_smears
    # 剔除 muon 專屬
    keep_scales = [x for x in keep_scales if x not in mu_scales]
    keep_smears = [x for x in keep_smears if x not in mu_smears]
  elif ch in ('mu', 'muon', 'm'):
    keep_scales = base_scales + mu_scales + user_scales
    keep_smears = base_smears + mu_smears + user_smears
    # 剔除 electron 專屬
    keep_scales = [x for x in keep_scales if x not in ele_scales]
    keep_smears = [x for x in keep_smears if x not in ele_smears]
  elif ch in ('leptons', 'emu', 'combined'):
    keep_scales = base_scales + ele_scales + mu_scales + user_scales
    keep_smears = base_smears + ele_smears + mu_smears + user_smears
  else:
    # 未知 channel：僅合併使用者輸入與 photon 預設
    keep_scales = base_scales + user_scales
    keep_smears = base_smears + user_smears

  keep_scales     = _uniq(keep_scales)
  keep_scalesCorr = _uniq(base_scalesCorr + user_scalesCorr)
  keep_scalesGlob = _uniq(base_scalesGlobal)  # 僅去重
  keep_smears     = _uniq(keep_smears)

  return {
    'scales'      : ','.join(keep_scales),
    'scalesCorr'  : ','.join(keep_scalesCorr),
    'scalesGlobal': ','.join(keep_scalesGlob),
    'smears'      : ','.join(keep_smears),
  }

def get_options():
  parser = OptionParser()
  parser.add_option('--mass_ALP', dest='mass_ALP', default=1, type='int', help="ALP mass") # PZ
  parser.add_option("--channel", dest='channel', default='', help="ele, mu, or leptons") # PZ

  parser.add_option("--xvar", dest='xvar', default='CMS_hza_mass', help="Observable to fit")
  parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  parser.add_option("--proc", dest='proc', default='GG2H', help="Signal process") # PZ
  parser.add_option("--cat", dest='cat', default='cat0', help="RECO category") # PZ
  parser.add_option("--year", dest='year', default='2016', help="Year")
  parser.add_option("--analysis", dest='analysis', default='HZa', help="Analysis handle: used to specify replacement map and XS*BR normalisations")
  parser.add_option('--massPoints', dest='massPoints', default='125', help="Mass points to fit")
  parser.add_option('--skipBeamspotReweigh', dest='skipBeamspotReweigh', default=False, action="store_true", help="Skip beamspot reweigh to match beamspot distribution in data") # PZ
  parser.add_option('--doPlots', dest='doPlots', default=True, action="store_true", help="Produce Signal Fitting plots") # PZ
  parser.add_option("--doVoigtian", dest='doVoigtian', default=False, action="store_true", help="Use Voigtians instead of Gaussians for signal models with Higgs width as parameter")
  parser.add_option("--useDCB", dest='useDCB', default=False, action="store_true", help="Use DCB in signal fitting")
  parser.add_option("--useDiagonalProcForShape", dest='useDiagonalProcForShape', default=False, action="store_true", help="Use shape of diagonal process, keeping normalisation (requires diagonal mapping produced by getDiagProc script)")
  parser.add_option('--skipVertexScenarioSplit', dest='skipVertexScenarioSplit', default=True, action="store_true", help="Skip vertex scenario split") # PZ
  parser.add_option('--skipZeroes', dest='skipZeroes', default=False, action="store_true", help="Skip proc x cat is numEntries = 0., or sumEntries < 0.")
  # For systematics
  parser.add_option('--skipSystematics', dest='skipSystematics', default=True, action="store_false", help="Skip shape systematics in signal model") # PZ
  # 新增一個直觀的別名旗標（出現即啟用系統誤差），與舊旗標相容
  parser.add_option('--doSystematics', dest='doSystematics', default=False, action="store_true", help="Enable shape systematics in signal model (alias)")
  parser.add_option('--useDiagonalProcForSyst', dest='useDiagonalProcForSyst', default=False, action="store_true", help="Use diagonal process for systematics (requires diagonal mapping produced by getDiagProc script)")
  parser.add_option("--scales", dest='scales', default='', help="Photon shape systematics: scales")
  parser.add_option("--scalesCorr", dest='scalesCorr', default='', help='Photon shape systematics: scalesCorr')
  parser.add_option("--scalesGlobal", dest='scalesGlobal', default='', help='Photon shape systematics: scalesGlobal')
  parser.add_option("--smears", dest='smears', default='', help='Photon shape systematics: smears')
  # Parameter values
  parser.add_option('--replacementThreshold', dest='replacementThreshold', default=20, type='int', help="Nevent threshold to trigger replacement dataset")
  parser.add_option('--beamspotWidthData', dest='beamspotWidthData', default=3.5, type='float', help="Width of beamspot in data [cm]")
  parser.add_option('--beamspotWidthMC', dest='beamspotWidthMC', default=5.14, type='float', help="Width of beamspot in MC [cm]") # PZ
  parser.add_option('--MHPolyOrder', dest='MHPolyOrder', default=1, type='int', help="Order of polynomial for MH dependence")
  parser.add_option('--nBins', dest='nBins', default=80, type='int', help="Number of bins for fit")
  # Minimizer options
  parser.add_option('--minimizerMethod', dest='minimizerMethod', default='TNC', help="(Scipy) Minimizer method")
  parser.add_option('--minimizerTolerance', dest='minimizerTolerance', default=1e-8, type='float', help="(Scipy) Minimizer toleranve")
  return parser.parse_args()
(opt,args) = get_options()

# Low ma: switch to DCB+Gaussian to handle the photon-merging / brem tail.
opt.useDCB = _force_dcb_for_low_ma(opt.mass_ALP, opt.useDCB)
print(f" [CFG] useDCB resolved to {opt.useDCB} for mA={opt.mass_ALP}")

# 解析與正規化最終的系統誤差開關，兼容 --skipSystematics 與 --doSystematics
# 規則：有 --doSystematics 或（有 --skipSystematics 旗標）=> 啟用系統
#      其他情況 => 依預設（目前為停用系統）
doSystematics = bool(getattr(opt, 'doSystematics', False) or (not opt.skipSystematics))
# 回填以維持 downstream 參考 opt.skipSystematics 的相容性
opt.skipSystematics = (not doSystematics)

# 清楚列印設定，避免混淆
print(f" [CFG] Systematics: {'ENABLED' if doSystematics else 'DISABLED'} (opt.skipSystematics={opt.skipSystematics})")

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

# 用 channel 計算系統誤差（若 skipSystematics 為 True 則不變）
_sys = build_channel_systematics(opt.channel, opt.skipSystematics, opt.scales, opt.scalesCorr, opt.scalesGlobal, opt.smears)
if doSystematics:
  print(" [INFO] 使用之系統誤差 (channel=%s)" % opt.channel)
  print("        scales      :", _sys['scales'] or "(none)")
  print("        scalesCorr  :", _sys['scalesCorr'] or "(none)")
  print("        scalesGlobal:", _sys['scalesGlobal'] or "(none)")
  print("        smears      :", _sys['smears'] or "(none)")

# 保留原本的 channel-based 覆寫，但改用 doSystematics 判斷，且僅印訊息
if doSystematics:
  if opt.channel == 'ele':
    _forced = {
      'scales'     : ['PhotonScale','ElectronScale'],
      'scalesCorr' : ['FNUF','Material'],
      'smears'     : ['PhotonSmear','ElectronSmear'],
    }
    print(" [INFO] channel=ele -> 覆寫系統列表:", _forced)
  elif opt.channel == 'mu':
    _forced = {
      'scales'     : ['PhotonScale','MuonScale'],
      'scalesCorr' : ['FNUF','Material'],
      'smears'     : ['PhotonSmear','MuonSmear'],
    }
    print(" [INFO] channel=mu -> 覆寫系統列表:", _forced)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETUP: signal fit
print(" --> Running fit for (proc,cat) = (%s,%s)"%(opt.proc,opt.cat))
if( len(opt.massPoints.split(",")) == 1 )&( opt.MHPolyOrder > 0 ):
  print(" --> [WARNING] Attempting to fit polynomials of O(MH^%g) for single mass point. Setting order to 0"%opt.MHPolyOrder)
  opt.MHPolyOrder=0

# Add stopwatch function

# Load replacement map
if opt.analysis not in globalReplacementMap:
  print(" --> [ERROR] replacement map does not exist for analysis (%s). Please add to tools/replacementMap.py"%opt.analysis)
  leave()
else: rMap = globalReplacementMap[opt.analysis]

# Load XSBR map
if opt.analysis not in globalXSBRMap:
  print(" --> [ERROR] XS * BR map does not exist for analysis (%s). Please add to tools/XSBRMap.py"%opt.analysis)
  leave()
else: xsbrMap = globalXSBRMap[opt.analysis]

# Load RooRealVars
print(f"{opt.inputWSDir}/ws_{opt.channel}_{opt.year}.root")
nominalWSFileName = glob.glob(f"{opt.inputWSDir}/ws_{opt.channel}_{opt.year}.root")[0] # PZ
f0 = ROOT.TFile(nominalWSFileName,"read")
inputWS0 = get_workspace(f0, inputWSName__)
xvar = inputWS0.var(opt.xvar)
# 若變數不存在則提示可能的變數名稱
if not xvar:
  print(f" --> [ERROR] Variable '{opt.xvar}' not found in workspace '{inputWS0.GetName()}'.")
  try:
    var_iter = inputWS0.allVars().createIterator()
    hints = []
    v = var_iter.Next()
    while v:
      n = v.GetName()
      if "mass" in n or "CMS" in n:
        hints.append(n)
      v = var_iter.Next()
    if hints:
      print("     Available candidates containing 'mass' or 'CMS': %s" % ", ".join(sorted(set(hints))[:15]))
  except Exception:
    pass
  leave()
xvarFit = xvar.Clone()
dZ = ROOT.RooRealVar("dZ", "dZ", 0)  # PZ
# dZ = inputWS0.var("dZ")
aset = ROOT.RooArgSet(xvar,dZ)
f0.Close()
# official sample structure
# inputWSName__ = "tagsDumper/cms_hgg_13TeV"
# procToData(proc.split("_")[0]) = ggh
# RooDataHist::ggh_120_13TeV_EEEB_highR9highR9_SmearingDown01sigma(CMS_hgg_mass)

# PZ sample structure
# inputWSName__ = "CMS_hza_workspace"
# RooDataSet::ggh_125_13TeV_cat0(CMS_hza_mass)

# Mass-dependent fit range: low ma suffers from low-side photon-merging tail,
# nGauss inflates sigma to cover the tail and σ_eff explodes.
# Tighten the fit window for low ma so the core resolution is recovered.
MHLow, MHHigh = _resolve_fit_range(opt.mass_ALP, MHLow, MHHigh)
print(f" [CFG] Fit range for mA={opt.mass_ALP}: [{MHLow}, {MHHigh}]")

# Create MH var
MH = ROOT.RooRealVar("MH","m_{H}", int(MHLow), int(MHHigh))
MH.setUnit("GeV")
MH.setConstant(True)
# Constrain the observable to the chosen fit window before reduceDataset cuts events.
try:
    xvar.setRange(int(MHLow), int(MHHigh))
    xvarFit.setRange(int(MHLow), int(MHHigh))
except Exception:
    pass

if opt.skipZeroes:
  # Extract nominal mass dataset and see if entries == 0
  WSFileName = glob.glob(f"{opt.inputWSDir}/ws_{opt.channel}_{opt.year}.root")[0]
  f = ROOT.TFile(WSFileName,"read")
  inputWS = get_workspace(f, inputWSName__)
  dname = "%s_%s_Za_%s_%s_%s"%(procToData(opt.proc.split("_")[0]),MHNominal,opt.channel,sqrts__,opt.cat)
  d = reduceDataset(get_dataset(inputWS, dname), aset)
  if (d.numEntries() == 0.) or (d.sumEntries() <= 0.):
    print(" --> (%s,%s) has zero events. Will not construct signal model"%(opt.proc,opt.cat))
    exit()
  inputWS.Delete()
  f.Close()
 
# Define proc x cat with which to extract shape: if skipVertexScenarioSplit label all events as "RV"
procRVFit, catRVFit = opt.proc, opt.cat
if opt.skipVertexScenarioSplit: 
  print(" --> Skipping vertex scenario split")
else:
  procWVFit, catWVFit = opt.proc, opt.cat

# Options for using diagonal process from getDiagProc output json
if opt.useDiagonalProcForShape:
  if not os.path.exists("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext)):
    print(" --> [ERROR] Diagonal process json from getDiagProc does not exist. Using nominal proc x cat for shape")
  else:
    with open("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext),"r") as jf: dproc = json.load(jf)
    procRVFit = dproc[opt.cat]
    print(" --> Using diagonal proc (%s,%s) for shape"%(procRVFit,opt.cat))
    if not opt.skipVertexScenarioSplit: procWVFit = dproc[opt.cat]

# Process for syst
procSyst = opt.proc
if opt.useDiagonalProcForSyst:
  if not os.path.exists("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext)):
    print(" --> [ERROR] Diagonal process json from getDiagProc does not exist. Using nominal proc x cat for systematics")
  else:
    with open("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext),"r") as jf: dproc = json.load(jf)
    procSyst = dproc[opt.cat]
    print(" --> Using diagonal proc (%s,%s) for systematics"%(procSyst,opt.cat))

# Define process with which to extract normalisation: nominal
procNorm, catNorm = opt.proc, opt.cat

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT DATASETS TO FIT (for each mass point)
nominalDatasets = od()
# For RV (or if skipping vertex scenario split)
datasetRVForFit = od()
for mp in opt.massPoints.split(","):
  WSFileName = glob.glob(f"{opt.inputWSDir}/ws_{opt.channel}_{opt.year}.root")[0] # PZ
  f = ROOT.TFile(WSFileName,"read")
  inputWS = get_workspace(f, inputWSName__)
  dname = "%s_%s_Za_%s_%s_%s"%(procToData(procRVFit.split("_")[0]),mp,opt.channel,sqrts__,opt.cat)  # PZ
  d = reduceDataset(get_dataset(inputWS, dname), aset)
  nominalDatasets[mp] = d.Clone()
  if opt.skipVertexScenarioSplit: datasetRVForFit[mp] = d
  else: datasetRVForFit[mp] = splitRVWV(d,aset,mode="RV")
  inputWS.Delete()
  f.Close()

# Check if nominal yield > threshold (or if +ve sum of weights). If not then use replacement proc x cat
if( datasetRVForFit[MHNominal].numEntries() < opt.replacementThreshold  )|( datasetRVForFit[MHNominal].sumEntries() < 0. ):
  nominal_numEntries = datasetRVForFit[MHNominal].numEntries()
  procReplacementFit, catReplacementFit = rMap['procRVMap'][opt.cat], rMap['catRVMap'][opt.cat]
  for mp in opt.massPoints.split(","):
    WSFileName = glob.glob(f"{opt.inputWSDir}/ws_{opt.channel}_{opt.year}.root")[0] # PZ
    f = ROOT.TFile(WSFileName,"read")
    inputWS = get_workspace(f, inputWSName__)
    dname = "%s_%s_Za_%s_%s_%s"%(procToData(procReplacementFit.split("_")[0]),mp,opt.channel,sqrts__,opt.cat)  # PZ
    d = reduceDataset(get_dataset(inputWS, dname), aset)
    if opt.skipVertexScenarioSplit: datasetRVForFit[mp] = d
    else: datasetRVForFit[mp] = splitRVWV(d,aset,mode="RV")
    inputWS.Delete()
    f.Close() 

  # Check if replacement dataset has too few entries: if so throw error
  if( datasetRVForFit[MHNominal].numEntries() < opt.replacementThreshold )|( datasetRVForFit[MHNominal].sumEntries() < 0. ):
    print(" --> [ERROR] replacement dataset (%s,%s) has too few entries (%g < %g)"%(procReplacementFit,catReplacementFit,datasetRVForFit[MHNominal].numEntries(),opt.replacementThreshold))
    sys.exit(1)

  else:
    procRVFit, catRVFit = procReplacementFit, catReplacementFit
    if opt.skipVertexScenarioSplit: 
      print(" --> Too few entries in nominal dataset (%g < %g). Using replacement (proc,cat) = (%s,%s) for extracting shape"%(nominal_numEntries,opt.replacementThreshold,procRVFit,catRVFit))
      for mp in opt.massPoints.split(","):
        print("     * MH = %s GeV: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit[mp].numEntries(),datasetRVForFit[mp].sumEntries()))
    else: 
      print(" --> RV: Too few entries in nominal dataset (%g < %g). Using replacement (proc,cat) = (%s,%s) for extracting shape"%(nominal_numEntries,opt.replacementThreshold,procRVFit,catRVFit))
      for mp in opt.massPoints.split(","):
        print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit[mp].numEntries(),datasetRVForFit[mp].sumEntries()))

else:
  if opt.skipVertexScenarioSplit: 
    print(" --> Using (proc,cat) = (%s,%s) for extracting shape"%(procRVFit,catRVFit))
    for mp in opt.massPoints.split(","):
      print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit[mp].numEntries(),datasetRVForFit[mp].sumEntries()))
  else: 
    print(" --> RV: Using (proc,cat) = (%s,%s) for extracting shape"%(procRVFit,catRVFit))
    for mp in opt.massPoints.split(","):
      print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit[mp].numEntries(),datasetRVForFit[mp].sumEntries()))

# Repeat for WV scenario
if not opt.skipVertexScenarioSplit:
  datasetWVForFit = od()
  for mp in opt.massPoints.split(","):
    WSFileName = glob.glob(f"{opt.inputWSDir}/ws_{opt.channel}_{opt.year}.root")[0] # PZ
    f = ROOT.TFile(WSFileName,"read")
    inputWS = get_workspace(f, inputWSName__)
    dname = "%s_%s_Za_%s_%s_%s"%(procToData(proc.split("_")[0]),mp,opt.channel,sqrts__,opt.cat)  # PZ
    d = reduceDataset(get_dataset(inputWS, dname), aset)
    datasetWVForFit[mp] = splitRVWV(d,aset,mode="WV")
    inputWS.Delete()
    f.Close()

  # Check nominal mass dataset
  if( datasetWVForFit[MHNominal].numEntries() < opt.replacementThreshold  )|( datasetWVForFit[MHNominal].sumEntries() < 0. ):
    nominal_numEntries = datasetWVForFit[MHNominal].numEntries()
    procReplacementFit, catReplacementFit = rMap['procWV'], rMap['catWV']
    for mp in opt.massPoints.split(","):
      WSFileName = glob.glob(f"{opt.inputWSDir}/ws_{opt.channel}_{opt.year}.root")[0] # PZ
      f = ROOT.TFile(WSFileName,"read")
      inputWS = get_workspace(f, inputWSName__)
      dname = "%s_%s_Za_%s_%s_%s"%(procToData(proc.split("_")[0]),MHNominal,opt.channel,mp,sqrts__,opt.cat)  # PZ
      d = reduceDataset(get_dataset(inputWS, dname), aset)
      datasetWVForFit[mp] = splitRVWV(d,aset,mode="WV")
      inputWS.Delete()
      f.Close()
    # Check if replacement dataset has too few entries: if so throw error
    if( datasetWVForFit[MHNominal].numEntries() < opt.replacementThreshold )|( datasetWVForFit[MHNominal].sumEntries() < 0. ):
      print(" --> [ERROR] replacement dataset (%s,%s) has too few entries (%g < %g)"%(procReplacementFit,catReplacementFit,datasetWVForFit[MHNominal].numEntries(),opt.replacementThreshold))
      sys.exit(1)
    else:
      procWVFit, catWVFit = procReplacementFit, catReplacementFit
      print(" --> WV: Too few entries in nominal dataset (%g < %g). Using replacement (proc,cat) = (%s,%s) for extracting shape"%(nominal_numEntries,opt.replacementThreshold,procWVFit,catWVFit))
      for mp in opt.massPoints.split(","):
        print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetWVForFit[mp].numEntries(),datasetWVForFit[mp].sumEntries()))
  else:
    print(" --> WV: Using (proc,cat) = (%s,%s) for extracting shape"%(procWVFit,catRVFit))
    for mp in opt.massPoints.split(","):
      print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetWVForFit[mp].numEntries(),datasetWVForFit[mp].sumEntries()))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEAMSPOT REWEIGH
if not opt.skipBeamspotReweigh:
  # Datasets for fit
  for mp,d in datasetRVForFit.items(): 
    drw = beamspotReweigh(datasetRVForFit[mp],opt.beamspotWidthData,opt.beamspotWidthMC,xvar,dZ,_x=opt.xvar)
    datasetRVForFit[mp] = drw
  if not opt.skipVertexScenarioSplit:
    for mp,d in datasetWVForFit.items(): 
      drw = beamspotReweigh(datasetWVForFit[mp],opt.beamspotWidthData,opt.beamspotWidthMC,xvar,dZ,_x=opt.xvar)
      datasetWVForFit[mp] = drw
    print(" --> Beamspot reweigh: RV(sumEntries) = %.6f, WV(sumEntries) = %.6f"%(datasetRVForFit[mp].sumEntries(),datasetWVForFit[mp].sumEntries()))
  else:
    print(" --> Beamspot reweigh: sumEntries = %.6f"%datasetRVForFit[mp].sumEntries())

  # Nominal datasets for saving to output Workspace: preserve norm for eff * acc calculation
  for mp,d in nominalDatasets.items():
    drw = beamspotReweigh(d,opt.beamspotWidthData,opt.beamspotWidthMC,xvar,dZ,_x=opt.xvar,preserveNorm=True)
    nominalDatasets[mp] = drw

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If using nGaussian fit then extract nGaussians from fTest json file
if not opt.useDCB:
  # with open("%s/outdir_%s/fTest/json/nGauss_%s.json"%(swd__,opt.ext,catRVFit)) as jf: ngauss = json.load(jf)
  with open(f"{swd__}/outdir_{opt.channel}/fTest/json/{opt.mass_ALP}_nGauss_{opt.year}_{opt.channel}_Hm125.json") as jf: ngauss = json.load(jf)
  nRV = int(ngauss["%s__%s"%(procRVFit,catRVFit)]['nRV'])
  if opt.skipVertexScenarioSplit: print(" --> Fitting function: convolution of nGaussians (%g)"%nRV)
  else: 
    with open(f"{swd__}/outdir_{opt.channel}/fTest/json/{opt.mass_ALP}_nGauss_{opt.year}_{opt.channel}_Hm125.json") as jf: ngauss = json.load(jf)
    nWV = int(ngauss["%s__%s"%(procWVFit,catWVFit)]['nWV'])
    print(" --> Fitting function: convolution of nGaussians (RV=%g,WV=%g)"%(nRV,nWV))
else:
  print(" --> Fitting function: DCB + 1 Gaussian")

if opt.doVoigtian:
  print(" --> Will add natural Higgs width as parameter in Pdf (Gaussians -> Voigtians)")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIT: simultaneous signal fit (ssf)
ssfMap = od()
name = "Total" if opt.skipVertexScenarioSplit else "RV"
ssfRV = SimultaneousFit(name,opt.proc,opt.cat,datasetRVForFit,xvar.Clone(),MH,MHLow,MHHigh,opt.massPoints,opt.nBins,opt.MHPolyOrder,opt.minimizerMethod,opt.minimizerTolerance)
if opt.useDCB: ssfRV.buildDCBplusGaussian()
else: ssfRV.buildNGaussians(nRV)
ssfRV.runFit()
ssfRV.buildSplines()
ssfMap[name] = ssfRV

if not opt.skipVertexScenarioSplit:
  name = "WV"
  ssfWV = SimultaneousFit(name,opt.proc,opt.cat,datasetWVForFit,xvar.Clone(),MH,MHLow,MHHigh,opt.massPoints,opt.nBins,opt.MHPolyOrder,opt.minimizerMethod,opt.minimizerTolerance)
  if opt.useDCB: ssfWV.buildDCBplusGaussian()
  else: ssfWV.buildNGaussians(nWV)
  ssfWV.runFit()
  ssfWV.buildSplines()
  ssfMap[name] = ssfWV

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FINAL MODEL: construction
print("\n --> Constructing final model")
fm = FinalModel(
  ssfMap,opt.proc,opt.cat,opt.ext,opt.year,sqrts__,
  nominalDatasets,xvar,MH,MHLow,MHHigh,opt.massPoints,
  xsbrMap,procSyst,
  _sys['scales'], _sys['scalesCorr'], _sys['scalesGlobal'], _sys['smears'],
  opt.doVoigtian,opt.useDCB,opt.skipVertexScenarioSplit,opt.skipSystematics,
  channel=opt.channel, mass_ALP=opt.mass_ALP
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE: to output workspace
foutDir = "%s/outdir_%s/signalFit/output"%(swd__,opt.channel)
os.makedirs(foutDir, exist_ok=True)
foutName = f"{swd__}/outdir_{opt.channel}/signalFit/output/{opt.mass_ALP}_CMS-HGG_sigfit_{opt.year}_{opt.channel}_Hm125.root"
print("\n --> Saving output workspace to file: %s"%foutName)
if not os.path.isdir(foutDir): os.system("mkdir %s"%foutDir)
fout = ROOT.TFile(foutName,"RECREATE")
outWS = ROOT.RooWorkspace("%s_%s"%(outputWSName__,sqrts__),"%s_%s"%(outputWSName__,sqrts__))
fm.save(outWS)
outWS.Write()
fout.Close()
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTTING
if opt.doPlots:
  print("\n --> Making plots...")
  if not os.path.isdir("%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel)): os.system("mkdir %s/outdir_%s/signalFit/Plots"%(swd__,opt.channel))
  # Plot failures are cosmetic; the workspace has already been written above.
  def _safe(fn, *a, **kw):
    try:
      fn(*a, **kw)
    except Exception as e:
      print(f"   [WARN] Plot step '{fn.__name__}' failed: {e}; continuing.")
  if opt.skipVertexScenarioSplit:
    _safe(plotPdfComponents, ssfRV,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel),_extension="total_",_proc=procRVFit,_cat=catRVFit, _Amass=opt.mass_ALP,_year=opt.year,_channel=opt.channel)
  if not opt.skipVertexScenarioSplit:
    _safe(plotPdfComponents, ssfRV,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel),_extension="RV_",_proc=procRVFit,_cat=catRVFit, _Amass=opt.mass_ALP,_year=opt.year,_channel=opt.channel)
    _safe(plotPdfComponents, ssfWV,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel),_extension="WV_",_proc=procWVFit,_cat=catRVFit, _Amass=opt.mass_ALP,_year=opt.year,_channel=opt.channel)
  # Plot interpolation
  _safe(plotInterpolation, fm,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel), _Amass=opt.mass_ALP,_year=opt.year,_channel=opt.channel)
  _safe(plotSplines, fm,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel),_nominalMass=MHNominal, _Amass=opt.mass_ALP,_year=opt.year,_channel=opt.channel)
