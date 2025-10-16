# Script to calculate photon systematics
# * Run script once per category, loops over signal processes
# * Output is pandas dataframe 

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZallgg SYST CALCULATOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
import ROOT
import pandas as pd
import pickle
import os, sys
from optparse import OptionParser
import glob
import re

# From tools
from plottingTools import * #getEffSigma function
from commonTools import *
from commonObjects import *

def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZallgg SYST CALCULATOR (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit(0)

def get_options():
  parser = OptionParser()
  parser.add_option('--mass_ALP', dest='mass_ALP', default=1, type='int', help="ALP mass") # PZ
  parser.add_option('--year', dest='year', default='16', help="year") # PZ
  parser.add_option("--channel", dest='channel', default='', help="ele, mu, or leptons") # PZ

  parser.add_option("--xvar", dest='xvar', default='CMS_hza_mass', help="Observable")
  parser.add_option("--cat", dest='cat', default='cat0', help="RECO category")
  parser.add_option("--procs", dest='procs', default='GG2H', help="Signal processes")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--scales", dest='scales', default='PhotonScale, ElectronScale, MuonPtScale', help="Photon shape systematics: scales")
  parser.add_option("--scalesCorr", dest='scalesCorr', default='FNUF, Material', help='Photon shape systematics: scalesCorr')
  parser.add_option("--scalesGlobal", dest='scalesGlobal', default='', help='Photon shape systematics: scalesGlobal')
  parser.add_option("--smears", dest='smears', default='PhotonSmear, ElectronSmear, MuonPtSmear', help='Photon shape systematics: smears')
  parser.add_option("--nBins", dest='nBins', default=80, type='int', help='Number of bins in histograms')
  parser.add_option("--thresholdMean", dest='thresholdMean', default=0.05, type='float', help='Reject mean variations if larger than thresholdMean')
  parser.add_option("--thresholdSigma", dest='thresholdSigma', default=0.5, type='float', help='Reject mean variations if larger than thresholdSigma')
  parser.add_option("--thresholdRate", dest='thresholdRate', default=0.05, type='float', help='Reject mean variations if larger than thresholdRate')
  parser.add_option("--reportMissing", dest='reportMissing', action='store_true', default=True, help="列出找不到的系統變化")
  parser.add_option("--debugCols", dest='debugCols', action='store_true', default=False, help="列印欄位與系統解析結果")
  parser.add_option("--flatNames", dest='flatNames', action='store_true', default=False, help="欄位名稱不加類型後綴")
  parser.add_option("--debugMissing", dest='debugMissing', action='store_true', default=False, help="列印找不到的系統所有嘗試名稱")
  return parser.parse_args()
(opt,args) = get_options()

# RooRealVar to fill histograms
mgg = ROOT.RooRealVar(opt.xvar,opt.xvar,125)

# --- helpers: robust file and workspace resolution ---
def resolve_ws_file(input_dir, channel, year):
  """Pick a ROOT file under input_dir in order of preference."""
  prefs = [
    f"{input_dir}/ws_{channel}_{year}.root",
    f"{input_dir}/ws_{channel}_*.root",
    f"{input_dir}/ws_*.root",
    f"{input_dir}/*.root",
  ]
  for pat in prefs:
    matches = sorted(glob.glob(pat))
    if matches:
      return matches[0]
  raise RuntimeError(f"No ROOT ws file found under {input_dir}")

def _find_workspace_recursive(dir_obj):
  """Recursively search first RooWorkspace inside a TFile/TDirectory."""
  for key in dir_obj.GetListOfKeys():
    name = key.GetName()
    obj = dir_obj.Get(name)
    if not obj:
      continue
    if obj.InheritsFrom("RooWorkspace"):
      return obj
    if obj.InheritsFrom("TDirectory"):
      ws = _find_workspace_recursive(obj)
      if ws:
        return ws
  return None

def get_workspace(rootfile_or_dir, name_hint=None):
  """Return RooWorkspace object. If name_hint points to a directory, search inside."""
  obj = None
  if name_hint:
    obj = rootfile_or_dir.Get(name_hint)
    if obj and obj.InheritsFrom("RooWorkspace"):
      return obj
    if obj and obj.InheritsFrom("TDirectory"):
      ws = _find_workspace_recursive(obj)
      if ws:
        return ws
  # fallback: first workspace anywhere
  return _find_workspace_recursive(rootfile_or_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to extract histograms from WS (精簡：只嘗試單一標準命名 <sysPrefix>_<systName>)
def getHistogramsFlexible(ws, nominalName, sysPrefix, sname):
  rds_nominal = ws.data(nominalName)
  if not rds_nominal:
    print(f" --> [WARN] Nominal RooDataSet 不存在: {nominalName}")
    return None
  up_name = f"{sysPrefix}_{sname}Up01sigma"
  down_name = f"{sysPrefix}_{sname}Down01sigma"
  rdh_up = ws.data(up_name)
  rdh_down = ws.data(down_name)
  if not (rdh_up and rdh_down):
    if opt.debugMissing:
      print(f"    [MISS] {up_name} 或 {down_name} 不存在")
    return None
  hists = {
    'nominal': ROOT.TH1F("nominal","nominal",opt.nBins,100,180),
    'up'     : ROOT.TH1F(f"{sname}_up",f"{sname}_up",opt.nBins,100,180),
    'down'   : ROOT.TH1F(f"{sname}_down",f"{sname}_down",opt.nBins,100,180),
  }
  rds_nominal.fillHistogram(hists['nominal'],ROOT.RooArgList(mgg))
  rdh_up.fillHistogram(hists['up'],ROOT.RooArgList(mgg))
  rdh_down.fillHistogram(hists['down'],ROOT.RooArgList(mgg))
  return hists

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to extract mean, sigma and rate variations
def getMeanVar(_hists):
  mu, muVar = {}, {}
  for htype,h in _hists.items(): mu[htype] = h.GetMean()
  if mu['nominal']==0: return 0
  for htype in ['up','down']: muVar[htype] = (mu[htype]-mu['nominal'])/mu['nominal']
  x = (abs(muVar['up'])+abs(muVar['down']))/2
  # Check for NaN
  if x!=x: return 0
  else: return min(x,opt.thresholdMean)

def getSigmaVar(_hists):
  sigma, sigmaVar = {}, {}
  for htype,h in _hists.items(): sigma[htype] = getEffSigma(h)
  if sigma['nominal']==0: return 0
  for htype in ['up','down']: sigmaVar[htype] = (sigma[htype]-sigma['nominal'])/sigma['nominal']
  x = (abs(sigmaVar['up'])+abs(sigmaVar['down']))/2
  if x!=x: return 0
  else: return min(x,opt.thresholdSigma)

def getRateVar(_hists):
  rate, rateVar = {}, {}
  for htype,h in _hists.items(): rate[htype] = h.Integral()
  # Shape variations can both be one sided therefore use midpoint as nominal
  rate['midpoint'] = 0.5*(rate['up']+rate['down'])
  if rate['midpoint']==0: return 0
  for htype in ['up','down']: rateVar[htype] = (rate[htype]-rate['midpoint'])/rate['midpoint']
  x = (abs(rateVar['up'])+abs(rateVar['down']))/2
  if x!=x: return 0
  else: return min(x,opt.thresholdRate)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------- 系統名稱最小清洗 (僅去空白與重複) ----------
def _sanitize_syst_list(raw):
  if not raw: return []
  out, seen = [], set()
  for s in raw.split(","):
    s2 = s.strip()
    if not s2 or s2 in seen: continue
    seen.add(s2); out.append(s2)
  return out

# 直接使用使用者輸入 (不做任何映射/白名單/自動補)
sanitized_systs = {
  'scales'     : _sanitize_syst_list(opt.scales),
  'scalesCorr' : _sanitize_syst_list(opt.scalesCorr),
  'smears'     : _sanitize_syst_list(opt.smears),
}
if opt.debugCols:
  print(" [DEBUG] 使用者系統列表 =", sanitized_systs)

# --------- 依 channel 覆寫 (需求：ele / mu) ---------
if opt.channel == 'ele':
  _forced = {
    'scales'     : ['PhotonScale','ElectronScale'],
    'scalesCorr' : ['FNUF','Material'],
    'smears'     : ['PhotonSmear','ElectronSmear'],
  }
  sanitized_systs = _forced
  print(" [INFO] channel=ele -> 覆寫系統列表:", sanitized_systs)
elif opt.channel == 'mu':
  _forced = {
    'scales'     : ['PhotonScale','MuonPtScale'],
    'scalesCorr' : ['FNUF','Material'],
    'smears'     : ['PhotonSmear','MuonPtSmear'],
  }
  sanitized_systs = _forced
  print(" [INFO] channel=mu -> 覆寫系統列表:", sanitized_systs)

# ---------- 欄位命名工具 ----------
def build_col_name(syst, stype, what):
  # what in {mean,sigma,rate}
  if opt.flatNames:
    return f"{syst}_{what}"
  ext = outputNuisanceExtMap[stype]
  ext_part = f"_{ext}" if ext else ""
  return f"{syst}{ext_part}_{what}"

# Define dataFrame (updated: add sysDataPrefix)  (修正：使用清洗後名稱，避免前導空白造成重複欄位)
columns_data = ['proc','cat','inputWSFile','nominalDataName','sysDataPrefix']
for stype, syst_list in sanitized_systs.items():
  for s_clean in syst_list:
    for x in ['mean','sigma','rate']:
      columns_data.append(build_col_name(s_clean, stype, x))
data = pd.DataFrame(columns=columns_data)
if opt.debugCols:
  print(" [DEBUG] Initial columns:", columns_data)

# Loop over processes and add row to dataframe (updated naming)
for _proc in opt.procs.split(","):
  _WSFileName = resolve_ws_file(opt.inputWSDir, opt.channel, opt.year)
  # nominal (無 Za_<channel>)
  _nominalDataName = "%s_125_Za_%s_%s_%s"%(procToData(_proc.split("_")[0]),opt.channel,sqrts__,opt.cat)
  # 系統 (有 Za_<channel>)
  _sysDataPrefix = "%s_125_Za_%s_%s_%s"%(procToData(_proc.split("_")[0]),opt.channel,sqrts__,opt.cat)
  data = pd.concat([data,pd.DataFrame([{
    'proc':_proc,
    'cat':opt.cat,
    'inputWSFile':_WSFileName,
    'nominalDataName':_nominalDataName,
    'sysDataPrefix':_sysDataPrefix
  }])], ignore_index=True, sort=False)

# ---------- 紀錄不存在的系統 ----------
_missing_syst_global = set()

# Loop over rows in dataFrame and open ws
for ir,r in data.iterrows():

  print(" --> Processing (%s,%s)"%(r['proc'],opt.cat))
  f = ROOT.TFile(r['inputWSFile'],"read")
  inputWS = get_workspace(f, inputWSName__)
  if not inputWS:
    print(" --> [ERROR] RooWorkspace not found (hint=%s) in file: %s. Skipping."%(inputWSName__, r['inputWSFile']))
    f.Close()
    continue

  for stype, syst_list in sanitized_systs.items():
    for s_clean in syst_list:

      hists = getHistogramsFlexible(inputWS, r['nominalDataName'], r['sysDataPrefix'], s_clean)

      if not hists:
        _meanVar = _sigmaVar = _rateVar = 0
        _missing_syst_global.add(f"{stype}:{s_clean}")
      else:
        if hists['nominal'].Integral() == 0:
          _meanVar = _sigmaVar = _rateVar = 0
        else:
          _meanVar = getMeanVar(hists)
          _sigmaVar = getSigmaVar(hists)
          _rateVar = getRateVar(hists)
        if opt.debugMissing:
          print(f"    [HIT] {r['sysDataPrefix']}_{s_clean} Up/Down01sigma")
        for h in hists.values(): h.Delete()

      data.at[ir, build_col_name(s_clean, stype, 'mean')]  = _meanVar
      data.at[ir, build_col_name(s_clean, stype, 'sigma')] = _sigmaVar
      data.at[ir, build_col_name(s_clean, stype, 'rate')]  = _rateVar

  inputWS.Delete()
  f.Close()

if opt.debugCols:
  print(" [DEBUG] Final columns:", list(data.columns))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Output dataFrame as pickle file to be read in by signalFit.py
if not os.path.isdir("%s/outdir_%s"%(swd__,opt.channel)): os.system("mkdir %s/outdir_%s"%(swd__,opt.channel))
if not os.path.isdir("%s/outdir_%s/calcPhotonSyst"%(swd__,opt.channel)): os.system("mkdir %s/outdir_%s/calcPhotonSyst"%(swd__,opt.channel))
if not os.path.isdir("%s/outdir_%s/calcPhotonSyst/pkl"%(swd__,opt.channel)): os.system("mkdir %s/outdir_%s/calcPhotonSyst/pkl"%(swd__,opt.channel))
with open("%s/outdir_%s/calcPhotonSyst/pkl/%s_%s.pkl"%(swd__,opt.channel,opt.mass_ALP,opt.year),"wb") as f: pickle.dump(data,f) 
print(" --> Successfully saved photon systematics as pkl file: %s/outdir_%s/calcPhotonSyst/pkl/%s_%s.pkl"%(swd__,opt.channel,opt.mass_ALP,opt.year))

# ---------- 額外：列出未找到的系統 ----------
if opt.reportMissing and _missing_syst_global:
  print(" --> 下列系統在任何資料集中未找到 (已設為 0)：")
  for tag in sorted(_missing_syst_global):
    print("     -", tag)
