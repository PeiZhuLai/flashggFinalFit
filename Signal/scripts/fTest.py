# Python script to perform signal modelling fTest: extract number of gaussians for final fit
# * run per category over single mass point

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG SIGNAL FTEST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
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
from simultaneousFit import *
from plottingTools import *

MHLow, MHHigh = '100', '180'

def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG SIGNAL FTEST (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit(0)

def get_options():
  parser = OptionParser()
  
  parser.add_option('--mass_ALP', dest='mass_ALP', default='1', type='string', help="ALP mass") # PZ
  parser.add_option('--year', dest='year', default='16', help="year") # PZ
  parser.add_option("--channel", dest='channel', default='', help="ele, mu, or leptons") # PZ

  parser.add_option("--xvar", dest='xvar', default='CMS_hza_mass', help="Observable to fit")
  parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  parser.add_option("--procs", dest='procs', default='GG2H', help="Signal processes") # PZ
  parser.add_option("--nProcsToFTest", dest='nProcsToFTest', default=1, type='int',help="Number of signal processes to fTest (ordered by sum entries), others are set to nRV=1,nWV=1. Set to -1 to run over all")
  parser.add_option("--cat", dest='cat', default='cat0', help="RECO category") # PZ
  parser.add_option('--mass', dest='mass', default='125', help="Mass point to fit") # PZ
  parser.add_option('--doPlots', dest='doPlots', default=True, action="store_true", help="Produce Signal fTest plots") # PZ
  parser.add_option('--nBins', dest='nBins', default=80, type='int', help="Number of bins for fit")
  parser.add_option('--threshold', dest='threshold', default=1, type='int', help="Threshold number of events")
  parser.add_option('--nGaussMax', dest='nGaussMax', default=5, type='int', help="Max number of gaussians to test")
  parser.add_option('--skipWV', dest='skipWV', default=True, action="store_true", help="Skip processing of WV case") # PZ
  # Minimizer options CMS_hza_workspace
  parser.add_option('--minimizerMethod', dest='minimizerMethod', default='TNC', help="(Scipy) Minimizer method")
  parser.add_option('--minimizerTolerance', dest='minimizerTolerance', default=1e-8, type='float', help="(Scipy) Minimizer toleranve")
  return parser.parse_args()
(opt,args) = get_options()

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

def get_workspace(rootfile, name_hint=None):
  """Return a RooWorkspace from rootfile; use name_hint if provided, else first RooWorkspace found."""
  obj = None
  if name_hint:
    obj = rootfile.Get(name_hint)
    # if hint points to a directory, search inside
    if obj and obj.InheritsFrom("TDirectory"):
      ws = _find_workspace_recursive(obj)
      if ws:
        return ws
  if obj and obj.InheritsFrom("RooWorkspace"):
    return obj
  # fallback: first RooWorkspace anywhere
  ws = _find_workspace_recursive(rootfile)
  return ws

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
if opt.doPlots: 
  # if not os.path.isdir("%s/outdir_%s/fTest/Plots"%(swd__,opt.ext)): os.system("mkdir -p %s/outdir_%s/fTest/Plots"%(swd__,opt.ext))
  if not os.path.isdir("%s/outdir_%s/fTest/Plots"%(swd__,opt.channel)): os.system("mkdir -p %s/outdir_%s/fTest/Plots"%(swd__,opt.channel))

# Load xvar to fit (robust)
print(f"Looking for files in: {opt.inputWSDir}")
ws_file_for_init = resolve_ws_file(opt.inputWSDir, opt.channel, opt.year)
f0 = ROOT.TFile(ws_file_for_init, "read")
inputWS0 = get_workspace(f0, inputWSName__)
if not inputWS0:
  print(f"Error: Could not find RooWorkspace (hint='{inputWSName__}') in file: {ws_file_for_init}")
  leave()
xvar = inputWS0.var(opt.xvar)
if not xvar:
  print(f"Error: Workspace does not contain variable '{opt.xvar}'")
  leave()
xvarFit = xvar.Clone()
# prefer workspace dZ if present, else create dummy
dZ_ws = inputWS0.var("dZ")
if dZ_ws:
  dZ = dZ_ws
  aset = ROOT.RooArgSet(xvar, dZ)
else:
  dZ = ROOT.RooRealVar("dZ", "dZ", 0)
  aset = ROOT.RooArgSet(xvar)  # do not include non-existing var
f0.Close()

# Create MH var
MH = ROOT.RooRealVar("MH","m_{H}", int(MHLow), int(MHHigh))
MH.setUnit("GeV")
MH.setConstant(True)

# Loop over processes: extract sum entries and fill dict. Default nRV,nWV = 1,1
df = pd.DataFrame(columns=['proc','sumEntries','nRV','nWV'])
procYields = od()
for proc in opt.procs.split(","):
  print(f"Looking for files in: {opt.inputWSDir}")
  WSFileName = resolve_ws_file(opt.inputWSDir, opt.channel, opt.year)
  f = ROOT.TFile(WSFileName,"read")
  inputWS = get_workspace(f, inputWSName__)
  if not inputWS:
    print(f"Error: Could not find RooWorkspace (hint='{inputWSName__}') in file: {WSFileName}")
    f.Close()
    continue
  d = reduceDataset(inputWS.data("%s_%s_Za_%s_%s_%s"%(procToData(proc.split("_")[0]),opt.mass,opt.channel,sqrts__,opt.cat)),aset)
  # d = d.reduce(aset, "abs(dZ) <= 1.")  # PZ
  df.loc[len(df)] = [proc,d.sumEntries(),1,1]
  inputWS.Delete()
  f.Close()
# official sample structure
# inputWSName__ = "tagsDumper/cms_hgg_13TeV"
# procToData(proc.split("_")[0]) = ggh
# RooDataHist::ggh_120_13TeV_EEEB_highR9highR9_SmearingDown01sigma(CMS_hgg_mass)

# PZ sample structure
# inputWSName__ = "CMS_hza_workspace"
# RooDataSet::ggh_125_13TeV_cat0(CMS_hza_mass)


# Extract processes to perform fTest (i.e. first nProcsToFTest):
if( opt.nProcsToFTest == -1)|( opt.nProcsToFTest > len(opt.procs.split(",")) ): procsToFTest = opt.procs.split(",")
else: procsToFTest = list(df.sort_values('sumEntries',ascending=False)[0:opt.nProcsToFTest].proc.values)
for pidx, proc in enumerate(procsToFTest): 

  print("\n --> Process (%g): %s"%(pidx,proc))

  # Split dataset to RV/WV: ssf requires input as dict (with mass point as key)
  datasets_RV, datasets_WV = od(), od()
  WSFileName = resolve_ws_file(opt.inputWSDir, opt.channel, opt.year) # PZ
  f = ROOT.TFile(WSFileName,"read")
  inputWS = get_workspace(f, inputWSName__)
  if not inputWS:
    print(f"Error: Could not find RooWorkspace (hint='{inputWSName__}') in file: {WSFileName}")
    f.Close()
    continue
  d = reduceDataset(inputWS.data("%s_%s_Za_%s_%s_%s"%(procToData(proc.split("_")[0]),opt.mass,opt.channel,sqrts__,opt.cat)),aset) # PZ
  # datasets_RV[opt.mass] = splitRVWV(d,aset,mode="RV")
  # datasets_WV[opt.mass] = splitRVWV(d,aset,mode="WV")
  datasets_RV[opt.mass] = d #PZ
  datasets_WV[opt.mass] = d #PZ

  # Run fTest: RV
  # If numEntries below threshold then keep as n = 1
  if datasets_RV[opt.mass].numEntries() < opt.threshold: continue  
  else:
    ssfs = od()
    min_reduced_chi2, nGauss_opt = 999, 1
    for nGauss in range(1,opt.nGaussMax+1):
      k = "nGauss_%g"%nGauss
      ssf = SimultaneousFit("fTest_RV_%g"%nGauss,proc,opt.cat,datasets_RV,xvar.Clone(),MH,MHLow,MHHigh,opt.mass,opt.nBins,0,opt.minimizerMethod,opt.minimizerTolerance,verbose=False)
      ssf.buildNGaussians(nGauss)
      ssf.runFit()
      ssf.buildSplines()
      if ssf.Ndof >= 1: 
        ssfs[k] = ssf
        if ssfs[k].getReducedChi2() < min_reduced_chi2: 
          min_reduced_chi2 = ssfs[k].getReducedChi2()
          nGauss_opt = nGauss
        print("   * (%s,%s,RV): nGauss = %g, chi^2/n(dof) = %.4f"%(proc,opt.cat,nGauss,ssfs[k].getReducedChi2()))
    # Set optimum
    df.loc[df['proc']==proc,'nRV'] = nGauss_opt
    # Make plots
    if( opt.doPlots )&( len(ssfs.keys())!=0 ):
      plotFTest(ssfs,_opt=nGauss_opt,_outdir="%s/outdir_%s/fTest/Plots"%(swd__,opt.channel),_extension="RV",_proc=proc,_cat=opt.cat,_mass=opt.mass, _Amass=opt.mass_ALP,_year=opt.year,_channel=opt.channel,_Hmass=opt.mass)
      plotFTestResults(ssfs,_opt=nGauss_opt,_outdir="%s/outdir_%s/fTest/Plots"%(swd__,opt.channel),_extension="RV",_proc=proc,_cat=opt.cat,_mass=opt.mass, _Amass=opt.mass_ALP,_year=opt.year,_channel=opt.channel,_Hmass=opt.mass)

  # Run fTest: WV
  # If numEntries below threshold then keep as n = 1
  if( datasets_WV[opt.mass].numEntries() < opt.threshold )|( opt.skipWV ): continue
  else:
    ssfs = od()
    min_reduced_chi2, nGauss_opt = 999, 1
    for nGauss in range(1,opt.nGaussMax+1):
      k = "nGauss_%g"%nGauss
      ssf = SimultaneousFit("fTest_WV_%g"%nGauss,proc,opt.cat,datasets_WV,xvar.Clone(),MH,MHLow,MHHigh,opt.mass,opt.nBins,0,opt.minimizerMethod,opt.minimizerTolerance,verbose=False)
      ssf.buildNGaussians(nGauss)
      ssf.runFit()
      ssf.buildSplines()
      if ssf.Ndof >= 1:
        ssfs[k] = ssf
        if ssfs[k].getReducedChi2() < min_reduced_chi2:
          min_reduced_chi2 = ssfs[k].getReducedChi2()
          nGauss_opt = nGauss
        print("   * (%s,%s,WV): nGauss = %g, chi^2/n(dof) = %.4f"%(proc,opt.cat,nGauss,ssfs[k].getReducedChi2()))
    # Set optimum
    df.loc[df['proc']==proc,'nWV'] = nGauss_opt
    # Make plots
    if( opt.doPlots )&( len(ssfs.keys())!=0 ):
      plotFTest(ssfs,_opt=nGauss_opt,_outdir="%s/outdir_%s/fTest/Plots"%(swd__,opt.channel),_extension="WV",_proc=proc,_cat=opt.cat,_mass=opt.mass)
      plotFTestResults(ssfs,_opt=nGauss_opt,_outdir="%s/outdir_%s/fTest/Plots"%(swd__,opt.channel),_extension="WV",_proc=proc,_cat=opt.cat,_mass=opt.mass)

  # Close ROOT file
  inputWS.Delete()
  f.Close()

# Make output
if not os.path.isdir("%s/outdir_%s/fTest/json"%(swd__,opt.channel)): os.system("mkdir %s/outdir_%s/fTest/json"%(swd__,opt.channel))
ff = open(f"{swd__}/outdir_{opt.channel}/fTest/json/{opt.mass_ALP}_nGauss_{opt.year}_{opt.channel}_Hm{opt.mass}.json","w")
ff.write("{\n")
# Iterate over rows in dataframe: sorted by sumEntries
pitr = 1
for ir,r in df.sort_values('sumEntries',ascending=False).iterrows():
  k = "\"%s__%s\""%(r['proc'],opt.cat)
  ff.write("    %-90s : {\"nRV\":%s,\"nWV\":%s}"%(k,r['nRV'],r['nWV']))
  # Drop comma for last proc
  if pitr == len(df): ff.write("\n")
  else: ff.write(",\n")
  pitr += 1
ff.write("}")
ff.close()
