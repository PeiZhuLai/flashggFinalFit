# Script for making signal model plot
import os, sys
import ROOT
import re, glob
import json
from optparse import OptionParser

from commonTools import *
from commonObjects import *
from tools.plottingTools import *

def get_options():
  parser = OptionParser()
  parser.add_option('--mass_ALP', dest='mass_ALP', default=1, type='int', help="ALP mass") # PZ
  parser.add_option("--channel", dest='channel', default='', help="ele, mu, or leptons") # PZ

  parser.add_option('--procs', dest='procs', default='GG2H', help="Comma separated list of processes to include. all = sum all signal procs")  
  parser.add_option('--years', dest='years', default='16,16APV,17,18', help="Comma separated list of years to include")  
  parser.add_option('--cats', dest='cats', default='cat0', help="Comma separated list of analysis categories to include. all = sum of all categories, wall = weighted sum of categories (requires S/S+B from ./Plots/getCatInfo.py)")
  parser.add_option('--loadCatWeights', dest='loadCatWeights', default='', help="Load S/S+B weights for analysis categories (path to weights json file)")
  parser.add_option('--ext', dest='ext', default='test', help="Extension: defines output dir where signal models are saved")
  parser.add_option("--xvar", dest="xvar", default='CMS_hza_mass:m_{ ll#gamma#gamma}:GeV', help="x-var (name:title:units)")
  parser.add_option("--mass", dest="mass", default='125', help="Mass of datasets")
  parser.add_option("--MH", dest="MH", default='125', help="Higgs mass (for pdf)")
  parser.add_option("--nBins", dest="nBins", default=80, type='int', help="Number of bins")
  parser.add_option("--pdf_nBins", dest="pdf_nBins", default=1600, type='int', help="Number of bins")
  parser.add_option("--threshold", dest="threshold", default=0.001, type='float', help="Threshold to prune process from plot default = 0.1% of total category norm")
  parser.add_option("--translateCats", dest="translateCats", default=None, help="JSON to store cat translations")
  parser.add_option("--translateProcs", dest="translateProcs", default=None, help="JSON to store proc translations")
  parser.add_option("--label", dest="label", default='Simulation', help="CMS Sub-label")
  parser.add_option("--doFWHM", dest="doFWHM", default=True, action='store_true', help="Do FWHM")
  parser.add_option("--input", dest="input", default="", help="Optional: explicit path to input ROOT file (overrides auto path building)")
  return parser.parse_args()
(opt,args) = get_options()

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

def _ws_var(ws, name, required=True):
  v = ws.var(name) if ws else None
  if required and not v:
    raise KeyError(f'Variable "{name}" not found in workspace')
  return v

def _ws_func(ws, name, required=False, context=""):
  f = ws.function(name) if ws else None
  if required and not f:
    msg = f'Function "{name}" not found in workspace'
    if context: msg += f" ({context})"
    raise KeyError(msg)
  if (not required) and (not f):
    if context:
      print(f'[WARN] Missing function "{name}" ({context}) -> skipping')
    else:
      print(f'[WARN] Missing function "{name}" -> skipping')
  return f

def _safe_getval(absreal):
  # RooFit pythonization in some setups prefers an explicit RooArgSet
  if not absreal:
    return 0.0
  try:
    return float(absreal.getVal(ROOT.RooArgSet()))
  except TypeError:
    return float(absreal.getVal())

# Extract input files: for first file extract xvar
inputFiles = od()
citr = 0
if opt.cats in ['all','wall']:
  fs = glob.glob("%s/outdir_%s/CMS-HGG_sigfit_%s_*.root"%(swd__,opt.ext,opt.ext))
  for f in fs:
    cat = re.sub(".root","",f.split("/")[-1].split("_%s_"%opt.ext)[-1])
    inputFiles[cat] = f
    if citr == 0:
      w = ROOT.TFile(f).Get("wsig_13p6TeV")
      xvar = w.var(opt.xvar.split(":")[0])
      xvar.setPlotLabel(opt.xvar.split(":")[1])
      xvar.setUnit(opt.xvar.split(":")[2])
      alist = ROOT.RooArgList(xvar)
    citr += 1
else:
  years_list = [y.strip() for y in opt.years.split(",") if y.strip()]
  years_tag = "_".join(years_list)  # IMPORTANT: avoid commas in filename
  for cat in opt.cats.split(","):
    tried = []

    if opt.input:
      f = opt.input
      tried.append(f)
      if not os.path.isfile(f):
        raise OSError("Failed to open input ROOT file. Tried:\n  - " + "\n  - ".join(tried))
      inputFiles[cat] = f
      if citr == 0:
        fin0 = ROOT.TFile.Open(f)
        if not fin0 or fin0.IsZombie():
          raise OSError(f"ROOT could not open file (zombie): {f}")
        w = fin0.Get("wsig_13p6TeV")
        if not w:
          raise KeyError(f'Workspace "wsig_13p6TeV" not found in file: {f}')
        xvar = w.var(opt.xvar.split(":")[0])
        if not xvar:
          raise KeyError(f'Variable "{opt.xvar.split(":")[0]}" not found in workspace in file: {f}')
        xvar.setPlotLabel(opt.xvar.split(":")[1])
        xvar.setUnit(opt.xvar.split(":")[2])
        alist = ROOT.RooArgList(xvar)
        fin0.Close()
      citr += 1
    else:
      # NEW: build per-year file mapping (year -> file) for this category
      yearFiles = od()
      for y in years_list:
        f_y = f"{swd__}/outdir_{opt.channel}/signalFit/output/{opt.mass_ALP}_CMS-HGG_sigfit_{y}_{opt.channel}_Hm125.root"
        tried.append(f_y)
        if os.path.isfile(f_y):
          yearFiles[y] = f_y

      # Optional fallback: combined-years file if no per-year files found
      if len(yearFiles) == 0:
        f_combined = f"{swd__}/outdir_{opt.channel}/signalFit/output/{opt.mass_ALP}_CMS-HGG_sigfit_{years_tag}_{opt.channel}_Hm125.root"
        tried.append(f_combined)
        if os.path.isfile(f_combined):
          yearFiles["combined"] = f_combined

      if len(yearFiles) == 0:
        raise OSError("Failed to open input ROOT file(s). Tried:\n  - " + "\n  - ".join(tried))

      inputFiles[cat] = yearFiles

      # init xvar from the first available file
      if citr == 0:
        f0 = next(iter(yearFiles.values()))
        fin0 = ROOT.TFile.Open(f0)
        if not fin0 or fin0.IsZombie():
          raise OSError(f"ROOT could not open file (zombie): {f0}")
        w = fin0.Get("wsig_13p6TeV")
        if not w:
          raise KeyError(f'Workspace "wsig_13p6TeV" not found in file: {f0}')
        xvar = w.var(opt.xvar.split(":")[0])
        if not xvar:
          raise KeyError(f'Variable "{opt.xvar.split(":")[0]}" not found in workspace in file: {f0}')
        xvar.setPlotLabel(opt.xvar.split(":")[1])
        xvar.setUnit(opt.xvar.split(":")[2])
        alist = ROOT.RooArgList(xvar)
        fin0.Close()
      citr += 1

# Load cat S/S+B weights
if opt.loadCatWeights != '':
  with open( opt.loadCatWeights ) as jsonfile: catsWeights = json.load(jsonfile)

# Define dict to store data histogram and inclusive + per-year pdf histograms
hists = od()
hists['data'] = xvar.createHistogram("h_data", ROOT.RooFit.Binning(opt.nBins))

# Pre-create per-year pdf histograms to guarantee they are real TH1 objects (never None)
_years_list = [y.strip() for y in opt.years.split(",") if y.strip()]
for year in _years_list:
  hists[f'pdf_{year}'] = xvar.createHistogram(f"h_pdf_{year}", ROOT.RooFit.Binning(opt.pdf_nBins))
  hists[f'pdf_{year}'].Reset()

# NEW: always create combined pdf hist (so plotSignalModel can always use it)
hists['pdf'] = xvar.createHistogram("h_pdf", ROOT.RooFit.Binning(opt.pdf_nBins))
hists['pdf'].Reset()

# Loop over files
for cat, f_or_map in inputFiles.items():
  print(" --> Processing %s: file = %s"%(cat, f_or_map))

  # Define cat weight
  wcat = catsWeights[cat] if opt.loadCatWeights != '' else 1.

  # Normalize to a per-cat year->file map
  if isinstance(f_or_map, (str, bytes)):
    yearFileMap = od()
    for y in _years_list:
      yearFileMap[y] = f_or_map
  else:
    yearFileMap = f_or_map  # already year -> file

  # Containers across years (for this category)
  norms = od()
  data_rwgt = od()
  hpdfs = od()

  # First pass: compute catNorm across years (from their own files)
  catNorm = 0.0
  for year, f in yearFileMap.items():
    fin = ROOT.TFile.Open(f)
    if not fin or fin.IsZombie():
      raise OSError(f"Failed to open file (zombie): {f}")
    w = fin.Get("wsig_13p6TeV")
    if not w:
      raise KeyError(f'Workspace "wsig_13p6TeV" not found in file: {f}')
    _ws_var(w, "MH").setVal(float(opt.MH))
    intLumiVar = _ws_var(w, "IntLumi", required=True)
    intLumiVar.setVal(lumiScaleFactor*lumiMap.get(year, lumiMap.get(str(year), 0.0)))

    if opt.procs == 'all':
      allNorms = w.allFunctions().selectByName(f"*{year}*normThisLumi")
      for norm in rooiter(allNorms):
        proc = norm.GetName().split("%s_"%outputWSObjectTitle__)[-1].split(f"_{year}")[0]
        k  =  "%s__%s"%(proc,year)
        _id = "%s_%s_%s_%s"%(proc,year,cat,sqrts__)
        fname = "%s_%s_normThisLumi"%(outputWSObjectTitle__,_id)
        fn = _ws_func(w, fname, required=False, context=f"cat={cat}, year={year}, proc={proc}")
        if fn:
          norms[k] = fn
          catNorm += _safe_getval(fn)
    else:
      for proc in opt.procs.split(","):
        k = "%s__%s"%(proc,year)
        _id = "%s_%s_%s_%s"%(proc,year,cat,sqrts__)
        fname = "%s_%s_normThisLumi"%(outputWSObjectTitle__,_id)
        fn = _ws_func(w, fname, required=False, context=f"cat={cat}, year={year}, proc={proc}")
        if fn:
          norms[k] = fn
          catNorm += _safe_getval(fn)

    fin.Close()

  if catNorm <= 0 and len(norms) == 0:
    raise RuntimeError(
      f"No valid normThisLumi functions found for cat={cat}. "
      f"Check procs/years and workspace contents in the provided year files."
    )

  # Second pass: fill datasets/pdfs per year from that year's file ONLY
  for year, f in yearFileMap.items():
    fin = ROOT.TFile.Open(f)
    if not fin or fin.IsZombie():
      raise OSError(f"Failed to open file (zombie): {f}")
    w = fin.Get("wsig_13p6TeV")
    if not w:
      raise KeyError(f'Workspace "wsig_13p6TeV" not found in file: {f}')
    _ws_var(w, "MH").setVal(float(opt.MH))
    intLumiVar = _ws_var(w, "IntLumi", required=True)
    intLumiVar.setVal(lumiScaleFactor*lumiMap.get(year, lumiMap.get(str(year), 0.0)))

    # gather procs for this year (only)
    if opt.procs == 'all':
      allNorms = w.allFunctions().selectByName(f"*{year}*normThisLumi")
      procs_this_year = []
      for norm in rooiter(allNorms):
        proc = norm.GetName().split("%s_"%outputWSObjectTitle__)[-1].split(f"_{year}")[0]
        procs_this_year.append(proc)
    else:
      procs_this_year = [p.strip() for p in opt.procs.split(",") if p.strip()]

    for proc in procs_this_year:
      _id = "%s_%s_%s_%s"%(proc,year,cat,sqrts__)
      fname = "%s_%s_normThisLumi"%(outputWSObjectTitle__,_id)
      norm = _ws_func(w, fname, required=False, context=f"cat={cat}, year={year}, proc={proc}")
      if not norm:
        continue

      nval = _safe_getval(norm)
      if catNorm > 0 and nval < opt.threshold*catNorm:
        continue

      d = w.data("sig_mass_m%s_%s"%(opt.mass,_id))
      if not d:
        print(f'[WARN] Missing dataset sig_mass_m{opt.mass}_{_id} -> skipping')
        continue
      d_rwgt = d.emptyClone(_id)

      nf = 0 if d.sumEntries() == 0 else nval/d.sumEntries()
      for i in range(d.numEntries()):
        p = d.get(i)
        rw, rwe = d.weight()*nf*wcat, d.weightError()*nf*wcat
        d_rwgt.add(p,rw,rwe)
      data_rwgt[_id] = d_rwgt

      pdf = w.pdf("extend%s_%sThisLumi"%(outputWSObjectTitle__,_id))
      if not pdf:
        print(f'[WARN] Missing pdf extend{outputWSObjectTitle__}_{_id}ThisLumi -> skipping')
        continue
      hpdfs[_id] = pdf.createHistogram("h_pdf_%s"%_id,xvar,ROOT.RooFit.Binning(opt.pdf_nBins))
      hpdfs[_id].Scale(wcat*float(opt.nBins)/80)

      # per-year accumulation is explicit (no string guessing)
      if f'pdf_{year}' in hists:
        hists[f'pdf_{year}'] += hpdfs[_id]

      # NEW: combined accumulation happens here too
      hists['pdf'] += hpdfs[_id]

    fin.Close()

  # Fill total data histogram (across years)
  for _id,d in data_rwgt.items(): 
    d.fillHistogram(hists['data'],alist)

# Make plot
if not os.path.isdir("%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel)): os.system("mkdir %s/outdir_%s/signalFit/Plots"%(swd__,opt.channel))
plotSignalModel(hists,opt,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.channel), _Amass=opt.mass_ALP,_year=opt.years,_channel=opt.channel)
