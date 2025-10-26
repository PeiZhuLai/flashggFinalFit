#!/usr/bin/env python3

import os, glob, shutil
from biasUtils import *

from optparse import OptionParser
parser = OptionParser()
parser.add_option('--mA', dest='mA', default=5, type='int', help="ALP mass") # PZ
parser.add_option("-d","--datacard",default="Datacard.root")
parser.add_option("-w","--workspace",default="w")
parser.add_option("-t","--toys",action="store_true", default=False)
parser.add_option("-n","--nToys",default=1000,type="int")
parser.add_option("-f","--fits",action="store_true", default=False)
parser.add_option("-p","--plots",action="store_true", default=False)
parser.add_option("-e","--expectSignal",default=1.,type="float")
parser.add_option("-m","--mH",default=125,type="float")
parser.add_option("-c","--combineOptions",default="")
parser.add_option("-s","--seed",default=-1,type="int")
parser.add_option("--dryRun",action="store_true", default=False)
parser.add_option("--poi",default="r")
parser.add_option("--split",default=500,type="int")
parser.add_option("--selectFunction",default=None)
parser.add_option("--gaussianFit",action="store_true", default=False)
# 新增選擇 expectSignal 來源的開關: arg | q50 | auto (預設 arg)
parser.add_option("--expectFrom", default="arg",help="Source of expectSignal: 'arg' (use -e), 'q50' (use median from AsymptoticLimits), 'auto' (prefer q50, fallback to -e)")
# 新增輸出與 limits 檔案絕對路徑選項
parser.add_option("--outdir", default="", help="Absolute base output directory for Condor (BiasToys, BiasFits, BiasJson, <mA>_BiasPlots will be created here)")
parser.add_option("--limitsFile", default="", help="Absolute path to AsymptoticLimits ROOT file (override default)")
(opts,args) = parser.parse_args()
print()
if opts.nToys>opts.split and not opts.nToys%opts.split==0: raise RuntimeError('The number of toys %g needs to be smaller than or divisible by the split number %g'%(opts.nToys, opts.split))

CombineDir = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results"
LimitFileDir = os.path.join(CombineDir, f"higgsCombine{opts.mA}.AsymptoticLimits.mH125.38.root")
BiasStudyDir = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal"

import ROOT as r
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(2211)
import json

# 將關鍵路徑正規化為絕對路徑（適配 condor）
opts.datacard = os.path.abspath(opts.datacard)
BASE_OUT = os.path.abspath(opts.outdir) if getattr(opts, "outdir", "") else os.path.abspath(BiasStudyDir)
ToyDir  = os.path.join(BASE_OUT, f"{opts.mA}_BiasToys")
FitDir  = os.path.join(BASE_OUT, f"{opts.mA}_BiasFits")
JsonDir = os.path.join(BASE_OUT, "BiasJson")
PlotDir = os.path.join(BASE_OUT, f"{opts.mA}_BiasPlots")
# 確保目錄存在
os.makedirs(ToyDir,  exist_ok=True)
os.makedirs(FitDir,  exist_ok=True)
os.makedirs(JsonDir, exist_ok=True)
os.makedirs(PlotDir, exist_ok=True)

# limits 檔案（允許由參數覆蓋）
LimitFilePath = os.path.abspath(opts.limitsFile) if getattr(opts, "limitsFile", "") else os.path.abspath(LimitFileDir)

# 包裝 biasUtils 的檔名產生器，確保回傳絕對路徑
def _ensure_abs(p, parent):
    return p if os.path.isabs(p) else os.path.join(parent, os.path.basename(p))

try:
    _toyName_orig
except NameError:
    try:
        _toyName_orig = toyName
    except NameError:
        _toyName_orig = None
    try:
        _fitName_orig = fitName
    except NameError:
        _fitName_orig = None
    try:
        _plotName_orig = plotName
    except NameError:
        _plotName_orig = None

if _toyName_orig:
    def toyName(name, split=None):
        p = _toyName_orig(name, split=split)
        return _ensure_abs(p, ToyDir)

if _fitName_orig:
    def fitName(name, split=None):
        p = _fitName_orig(name, split=split)
        return _ensure_abs(p, FitDir)

if _plotName_orig:
    def plotName(name):
        p = _plotName_orig(name)
        base = _ensure_abs(p, PlotDir)
        # 避免重複副檔名
        for ext in (".pdf", ".png", ".jpg", ".jpeg"):
            if base.lower().endswith(ext):
                base = base[: -len(ext)]
                break
        return base

def read_q50_from_limits(fname):
    f = r.TFile.Open(fname, "READ")
    if not f or f.IsZombie():
        print(f"[RunBiasStudy] Cannot open limits file: {fname}")
        return None
    t = f.Get("limit")
    if not t:
        print(f"[RunBiasStudy] No TTree 'limit' in: {fname}")
        f.Close()
        return None
    has_quant = bool(t.GetBranch("quantileExpected"))
    q50 = None
    eps = 1e-3
    entries = []
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        try:
            val = float(getattr(t, "limit"))
        except Exception:
            continue
        qv = None
        if has_quant:
            try:
                qv = float(getattr(t, "quantileExpected"))
            except Exception:
                qv = None
        entries.append((qv, val))
        if has_quant and qv is not None and abs(qv - 0.5) < eps:
            q50 = val
    if q50 is None and len(entries) >= 5:
        offset = 1 if (entries[0][0] is not None and entries[0][0] < 0) else 0
        try:
            q50 = entries[offset + 2][1]
        except Exception:
            q50 = None
    f.Close()
    return q50

# 依 --expectFrom 決定是否使用 q50 覆蓋 expectSignal
_exp_q50 = read_q50_from_limits(LimitFilePath)
_src = (opts.expectFrom or "auto").lower()
if _src == "q50":
    if _exp_q50 is not None:
        opts.expectSignal = round(float(_exp_q50), 4)
        print(f"[RunBiasStudy] Using expectSignal from limits q50: {opts.expectSignal} (--expectFrom=q50)")
    else:
        print("[RunBiasStudy] --expectFrom=q50 but failed to read q50; keep --expectSignal as given.")
elif _src == "arg":
    print(f"[RunBiasStudy] Using expectSignal from arguments: {opts.expectSignal} (--expectFrom=arg)")
else:
    if _exp_q50 is not None:
        opts.expectSignal = round(float(_exp_q50), 4)
        print(f"[RunBiasStudy] Using expectSignal from limits q50: {opts.expectSignal} (--expectFrom=auto)")
    else:
        print("[RunBiasStudy] Failed to read q50 from limits file, keep --expectSignal as given. (--expectFrom=auto)")

ws = r.TFile(opts.datacard).Get(opts.workspace)

pdfs = rooArgSetToList(ws.allPdfs())
multipdfName = None
for pdf in pdfs:
    if pdf.InheritsFrom("RooMultiPdf"):
        if multipdfName is not None: raiseMultiError() 
        multipdfName = pdf.GetName()
        print('Conduct bias study for multipdf called %s'%multipdfName)
multipdf = ws.pdf(multipdfName)
print()

varlist = rooArgSetToList(ws.allCats())
indexName = None
for var in varlist:
    if var.GetName().startswith('pdfindex'):
        if indexName is not None: raiseMultiError()
        indexName = var.GetName()
        print('Found index called %s'%indexName)
print()

from collections import OrderedDict as od
indexNameMap = od()
for ipdf in range(multipdf.getNumPdfs()):
    if opts.selectFunction is not None:
        if not multipdf.getPdf(ipdf).GetName().count(opts.selectFunction): continue
    indexNameMap[ipdf] = multipdf.getPdf(ipdf).GetName()

# 小工具：挑選並搬移單一輸出檔案（避免 shell mv 對檔名目標報錯）
def move_artifact(src_glob, dest_path, expect_tag=None, dry=False):
    files = sorted(glob.glob(src_glob), key=os.path.getmtime)
    if not files:
        print(f"[RunBiasStudy] No files match: {src_glob}")
        return False
    cand = None
    if expect_tag:
        for f in reversed(files):
            if expect_tag in os.path.basename(f):
                cand = f
                break
    if cand is None:
        cand = files[-1]
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    if dry:
        print(f"[RunBiasStudy][dry] mv {cand} -> {dest_path}")
        return True
    try:
        if os.path.exists(dest_path):
            os.remove(dest_path)
    except Exception:
        pass
    shutil.move(cand, dest_path)
    return True

if opts.toys:
    # 使用 mA 前綴的絕對 ToyDir
    if not os.path.isdir(ToyDir): os.makedirs(ToyDir, exist_ok=True)
    # 修正：不要把 combineOptions 當成 --saveToys 的值
    toyCmdBase = 'combine -m %.4f -d %s -M GenerateOnly --expectSignal %.4f -s %g --saveToys'%(opts.mH, opts.datacard, opts.expectSignal, opts.seed)
    for ipdf,pdfName in indexNameMap.items():
        name = shortName(pdfName)
        if opts.nToys > opts.split:
            for isplit in range(opts.nToys//opts.split):
                toyCmd = toyCmdBase + ' -t %g -n _%s_split%g --setParameters %s=%g --freezeParameters %s %s'%(opts.split, name, isplit, indexName, ipdf, indexName, opts.combineOptions)
                run(toyCmd, dry=opts.dryRun)
                # 以 GenerateOnly 標籤挑檔並改名存到 toyName(...)
                src = f'higgsCombine_{name}*.GenerateOnly*.root'
                move_artifact(src, toyName(name,split=isplit), expect_tag='GenerateOnly', dry=opts.dryRun)
        else:
            toyCmd = toyCmdBase + ' -t %g -n _%s --setParameters %s=%g --freezeParameters %s %s'%(opts.nToys, name, indexName, ipdf, indexName, opts.combineOptions)
            run(toyCmd, dry=opts.dryRun)
            src = f'higgsCombine_{name}*.GenerateOnly*.root'
            move_artifact(src, toyName(name), expect_tag='GenerateOnly', dry=opts.dryRun)
print()

if opts.fits:
    # 使用 mA 前綴的絕對 FitDir
    if not os.path.isdir(FitDir): os.makedirs(FitDir, exist_ok=True)
    fitCmdBase = 'combine -m %.4f -d %s -M MultiDimFit -P %s --algo singles %s '%(opts.mH, opts.datacard, opts.poi, opts.combineOptions)
    for ipdf,pdfName in indexNameMap.items():
        name = shortName(pdfName)
        if opts.nToys > opts.split:
            for isplit in range(opts.nToys//opts.split):
                fitCmd = fitCmdBase + ' -t %g -n _%s_split%g --toysFile=%s'%(opts.split, name, isplit, toyName(name,split=isplit))
                run(fitCmd, dry=opts.dryRun)
                # 以 MultiDimFit 標籤挑檔並改名存到 fitName(...)
                src = f'higgsCombine_{name}*.MultiDimFit*.root'
                move_artifact(src, fitName(name,split=isplit), expect_tag='MultiDimFit', dry=opts.dryRun)
            # 使用 FitDir 的絕對路徑樣式合併 split 結果，強制覆蓋 (-f)
            pattern = os.path.join(FitDir, f"*{name}*split*_fits.root")
            run(f'hadd -f {fitName(name)} {pattern}', dry=opts.dryRun)
        else:
            fitCmd = fitCmdBase + ' -t %g -n _%s --toysFile=%s'%(opts.nToys, name, toyName(name))
            run(fitCmd, dry=opts.dryRun)
            src = f'higgsCombine_{name}*.MultiDimFit*.root'
            move_artifact(src, fitName(name), expect_tag='MultiDimFit', dry=opts.dryRun)

if opts.plots:
    # 使用 mA 前綴的絕對 PlotDir
    if not os.path.isdir(PlotDir): os.makedirs(PlotDir, exist_ok=True)
    fit_results = {} if opts.gaussianFit else None
    for ipdf,pdfName in indexNameMap.items():
        name = shortName(pdfName)
        # 更嚴謹的打開與檢查
        tfile = r.TFile.Open(fitName(name), "READ")
        if not tfile or tfile.IsZombie():
            print(f"[RunBiasStudy] Cannot open fit file: {fitName(name)}")
            if tfile: tfile.Close()
            continue
        tree = tfile.Get('limit')
        # 防止 'TObject' 無 GetEntry 的錯誤
        if not tree or not isinstance(tree, r.TTree):
            print(f"[RunBiasStudy] No TTree 'limit' in file: {fitName(name)}")
            tfile.Close()
            continue

        pullHist = r.TH1F('pullsForTruth_%s'%name, 'Pull distribution using the envelope to fit %s'%name, 80, -4., 4.)
        pullHist.GetXaxis().SetTitle('Pull')
        pullHist.GetYaxis().SetTitle('Entries')

        nentries = int(tree.GetEntries())
        ntoys_avail = nentries // 3
        ntoys_loop = min(int(opts.nToys), ntoys_avail)
        tol = 1e-3

        for itoy in range(ntoys_loop):
            # 0: bestfit (quantileExpected == -1)
            if tree.GetEntry(3*itoy) <= 0:
                continue
            if not (hasattr(tree,'quantileExpected') and abs(float(getattr(tree,'quantileExpected')) - (-1.0)) < tol):
                # 不是 bestfit，跳過這個 toy
                continue
            bf = float(getattr(tree, 'r'))
            # 1: lower (quantileExpected == -0.32)
            if tree.GetEntry(3*itoy+1) <= 0:
                continue
            if not (hasattr(tree,'quantileExpected') and abs(float(getattr(tree,'quantileExpected')) - (-0.32)) < tol):
                continue
            lo = float(getattr(tree, 'r'))
            # 2: upper (quantileExpected == +0.32)
            if tree.GetEntry(3*itoy+2) <= 0:
                continue
            if not (hasattr(tree,'quantileExpected') and abs(float(getattr(tree,'quantileExpected')) - (0.32)) < tol):
                continue
            hi = float(getattr(tree, 'r'))

            diff = bf - float(opts.expectSignal)
            unc = 0.5 * (hi - lo)
            if unc > 0.:
                pullHist.Fill(diff/unc)

        canv = r.TCanvas()
        pullHist.Draw()
        if opts.gaussianFit:
           r.gStyle.SetOptFit(111)
           pullHist.Fit('gaus')
           f = pullHist.GetFunction('gaus')
           if f:
               fit_results[name] = {
                   'mean': round(float(f.GetParameter(1)),4),
                   'sigma': round(float(f.GetParameter(2)),4)
               }
        canv.SaveAs('%s.pdf'%plotName(name))
        canv.SaveAs('%s.png'%plotName(name))
        tfile.Close()

    # 使用 mA 前綴的絕對 JsonDir
    if not os.path.isdir(JsonDir): os.makedirs(JsonDir, exist_ok=True)
    if opts.gaussianFit and fit_results is not None:
        out_json = os.path.join(JsonDir, f"{opts.mA}_gaussfit.json")
        payload = {
            'mA': int(opts.mA),
            'exp': round(float(opts.expectSignal),4) if opts.expectSignal is not None else None,
            'fit_results': fit_results
        }
        with open(out_json, 'w') as jf:
            json.dump(payload, jf, indent=2, sort_keys=True)
    else:
        out_json = os.path.join(JsonDir, f"{opts.mA}_exp.json")
        payload = {
            'mA': int(opts.mA),
            'exp': round(float(opts.expectSignal), 4) if opts.expectSignal is not None else None
        }
        with open(out_json, 'w') as jf:
            json.dump(payload, jf, indent=2, sort_keys=True)
