#!/usr/bin/env python3

from biasUtils import *

from optparse import OptionParser
parser = OptionParser()
parser.add_option('--mA', dest='mA', default=5, type='int', help="ALP mass") # PZ
parser.add_option("-d","--datacard",default="Datacard.root")
parser.add_option("-w","--workspace",default="w")
parser.add_option("-t","--toys",action="store_true", default=False)
parser.add_option("-n","--nToys",default=2000,type="int")
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
parser.add_option("--expectFrom", default="q50",help="Source of expectSignal: 'arg' (use -e), 'q50' (use median from AsymptoticLimits), 'auto' (prefer q50, fallback to -e)")
(opts,args) = parser.parse_args()
print()
if opts.nToys>opts.split and not opts.nToys%opts.split==0: raise RuntimeError('The number of toys %g needs to be smaller than or divisible by the split number %g'%(opts.nToys, opts.split))

COMBINE_BASW = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results"
limits_path = os.path.join(COMBINE_BASW, f"higgsCombine{opts.mA}.AsymptoticLimits.mH125.38.root")

import ROOT as r
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(2211)
import os
import json

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
_exp_q50 = read_q50_from_limits(limits_path)
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

if opts.toys:
    if not os.path.isdir('BiasToys'): os.makedirs('BiasToys', exist_ok=True)
    toyCmdBase = 'combine -m %.4f -d %s -M GenerateOnly --expectSignal %.4f -s %g --saveToys %s '%(opts.mH, opts.datacard, opts.expectSignal, opts.seed, opts.combineOptions)
    for ipdf,pdfName in indexNameMap.items():
        name = shortName(pdfName)
        if opts.nToys > opts.split:
            for isplit in range(opts.nToys//opts.split):
                toyCmd = toyCmdBase + ' -t %g -n _%s_split%g --setParameters %s=%g --freezeParameters %s'%(opts.split, name, isplit, indexName, ipdf, indexName)
                run(toyCmd, dry=opts.dryRun)
                os.system('mv higgsCombine_%s* %s'%(name, toyName(name,split=isplit)))
        else: 
            toyCmd = toyCmdBase + ' -t %g -n _%s --setParameters %s=%g --freezeParameters %s'%(opts.nToys, name, indexName, ipdf, indexName)
            run(toyCmd, dry=opts.dryRun)
            os.system('mv higgsCombine_%s* %s'%(name, toyName(name)))
print()

if opts.fits:
    if not os.path.isdir('BiasFits'): os.makedirs('BiasFits', exist_ok=True)
    fitCmdBase = 'combine -m %.4f -d %s -M MultiDimFit -P %s --algo singles %s '%(opts.mH, opts.datacard, opts.poi, opts.combineOptions)
    for ipdf,pdfName in indexNameMap.items():
        name = shortName(pdfName)
        if opts.nToys > opts.split:
            for isplit in range(opts.nToys//opts.split):
                fitCmd = fitCmdBase + ' -t %g -n _%s_split%g --toysFile=%s'%(opts.split, name, isplit, toyName(name,split=isplit))
                run(fitCmd, dry=opts.dryRun)
                os.system('mv higgsCombine_%s* %s'%(name, fitName(name,split=isplit)))
            run('hadd %s BiasFits/*%s*split*.root'%(fitName(name),name), dry=opts.dryRun)
        else:
            fitCmd = fitCmdBase + ' -t %g -n _%s --toysFile=%s'%(opts.nToys, name, toyName(name))
            run(fitCmd, dry=opts.dryRun)
            os.system('mv higgsCombine_%s* %s'%(name, fitName(name)))

if opts.plots:
    if not os.path.isdir('BiasPlots'): os.makedirs('BiasPlots', exist_ok=True)
    fit_results = {} if opts.gaussianFit else None
    for ipdf,pdfName in indexNameMap.items():
        name = shortName(pdfName)
        tfile = r.TFile(fitName(name))
        tree = tfile.Get('limit')
        pullHist = r.TH1F('pullsForTruth_%s'%name, 'Pull distribution using the envelope to fit %s'%name, 80, -4., 4.)
        pullHist.GetXaxis().SetTitle('Pull')
        pullHist.GetYaxis().SetTitle('Entries')
        for itoy in range(opts.nToys):
            tree.GetEntry(3*itoy)
            if not getattr(tree,'quantileExpected')==-1: 
                raiseFailError(itoy,True) 
                continue
            bf = getattr(tree, 'r')
            tree.GetEntry(3*itoy+1)
            if not abs(getattr(tree,'quantileExpected')--0.32)<0.001: 
                raiseFailError(itoy,True) 
                continue
            lo = getattr(tree, 'r')
            tree.GetEntry(3*itoy+2)
            if not abs(getattr(tree,'quantileExpected')-0.32)<0.001: 
                raiseFailError(itoy,True) 
                continue
            hi = getattr(tree, 'r')
            diff = bf - opts.expectSignal
            unc = 0.5 * (hi-lo)
            # mu_true = float(opts.expectSignal)
            # if bf > mu_true:
            #     unc = hi - bf           # 上誤差（應為正）
            # else:
            #     unc = bf - lo           # 下誤差（應為正）
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
                   'mean_err': round(float(f.GetParError(1)),4),
                   'sigma': round(float(f.GetParameter(2)),4),
                   'sigma_err': round(float(f.GetParError(2)),4)
               }
        canv.SaveAs('%s.pdf'%plotName(name))
        canv.SaveAs('%s.png'%plotName(name))
    if opts.gaussianFit and fit_results is not None:
        if not os.path.isdir('BiasJson'): os.makedirs('BiasJson', exist_ok=True)
        out_json = os.path.join('BiasJson', f"{opts.mA}_gaussfit.json")
        payload = {
            'mA': int(opts.mA),
            'exp': round(float(opts.expectSignal),4) if opts.expectSignal is not None else None,
            'fit_results': fit_results
        }
        with open(out_json, 'w') as jf:
            json.dump(payload, jf, indent=2, sort_keys=True)
    else:
        out_json = os.path.join('BiasJson', f"{opts.mA}_exp.json")
        payload = {
            'mA': int(opts.mA),
            'exp': round(float(opts.expectSignal), 4) if opts.expectSignal is not None else None
        }
        with open(out_json, 'w') as jf:
            json.dump(payload, jf, indent=2, sort_keys=True)
