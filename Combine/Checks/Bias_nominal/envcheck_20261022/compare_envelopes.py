#!/usr/bin/env python3
"""Compare the rebuilt envelope of each mass point with the production one.

Prints the members of both multipdfs and, for members present in both, the
largest relative difference of their fitted parameters. The rebuilt envelope
should equal the production one plus the restored member(s); any other
difference means the rerun is not a like-for-like reconstruction.

Usage: python3 compare_envelopes.py [MASS ...]
"""
import os
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

EC = os.path.dirname(os.path.abspath(__file__))
PROD = ("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/"
        "Background/ALP_BkgModel_ReReco/fit_results_run3")
NAME = "CMS-HGG_mva_13p6TeV_multipdf.root"


def members(path):
    f = ROOT.TFile.Open(path)
    w = f.Get("multipdf")
    mp = w.pdf("CMS_hgg_cat0_13p6TeV_bkgshape")
    out = {}
    for i in range(mp.getNumPdfs()):
        pdf = mp.getPdf(i)
        params = {}
        for v in pdf.getVariables():
            if v.GetName() != "CMS_hgg_mass" and v.InheritsFrom("RooRealVar"):
                params[v.GetName()] = v.getVal()
        out[pdf.GetName()] = params
    return f, out


masses = [int(m) for m in sys.argv[1:]] or [2, 3, 4, 20, 24, 25]
for m in masses:
    _fp, prod = members(os.path.join(PROD, str(m), NAME))
    _fn, new = members(os.path.join(EC, "bkg", str(m), NAME))
    print("mA=%d  production: %s" % (m, ", ".join(prod)))
    print("       rebuilt   : %s" % ", ".join(new))
    print("       added %s, dropped %s" % (sorted(set(new) - set(prod)), sorted(set(prod) - set(new))))
    for name in sorted(set(new) & set(prod)):
        diffs = [abs(new[name][k] - v) / max(abs(v), 1e-9)
                 for k, v in prod[name].items() if k in new[name]]
        print("       %-45s max rel. param diff %.2e" % (name, max(diffs) if diffs else float("nan")))
