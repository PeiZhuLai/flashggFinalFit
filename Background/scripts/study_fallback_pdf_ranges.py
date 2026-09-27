#!/usr/bin/env python3
"""Test whether the fallback background functions are limited by their parameter ranges.

Every function that only entered the envelope through the fallback rule has one
or more shape parameters sitting on a range boundary (and, at m_a = 22 and 28,
the Laurent turn-on even sits outside its declared range).  This script reloads
each of those fits from the multipdf workspace, refits them as they are, then
refits them with the turn-on window, the Gaussian smearing and the step width
opened up, and compares the chi2.

Nothing is written back: the workspaces are opened read-only and only in-memory
copies of the parameters are touched.  The purpose is to decide whether a
re-tuning of the ranges is worth doing, not to produce the final numbers.

Input
  <bkgmodel>/fit_results_run3/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root

Usage
  cmsenv && python3 scripts/study_fallback_pdf_ranges.py
"""

import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BKG_DIR = os.path.dirname(SCRIPT_DIR)
FITRES = os.path.join(BKG_DIR, "ALP_BkgModel_ReReco", "fit_results_run3")

# Functions that only entered the envelope through the fallback rule.
FALLBACK = {
    2: ["Exp1", "Pow3", "Lau1", "Bern2"],
    22: ["Pow3", "Lau1"],
    28: ["Exp1", "Pow1", "Lau1"],
    29: ["Pow1", "Lau1"],
    30: ["Lau1"],
}

# How far each kind of shape parameter is allowed to move once opened up.
WIDE = {
    "turnon": (98.0, 125.0),
    "step": (98.0, 125.0),
    "sigma": (0.5, 25.0),
    "width": (0.2, 25.0),
    "gsigma": (0.5, 25.0),
    "stepWidth": (0.2, 25.0),
}


def classify(name):
    low = name.lower()
    for key in ("turnon", "stepwidth", "gsigma", "sigma", "width", "step"):
        if key in low:
            return {"stepwidth": "stepWidth", "gsigma": "gsigma"}.get(key, key)
    return None


def float_params(pdf, obs):
    pars = pdf.getParameters(ROOT.RooArgSet(obs))
    out = []
    it = pars.createIterator()
    v = it.Next()
    while v:
        if v.InheritsFrom("RooRealVar") and not v.isConstant():
            out.append(v)
        v = it.Next()
    return out


def fit_and_chi2(pdf, data, mass, nbins):
    nll = pdf.createNLL(data, ROOT.RooFit.Offset(True))
    minim = ROOT.RooMinimizer(nll)
    minim.setPrintLevel(-1)
    minim.setStrategy(1)
    minim.setOffsetting(True)
    status = minim.minimize("Minuit2", "migrad")

    frame = mass.frame()
    data.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(nbins),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    pdf.plotOn(frame, ROOT.RooFit.Name("pdf"))
    npar = len(float_params(pdf, mass))
    chi2ndf = frame.chiSquare("pdf", "data", npar)
    ndf = max(nbins - npar, 1)
    prob = ROOT.TMath.Prob(chi2ndf * ndf, ndf)
    return status, chi2ndf, prob, npar


def main():
    nbins = 85  # 1 GeV bins over 95-180 GeV

    print("%-4s %-6s | %-28s | %-28s" % ("mA", "pdf", "as-is (chi2/ndf, p, status)",
                                         "ranges opened (chi2/ndf, p, status)"))
    print("-" * 96)

    for mA in sorted(FALLBACK):
        path = os.path.join(FITRES, str(mA), "CMS-HGG_mva_13p6TeV_multipdf.root")
        fin = ROOT.TFile.Open(path)
        if not fin or fin.IsZombie():
            print("[skip] cannot open %s" % path)
            continue

        ws = None
        for key in fin.GetListOfKeys():
            obj = key.ReadObj()
            if obj.InheritsFrom("RooWorkspace"):
                ws = obj
                break
        if ws is None:
            print("[skip] no workspace in %s" % path)
            fin.Close()
            continue

        data = None
        for d in ws.allData():
            data = d
        mass = ws.var("CMS_hza_mass")
        if data is None or mass is None:
            print("[skip] mA=%d: no data/observable" % mA)
            fin.Close()
            continue
        mass.setRange(95.0, 180.0)

        multipdf = None
        it = ws.allPdfs().createIterator()
        p = it.Next()
        while p:
            if p.ClassName() == "RooMultiPdf":
                multipdf = p
            p = it.Next()

        for i in range(multipdf.getNumPdfs()):
            pdf = multipdf.getPdf(i)
            if pdf.GetName() not in FALLBACK[mA]:
                continue

            pars = float_params(pdf, mass)
            snapshot = [(v, v.getVal(), v.getMin(), v.getMax()) for v in pars]

            st0, c0, p0, _ = fit_and_chi2(pdf, data, mass, nbins)

            # Reset, then widen every step/Gaussian parameter and let the
            # falling-component parameters move twice as far as before.
            for v, val, lo, hi in snapshot:
                v.setRange(lo, hi)
                v.setVal(val)
            opened = []
            for v, val, lo, hi in snapshot:
                kind = classify(v.GetName())
                if kind is not None:
                    nlo, nhi = WIDE[kind]
                    v.setRange(min(nlo, lo), max(nhi, hi))
                    opened.append(v.GetName())
                else:
                    span = hi - lo
                    v.setRange(lo - 0.5 * span, hi + 0.5 * span)
                v.setVal(val)

            st1, c1, p1, _ = fit_and_chi2(pdf, data, mass, nbins)

            print("%-4d %-6s | %8.3f  %9.3g  st=%-3d | %8.3f  %9.3g  st=%-3d %s"
                  % (mA, pdf.GetName(), c0, p0, st0, c1, p1, st1,
                     "IMPROVED" if c1 < c0 - 1e-6 else ""))

            for v, val, lo, hi in snapshot:
                v.setRange(lo, hi)
                v.setVal(val)

        fin.Close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
