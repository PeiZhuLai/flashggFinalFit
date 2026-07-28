#!/usr/bin/env python3
"""Build DY-MC pseudo-data workspaces for the HZa sculpting closure test.

For each mA, select the combined DY background MC with the SAME adopted BDT working
point the analysis uses (validated in validate_cuts.py), histogram the Higgs-candidate
mass at 1 GeV (the roohist binning), floor negative bins to zero (NLO negative-weight
undershoots in sparse sidebands), and round to integer counts. This deterministic
"Asimov-from-MC" template is materialized as an UNWEIGHTED RooDataSet 'Data_13p6TeV'
(unit-weight entries at bin centers) inside RooWorkspace 'CMS_hza_workspace', exactly
mimicking the real-data trees2ws output (run3.root). Downstream (fTest envelope fit,
datacard, combine) then runs bit-for-bit like the real analysis, but on MC-as-data.

Why deterministic rounded template (not a Poisson toy): it isolates the SHAPE question
(does the BDT-sculpted DY spectrum fake a 125 peak the smooth envelope cannot absorb?)
without single-toy luck. Negative-weight MC forbids feeding weighted data straight to
the (unbinned) envelope fit; the floored integer template is non-negative and data-like.

Output: <OUTBASE>/mA_M<mA>/ws/run3.root   (ws 'CMS_hza_workspace', RooDataSet 'Data_13p6TeV')

Env: higgs-alp-ana (uproot + ROOT). Writing a plain RooDataSet workspace needs no
combine libraries (unlike READING a multipdf file).
"""
import argparse
import os
import sys

import numpy as np
import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_lib import (ADOPTED_CUTS, MASS_HIGH, load_combined, mass_low, select)

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

WSNAME = "CMS_hza_workspace"          # inputWSName__ (tools/commonObjects.py)
DNAME = "Data_13p6TeV"                # sqrts__ = 13p6TeV
OUTBASE_DEFAULT = ("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/"
                   "flashggFinalFit/pseudodata_closure/pseudodata")


def build_one(mA, outbase, verbose=True, cut=None):
    cut = ADOPTED_CUTS[int(mA)] if cut is None else float(cut)
    lo = mass_low(mA)
    nbin_1gev = int(round(MASS_HIGH - lo))          # 1 GeV bins (roohist binning)
    m, w, s = load_combined(mA)
    msk = select(mA, m, w, s, cut)
    mm, ww = m[msk], w[msk]
    hist, edges = np.histogram(mm, bins=nbin_1gev, range=(lo, MASS_HIGH), weights=ww)
    centers = 0.5 * (edges[:-1] + edges[1:])
    counts = np.rint(np.clip(hist, 0.0, None)).astype(int)   # floor<0, round to int

    # --- workspace mirroring trees2ws_data.py ---
    ws = ROOT.RooWorkspace(WSNAME, WSNAME)
    intLumi = ROOT.RooRealVar("intLumi", "intLumi", 1000., 0., 999999999.)
    intLumi.setConstant(True)
    getattr(ws, "import")(intLumi)
    mass = ROOT.RooRealVar("CMS_hza_mass", "CMS_hza_mass", 125., float(lo), MASS_HIGH)
    mass.setBins(int(round((MASS_HIGH - lo) / 0.5)))         # 0.5 GeV, as trees2ws
    weight = ROOT.RooRealVar("weight", "weight", 0.)
    getattr(ws, "import")(mass, ROOT.RooFit.Silence())
    getattr(ws, "import")(weight, ROOT.RooFit.Silence())
    aset = ROOT.RooArgSet(mass, weight)
    d = ROOT.RooDataSet(DNAME, DNAME, aset, "weight")
    for c, n in zip(centers, counts):
        if n <= 0:
            continue
        mass.setVal(float(c))
        for _ in range(int(n)):
            d.add(aset, 1.0)
    getattr(ws, "import")(d)

    outdir = os.path.join(outbase, f"mA_M{mA}", "ws")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "run3.root")
    ws.writeToFile(outpath)
    if verbose:
        print(f"mA{mA:>2}: cut={cut:.4f} range=[{lo:.0f},{MASS_HIGH:.0f}] "
              f"sumW={ww.sum():.1f} -> pseudo-data N={int(counts.sum())} "
              f"(neg_bins_floored={int((hist<0).sum())})  -> {outpath}")
    return int(counts.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mA", type=str, default="1-30",
                    help="comma list or A-B range, e.g. '3' or '1,2,3' or '1-30'")
    ap.add_argument("--outbase", default=OUTBASE_DEFAULT)
    ap.add_argument("--cut", type=float, default=None,
                    help="override the adopted BDT cut (for WP scans)")
    a = ap.parse_args()
    if "-" in a.mA and "," not in a.mA:
        lo, hi = a.mA.split("-"); mAs = list(range(int(lo), int(hi) + 1))
    else:
        mAs = [int(x) for x in a.mA.split(",")]
    for mA in mAs:
        build_one(mA, a.outbase, cut=a.cut)


if __name__ == "__main__":
    main()
