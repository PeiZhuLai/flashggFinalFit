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
import uproot
from select_lib import (ADOPTED_CUTS, MASS_HIGH, load_combined, mass_low, select,
                        SCORED_BASE, ERAS, TREE, MASS_BR, SCORE_FMT)

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

WSNAME = "CMS_hza_workspace"          # inputWSName__ (tools/commonObjects.py)
DNAME = "Data_13p6TeV"                # sqrts__ = 13p6TeV
OUTBASE_DEFAULT = ("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/"
                   "flashggFinalFit/pseudodata_closure/pseudodata")


BLIND_LOW, BLIND_HIGH = 115.0, 135.0   # blinded Higgs window (same as bkgmcSculptingCheck)


def _load_data_mass(mA):
    """Real Data (unweighted) H_mass + BDT score from run3_bdt_scored_nominal for mass mA."""
    sb = SCORE_FMT.format(int(mA))
    ms, ss = [], []
    for era in ERAS:
        path = os.path.join(SCORED_BASE, "Data", "%s.root" % era)
        if not os.path.exists(path):
            continue
        arr = uproot.open(path)[TREE].arrays([MASS_BR, sb], library="np")
        ms.append(arr[MASS_BR])
        ss.append(arr[sb])
    if not ms:
        return np.array([]), np.array([])
    return np.concatenate(ms), np.concatenate(ss)


def _sideband_scale(mA, cut, mm, ww):
    """data/MC yield ratio in the sidebands (fit window excluding the blind [115,135]).

    Identical normalization to bkgmcSculptingCheck._scale_bkg_to_data_sideband, so the
    closure pseudo-data carries the SAME background-MC yield as the data/MC sculpting-check
    plots (the raw DY MC over-predicts the data by ~1/scale; only shapes are meaningful).
    """
    def _sb(x):
        return (x < BLIND_LOW) | (x > BLIND_HIGH)
    mc_sb = float(ww[_sb(mm)].sum())
    dm, ds = _load_data_mass(mA)
    if dm.size == 0 or mc_sb <= 0.0:
        return 1.0
    dmsk = select(mA, dm, np.ones_like(dm), ds, cut)   # Data: factor == 1
    data_sb = float(_sb(dm[dmsk]).sum())
    return data_sb / mc_sb if mc_sb > 0.0 else 1.0


def build_one(mA, outbase, verbose=True, cut=None, bootstrap=None):
    cut = ADOPTED_CUTS[int(mA)] if cut is None else float(cut)
    lo = mass_low(mA)
    nbin_1gev = int(round(MASS_HIGH - lo))          # 1 GeV bins (roohist binning)
    m, w, s = load_combined(mA)
    msk = select(mA, m, w, s, cut)
    mm, ww = m[msk], w[msk]
    if bootstrap is not None:
        # MC-statistics replica (2026-09-27): Poisson(1) weight per selected MC event, so the
        # spread of the closure r_hat over replicas measures how much of it is MC fluctuation.
        ww = ww * np.random.default_rng(int(bootstrap)).poisson(1.0, ww.size)
    sf = _sideband_scale(mA, cut, mm, ww)           # normalize MC to data in sidebands
    ww = ww * sf                                      # (match bkgmcSculptingCheck)
    if verbose:
        print("[mA=%s] sideband data/MC scale = %.4f" % (mA, sf))
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
    ap.add_argument("--bootstrap", type=int, default=None,
                    help="seed for a Poisson(1) MC-statistics bootstrap replica (default: off)")
    a = ap.parse_args()
    if "-" in a.mA and "," not in a.mA:
        lo, hi = a.mA.split("-"); mAs = list(range(int(lo), int(hi) + 1))
    else:
        mAs = [int(x) for x in a.mA.split(",")]
    for mA in mAs:
        build_one(mA, a.outbase, cut=a.cut, bootstrap=a.bootstrap)


if __name__ == "__main__":
    main()
