#!/usr/bin/env python3
"""Validate the adopted per-mA BDT cut table + selection machinery against the REAL
analysis, by reproducing the observed-data yield.

For each mA 1..30:
  - real yield  = sumEntries of roohist_data_mass_cat0 inside the existing
                  fit_results_run3/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root
  - my   yield  = # scored Data events passing the adopted BDT cut in [massLow,180]

They must agree (Data is unweighted, factor==1). Agreement confirms the cut table
and the whole selection recipe, so the SAME recipe can be applied to DY MC to build
faithful pseudo-data. mA1 is additionally tested against the R=1 alternative cut.

Env: higgs-alp-ana (uproot + ROOT/RooFit).
"""
import os
import sys
import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_lib import (ADOPTED_CUTS, ERAS, MASS_HIGH, SCORED_BASE, TREE,
                        MASS_BR, WEIGHT_BR, SCORE_FMT, mass_low, select)
import numpy as np
import uproot

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

FFIT = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
FITDIR = f"{FFIT}/Background/ALP_BkgModel_ReReco/fit_results_run3"
WSNAME = "multipdf"
ROOHIST = "roohist_data_mass_cat0"

DATA_SAMPLES = {era: ["Data"] for era in ERAS}


def real_yield(mA):
    """sumEntries of roohist_data_mass_cat0 in the real multipdf file for this mA."""
    path = f"{FITDIR}/{mA}/CMS-HGG_mva_13p6TeV_multipdf.root"
    if not os.path.exists(path):
        return None, None
    f = ROOT.TFile.Open(path)
    w = f.Get(WSNAME)
    dh = w.data(ROOHIST) if w else None
    if dh is None:
        f.Close(); return None, None
    n = dh.sumEntries()
    var = w.var("CMS_hza_mass")
    rng = (var.getMin(), var.getMax()) if var else (None, None)
    f.Close()
    return n, rng


def load_data(mA):
    br = SCORE_FMT.format(int(mA))
    m, w, s = [], [], []
    for era in ERAS:
        path = os.path.join(SCORED_BASE, "Data", f"{era}.root")
        a = uproot.open(path)[TREE].arrays([MASS_BR, WEIGHT_BR, br], library="np")
        m.append(a[MASS_BR]); w.append(a[WEIGHT_BR]); s.append(a[br])
    return np.concatenate(m), np.concatenate(w), np.concatenate(s)


def main():
    print(f"{'mA':>3} {'cut':>7} {'real_roohist':>12} {'my_dataYield':>12} "
          f"{'diff':>6} {'range':>12}  status")
    nbad = 0
    for mA in range(1, 31):
        cut = ADOPTED_CUTS[mA]
        rn, rng = real_yield(mA)
        m, w, s = load_data(mA)
        mask = select(mA, m, w, s, cut)
        my = int(mask.sum())
        if rn is None:
            print(f"{mA:>3} {cut:>7.4f} {'MISSING':>12} {my:>12} "
                  f"{'--':>6} {'--':>12}  no real file")
            continue
        diff = my - rn
        ok = abs(diff) <= max(2, 0.01 * rn)  # allow tiny rounding/boundary
        tag = "OK" if ok else "*** MISMATCH"
        if not ok:
            nbad += 1
        rngs = f"[{rng[0]:.0f},{rng[1]:.0f}]" if rng[0] is not None else "?"
        print(f"{mA:>3} {cut:>7.4f} {rn:>12.1f} {my:>12} "
              f"{diff:>6.0f} {rngs:>12}  {tag}")
        # mA1: also probe R=1 alternative cuts if adopted mismatches
        if mA == 1 and not ok:
            for alt in (0.955, 0.945, 0.905, 0.9):
                ym = int(select(mA, m, w, s, alt).sum())
                print(f"       mA1 alt cut {alt:.4f} -> Data yield {ym} "
                      f"(target {rn:.0f})")
    print(f"\n{'ALL MATCH' if nbad == 0 else str(nbad)+' MISMATCH(es)'}")
    return nbad


if __name__ == "__main__":
    sys.exit(main())
