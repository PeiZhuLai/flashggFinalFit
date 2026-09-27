#!/usr/bin/env python3
"""Build a real-data workspace for one mass point at an arbitrary BDT cut (2026-09-27).

Used by the mA3 working-point re-optimization to get the DATA background-envelope GOF at each
candidate cut without rerunning the MVA-cut + Tree2WS pipeline. Same selection as the closure
(select_lib: scored P2Root, H_mass in [mass_low, 180], BDT score > cut) and the same workspace
layout as trees2ws_data.py (ws 'CMS_hza_workspace', RooDataSet 'Data_13p6TeV', CMS_hza_mass),
filled with the exact event masses (weight 1).
Validated at mA3 / 0.988: 394 selected events == the official MVA-cut data tree.

Run with CLOSURE_FSRFIX=1 to read the FSR-fix re-production.
Output: <outbase>/mA_M<mA>/ws/run3.root
"""
import argparse
import os

import numpy as np
import ROOT

from select_lib import mass_low, MASS_HIGH
from build_pseudodata_ws import _load_data_mass, WSNAME, DNAME

ROOT.gROOT.SetBatch(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mA", type=int, required=True)
    ap.add_argument("--cut", type=float, required=True)
    ap.add_argument("--outbase", required=True)
    a = ap.parse_args()
    lo = mass_low(a.mA)
    m, s = _load_data_mass(a.mA)
    sel = (s > a.cut) & (m > lo) & (m < MASS_HIGH)
    mm = m[sel]
    ws = ROOT.RooWorkspace(WSNAME, WSNAME)
    intLumi = ROOT.RooRealVar("intLumi", "intLumi", 1000., 0., 999999999.)
    intLumi.setConstant(True)
    getattr(ws, "import")(intLumi)
    mass = ROOT.RooRealVar("CMS_hza_mass", "CMS_hza_mass", 125., float(lo), MASS_HIGH)
    mass.setBins(int(round((MASS_HIGH - lo) / 0.5)))
    weight = ROOT.RooRealVar("weight", "weight", 0.)
    getattr(ws, "import")(mass, ROOT.RooFit.Silence())
    getattr(ws, "import")(weight, ROOT.RooFit.Silence())
    aset = ROOT.RooArgSet(mass, weight)
    d = ROOT.RooDataSet(DNAME, DNAME, aset, "weight")
    for x in mm:
        mass.setVal(float(x))
        d.add(aset, 1.0)
    getattr(ws, "import")(d)
    outdir = os.path.join(a.outbase, "mA_M%d" % a.mA, "ws")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "run3.root")
    ws.writeToFile(out)
    print("mA%d: cut=%.4f data N=%d -> %s" % (a.mA, a.cut, mm.size, out))


if __name__ == "__main__":
    main()
