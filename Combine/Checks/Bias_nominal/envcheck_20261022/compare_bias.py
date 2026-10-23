#!/usr/bin/env python3
"""Bias of the restored pre-change envelopes vs. the production envelopes.

Restored: outputs/mA_<m>/merged/BiasJson/<m>_gaussfit.json (this study).
Production: ../bias_outputs_highstat/mA_<m>/merged/BiasJson/<m>_gaussfit.json.
For each mass point prints every truth function of the restored envelope, its
production value where the function is shared, the value that triggered the
envelope change, the number of converged toys, and flags |bias| >= 0.2.
"""
import json
import os
import ROOT

ROOT.gROOT.SetBatch(True)
HERE = os.path.dirname(os.path.abspath(__file__))
PROD = os.path.join(os.path.dirname(HERE), "bias_outputs_highstat")
TRIGGER = {(2, "Bern2"): -0.222, (3, "Bern2"): -0.211, (4, "Bern4"): -0.221,
           (20, "Bern6"): +0.208, (24, "Lau1"): -0.232, (25, "Bern6"): +0.212, (25, "Lau1"): -0.215}
MASSES = [int(x) for x in os.environ.get("MASSES", "2 3 4 20 24 25").split()]


def gauss(base, m):
    return json.load(open(os.path.join(base, "mA_%d" % m, "merged", "BiasJson", "%d_gaussfit.json" % m)))


def ntoys(base, m, f):
    tf = ROOT.TFile.Open(os.path.join(base, "mA_%d" % m, "merged", "BiasFits", "biasStudy_%s_fits.root" % f))
    t = tf.Get("limit") if tf and not tf.IsZombie() else None
    n = t.GetEntries() // 3 if t else -1
    if tf:
        tf.Close()
    return n


print("%4s %-6s %16s %16s %9s %7s" % ("mA", "func", "restored", "production", "trigger", "toys"))
for m in MASSES:
    new = gauss(os.path.join(HERE, "outputs"), m)
    prod = gauss(PROD, m)["fit_results"]
    for f, r in sorted(new["fit_results"].items()):
        p = prod.get(f)
        ps = "%+.3f+-%.3f" % (p["mean"], p["mean_err"]) if p else "-"
        tr = TRIGGER.get((m, f))
        flag = "  <-- |bias|>=0.2" if abs(r["mean"]) >= 0.2 else ""
        print("%4d %-6s %+.3f+-%.3f %16s %9s %7d%s"
              % (m, f, r["mean"], r["mean_err"], ps, "%+.3f" % tr if tr is not None else "",
                 ntoys(os.path.join(HERE, "outputs"), m, f), flag))
    print("     mu injected %.4f (production %.4f)" % (new["exp"], gauss(PROD, m)["exp"]))
