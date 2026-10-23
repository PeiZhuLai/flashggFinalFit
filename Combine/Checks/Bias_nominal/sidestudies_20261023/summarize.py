#!/usr/bin/env python3
"""Numbers of the bias side studies (chain_side.sh); no interpretation.

For every merged study: mean pull +- error per truth function (gaussfit json), the number of toys with a
complete best fit + interval, and how many of those hit the r range (-5 or +5) on either side.
Production values (bias_outputs_highstat) are printed next to them where the truth function is shared.
"""
import json
import os
import ROOT

ROOT.gROOT.SetBatch(True)
HERE = os.path.dirname(os.path.abspath(__file__))
BN = os.path.dirname(HERE)
PROD = os.path.join(BN, "bias_outputs_highstat")
STUDIES = [
    ("truthout: Bern6 truth, fitted with the production envelope", os.path.join(HERE, "truthout", "outputs"), [20, 25]),
    ("ma29 fallback", os.path.join(HERE, "ma29", "fallback"), [29]),
    ("ma29 fallback + robustFit", os.path.join(HERE, "ma29", "fbrobust"), [29]),
    ("r2 v11full (m_a = 11 without the low-order restriction)", os.path.join(HERE, "r2", "v11full", "outputs"), [11]),
    ("r2 v11b3 (m_a = 11 low-order, Bernstein up to 3)", os.path.join(HERE, "r2", "v11b3", "outputs"), [11]),
    ("r2 v27b6 (m_a = 27 with Bernstein, cap 6)", os.path.join(HERE, "r2", "v27b6", "outputs"), [27]),
]


def gauss(base, m):
    p = os.path.join(base, "mA_%d" % m, "merged", "BiasJson", "%d_gaussfit.json" % m)
    return json.load(open(p)) if os.path.exists(p) else None


def toys(base, m, f):
    p = os.path.join(base, "mA_%d" % m, "merged", "BiasFits", "biasStudy_%s_fits.root" % f)
    tf = ROOT.TFile.Open(p) if os.path.exists(p) else None
    t = tf.Get("limit") if tf and not tf.IsZombie() else None
    if not t:
        return -1, -1, -1
    n = lo_rail = hi_rail = 0
    for i in range(t.GetEntries() // 3):
        q = []
        for k in range(3):
            t.GetEntry(3 * i + k)
            q.append((t.quantileExpected, t.r))
        if abs(q[0][0] + 1) > 1e-3 or abs(q[1][0] + 0.32) > 1e-3 or abs(q[2][0] - 0.32) > 1e-3:
            continue
        n += 1
        lo_rail += q[1][1] <= -4.99
        hi_rail += q[2][1] >= 4.99
    tf.Close()
    return n, lo_rail, hi_rail


for title, base, masses in STUDIES:
    print("== " + title)
    for m in masses:
        g = gauss(base, m)
        if g is None:
            print("   mA%-2d no merged result" % m)
            continue
        pg = gauss(PROD, m) or {"fit_results": {}, "exp": None}
        print("   mA%-2d mu injected %.4f (production %s)" % (m, g["exp"], pg["exp"]))
        for f, r in sorted(g["fit_results"].items()):
            p = pg["fit_results"].get(f)
            n, lo, hi = toys(base, m, f)
            pn, plo, phi = toys(PROD, m, f) if p else (-1, -1, -1)
            print("     %-6s %+.3f+-%.3f  toys %5d (r=-5 %4d, r=+5 %4d)   production %s  toys %5d (r=-5 %4d, r=+5 %4d)%s"
                  % (f, r["mean"], r["mean_err"], n, lo, hi,
                     "%+.3f+-%.3f" % (p["mean"], p["mean_err"]) if p else "     -     ", pn, plo, phi,
                     "  <-- |bias|>=0.2" if abs(r["mean"]) >= 0.2 else ""))
