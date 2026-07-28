#!/usr/bin/env python3
"""Collect the sculpting-closure combine outputs into a single JSON.

Per mA it reads:
  output_limits/higgsCombinePS<mA>.AsymptoticLimits.mH125.38.root  (tree 'limit')
     quantileExpected: -1=observed(on MC pseudo-data), .025/.16/.5/.84/.975=expected band
  output_significance/higgsCombinePS<mA>.Significance.mH125.38.root (tree 'limit', limit=Z)
  output_fitdiag/fitDiagnosticsPS<mA>.root                          (tree 'tree_fit_sb': r,rLoErr,rHiErr)

Closure logic: on background-only MC pseudo-data the observed 95% CL limit should sit
inside the expected band and r-hat/significance should be ~0 across mA (esp. low-mA),
i.e. the BDT-sculpted DY spectrum does not fake a Higgs-mass peak.

Env: cmsenv (ROOT). Writes pseudodata_closure/closure_results.json.
"""
import json
import os
import ROOT

ROOT.gROOT.SetBatch(True)
PC = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/pseudodata_closure"


def read_limit(mA):
    p = f"{PC}/output_limits/higgsCombinePS{mA}.AsymptoticLimits.mH125.38.root"
    if not os.path.exists(p):
        return None
    f = ROOT.TFile.Open(p); t = f.Get("limit")
    q = {}
    for e in t:
        q[round(e.quantileExpected, 3)] = e.limit
    f.Close()
    return {
        "obs": q.get(-1.0),
        "exp_m2": q.get(0.025), "exp_m1": q.get(0.16), "exp_med": q.get(0.5),
        "exp_p1": q.get(0.84), "exp_p2": q.get(0.975),
    }


def read_signif(mA):
    p = f"{PC}/output_significance/higgsCombinePS{mA}.Significance.mH125.38.root"
    if not os.path.exists(p):
        return None
    f = ROOT.TFile.Open(p); t = f.Get("limit")
    z = None
    for e in t:
        z = e.limit
    f.Close()
    return z


def read_rhat(mA):
    p = f"{PC}/output_fitdiag/fitDiagnosticsPS{mA}.root"
    if not os.path.exists(p):
        return None
    f = ROOT.TFile.Open(p); t = f.Get("tree_fit_sb")
    if not t:
        f.Close(); return None
    r = rlo = rhi = None
    for e in t:
        r, rlo, rhi = e.r, e.rLoErr, e.rHiErr
    f.Close()
    return {"r": r, "rLoErr": rlo, "rHiErr": rhi}


def main():
    out = {}
    for mA in range(1, 31):
        lim = read_limit(mA)
        if lim is None:
            continue
        out[mA] = {"limit": lim, "signif": read_signif(mA), "rhat": read_rhat(mA)}
    with open(f"{PC}/closure_results.json", "w") as fo:
        json.dump(out, fo, indent=1)
    # console summary
    print(f"{'mA':>3} {'obs':>8} {'exp_med':>8} {'[-1,+1]band':>18} "
          f"{'in?':>4} {'rhat':>16} {'Z':>6}")
    for mA in sorted(out):
        L = out[mA]["limit"]; rh = out[mA]["rhat"]; z = out[mA]["signif"]
        o, em, m1, p1 = L["obs"], L["exp_med"], L["exp_m1"], L["exp_p1"]
        inband = (m1 is not None and o is not None and m1 <= o <= p1)
        rs = f"{rh['r']:+.3f}-{rh['rLoErr']:.3f}+{rh['rHiErr']:.3f}" if rh else "NA"
        print(f"{mA:>3} {o:>8.4f} {em:>8.4f} [{m1:>7.4f},{p1:>7.4f}] "
              f"{'Y' if inband else 'n':>4} {rs:>16} {z if z is not None else -9:>6.2f}")


if __name__ == "__main__":
    main()
