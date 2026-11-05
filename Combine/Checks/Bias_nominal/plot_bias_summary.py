#!/usr/bin/env python3
"""Bias-study summary versus m_a: mean fitted signal strength (top) and mean pull (bottom).

One marker per truth function of each mass point's envelope, grouped by family
(Bernstein, exponential, power law, Laurent); the order of each family differs between
mass points and is listed in the AN bias tables.

Inputs (high-statistics run, the one quoted in the AN):
  bias_outputs_highstat/mA_<m>/merged/BiasJson/<m>_gaussfit.json
      exp            injected signal strength mu_inj (expected-limit median)
      fit_results.<f>.mean / mean_err : Gaussian-fit mean of the pull
                     (r_hat - mu_inj) / sigma, sigma on the side facing mu_inj
  bias_outputs_highstat/mA_<m>/merged/BiasFits/biasStudy_<f>_fits.root
      tree 'limit', toy triplets (quantileExpected -1, -0.32, 0.32); the best fit
      r_hat of every valid triplet gives the mean of r_hat / mu_inj
Output: <outdir>/bias_summary_vs_mA.{pdf,png}
"""
import argparse, glob, json, os, re, sys
import numpy as np
import ROOT

if ROOT.gROOT.GetVersion().startswith("6.34"):
    sys.exit("ROOT 6.34 drops histograms from PDF output; use LCG_104 (ROOT 6.28)")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

HERE = os.path.dirname(os.path.abspath(__file__))
CONSTANTS = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"
# CMS Petroff-10, same order as the HZgamma style helpers
PETROFF = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2"]
FAMILIES = [("Bern", "Bernstein", 20), ("Exp", "Exponential", 21), ("Pow", "Power law", 22), ("Lau", "Laurent", 23)]
OFFSET = {"Bern": -0.24, "Exp": -0.08, "Pow": 0.08, "Lau": 0.24}
THRESHOLD = 0.2


def run3_lumi():
    # the module cannot be imported standalone, so read the three LUMI entries directly
    txt = open(CONSTANTS).read()
    block = txt[txt.index("LUMI = {"):]
    return sum(float(re.search(r'"%s"\s*:\s*([0-9.]+)' % e, block).group(1)) for e in ("2022", "2023", "2024"))


def mean_rhat(path):
    rdf = ROOT.RDataFrame("limit", path)
    a = rdf.AsNumpy(["r", "quantileExpected"])
    r, q = a["r"], a["quantileExpected"]
    n = len(r) // 3
    r, q = r[:3 * n].reshape(n, 3), q[:3 * n].reshape(n, 3)
    ok = (np.abs(q[:, 0] + 1) < 1e-3) & (np.abs(q[:, 1] + 0.32) < 1e-3) & (np.abs(q[:, 2] - 0.32) < 1e-3) \
        & (r[:, 2] > r[:, 1]) & np.isfinite(r[:, 0])
    bf = r[ok, 0]
    return bf.mean(), bf.std(ddof=1) / np.sqrt(len(bf)), len(bf)


def collect(base):
    rows = []
    for d in sorted(glob.glob(os.path.join(base, "mA_*", "merged")), key=lambda x: int(re.search(r"mA_(\d+)", x).group(1))):
        m = int(re.search(r"mA_(\d+)", d).group(1))
        js = json.load(open(os.path.join(d, "BiasJson", "%d_gaussfit.json" % m)))
        mu_inj = float(js["exp"])
        for f, v in js["fit_results"].items():
            fam = re.match(r"[A-Za-z]+", f).group(0)
            rh, rh_err, ntoy = mean_rhat(os.path.join(d, "BiasFits", "biasStudy_%s_fits.root" % f))
            rows.append(dict(mA=m, func=f, fam=fam, ratio=rh / mu_inj, ratio_err=rh_err / mu_inj,
                             pull=v["mean"], pull_err=v["mean_err"], ntoy=ntoy))
    return rows


def style_axis(h, pad, xtitle, ytitle):
    s = 1.0 / pad.GetAbsHNDC()
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(0.055 * s)
        ax.SetLabelSize(0.050 * s)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().SetTitleOffset(1.10)
    h.GetXaxis().SetTitleOffset(1.05)
    h.GetXaxis().SetLabelOffset(0.010)
    h.GetYaxis().SetNdivisions(505)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=os.path.join(HERE, "bias_outputs_highstat"))
    ap.add_argument("--outdir", default="/eos/home-p/pelai/HZa/output_plots/plots_bias")
    args = ap.parse_args()

    rows = collect(args.base)
    for r in rows:
        print("mA=%2d %-6s ntoy=%5d  r_hat/mu_inj=%.3f+-%.3f  pull=%+.3f+-%.3f" % (
            r["mA"], r["func"], r["ntoy"], r["ratio"], r["ratio_err"], r["pull"], r["pull_err"]))
    lumi = run3_lumi()

    c = ROOT.TCanvas("c", "", 900, 950)
    split = 0.50
    p1 = ROOT.TPad("p1", "", 0, split, 1, 1); p2 = ROOT.TPad("p2", "", 0, 0, 1, split)
    lm = (2.1 * 1.10 + 1.3) * 0.055 * 950 / 900 + 0.012  # house left-margin formula
    for p in (p1, p2):
        p.SetLeftMargin(lm); p.SetRightMargin(0.05); p.Draw()
    p1.SetTopMargin(0.09 / (1 - split)); p1.SetBottomMargin(0.02)
    p2.SetTopMargin(0.02); p2.SetBottomMargin(0.15 / split)

    ratios = [r["ratio"] for r in rows]
    span = max(abs(min(ratios) - 1), abs(max(ratios) - 1)) + 0.05
    p1.cd()
    f1 = ROOT.TH1F("f1", "", 1, 0.3, 30.7)
    f1.SetMinimum(1 - 1.6 * span); f1.SetMaximum(1 + 2.6 * span)
    style_axis(f1, p1, "", "Mean #hat{#mu}/#mu_{inj}")
    f1.GetXaxis().SetLabelSize(0)
    f1.Draw("axis")
    one = ROOT.TLine(0.3, 1, 30.7, 1); one.SetLineStyle(2); one.Draw()

    p2.cd()
    f2 = ROOT.TH1F("f2", "", 1, 0.3, 30.7)
    f2.SetMinimum(-0.32); f2.SetMaximum(0.32)
    style_axis(f2, p2, "m_{a} [GeV]", "Mean pull")
    f2.Draw("axis")
    lines = []
    for y, st in ((0, 2), (THRESHOLD, 7), (-THRESHOLD, 7)):
        l = ROOT.TLine(0.3, y, 30.7, y); l.SetLineStyle(st); l.SetLineColor(ROOT.kGray + 2 if y else ROOT.kBlack); l.Draw(); lines.append(l)

    graphs = []
    ftop = 1 - p1.GetTopMargin()
    leg = ROOT.TLegend(0.45, ftop - 0.02 - 2 * 0.075 / (1 - split), 0.93, ftop - 0.02)
    leg.SetNColumns(2); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040 / (1 - split))
    for i, (fam, name, mk) in enumerate(FAMILIES):
        sel = [r for r in rows if r["fam"] == fam]
        if not sel:
            continue
        col = ROOT.TColor.GetColor(PETROFF[i])
        for pad, key in ((p1, "ratio"), (p2, "pull")):
            g = ROOT.TGraphErrors(len(sel))
            for k, r in enumerate(sel):
                g.SetPoint(k, r["mA"] + OFFSET[fam], r[key]); g.SetPointError(k, 0, r[key + "_err"])
            g.SetMarkerStyle(mk); g.SetMarkerSize(1.4); g.SetMarkerColor(col); g.SetLineColor(col); g.SetLineWidth(2)
            pad.cd(); g.Draw("P same"); graphs.append(g)
            if pad is p1:
                leg.AddEntry(g, name, "p")
    p1.cd(); leg.Draw()

    c.cd()
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042)
    ytop = 1 - 0.09 + 0.020
    t.SetTextAlign(11); t.DrawLatex(lm + 0.005, ytop, "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.DrawLatex(0.95, ytop, "%.2f fb^{-1} (13.6 TeV)" % lumi)

    os.makedirs(args.outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        out = os.path.join(args.outdir, "bias_summary_vs_mA.%s" % ext)
        c.SaveAs(out)
    print("[plot_bias_summary] wrote %s/bias_summary_vs_mA.{pdf,png}; %d functions at %d mass points" % (
        args.outdir, len(rows), len({r['mA'] for r in rows})))


if __name__ == "__main__":
    main()
