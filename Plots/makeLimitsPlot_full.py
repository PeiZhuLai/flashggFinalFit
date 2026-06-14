#!/usr/bin/env python3
"""Full-range 95% CL XS limit plot: m_a = 0.1 .. 30 GeV on ONE log-x canvas.

Stitches the two H->Za regimes that use different reconstruction:
  * m_a < 1 GeV  : MERGED single-photon analysis (this work, Run-3 172.13 fb^-1).
                   combine roots = Combine/merged_limits/higgsCombine_merged_flashgg_M0p{i}...
  * m_a >= 1 GeV : RESOLVED di-photon analysis.
                   combine roots = Combine/output_combine_results/higgsCombine{m}...

Both pipelines use the same flashgg signal normalisation, so the limit r maps to
sigma x B identically (sigma x B [fb] = r * assume_xs, assume_xs = 100 fb). A
vertical dashed line at m_a = 1 GeV marks the merged|resolved boundary.

Reuses read_limits_from_file() from makeLimitsPlot.py. Run in env higgs-alp-ana.
"""
import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError + 1   # silence the harmless TList teardown spam
ROOT.gStyle.SetOptStat(0)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from makeLimitsPlot import read_limits_from_file  # noqa: E402
try:
    from makeLimitsPlot import _draw_cms_prelim_and_lumi
except Exception:
    _draw_cms_prelim_and_lumi = None

FLASHGG = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"


def cms_label(lumi):
    if _draw_cms_prelim_and_lumi is not None:
        _draw_cms_prelim_and_lumi(ROOT.gPad, textSize=0.05, lumi_fb=lumi, sqrts_tev=13.6,
                                  left_margin=0.13, right_margin=0.05, top_margin=0.08,
                                  drawLumi=True)
        return []
    keep = []
    l = ROOT.TLatex(); l.SetNDC(); l.SetTextFont(42)
    l.SetTextSize(0.05); l.DrawLatex(0.13, 0.945, "#bf{CMS} #it{Preliminary}")
    l.SetTextSize(0.045); l.DrawLatex(0.62, 0.945, f"{lumi:.2f} fb^{{-1}} (13.6 TeV)")
    keep.append(l)
    return keep


def collect(masses, fname_fn, assume_xs):
    """Return list of (ma, exp, e1lo, e1hi, e2lo, e2hi) for valid points."""
    pts, missing = [], []
    for ma in masses:
        ok, q025, q16, q50, q84, q975, obs, hasObs = read_limits_from_file(fname_fn(ma))
        if not ok:
            missing.append(ma); continue
        exp = q50 * assume_xs
        pts.append((ma, exp,
                    (q50 - q16) * assume_xs, (q84 - q50) * assume_xs,
                    (q50 - q025) * assume_xs, (q975 - q50) * assume_xs))
    return pts, missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=f"{FLASHGG}/Plots/plot_limits")
    ap.add_argument("--merged-dir", default=f"{FLASHGG}/Combine/merged_limits")
    ap.add_argument("--merged-pattern",
                    default="higgsCombine_merged_flashgg_M0p{i}.AsymptoticLimits.mH125.root")
    ap.add_argument("--resolved-dir", default=f"{FLASHGG}/Combine/output_combine_results")
    ap.add_argument("--resolved-pattern",
                    default="higgsCombine{m}.AsymptoticLimits.mH125.38.root")
    ap.add_argument("--merged-masses", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9")
    ap.add_argument("--resolved-masses",
                    default=",".join(str(m) for m in range(1, 31)))
    ap.add_argument("--assume-xs", type=float, default=100.0)
    ap.add_argument("--lumi", type=float, default=172.13)
    ap.add_argument("--formats", default="pdf,png")
    ap.add_argument("--tag", default="full_0p1_30")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    m_masses = [float(x) for x in args.merged_masses.split(",") if x]
    r_masses = [float(x) for x in args.resolved_masses.split(",") if x]

    m_pts, m_missing = collect(
        m_masses,
        lambda ma: os.path.join(args.merged_dir,
                                args.merged_pattern.format(i=int(round(ma * 10)))),
        args.assume_xs)
    r_pts, r_missing = collect(
        r_masses,
        lambda ma: os.path.join(args.resolved_dir,
                                args.resolved_pattern.format(m=int(round(ma)))),
        args.assume_xs)

    if not m_pts and not r_pts:
        print("[full] no valid combine files found", file=sys.stderr); sys.exit(1)
    if m_missing:
        print(f"[full] missing merged m_a: {m_missing}")
    if r_missing:
        print(f"[full] missing resolved m_a: {r_missing}")

    all_pts = sorted(m_pts + r_pts, key=lambda p: p[0])

    g_exp = ROOT.TGraph()
    g_1s = ROOT.TGraphAsymmErrors()
    g_2s = ROOT.TGraphAsymmErrors()
    for j, (ma, exp, e1lo, e1hi, e2lo, e2hi) in enumerate(all_pts):
        g_exp.SetPoint(j, ma, exp)
        g_1s.SetPoint(j, ma, exp); g_2s.SetPoint(j, ma, exp)
        g_1s.SetPointError(j, 0, 0, e1lo, e1hi)
        g_2s.SetPointError(j, 0, 0, e2lo, e2hi)

    c = ROOT.TCanvas("cLimits_full", "", 900, 600)
    c.SetLeftMargin(0.13); c.SetRightMargin(0.05); c.SetBottomMargin(0.13); c.SetTopMargin(0.08)
    c.SetLogx(); c.SetLogy()

    xlo, xhi = 0.08, 35.0
    ys = [p[1] for p in all_pts]
    e2 = [p[5] for p in all_pts]
    ymax = max(y + e for y, e in zip(ys, e2)) * 2.5
    ymin = min(ys) * 0.25

    ytitle = "#sigma(pp #rightarrow H) #times B(#rightarrow Za #rightarrow 2l + 2#gamma) [fb]"
    frame = ROOT.TH1F("frame", f";m_{{a}} (GeV);{ytitle}", 100, xlo, xhi)
    frame.SetMinimum(ymin); frame.SetMaximum(ymax)
    frame.GetXaxis().SetTitleSize(0.05); frame.GetXaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetMoreLogLabels(); frame.GetXaxis().SetNoExponent()
    frame.GetYaxis().SetTitleSize(0.048); frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw()

    # colour scheme matched to Limits_XS.pdf: 2sigma = kYellow(5), 1sigma = kGreen(3),
    # median = black dashed (no markers).
    g_2s.SetFillColor(5); g_2s.SetLineColor(5); g_2s.SetFillStyle(1001)
    g_1s.SetFillColor(3); g_1s.SetLineColor(3); g_1s.SetFillStyle(1001)
    g_exp.SetLineColor(ROOT.kBlack); g_exp.SetLineWidth(2); g_exp.SetLineStyle(2)
    g_2s.Draw("3 same"); g_1s.Draw("3 same"); g_exp.Draw("L same")

    # boundary marker at m_a = 1 GeV + regime labels in the empty lower band
    bnd = ROOT.TLine(1.0, ymin, 1.0, ymax)
    bnd.SetLineColor(ROOT.kGray + 2); bnd.SetLineStyle(3); bnd.SetLineWidth(2); bnd.Draw()
    txt = ROOT.TLatex(); txt.SetTextFont(42); txt.SetTextSize(0.034)
    txt.SetTextColor(ROOT.kGray + 3)
    ylab = ymin * 2.0
    txt.SetTextAlign(32); txt.DrawLatex(0.92, ylab, "merged #gamma #leftarrow")
    txt.SetTextAlign(12); txt.DrawLatex(1.1, ylab, "#rightarrow resolved #gamma#gamma")

    leg = ROOT.TLegend(0.42, 0.69, 0.93, 0.90)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.038)
    leg.AddEntry(g_exp, "Median expected", "l")
    leg.AddEntry(g_1s, "68% expected", "f")
    leg.AddEntry(g_2s, "95% expected", "f")
    leg.AddEntry(0, "m_{a} < 1 GeV merged / #geq 1 GeV resolved", "")
    leg.Draw()

    ROOT.gPad.SetTicks(1, 1); ROOT.gPad.RedrawAxis()
    keep = cms_label(args.lumi)

    for ext in args.formats.split(","):
        out = os.path.join(args.outdir, f"Limits_XS_{args.tag}.{ext}")
        c.SaveAs(out)
        print(f"[full] wrote {out}")


if __name__ == "__main__":
    main()
