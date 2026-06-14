#!/usr/bin/env python3
"""Brazilian-band 95% CL limit plot for the LOW-mA merged points (m_a = 0.1-0.9 GeV).

Companion to makeLimitsPlot.py (which covers the resolved m_a = 1-30 GeV). The
sub-GeV points span a decade, so the **x-axis is log scale**. Reuses
read_limits_from_file() from makeLimitsPlot.py to read the AsymptoticLimits
`limit` TTree of each combine output.

Default inputs are the merged real-data combine roots
(`MergedAna/output/limits/higgsCombine_merged_M0p<i>.AsymptoticLimits.mH125.root`);
point --combine-dir / --pattern elsewhere once the merged flashgg combine is run.
"""
import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# reuse the limit reader (+ CMS label drawer if present) from the resolved script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from makeLimitsPlot import read_limits_from_file  # noqa: E402
try:
    from makeLimitsPlot import _draw_cms_prelim_and_lumi
except Exception:
    _draw_cms_prelim_and_lumi = None


def cms_label(lumi, left=0.15, top=0.08):
    if _draw_cms_prelim_and_lumi is not None:
        _draw_cms_prelim_and_lumi(ROOT.gPad, textSize=0.05, lumi_fb=lumi,
                                  sqrts_tev=13.6, left_margin=left, right_margin=0.05,
                                  top_margin=top, drawLumi=True)
        return []
    keep = []
    l = ROOT.TLatex(); l.SetNDC(); l.SetTextFont(42)
    l.SetTextSize(0.05); l.DrawLatex(left, 0.945, "#bf{CMS} #it{Preliminary}")
    l.SetTextSize(0.045); l.DrawLatex(0.66, 0.945, f"{lumi:.1f} fb^{{-1}} (13.6 TeV)")
    keep.append(l)
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="plot_limits")
    ap.add_argument("--combine-dir",
                    default="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/output/limits")
    ap.add_argument("--pattern", default="higgsCombine_merged_M0p{i}.AsymptoticLimits.mH125.root",
                    help="filename pattern with {i} = sub-GeV index 1..9 (m_a = 0.i GeV)")
    ap.add_argument("--masses", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9")
    ap.add_argument("--assume-xs", type=float, default=100.0, help="reference signal xs [fb] (r unit)")
    ap.add_argument("--lumi", type=float, default=109.82, help="lumi [fb^-1] for label")
    ap.add_argument("--formats", default="pdf,png")
    ap.add_argument("--tag", default="lowmA")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    masses = [float(x) for x in args.masses.split(",") if x]

    g_exp = ROOT.TGraph()
    g_1s = ROOT.TGraphAsymmErrors()
    g_2s = ROOT.TGraphAsymmErrors()
    j = 0
    failed = []
    for ma in masses:
        i = int(round(ma * 10))  # 0.5 -> 5  -> M0p5
        fname = os.path.join(args.combine_dir, args.pattern.format(i=i))
        ok, q025, q16, q50, q84, q975, obs, hasObs = read_limits_from_file(fname)
        if not ok:
            failed.append(ma)
            continue
        exp = q50 * args.assume_xs
        g_exp.SetPoint(j, ma, exp)
        g_1s.SetPoint(j, ma, exp)
        g_2s.SetPoint(j, ma, exp)
        g_1s.SetPointError(j, 0, 0, (q50 - q16) * args.assume_xs, (q84 - q50) * args.assume_xs)
        g_2s.SetPointError(j, 0, 0, (q50 - q025) * args.assume_xs, (q975 - q50) * args.assume_xs)
        j += 1
    if j == 0:
        print(f"[lowmA] no valid combine files in {args.combine_dir} (pattern {args.pattern})",
              file=sys.stderr)
        if failed:
            print("  missing m_a:", failed, file=sys.stderr)
        sys.exit(1)
    if failed:
        print(f"[lowmA] warning: missing m_a points {failed}")

    c = ROOT.TCanvas("cLimits_lowmA", "", 800, 600)
    left = 0.15
    c.SetLeftMargin(left); c.SetRightMargin(0.05); c.SetBottomMargin(0.13); c.SetTopMargin(0.08)
    c.SetLogx()   # <-- log x-axis: sub-GeV masses span a decade
    c.SetLogy()

    xlo, xhi = min(masses) * 0.8, max(masses) * 1.15
    ys = [g_exp.GetPointY(k) for k in range(g_exp.GetN())]
    e2 = [g_2s.GetErrorYhigh(k) for k in range(g_2s.GetN())]
    ymax = max(y + e for y, e in zip(ys, e2)) * 2.0
    ymin = min(ys) * 0.3

    ytitle = "#sigma(pp #rightarrow H) #times B(#rightarrow Za #rightarrow 2l + 2#gamma) [fb]"
    frame = ROOT.TH1F("frame", f";m_{{a}} (GeV);{ytitle}", 100, xlo, xhi)
    frame.SetMinimum(ymin); frame.SetMaximum(ymax)
    frame.GetXaxis().SetTitleSize(0.055); frame.GetXaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetMoreLogLabels(); frame.GetXaxis().SetNoExponent()
    frame.GetYaxis().SetTitleSize(0.05); frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.Draw()

    g_2s.SetFillColor(ROOT.kOrange); g_2s.SetLineColor(ROOT.kOrange); g_2s.SetFillStyle(1001)
    g_1s.SetFillColor(ROOT.kGreen + 1); g_1s.SetLineColor(ROOT.kGreen + 1); g_1s.SetFillStyle(1001)
    g_exp.SetLineColor(ROOT.kBlack); g_exp.SetLineWidth(3); g_exp.SetLineStyle(2)
    g_2s.Draw("3 same"); g_1s.Draw("3 same"); g_exp.Draw("L same")

    leg = ROOT.TLegend(0.55, 0.68, 0.94, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.042)
    leg.AddEntry(g_exp, "Median expected", "l")
    leg.AddEntry(g_1s, "68% expected", "f")
    leg.AddEntry(g_2s, "95% expected", "f")
    leg.AddEntry(0, "merged photon, m_{a} < 1 GeV", "")
    leg.Draw()

    ROOT.gPad.SetTicks(1, 1); ROOT.gPad.RedrawAxis()
    keep = cms_label(args.lumi, left=left)

    for ext in args.formats.split(","):
        out = os.path.join(args.outdir, f"Limits_XS_{args.tag}.{ext}")
        c.SaveAs(out)
        print(f"[lowmA] wrote {out}")


if __name__ == "__main__":
    main()
