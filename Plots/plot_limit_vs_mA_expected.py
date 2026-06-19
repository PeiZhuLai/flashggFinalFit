#!/usr/bin/env python3
"""Expected (blind) 95% CL limit on mu vs mA (run3), styled like the AN plot.

Reads AsymptoticLimits --run blind outputs from output_combine_results/ for mA = 1..30
and draws: Median expected (dashed), 68% (green), 95% (yellow). No observed line.

mA1 uses the R=1 working point; mA2..30 are the production cut. Output -> doc/HZa.
Run in higgs-alp-ana (ROOT 6.24).
"""
import os, array, ROOT

EXPDIR = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results"
OUTDIR = "/afs/cern.ch/user/p/pelai/doc/HZa"
MASSES = list(range(1, 31))
LUMI   = "172.13 fb^{-1} (13.6 TeV)"


def read_limits(m):
    f = ROOT.TFile(f"{EXPDIR}/higgsCombine{m}.AsymptoticLimits.mH125.38.root")
    t = f.Get("limit")
    q = {}
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        q[round(t.quantileExpected, 3)] = t.limit
    f.Close()
    return q


def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    med = array.array('d'); x = array.array('d')
    g68_lo, g68_hi, g95_lo, g95_hi = (array.array('d') for _ in range(4))
    for m in MASSES:
        q = read_limits(m)
        x.append(float(m))
        med.append(q.get(0.5))
        g68_lo.append(q.get(0.16)); g68_hi.append(q.get(0.84))
        g95_lo.append(q.get(0.025)); g95_hi.append(q.get(0.975))

    n = len(x)
    ex = array.array('d', [0.0] * n)
    e68lo = array.array('d', [med[i] - g68_lo[i] for i in range(n)])
    e68hi = array.array('d', [g68_hi[i] - med[i] for i in range(n)])
    e95lo = array.array('d', [med[i] - g95_lo[i] for i in range(n)])
    e95hi = array.array('d', [g95_hi[i] - med[i] for i in range(n)])

    g95 = ROOT.TGraphAsymmErrors(n, x, med, ex, ex, e95lo, e95hi)
    g68 = ROOT.TGraphAsymmErrors(n, x, med, ex, ex, e68lo, e68hi)
    gmed = ROOT.TGraph(n, x, med)

    g95.SetFillColor(ROOT.kOrange); g95.SetLineColor(ROOT.kOrange)
    g68.SetFillColor(ROOT.kGreen + 1); g68.SetLineColor(ROOT.kGreen + 1)
    gmed.SetLineStyle(2); gmed.SetLineWidth(2); gmed.SetLineColor(ROOT.kBlack)

    c = ROOT.TCanvas("c", "", 800, 600)
    c.SetLogy(); c.SetTickx(1); c.SetTicky(1)
    c.SetLeftMargin(0.12); c.SetRightMargin(0.04); c.SetTopMargin(0.08); c.SetBottomMargin(0.12)

    ymin = min(g95_lo) * 0.6
    ymax = max(g95_hi) * 1.6
    frame = c.DrawFrame(0.0, ymin, 32.0, ymax)
    frame.GetXaxis().SetTitle("m_{a} [GeV]")
    frame.GetYaxis().SetTitle("95% CL expected limit on  #mu")
    frame.GetXaxis().SetTitleSize(0.045); frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0.04); frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.25)

    g95.Draw("3 SAME")
    g68.Draw("3 SAME")
    gmed.Draw("L SAME")
    ROOT.gPad.RedrawAxis()

    leg = ROOT.TLegend(0.52, 0.70, 0.88, 0.90)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.038)
    leg.AddEntry(gmed, "Median expected", "l")
    leg.AddEntry(g68, "68% expected", "f")
    leg.AddEntry(g95, "95% expected", "f")
    leg.Draw()

    tl = ROOT.TLatex(); tl.SetNDC(); tl.SetTextFont(42)
    tl.SetTextSize(0.05); tl.DrawLatex(0.12, 0.935, "#bf{CMS} #it{Preliminary}")
    tl.SetTextSize(0.042); tl.SetTextAlign(31); tl.DrawLatex(0.96, 0.935, LUMI)

    os.makedirs(OUTDIR, exist_ok=True)
    out = f"{OUTDIR}/limit_vs_mA_run3"
    c.SaveAs(out + ".pdf"); c.SaveAs(out + ".png")
    print("saved", out + ".pdf/.png")
    print(f"  mA1 expected median = {med[0]:.4f}")


if __name__ == "__main__":
    main()
