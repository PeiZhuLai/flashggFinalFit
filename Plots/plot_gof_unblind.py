#!/usr/bin/env python3
"""Per-mA goodness-of-fit plots on the unblinded data (after 1d_gof_unblind.sh).

For each mA reads the saturated-GOF combine outputs:
  unblind_GOF/roots/higgsCombine{mA}_obs.GoodnessOfFit.mH125.38.root        (observed t)
  unblind_GOF/roots/higgsCombine{mA}_toys.GoodnessOfFit.mH125.38.12345.root (toy t distribution)
draws the toy distribution + observed arrow, computes the GOF p-value
(= fraction of toys with t_toy >= t_obs), and saves:
  unblind_GOF/gof_mA{mA}.{pdf,png}     per-mA
  unblind_GOF/gof_pvalue_vs_mA.{pdf,png}  summary p-value vs mA
Run in higgs-alp-ana (ROOT 6.24).
"""
import os, glob, array, ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError + 1
ROOT.gStyle.SetOptStat(0)

BASE = os.path.dirname(os.path.abspath(__file__))
ROOTS = f"{BASE}/unblind_GOF/roots"
OUT = f"{BASE}/unblind_GOF"
MASSES = list(range(1, 31))
LUMI = "172.13 fb^{-1} (13.6 TeV)"


def read_t(path):
    if not os.path.exists(path):
        return None
    f = ROOT.TFile(path)
    t = f.Get("limit")
    vals = []
    for i in range(t.GetEntries()):
        t.GetEntry(i); vals.append(t.limit)
    f.Close()
    return vals


def main():
    os.makedirs(OUT, exist_ok=True)
    pmass, pval = array.array('d'), array.array('d')
    for m in MASSES:
        obs = read_t(f"{ROOTS}/higgsCombine{m}_obs.GoodnessOfFit.mH125.38.root")
        toyfiles = glob.glob(f"{ROOTS}/higgsCombine{m}_toys.GoodnessOfFit.mH125.38.*.root")
        if not obs or not toyfiles:
            print(f"[gof] mA{m}: missing inputs, skip"); continue
        t_obs = obs[0]
        toys = []
        for tf in toyfiles:
            toys += read_t(tf) or []
        toys = [t for t in toys if t == t]  # drop NaN
        if not toys:
            print(f"[gof] mA{m}: no toys, skip"); continue
        n_ge = sum(1 for t in toys if t >= t_obs)
        p = n_ge / float(len(toys))
        pmass.append(float(m)); pval.append(p)

        # per-mA distribution plot
        tmin = min(min(toys), t_obs) * 0.9
        tmax = max(max(toys), t_obs) * 1.1
        h = ROOT.TH1F(f"h{m}", ";saturated GOF test statistic;toys", 40, tmin, tmax)
        for t in toys:
            h.Fill(t)
        c = ROOT.TCanvas(f"c{m}", "", 800, 600)
        c.SetLeftMargin(0.12); c.SetRightMargin(0.04); c.SetTopMargin(0.08); c.SetBottomMargin(0.12)
        h.SetFillColor(ROOT.kAzure - 9); h.SetLineColor(ROOT.kAzure + 2)
        h.GetYaxis().SetRangeUser(0, h.GetMaximum() * 1.35)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.Draw("HIST")
        ln = ROOT.TLine(t_obs, 0, t_obs, h.GetMaximum() * 1.1)
        ln.SetLineColor(ROOT.kRed + 1); ln.SetLineWidth(3); ln.Draw()
        tl = ROOT.TLatex(); tl.SetNDC(); tl.SetTextFont(42)
        tl.SetTextSize(0.05); tl.DrawLatex(0.12, 0.935, "#bf{CMS} #it{Preliminary}")
        tl.SetTextSize(0.042); tl.SetTextAlign(31); tl.DrawLatex(0.96, 0.935, LUMI)
        tl.SetTextAlign(11); tl.SetTextSize(0.044)
        tl.DrawLatex(0.16, 0.86, f"m_{{a}} = {m} GeV")
        tl.SetTextColor(ROOT.kRed + 1); tl.SetTextSize(0.04)
        tl.DrawLatex(0.16, 0.80, f"observed t = {t_obs:.1f}")
        tl.DrawLatex(0.16, 0.75, f"GOF p-value = {p:.2f}")
        for ext in ("pdf", "png"):
            c.SaveAs(f"{OUT}/gof_mA{m:02d}.{ext}")
        print(f"[gof] mA{m}: t_obs={t_obs:.1f}  p={p:.3f}  ({len(toys)} toys)")

    # summary: p-value vs mA
    if len(pmass):
        g = ROOT.TGraph(len(pmass), pmass, pval)
        c = ROOT.TCanvas("csum", "", 850, 550)
        c.SetLeftMargin(0.12); c.SetRightMargin(0.04); c.SetTopMargin(0.08); c.SetBottomMargin(0.13)
        c.SetGridy()
        frame = c.DrawFrame(0, 0, 31, 1.05)
        frame.GetXaxis().SetTitle("m_{a} [GeV]"); frame.GetYaxis().SetTitle("GOF p-value (saturated)")
        frame.GetXaxis().SetTitleSize(0.05); frame.GetYaxis().SetTitleSize(0.05)
        g.SetMarkerStyle(20); g.SetMarkerSize(1.1); g.SetMarkerColor(ROOT.kBlack); g.SetLineColor(ROOT.kBlack)
        g.Draw("P SAME")
        l05 = ROOT.TLine(0, 0.05, 31, 0.05); l05.SetLineColor(ROOT.kRed + 1); l05.SetLineStyle(2); l05.SetLineWidth(2); l05.Draw()
        tl = ROOT.TLatex(); tl.SetNDC(); tl.SetTextFont(42)
        tl.SetTextSize(0.05); tl.DrawLatex(0.12, 0.935, "#bf{CMS} #it{Preliminary}")
        tl.SetTextSize(0.042); tl.SetTextAlign(31); tl.DrawLatex(0.96, 0.935, LUMI)
        for ext in ("pdf", "png"):
            c.SaveAs(f"{OUT}/gof_pvalue_vs_mA.{ext}")
        print(f"[gof] summary -> {OUT}/gof_pvalue_vs_mA.pdf")


if __name__ == "__main__":
    main()
