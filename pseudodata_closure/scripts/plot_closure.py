#!/usr/bin/env python3
"""Sculpting-closure figure: 95% CL limits on background-only MC pseudo-data vs the
Asimov expectation, plus the best-fit signal strength r-hat, per ALP mass.

Reads pseudodata_closure/closure_results.json (extract_results.py).
Top pad : expected median (dashed) with +/-1sigma (green) / +/-2sigma (yellow) band,
          and the OBSERVED-on-MC-pseudo-data limit (black markers). If the observed
          markers sit inside the band across mA, the BDT-sculpted DY spectrum does not
          fake a Higgs-mass peak -> no sculpting bias.
Bottom  : best-fit r-hat with asymmetric errors; a line at r=0. r-hat ~ 0 (esp. low-mA)
          confirms no spurious signal is pulled out of the sculpted background.

PyROOT, CMS style (CLAUDE.md): '#bf{CMS} #it{Preliminary}' top-left, lumi top-right.
Env: cmsenv (or any ROOT). Writes plots/ under pseudodata_closure.
"""
import json
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

PC = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/pseudodata_closure"
LUMI_LABEL = "172 fb^{-1} (13.6 TeV)"
LEFT, RIGHT = 0.13, 0.04
TSIZE, LSIZE = 0.055, 0.050


def cms_label(pad, lumi=LUMI_LABEL):
    pad.cd()
    lt = ROOT.TLatex(); lt.SetNDC(); lt.SetTextFont(42)
    lt.SetTextSize(0.05); lt.SetTextAlign(11)
    lt.DrawLatex(LEFT, 0.915, "#bf{CMS} #it{Preliminary}")
    lt.SetTextAlign(31)
    lt.DrawLatex(1 - RIGHT, 0.915, lumi)


def main():
    res = json.load(open(f"{PC}/closure_results.json"))
    mAs = sorted(int(k) for k in res)
    n = len(mAs)
    from array import array
    x = array('d', [float(m) for m in mAs]); ex = array('d', [0.0] * n)
    obs = array('d'); med = array('d')
    e1lo = array('d'); e1hi = array('d'); e2lo = array('d'); e2hi = array('d')
    rh = array('d'); rlo = array('d'); rhi = array('d'); zz = array('d')
    for m in mAs:
        L = res[str(m)]["limit"]
        obs.append(L["obs"]); med.append(L["exp_med"])
        e1lo.append(L["exp_med"] - L["exp_m1"]); e1hi.append(L["exp_p1"] - L["exp_med"])
        e2lo.append(L["exp_med"] - L["exp_m2"]); e2hi.append(L["exp_p2"] - L["exp_med"])
        rr = res[str(m)]["rhat"] or {"r": 0, "rLoErr": 0, "rHiErr": 0}
        rh.append(rr["r"]); rlo.append(abs(rr["rLoErr"])); rhi.append(abs(rr["rHiErr"]))
        zz.append(res[str(m)]["signif"] if res[str(m)]["signif"] is not None else 0.0)

    c = ROOT.TCanvas("c", "c", 900, 900)
    p1 = ROOT.TPad("p1", "p1", 0, 0.34, 1, 1); p1.SetBottomMargin(0.02)
    p1.SetLeftMargin(LEFT); p1.SetRightMargin(RIGHT); p1.SetTopMargin(0.09); p1.Draw()
    p2 = ROOT.TPad("p2", "p2", 0, 0, 1, 0.34); p2.SetTopMargin(0.03)
    p2.SetBottomMargin(0.30); p2.SetLeftMargin(LEFT); p2.SetRightMargin(RIGHT); p2.Draw()

    # ---- top: limits ----
    p1.cd()
    g2 = ROOT.TGraphAsymmErrors(n, x, med, ex, ex, e2lo, e2hi)
    g2.SetFillColor(ROOT.kOrange);
    g1 = ROOT.TGraphAsymmErrors(n, x, med, ex, ex, e1lo, e1hi)
    g1.SetFillColor(ROOT.kGreen + 1)
    gm = ROOT.TGraph(n, x, med); gm.SetLineStyle(2); gm.SetLineWidth(2); gm.SetLineColor(ROOT.kBlack)
    go = ROOT.TGraph(n, x, obs); go.SetMarkerStyle(20); go.SetMarkerSize(1.0); go.SetLineWidth(2)
    frame = p1.DrawFrame(0.5, 0.0, 30.5, max(obs.tolist() + [max(med)]) * 1.6)
    frame.GetYaxis().SetTitle("95% CL upper limit on #mu")
    frame.GetYaxis().SetTitleSize(TSIZE); frame.GetYaxis().SetLabelSize(LSIZE)
    frame.GetYaxis().SetTitleOffset(1.15)
    frame.GetXaxis().SetLabelSize(0)
    g2.Draw("3 same"); g1.Draw("3 same"); gm.Draw("L same"); go.Draw("LP same")
    leg = ROOT.TLegend(0.45, 0.60, 1 - RIGHT - 0.01, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040)
    leg.AddEntry(go, "Observed (bkg-only MC pseudo-data)", "LP")
    leg.AddEntry(gm, "Expected (Asimov)", "L")
    leg.AddEntry(g1, "Expected #pm1#sigma", "F")
    leg.AddEntry(g2, "Expected #pm2#sigma", "F")
    leg.Draw()
    cms_label(p1)

    # ---- bottom: r-hat ----
    p2.cd()
    gr = ROOT.TGraphAsymmErrors(n, x, rh, ex, ex, rlo, rhi)
    gr.SetMarkerStyle(20); gr.SetMarkerSize(0.9); gr.SetLineWidth(2)
    rmax = max(1.0, max(abs(rh[i]) + max(rlo[i], rhi[i]) for i in range(n)) * 1.2)
    fr2 = p2.DrawFrame(0.5, -rmax, 30.5, rmax)
    fr2.GetXaxis().SetTitle("m_{a} [GeV]")
    fr2.GetYaxis().SetTitle("#hat{#mu}")
    fr2.GetXaxis().SetTitleSize(0.11); fr2.GetXaxis().SetLabelSize(0.095)
    fr2.GetYaxis().SetTitleSize(0.11); fr2.GetYaxis().SetLabelSize(0.095)
    fr2.GetYaxis().SetTitleOffset(0.55); fr2.GetYaxis().SetNdivisions(505)
    fr2.GetXaxis().SetTitleOffset(1.05)
    l0 = ROOT.TLine(0.5, 0, 30.5, 0); l0.SetLineStyle(2); l0.SetLineColor(ROOT.kRed + 1); l0.Draw()
    gr.Draw("P same")

    os.makedirs(f"{PC}/plots", exist_ok=True)
    for ext in ("pdf", "png"):
        c.SaveAs(f"{PC}/plots/sculpting_closure_pseudodata.{ext}")
    print(f"[plot] -> {PC}/plots/sculpting_closure_pseudodata.pdf")

    # text summary
    print(f"\n{'mA':>3} {'obs':>8} {'exp':>8} {'in1s':>5} {'in2s':>5} {'rhat':>18} {'Z':>6}")
    n_out1 = n_out2 = 0
    for i, m in enumerate(mAs):
        L = res[str(m)]["limit"]
        in1 = L["exp_m1"] <= L["obs"] <= L["exp_p1"]
        in2 = L["exp_m2"] <= L["obs"] <= L["exp_p2"]
        n_out1 += (not in1); n_out2 += (not in2)
        print(f"{m:>3} {L['obs']:>8.4f} {L['exp_med']:>8.4f} {'Y' if in1 else '.':>5} "
              f"{'Y' if in2 else '.':>5} {rh[i]:>+8.3f}-{rlo[i]:.3f}+{rhi[i]:.3f} {zz[i]:>6.2f}")
    print(f"\nInside +/-1sigma: {n-n_out1}/{n}   inside +/-2sigma: {n-n_out2}/{n}")
    print(f"max |r-hat| = {max(abs(v) for v in rh):.3f}   max Z = {max(zz):.2f}")


if __name__ == "__main__":
    main()
