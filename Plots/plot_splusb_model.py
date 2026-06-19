#!/usr/bin/env python3
# S+B post-fit model plot for one mA, from a FitDiagnostics output (shapes_fit_s).
# Draws: observed data points, the post-fit S+B curve (total) and the B component
# (total_background) on the m_llgg axis, with the CMS Preliminary + lumi label.
# Usage: python3 plot_splusb_model.py <mA> <fitDiagnostics.root> <out.pdf> "<lumi fb^-1>" [blind=0] [blo=115] [bhi=135]
#   blind=1 blanks the observed data points in the (blo, bhi) GeV window (S+B/B curves stay full range).
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)

mA, infile, outpdf = sys.argv[1], sys.argv[2], sys.argv[3]
lumi = sys.argv[4] if len(sys.argv) > 4 else "172.13 fb^{-1}"
blind = (len(sys.argv) > 5 and sys.argv[5] not in ("0", "", "false", "False"))
blo = float(sys.argv[6]) if len(sys.argv) > 6 else 115.0
bhi = float(sys.argv[7]) if len(sys.argv) > 7 else 135.0

f = ROOT.TFile.Open(infile)
d = f.Get("shapes_fit_s/cat0")
if not d:
    print(f"[ERROR] shapes_fit_s/cat0 not in {infile}"); sys.exit(1)
data = d.Get("data"); bkg = d.Get("total_background"); sb = d.Get("total")

# Blinding: drop observed data points inside (blo, bhi). The S+B/B curves are unaffected.
if blind:
    g = ROOT.TGraphAsymmErrors(); j = 0
    for i in range(data.GetN()):
        x = data.GetPointX(i)
        if blo < x < bhi:
            continue
        g.SetPoint(j, x, data.GetPointY(i))
        g.SetPointError(j, data.GetErrorXlow(i), data.GetErrorXhigh(i),
                        data.GetErrorYlow(i), data.GetErrorYhigh(i))
        j += 1
    data = g

c = ROOT.TCanvas("c", "", 800, 700)
c.SetLeftMargin(0.13); c.SetRightMargin(0.05); c.SetTopMargin(0.07); c.SetBottomMargin(0.13)

bkg.SetLineColor(ROOT.kAzure+2); bkg.SetLineWidth(3); bkg.SetLineStyle(2); bkg.SetFillStyle(0)
sb.SetLineColor(ROOT.kRed+1);   sb.SetLineWidth(3);  sb.SetFillStyle(0)
bkg.SetTitle("")
bkg.GetXaxis().SetTitle("m_{ll#gamma#gamma} (GeV)"); bkg.GetXaxis().SetTitleSize(0.05); bkg.GetXaxis().SetLabelSize(0.04)
bkg.GetYaxis().SetTitle("Events / 1 GeV");           bkg.GetYaxis().SetTitleSize(0.05); bkg.GetYaxis().SetLabelSize(0.04)

ymax = max(bkg.GetMaximum(), sb.GetMaximum())
for i in range(data.GetN()):
    ymax = max(ymax, data.GetPointY(i) + data.GetErrorYhigh(i))
bkg.SetMaximum(ymax*1.4); bkg.SetMinimum(0.0)

bkg.Draw("HIST"); sb.Draw("HIST SAME")
data.SetMarkerStyle(20); data.SetMarkerSize(0.9); data.SetLineColor(ROOT.kBlack); data.Draw("PE SAME")

leg = ROOT.TLegend(0.58, 0.68, 0.93, 0.88)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.04)
leg.AddEntry(data, "Data", "PE"); leg.AddEntry(sb, "S+B fit", "L"); leg.AddEntry(bkg, "B component", "L")
leg.Draw()

lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextFont(42)
lat.SetTextSize(0.045); lat.SetTextAlign(11); lat.DrawLatex(0.13, 0.94, "#bf{CMS} #it{Preliminary}")
lat.SetTextSize(0.040); lat.SetTextAlign(31); lat.DrawLatex(0.95, 0.94, f"{lumi} (13.6 TeV)")
lat.SetTextSize(0.045); lat.SetTextAlign(11); lat.DrawLatex(0.17, 0.85, f"m_{{a}} = {mA} GeV")

c.SaveAs(outpdf)
print(f"[plot] {outpdf}")
