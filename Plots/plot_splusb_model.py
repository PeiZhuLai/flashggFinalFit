#!/usr/bin/env python3
# S+B post-fit model plot for one mA, drawn as CONTINUOUS post-fit pdf curves.
# Instead of the binned shapes_fit_s histograms (which look like a staircase), this
# loads the combine workspace, applies the post-fit parameter values + the post-fit
# discrete background pdf index from FitDiagnostics (fit_s), and plots the S+B and
# B-only pdfs as smooth RooCurves on the m_llgg axis, with the observed data points.
# Usage:
#   plot_splusb_model.py <mA> <workspace.root> <fitDiagnostics.root> <out.pdf> "<lumi fb^-1>" [blind=0] [blo=115] [bhi=135]
#   blind=1 blanks the observed data points in the (blo, bhi) GeV window (the curves stay full range).
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
RF = ROOT.RooFit

mA, wsfile, infile, outpdf = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
lumi  = sys.argv[5] if len(sys.argv) > 5 else "172.13 fb^{-1}"
blind = (len(sys.argv) > 6 and sys.argv[6] not in ("0", "", "false", "False"))
blo   = float(sys.argv[7]) if len(sys.argv) > 7 else 115.0
bhi   = float(sys.argv[8]) if len(sys.argv) > 8 else 135.0

# --- workspace: observable + data + S+B / B-only channel pdfs ---
wf = ROOT.TFile.Open(wsfile)
w = wf.Get("w")
if not w:
    print(f"[ERROR] no RooWorkspace 'w' in {wsfile}"); sys.exit(1)
x = w.var("CMS_hza_mass")
data = w.data("data_obs")
sb = w.pdf("pdf_bincat0_nuis")           # S+B channel pdf (RooAddPdf: signal + background)
b  = w.pdf("pdf_bincat0_bonly_nuis")     # B-only channel pdf
if not (x and data and sb and b):
    print(f"[ERROR] missing observable/data/pdf in {wsfile}"); sys.exit(1)

# --- apply the post-fit state from FitDiagnostics (fit_s) ---
ff = ROOT.TFile.Open(infile)
fit = ff.Get("fit_s")
if not fit:
    print(f"[ERROR] fit_s not in {infile}"); sys.exit(1)
w.allVars().assign(fit.floatParsFinal())                 # continuous post-fit params (incl r)
idxpar = fit.constPars().find("pdfindex_cat0_13p6TeV")   # discrete multipdf index (post-fit, kept const, RooCategory)
cat = w.cat("pdfindex_cat0_13p6TeV")
if idxpar and cat:
    cat.setIndex(idxpar.getIndex())

obs = ROOT.RooArgSet(x)
Nsb = sb.expectedEvents(obs)             # post-fit total S+B yield (reflects r-hat)
Nb  = b.expectedEvents(obs)             # post-fit background yield

# --- frame: continuous curves + data points (1 GeV bins to match "Events / 1 GeV") ---
xlo, xhi = x.getMin(), x.getMax()
nbin = int(round(xhi - xlo))
frame = x.frame(RF.Range(xlo, xhi), RF.Bins(nbin))
frame.SetTitle("")

# B (dashed azure) and S+B (solid red) as smooth post-fit pdf curves
b.plotOn(frame,  RF.Normalization(Nb,  ROOT.RooAbsReal.NumEvent),
         RF.LineColor(ROOT.kAzure + 2), RF.LineStyle(2), RF.LineWidth(3), RF.Name("bkg"))
sb.plotOn(frame, RF.Normalization(Nsb, ROOT.RooAbsReal.NumEvent),
          RF.LineColor(ROOT.kRed + 1), RF.LineWidth(3), RF.Name("sb"))

# observed data points (Poisson errors), drawn last so they sit on top
data.plotOn(frame, RF.Binning(nbin), RF.Name("dh"),
            RF.MarkerStyle(20), RF.MarkerSize(0.9), RF.LineColor(ROOT.kBlack))
dh = frame.getHist("dh")
if blind:                                # drop observed points inside (blo, bhi)
    i = 0
    while i < dh.GetN():
        if blo < dh.GetPointX(i) < bhi:
            dh.RemovePoint(i)
        else:
            i += 1

# y-range: cover data (incl. error bars) and the curve peak
ymax = 0.0
for i in range(dh.GetN()):
    ymax = max(ymax, dh.GetPointY(i) + dh.GetErrorYhigh(i))
csb = frame.getCurve("sb")
for i in range(csb.GetN()):
    ymax = max(ymax, csb.GetPointY(i))
frame.SetMaximum(ymax * 1.4); frame.SetMinimum(0.0)

frame.GetXaxis().SetTitle("m_{ll#gamma#gamma} (GeV)"); frame.GetXaxis().SetTitleSize(0.05); frame.GetXaxis().SetLabelSize(0.04)
frame.GetYaxis().SetTitle("Events / 1 GeV");           frame.GetYaxis().SetTitleSize(0.05); frame.GetYaxis().SetLabelSize(0.04)

c = ROOT.TCanvas("c", "", 800, 700)
c.SetLeftMargin(0.13); c.SetRightMargin(0.05); c.SetTopMargin(0.07); c.SetBottomMargin(0.13)
frame.Draw()

leg = ROOT.TLegend(0.58, 0.68, 0.93, 0.88)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.04)
leg.AddEntry(dh, "Data", "PE")
leg.AddEntry(frame.getCurve("sb"),  "S+B fit", "L")
leg.AddEntry(frame.getCurve("bkg"), "B component", "L")
leg.Draw()

lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextFont(42)
lat.SetTextSize(0.045); lat.SetTextAlign(11); lat.DrawLatex(0.13, 0.94, "#bf{CMS} #it{Preliminary}")
lat.SetTextSize(0.040); lat.SetTextAlign(31); lat.DrawLatex(0.95, 0.94, f"{lumi} (13.6 TeV)")
lat.SetTextSize(0.045); lat.SetTextAlign(11); lat.DrawLatex(0.17, 0.85, f"m_{{a}} = {mA} GeV")

c.SaveAs(outpdf)
print(f"[plot] {outpdf}  (S+B post-fit pdf curve; Nb={Nb:.1f}, Nsb={Nsb:.1f}, bkg index={cat.getIndex() if cat else 'NA'})")
