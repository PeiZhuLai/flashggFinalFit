import os
import sys
import math
from array import array
from typing import Tuple, List

try:
    from cmsLumi import CMS_lumi
    _has_cms_lumi = True
except Exception:
    CMS_lumi = None
    _has_cms_lumi = False

from ROOT import (
    TFile, TTree, TGraph, TGraphAsymmErrors, TH1F, TCanvas, TLegend,
    gSystem, gPad, gROOT
)
import json
from ROOT import TPad

run2XsLimitsJSON = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Plots/run2/xs_limits.json"
run2WilsonLimitsJSON = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Plots/run2/wilson_limits.json"

import ROOT
from ROOT import TLatex

# Global TLatex used by fallback drawer
latex = TLatex()
latex.SetNDC(True)
latex.SetTextFont(42)

# mh = 125.38
mh = 125.18  # H 質量 (GeV) Zebing
mz = 91.1876
# gamma_HToSM = 4.1e-3 # SM
gamma_HToSM = 3.2e-3  # H 總寬度 (GeV) Zebing

decoupling_energy_scale = 1000.0
ZToll_br = 0.06729

def lamda_formula(x, y):
    return (1 - x - y) ** 2 - 4 * x * y

def calc_Wilson_coupling(br, ma):
    gamma_HToZa = (br * gamma_HToSM) / (1.0 - br)
    ratio = mh**3 / (16.0 * math.pi)
    Wilson_Zh = ((gamma_HToZa * decoupling_energy_scale**2) / (ratio * lamda_formula((mz/mh)**2, (ma/mh)**2) ** 1.5)) ** 0.5
    return Wilson_Zh

def calc_indirect_Wilson(ma):
    br = 0.08
    gamma_HToZa = (br * gamma_HToSM) / (1.0 - br)
    Wilson_Zh = (gamma_HToZa / ((mh**3 / (16.0 * math.pi)) * lamda_formula((mz/mh)**2, (ma/mh)**2) ** 1.5) * decoupling_energy_scale**2) ** 0.5
    return Wilson_Zh

def read_limits_from_file(fname: str) -> Tuple[bool, float, float, float, float, float, float, bool]:
    """
    Return (ok, q025, q16, q50, q84, q975, obs, hasObs)
    """
    q025 = q16 = q50 = q84 = q975 = -1.0
    obs = -1.0
    hasObs = False

    if not os.path.exists(fname):
        print(f"[read_limits_from_file] File not found: {fname}", file=sys.stderr)
        return False, q025, q16, q50, q84, q975, obs, hasObs

    f = TFile.Open(fname, "READ")
    if not f or f.IsZombie():
        print(f"[read_limits_from_file] Cannot open file: {fname}", file=sys.stderr)
        return False, q025, q16, q50, q84, q975, obs, hasObs

    t = f.Get("limit")
    if not t:
        print(f"[read_limits_from_file] No TTree 'limit' in file: {fname}", file=sys.stderr)
        f.Close()
        return False, q025, q16, q50, q84, q975, obs, hasObs

    has_quant = bool(t.GetBranch("quantileExpected"))
    n = t.GetEntries()
    eps = 1e-3
    entries = []
    debug = os.environ.get("LIMITS_DEBUG", "") != ""

    for i in range(n):
        t.GetEntry(i)
        try:
            val = float(t.limit)
        except Exception:
            continue

        qv = None
        if has_quant:
            try:
                qv = float(t.quantileExpected)
            except Exception:
                qv = None

        entries.append((qv, val, i))

        if has_quant and qv is not None:
            if qv < 0:
                hasObs = True
                obs = val
            elif abs(qv - 0.500) < eps:
                q50 = val
            elif abs(qv - 0.160) < eps:
                q16 = val
            elif abs(qv - 0.840) < eps:
                q84 = val
            elif abs(qv - 0.025) < eps:
                q025 = val
            elif abs(qv - 0.975) < eps:
                q975 = val

    if q50 < 0 and has_quant and len(entries) >= 5:
        if debug:
            print(f"[read_limits_from_file][DEBUG] Fallback triggered for {fname}")
        offset = 1 if entries[0][0] is not None and entries[0][0] < 0 else 0
        try:
            q025 = entries[offset + 0][1]
            q16  = entries[offset + 1][1]
            q50  = entries[offset + 2][1]
            q84  = entries[offset + 3][1]
            q975 = entries[offset + 4][1]
            if entries[0][0] is not None and entries[0][0] < 0:
                hasObs = True
                obs = entries[0][1]
        except Exception:
            pass

    f.Close()

    if q50 < 0:
        if debug:
            print("[read_limits_from_file][DEBUG] Entries summary (no median found):")
            for qv, lv, idx in entries:
                print(f"  {idx}: quant={qv} limit={lv}")
        print(f"[read_limits_from_file] Median (q50) not found in file: {fname}", file=sys.stderr)
        return False, q025, q16, q50, q84, q975, obs, hasObs

    return True, q025, q16, q50, q84, q975, obs, hasObs

def _draw_cms_prelim_and_lumi(pad,
                             textSize: float = 0.05,
                             lumi_fb: float = 170.84,
                             sqrts_tev: float = 13.6,
                             left_margin: float = 0.05,
                             right_margin: float = 0.05,
                             top_margin: float = 0.08,
                             cms_text: str = "#bf{CMS} #it{Preliminary}",
                             drawLumi: bool = True):
    """
    Draw CMS label + lumi. Uses cmsLumi if available, otherwise TLatex fallback.
    Must be called after pad is created and cd()'d.
    """
    pad.cd()
    if _has_cms_lumi and CMS_lumi is not None:
        try:
            # Keep it simple/robust: delegate to CMS_lumi when present.
            # iPosX=0 usually means "out of frame"; adjust if you prefer.
            CMS_lumi(pad, 0, 0)
            return
        except Exception:
            pass

    # Fallback TLatex-based header
    latex.SetTextSize(textSize)
    latex.SetTextAlign(11)  # left-top
    latex.DrawLatex(left_margin, 1.0 - top_margin + 0.01, cms_text)

    if drawLumi:
        latex.SetTextAlign(31)  # right-top
        latex.DrawLatex(1.0 - right_margin, 1.0 - top_margin + 0.01,
                        f"{lumi_fb:.2f} fb^{{-1}} ({sqrts_tev:.1f} TeV)")
    else:
        latex.SetTextAlign(31)  # right-top
        latex.DrawLatex(1.0 - right_margin, 1.0 - top_margin + 0.01,
                        f"({sqrts_tev:.1f} TeV)")

def build_graph_from_hepdata_values(jmap: dict, masses: List[int], which: str = "exp") -> TGraph:
    """
    Build TGraph from the provided hepdata-like schema:
      - jmap["values"] is a list
      - each row contains:
          row["x"][0]["value"] -> mass
          row["y"] list with entries:
              entry["group"] == 0 (Observed) or 1 (Expected)
              entry["value"] -> y
    which: "exp" (group=1) or "obs" (group=0)
    """
    which = (which or "exp").lower()
    group = 1 if which in ("exp", "expected") else 0

    mass_map = {}
    for row in jmap.get("values", []):
        try:
            xarr = row.get("x", [])
            if not xarr:
                continue
            m = float(xarr[0].get("value"))
            yarr = row.get("y", [])
            yv = None
            for ent in yarr:
                if int(ent.get("group", -999)) == group:
                    yv = float(ent.get("value"))
                    break
            if yv is None:
                continue
            mass_map[m] = yv
        except Exception:
            continue

    g = TGraph()
    j = 0
    for m in masses:
        mv = float(m)
        if mv not in mass_map:
            continue
        g.SetPoint(j, mv, float(mass_map[mv]))
        j += 1
    return g

def _infer_masses_from_hepdata(jmap: dict) -> List[int]:
    ms = []
    for row in jmap.get("values", []):
        try:
            xarr = row.get("x", [])
            if not xarr:
                continue
            ms.append(int(round(float(xarr[0].get("value")))))
        except Exception:
            pass
    ms = sorted(list(dict.fromkeys(ms)))
    return ms

def build_graph_from_combine(masses: List[int], mode: str = "xs",
                            assume_xs: float = 100.0, ggF_xs: float = 52170.0):
    mode = (mode or "xs").lower()
    g = TGraph()
    failed = []
    j = 0

    def _fname(m):
        return f"/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results/higgsCombine{m}.AsymptoticLimits.mH125.38.root"

    for m in masses:
        ok, q025, q16, q50, q84, q975, obs, hasObs = read_limits_from_file(_fname(m))
        if not ok:
            failed.append(m)
            continue

        if mode == "xs":
            y = float(q50) * float(assume_xs)
        elif mode == "wilson":
            br = float(q50) * float(assume_xs) / (float(ggF_xs) * float(ZToll_br))
            y = float(calc_Wilson_coupling(br, m))
        else:
            raise ValueError(f"Unknown mode: {mode}")

        g.SetPoint(j, float(m), y)
        j += 1

    return g, failed

def make_comparison_plot(
    outpath_base: str,
    masses: List[int],
    g_this: TGraph,
    g_ref: TGraph,
    title_y: str,
    label_this: str,
    label_ref: str,
    formats: List[str] = None,
    logy: bool = True,
    y_max: float = None,
    y_min: float = None,
    title_size: float = 0.05,
    lumi_fb: float = 61.89,
    sqrts_tev: float = 13.6,
    header_left: str = "#bf{CMS} #it{Preliminary}",
    ratio_title: str = "(this-ref)/ref",
):
    """
    兩 pad：
      上 pad：this vs ref
      下 pad： (this - ref)/ref
    """
    if formats is None:
        formats = ["pdf"]

    def to_map(g):
        d = {}
        for i in range(g.GetN()):
            x = g.GetPointX(i)
            y = g.GetPointY(i)
            d[float(x)] = float(y)
        return d

    this_map = to_map(g_this)
    ref_map = to_map(g_ref)

    g_ratio = TGraph()
    j = 0
    for m in masses:
        mv = float(m)
        if mv not in this_map or mv not in ref_map:
            continue
        denom = ref_map[mv]
        if denom == 0:
            continue
        r = (this_map[mv] - denom) / denom
        g_ratio.SetPoint(j, mv, r)
        j += 1

    c = TCanvas(f"cComp_{os.path.basename(outpath_base)}", "", 800, 800)

    rightMargin = 0.05
    topMargin = 0.08
    leftMargin = 0.16

    c.SetLeftMargin(leftMargin)
    c.SetRightMargin(rightMargin)
    c.SetBottomMargin(0.135)
    c.SetTopMargin(topMargin)

    splitY = 0.4
    pad1 = TPad("pad1", "pad1", 0.0, splitY, 1.0, 0.98)
    pad2 = TPad("pad2", "pad2", 0.0, 0.02, 1.0, splitY)

    pad1.SetLeftMargin(leftMargin)
    pad1.SetRightMargin(rightMargin)
    pad1.SetBottomMargin(0.04)
    pad1.SetTopMargin(topMargin)
    pad1.SetTicks(1, 1)

    pad2.SetLeftMargin(leftMargin)
    pad2.SetRightMargin(rightMargin)
    pad2.SetBottomMargin(0.28)
    pad2.SetTopMargin(0.045)
    pad2.SetTicks(1, 1)

    pad2.Draw()
    pad1.Draw()

    xmin = min(masses) - 1
    xmax = max(masses) + 1

    pad1.cd()
    frame1 = TH1F("frame1", f";;{title_y}", 100, xmin, xmax)
    frame1.SetStats(0)
    frame1.GetYaxis().SetTitleSize(title_size)
    frame1.GetYaxis().SetLabelSize(0.075)
    frame1.GetYaxis().SetTitleOffset(1.10)
    frame1.GetXaxis().SetLabelSize(0)

    # vals = []
    # for m in masses:
    #     mv = float(m)
    #     if mv in this_map: vals.append(this_map[mv])
    #     if mv in ref_map:  vals.append(ref_map[mv])
    # if vals:
    #     vmin = min(vals)
    #     vmax = max(vals)
    #     if logy:
    #         ymin = max(1e-12, vmin * 0.8)
    #         ymax = vmax * 2.0
    #     else:
    #         span = (vmax - vmin) if vmax > vmin else (abs(vmax) if vmax != 0 else 1.0)
    #         ymin = vmin - 0.2 * span
    #         ymax = vmax + 0.3 * span
    #     frame1.SetMinimum(ymin)
    #     frame1.SetMaximum(ymax)

    frame1.SetMinimum(y_min)
    frame1.SetMaximum(y_max)
    frame1.Draw()

    g_this.SetLineColor(ROOT.TColor.GetColor("#E31A1C"))
    g_this.SetLineWidth(3)
    g_this.SetLineStyle(1)
    g_this.SetMarkerStyle(20)
    g_this.SetMarkerSize(1.3)
    g_this.SetMarkerColor(ROOT.TColor.GetColor("#E31A1C"))

    g_ref.SetLineColor(ROOT.TColor.GetColor("#1F78B4"))
    g_ref.SetLineWidth(3)
    g_ref.SetLineStyle(1)
    g_ref.SetMarkerStyle(21)
    g_ref.SetMarkerSize(1.3)
    g_ref.SetMarkerColor(ROOT.TColor.GetColor("#1F78B4"))

    g_ref.Draw("LP same")
    g_this.Draw("LP same")

    leg = TLegend(0.42, 0.55, 0.85, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.058)
    leg.AddEntry(g_ref, label_ref, "lp")
    leg.AddEntry(g_this, label_this, "lp")
    leg.Draw()

    if logy:
        pad1.SetLogy()

    _draw_cms_prelim_and_lumi(pad1, textSize=0.08, lumi_fb=lumi_fb, sqrts_tev=sqrts_tev,
                              left_margin=leftMargin, right_margin=rightMargin, top_margin=topMargin,
                              cms_text=header_left, drawLumi=False)

    pad2.cd()
    frame2 = TH1F("frame2", f";m_{{a}} (GeV);{ratio_title}", 100, xmin, xmax)
    frame2.SetStats(0)
    frame2.SetMinimum(-0.71)
    frame2.SetMaximum(+0.71)
    frame2.GetXaxis().SetTitleSize(0.132)
    frame2.GetXaxis().SetTitleOffset(1.05)
    frame2.GetXaxis().SetLabelSize(0.122)
    frame2.GetXaxis().SetLabelOffset(0.007)

    frame2.GetYaxis().SetTitleSize(0.10)
    frame2.GetYaxis().SetLabelSize(0.12)
    frame2.GetYaxis().SetTitleOffset(0.62)
    frame2.GetYaxis().SetNdivisions(505)
    frame2.GetYaxis().CenterTitle(True)
    frame2.Draw()

    g_ratio.SetLineColor(ROOT.TColor.GetColor("#E31A1C"))
    g_ratio.SetLineWidth(3)
    g_ratio.SetMarkerStyle(20)
    g_ratio.SetMarkerSize(1.3)
    g_ratio.SetMarkerColor(ROOT.TColor.GetColor("#E31A1C"))
    g_ratio.Draw("LP same")

    gPad.SetTicks(1, 1)
    gPad.RedrawAxis()

    for ext in formats:
        c.SaveAs(f"{outpath_base}.{ext}")

def read_hepdata_limits_json(path: str) -> dict:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with open(path, "r") as f:
        return json.load(f)

def BrazilianPlots(sample: int = 0,
                   isInt: bool = True,
                   year: int = 0,
                   APV: bool = False,
                   drawObs: bool = True,
                   setLimitsOnBR: bool = False,
                   setLimitsOnWilsonCoefficient: bool = False,
                   masses: List[int] = None,
                   outdir: str = "output_plots",
                   assume_xs: float = 100.0,
                   ggF_xs: float = 52170.0,
                   lumi_fb: float = 61.89,
                   formats: List[str] = None,
                   logy: bool = True,
                   save_root: bool = False,
                   tag_suffix: str = "",
                   indirect_min: float = None,
                   indirect_max: float = None,
                   indirect_step: float = None):

    def make_file_name(m):
        return f"/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results/higgsCombine{m}.AsymptoticLimits.mH125.38.root"

    if masses is None:
        masses = list(range(5, 31))
    if formats is None:
        formats = ["pdf"]

    os.makedirs(outdir, exist_ok=True)

    g_exp = TGraph()
    g_exp_1s = TGraphAsymmErrors()
    g_exp_2s = TGraphAsymmErrors()
    g_obs = TGraph()
    g_indirect = TGraph() if setLimitsOnWilsonCoefficient else None

    idx = 0
    failed = []

    for m in masses:
        ok, q025, q16, q50, q84, q975, obs, hasObs = read_limits_from_file(make_file_name(m))
        if not ok:
            failed.append(m)
            continue

        if not setLimitsOnBR:
            exp = q50 * assume_xs
            p1s = (q84 - q50) * assume_xs
            m1s = (q50 - q16) * assume_xs
            p2s = (q975 - q50) * assume_xs
            m2s = (q50 - q025) * assume_xs
        else:
            exp = q50 * assume_xs / ggF_xs
            p1s = (q84 - q50) * assume_xs / ggF_xs
            m1s = (q50 - q16) * assume_xs / ggF_xs
            p2s = (q975 - q50) * assume_xs / ggF_xs
            m2s = (q50 - q025) * assume_xs / ggF_xs

        if setLimitsOnWilsonCoefficient:
            mid_br = q50 * assume_xs / (ggF_xs * ZToll_br)
            up1_br = q84 * assume_xs / (ggF_xs * ZToll_br)
            dn1_br = q16 * assume_xs / (ggF_xs * ZToll_br)
            up2_br = q975 * assume_xs / (ggF_xs * ZToll_br)
            dn2_br = q025 * assume_xs / (ggF_xs * ZToll_br)

            mid_c = calc_Wilson_coupling(mid_br, m)
            up1_c = calc_Wilson_coupling(up1_br, m)
            dn1_c = calc_Wilson_coupling(dn1_br, m)
            up2_c = calc_Wilson_coupling(up2_br, m)
            dn2_c = calc_Wilson_coupling(dn2_br, m)

            exp = mid_c
            p1s = max(0.0, up1_c - mid_c)
            m1s = max(0.0, mid_c - dn1_c)
            p2s = max(0.0, up2_c - mid_c)
            m2s = max(0.0, mid_c - dn2_c)

        g_exp.SetPoint(idx, m, exp)
        g_exp_1s.SetPoint(idx, m, exp)
        g_exp_2s.SetPoint(idx, m, exp)
        g_exp_1s.SetPointError(idx, 0.0, 0.0, m1s, p1s)
        g_exp_2s.SetPointError(idx, 0.0, 0.0, m2s, p2s)

        if drawObs and hasObs and (not setLimitsOnWilsonCoefficient):
            yobs = obs * assume_xs if not setLimitsOnBR else obs * assume_xs / ggF_xs
            g_obs.SetPoint(idx, m, yobs)

        idx += 1

    if idx == 0:
        print("[BrazilianPlots] No valid points. Abort.", file=sys.stderr)
        if failed:
            print("  Missing:", failed, file=sys.stderr)
        return

    tag = "BR" if setLimitsOnBR else "XS"
    if setLimitsOnWilsonCoefficient:
        tag = "Wilson"
    if tag_suffix:
        tag = f"{tag}_{tag_suffix}"

    c = TCanvas(f"cLimits_{tag}", "", 800, 600)
    leftMargin = 0.15
    c.SetBottomMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(leftMargin)
    c.SetTopMargin(0.08)

    xmin = min(masses) - 1
    xmax = max(masses) + 1

    if setLimitsOnWilsonCoefficient:
        ytitle = "|C^{eff}_{ZH}| [#frac{#Lambda}{1 TeV}]"
        ymin, ymax = 1e-2, 10.0
    else:
        if not setLimitsOnBR:
            ytitle = "#sigma(pp #rightarrow H) #times B(#rightarrow Za #rightarrow 2l + 2#gamma) [fb]"
            ymin, ymax = 4e-1, 100.0
        else:
            ytitle = "Br(H #rightarrow Za #rightarrow 2l + 2#gamma)"
            ymin, ymax = 1e-6, 2e-2

    frame = TH1F("frame", f";m_{{a}} (GeV);{ytitle}", 100, xmin, xmax)
    frame.SetStats(0)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.Draw()

    g_exp_2s.SetFillColor(5)
    g_exp_2s.SetLineColor(5)
    g_exp_2s.SetFillStyle(1001)

    g_exp_1s.SetFillColor(3)
    g_exp_1s.SetLineColor(3)
    g_exp_1s.SetFillStyle(1001)

    g_exp.SetLineColor(1)
    g_exp.SetLineWidth(3)
    g_exp.SetLineStyle(2)

    g_exp_2s.Draw("3 same")
    g_exp_1s.Draw("3 same")
    g_exp.Draw("L same")

    if drawObs and g_obs.GetN() > 0 and (not setLimitsOnWilsonCoefficient):
        g_obs.SetMarkerStyle(20)
        g_obs.SetLineColor(1)
        g_obs.SetLineWidth(3)
        g_obs.Draw("LP same")

    if setLimitsOnWilsonCoefficient and g_indirect is not None:
        g_indirect.SetLineColor(2)
        g_indirect.SetLineWidth(3)

        _imin = indirect_min if indirect_min is not None else float(min(masses))
        _imax = indirect_max if indirect_max is not None else float(max(masses))
        _istep = indirect_step if indirect_step is not None else 0.25
        if _istep <= 0 or _imax <= _imin:
            _imin, _imax, _istep = float(min(masses)), float(max(masses)), 0.25

        j = 0
        x = _imin
        while x <= _imax + 1e-9:
            g_indirect.SetPoint(j, x, calc_indirect_Wilson(x))
            j += 1
            x += _istep
        g_indirect.Draw("L same")

    leg = TLegend(0.58, 0.68, 0.95, 0.87)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry(g_exp, "Median expected", "l")
    leg.AddEntry(g_exp_1s, "68% expected", "f")
    leg.AddEntry(g_exp_2s, "95% expected", "f")
    if setLimitsOnWilsonCoefficient and g_indirect is not None:
        leg.AddEntry(g_indirect, "Indirect search", "l")
    leg.Draw()

    if logy:
        c.SetLogy()
    gPad.SetTicks(1, 1)
    gPad.RedrawAxis()

    _draw_cms_prelim_and_lumi(c, textSize=0.05, lumi_fb=lumi_fb, sqrts_tev=13.6, left_margin=leftMargin, right_margin=0.05, top_margin=0.08, drawLumi=True)

    for ext in formats:
        c.SaveAs(os.path.join(outdir, f"Limits_{tag}.{ext}"))

    if save_root:
        rpath = os.path.join(outdir, f"Limits_{tag}.root")
        rf = TFile(rpath, "RECREATE")
        g_exp.Write("g_exp")
        g_exp_1s.Write("g_exp_1s")
        g_exp_2s.Write("g_exp_2s")
        if drawObs and g_obs.GetN() > 0:
            g_obs.Write("g_obs")
        if setLimitsOnWilsonCoefficient and g_indirect is not None:
            g_indirect.Write("g_indirect")
        rf.Close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Make Brazilian limit plots + optional Run2 comparison (輸出 5 張圖)")
    parser.add_argument("-y", "--year", type=int, default=2018, help="Year")
    parser.add_argument("--apv", action="store_true", help="APV flag")
    parser.add_argument("--no-obs", action="store_true", help="Do not draw observed")
    parser.add_argument("--outdir", default="output_plots", help="Output directory")
    parser.add_argument("--assume-xs", type=float, default=100.0, help="Assumed signal cross section (fb)")
    parser.add_argument("--ggf-xs", type=float, default=52170.0, help="ggF total cross section (fb)")
    parser.add_argument("--lumi", type=float, default=61.89, help="Luminosity in fb^-1 for label")
    parser.add_argument("--formats", default="pdf", help="Output formats, e.g. png,pdf,root")
    parser.add_argument("--linear-y", action="store_true", help="Use linear y-axis (default log)")
    parser.add_argument("--save-root", action="store_true", help="Also save TGraphs to ROOT file")
    parser.add_argument("--tag-suffix", default="", help="Extra tag suffix for output filenames")
    parser.add_argument("--masses", default="", help="逗號分隔質量點 (例: 5,15,30) 留空使用內建")
    parser.add_argument("--only", default="",
                        help="只輸出哪些: xs,br,wilson,compare (逗號分隔), 留空=全部(=5張)")
    parser.add_argument("--indirect-min", type=float, default=None, help="Wilson 紅線最小 ma")
    parser.add_argument("--indirect-max", type=float, default=None, help="Wilson 紅線最大 ma")
    parser.add_argument("--indirect-step", type=float, default=0.25, help="Wilson 紅線步長")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    formats = [x.strip() for x in (args.formats or "pdf").split(",") if x.strip()]
    logy = not args.linear_y

    if args.masses.strip():
        masses = [int(x) for x in args.masses.split(",") if x.strip()]
        masses = sorted(list(dict.fromkeys(masses)))
    else:
        masses = []
        if args.compare_run2:
            try:
                jxs = read_hepdata_limits_json(run2XsLimitsJSON)
                masses = _infer_masses_from_hepdata(jxs)
            except Exception:
                masses = []
        if not masses:
            masses = list(range(5, 31))

    only_req = [s.strip().lower() for s in (args.only or "").split(",") if s.strip()]
    valid = {"xs", "br", "wilson", "compare"}
    if only_req:
        req = [x for x in only_req if x in valid]
        if not req:
            req = ["xs", "br", "wilson", "compare"]
    else:
        req = ["xs", "br", "wilson", "compare"]

    save_root = args.save_root or ("root" in formats)
    if "root" in formats:
        formats = [f for f in formats if f != "root"]

    if "xs" in req:
        BrazilianPlots(
            masses=masses, outdir=args.outdir, assume_xs=args.assume_xs, ggF_xs=args.ggf_xs,
            lumi_fb=args.lumi, formats=formats, logy=logy, save_root=save_root, tag_suffix=args.tag_suffix,
            drawObs=not args.no_obs, setLimitsOnBR=False, setLimitsOnWilsonCoefficient=False,
        )
    if "br" in req:
        BrazilianPlots(
            masses=masses, outdir=args.outdir, assume_xs=args.assume_xs, ggF_xs=args.ggf_xs,
            lumi_fb=args.lumi, formats=formats, logy=logy, save_root=save_root, tag_suffix=args.tag_suffix,
            drawObs=not args.no_obs, setLimitsOnBR=True, setLimitsOnWilsonCoefficient=False,
        )
    if "wilson" in req:
        BrazilianPlots(
            masses=masses, outdir=args.outdir, assume_xs=args.assume_xs, ggF_xs=args.ggf_xs,
            lumi_fb=args.lumi, formats=formats, logy=logy, save_root=save_root, tag_suffix=args.tag_suffix,
            drawObs=not args.no_obs, setLimitsOnBR=True, setLimitsOnWilsonCoefficient=True,
            indirect_min=args.indirect_min, indirect_max=args.indirect_max, indirect_step=args.indirect_step,
        )

    if "compare" in req:
        args.compare_run2 = True

        j_xs = read_hepdata_limits_json(run2XsLimitsJSON)
        g_run2_xs_exp = build_graph_from_hepdata_values(j_xs, masses, which="exp")
        g_run2_xs_obs = build_graph_from_hepdata_values(j_xs, masses, which="obs")

        g_this_xs_exp, failed_xs = build_graph_from_combine(
            masses, mode="xs", assume_xs=args.assume_xs, ggF_xs=args.ggf_xs
        )

        outbase_xs = os.path.join(args.outdir, f"compare_xs_run2{args.tag_suffix}")
        make_comparison_plot(
            outpath_base=outbase_xs,
            masses=masses,
            g_this=g_this_xs_exp,
            g_ref=g_run2_xs_exp,
            title_y="#sigma(pp #rightarrow H) #times B(#rightarrow Za #rightarrow 2l + 2#gamma) [fb]",
            label_this="#splitline{Median expected}{(171 fb^{-1} Run3 2022+23+24)}",
            label_ref="#splitline{Median expected}{(138 fb^{-1} Full Run2)}",
            formats=formats,
            logy=logy,
            lumi_fb=args.lumi,
            y_max=100,
            y_min=0.4,
            title_size=0.06,
            sqrts_tev=13.6,
            header_left="#bf{CMS} #it{Preliminary}",
            ratio_title="#frac{Run3 - Run2}{Run2}",
        )

        j_w = read_hepdata_limits_json(run2WilsonLimitsJSON)
        g_run2_w_exp = build_graph_from_hepdata_values(j_w, masses, which="exp")
        g_run2_w_obs = build_graph_from_hepdata_values(j_w, masses, which="obs")

        g_this_w_exp, failed_w = build_graph_from_combine(
            masses, mode="wilson", assume_xs=args.assume_xs, ggF_xs=args.ggf_xs
        )

        outbase_w = os.path.join(args.outdir, f"compare_wilson_run2{args.tag_suffix}")
        make_comparison_plot(
            outpath_base=outbase_w,
            masses=masses,
            g_this=g_this_w_exp,
            g_ref=g_run2_w_exp,
            title_y="|C^{eff}_{ZH}| [#frac{#Lambda}{1 TeV}]",
            label_this="#splitline{Median expected}{(171 fb^{-1} Run3 2022+23+24)}",
            label_ref="#splitline{Median expected}{(138 fb^{-1} Full Run2)}",
            formats=formats,
            logy=logy,
            lumi_fb=args.lumi,
            y_max=10,
            y_min=1E-2,
            title_size=0.07,
            sqrts_tev=13.6,
            header_left="#bf{CMS} #it{Preliminary}",
            ratio_title="#frac{Run3 - Run2}{Run2}",
        )

        if failed_xs or failed_w:
            sys.stderr.write(
                "[compare-run2] missing combine points:"
                f" xs_failed={failed_xs} wilson_failed={failed_w}\n"
            )

if __name__ == "__main__":
    main()
