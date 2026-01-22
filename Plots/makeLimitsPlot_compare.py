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
from ROOT import TLatex  # <-- add

# Global TLatex used by fallback drawer
latex = TLatex()
latex.SetNDC(True)
latex.SetTextFont(42)

def _draw_cms_prelim_and_lumi(pad,
                             lumi_fb: float = 61.89,
                             sqrts_tev: float = 13.6,
                             right_margin: float = 0.05,
                             top_margin: float = 0.08,
                             cms_text: str = "#bf{CMS} #it{Preliminary}"):
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
    latex.SetTextSize(0.05)
    latex.SetTextAlign(11)  # left-top
    latex.DrawLatex(0.14, 1.0 - top_margin + 0.01, cms_text)

    latex.SetTextAlign(31)  # right-top
    latex.DrawLatex(1.0 - right_margin, 1.0 - top_margin + 0.01,
                    f"{lumi_fb:.2f} fb^{{-1}} ({sqrts_tev:.1f} TeV)")

# --- FIX: load BrazilianPlots reliably (avoid broken import stubs) ---
import importlib.util

def _load_brazilianplots_from_copy():
    """
    載入同目錄的 'makeLimitsPlot copy.py' 內的 BrazilianPlots。
    這樣可避免目前檔案內部曾被破壞的 import stub/片段導致 NameError。
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))
    copy_path = os.path.join(this_dir, "makeLimitsPlot copy.py")
    if not os.path.exists(copy_path):
        raise FileNotFoundError(f"Cannot find: {copy_path}")

    spec = importlib.util.spec_from_file_location("makeLimitsPlot_copy", copy_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot create module spec for: {copy_path}")

    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, "BrazilianPlots"):
        raise AttributeError(f"'BrazilianPlots' not found in: {copy_path}")
    return mod.BrazilianPlots

try:
    BrazilianPlots = _load_brazilianplots_from_copy()
except Exception as e:
    # 讓錯誤在啟動時就顯示，而不是跑到 main() 才 NameError
    raise RuntimeError(
        "Failed to load BrazilianPlots. "
        "Expected it in 'makeLimitsPlot copy.py' (same directory). "
        f"Original error: {e}"
    )

# --- add: helpers for Run2 JSON + combine extraction (used by --compare) ---
def read_hepdata_limits_json(path: str) -> dict:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with open(path, "r") as f:
        return json.load(f)

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

# --- keep old name for compatibility, but route to schema-aware builder ---
def build_graph_from_json(jmap: dict, masses: List[int], which: str = "exp") -> TGraph:
    return build_graph_from_hepdata_values(jmap, masses, which=which)

def build_graph_from_combine(masses: List[int], mode: str = "xs",
                            assume_xs: float = 100.0, ggF_xs: float = 52170.0):
    """
    Returns (TGraph median_expected, failed_masses)
    mode: "xs" -> q50*assume_xs
          "wilson" -> use BR->Wilson conversion from copy module
    """
    mode = (mode or "xs").lower()
    g = TGraph()
    failed = []
    j = 0

    # reuse functions from the loaded copy-module via BrazilianPlots.__globals__
    gbl = getattr(BrazilianPlots, "__globals__", {})
    read_limits_from_file = gbl.get("read_limits_from_file", None)
    calc_Wilson_coupling = gbl.get("calc_Wilson_coupling", None)
    ZToll_br = gbl.get("ZToll_br", 0.06729)

    if read_limits_from_file is None:
        raise RuntimeError("Cannot access read_limits_from_file from 'makeLimitsPlot copy.py' globals.")
    if mode == "wilson" and calc_Wilson_coupling is None:
        raise RuntimeError("Cannot access calc_Wilson_coupling from 'makeLimitsPlot copy.py' globals.")

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
            # replicate copy.py convention: BR = q50*assume_xs/(ggF_xs*ZToll_br)
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

    c = TCanvas(f"cComp_{os.path.basename(outpath_base)}", "", 900, 600)

    rightMargin = 0.05
    topMargin = 0.08
    leftMargin = 0.14

    c.SetLeftMargin(leftMargin)
    c.SetRightMargin(rightMargin)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(topMargin)

    splitY = 0.36
    pad1 = TPad("pad1", "pad1", 0.0, splitY, 1.0, 0.98)
    pad2 = TPad("pad2", "pad2", 0.0, 0.02, 1.0, splitY)

    pad1.SetLeftMargin(leftMargin)
    pad1.SetRightMargin(rightMargin)
    pad1.SetBottomMargin(0.05)
    pad1.SetTopMargin(topMargin)
    pad1.SetTicks(1, 1)

    pad2.SetLeftMargin(leftMargin)
    pad2.SetRightMargin(rightMargin)
    pad2.SetBottomMargin(0.40)
    pad2.SetTopMargin(0.045)
    pad2.SetTicks(1, 1)

    pad2.Draw()
    pad1.Draw()

    xmin = min(masses) - 1
    xmax = max(masses) + 1

    pad1.cd()
    frame1 = TH1F("frame1", f";;{title_y}", 100, xmin, xmax)
    frame1.SetStats(0)
    frame1.GetYaxis().SetTitleSize(0.055)
    frame1.GetYaxis().SetLabelSize(0.05)
    frame1.GetYaxis().SetTitleOffset(1.2)
    frame1.GetXaxis().SetLabelSize(0)

    vals = []
    for m in masses:
        mv = float(m)
        if mv in this_map: vals.append(this_map[mv])
        if mv in ref_map:  vals.append(ref_map[mv])
    if vals:
        vmin = min(vals)
        vmax = max(vals)
        if logy:
            ymin = max(1e-12, vmin * 0.5)
            ymax = vmax * 2.0
        else:
            span = (vmax - vmin) if vmax > vmin else (abs(vmax) if vmax != 0 else 1.0)
            ymin = vmin - 0.2 * span
            ymax = vmax + 0.3 * span
        frame1.SetMinimum(ymin)
        frame1.SetMaximum(ymax)

    frame1.Draw()

    g_this.SetLineColor(1)
    g_this.SetLineWidth(3)
    g_this.SetLineStyle(2)
    g_this.SetMarkerStyle(24)
    g_this.SetMarkerSize(1.0)
    g_this.SetMarkerColor(1)

    g_ref.SetLineColor(4)
    g_ref.SetLineWidth(3)
    g_ref.SetLineStyle(1)
    g_ref.SetMarkerStyle(20)
    g_ref.SetMarkerSize(1.0)
    g_ref.SetMarkerColor(4)

    g_ref.Draw("LP same")
    g_this.Draw("LP same")

    leg = TLegend(0.55, 0.72, 0.93, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry(g_this, label_this, "lp")
    leg.AddEntry(g_ref, label_ref, "lp")
    leg.Draw()

    if logy:
        pad1.SetLogy()

    _draw_cms_prelim_and_lumi(pad1, lumi_fb=lumi_fb, sqrts_tev=sqrts_tev,
                              right_margin=rightMargin, top_margin=topMargin,
                              cms_text=header_left)

    pad2.cd()
    frame2 = TH1F("frame2", f";m_{{a}} (GeV);{ratio_title}", 100, xmin, xmax)
    frame2.SetStats(0)
    frame2.SetMinimum(-1.0)
    frame2.SetMaximum(+1.0)
    frame2.GetXaxis().SetTitleSize(0.17)
    frame2.GetXaxis().SetLabelSize(0.15)
    frame2.GetYaxis().SetTitleSize(0.127)
    frame2.GetYaxis().SetLabelSize(0.12)
    frame2.GetYaxis().SetTitleOffset(0.55)
    frame2.GetYaxis().SetNdivisions(505)
    frame2.Draw()

    g_ratio.SetLineColor(2)
    g_ratio.SetLineWidth(3)
    g_ratio.SetMarkerStyle(20)
    g_ratio.SetMarkerSize(1.0)
    g_ratio.SetMarkerColor(2)
    g_ratio.Draw("LP same")

    gPad.SetTicks(1, 1)
    gPad.RedrawAxis()

    for ext in formats:
        c.SaveAs(f"{outpath_base}.{ext}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Make Brazilian limit plots (自動輸出 XS / BR / Wilson 三種)")
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

    parser.add_argument("--compare-run2", action="store_true",
                        help="另外輸出比較圖：this(combine) vs Run2 JSON (xs_limits / wilson_limits)")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    formats = [x.strip() for x in (args.formats or "pdf").split(",") if x.strip()]
    logy = not args.linear_y

    # ...existing code... (你原本的 BrazilianPlots 畫圖流程)

    # --- add: parse mass points ---
    if args.masses.strip():
        masses = [int(x) for x in args.masses.split(",") if x.strip()]
        masses = sorted(list(dict.fromkeys(masses)))
    else:
        # 若未指定，優先用 Run2 JSON 內建 mass grid，否則退回 copy.py 內建邏輯(你原本的)
        masses = []
        if args.compare_run2:
            try:
                jxs = read_hepdata_limits_json(run2XsLimitsJSON)
                masses = _infer_masses_from_hepdata(jxs)
            except Exception:
                masses = []
        if not masses:
            # 你可改成你原本使用的 mass 列表來源
            masses = list(range(5, 31))  # minimal fallback

    # --- add: make compare plots ---
    if args.compare_run2:
        # 1) XS compare: this (combine xs) vs run2 xs_limits.json
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
            title_y="95% CL UL on #sigma#timesBR (fb)",
            label_this="This (combine) expected",
            label_ref="Run2 JSON expected",
            formats=formats,
            logy=logy,
            lumi_fb=args.lumi,
            sqrts_tev=13.6,
            header_left="#bf{CMS} #it{Preliminary}",
            ratio_title="(this-run2)/run2",
        )

        # 2) Wilson compare: this (combine->wilson) vs run2 wilson_limits.json
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
            title_y="95% CL UL on C_{Zh}^{eff} (#Lambda/1 TeV)",
            label_this="This (combine) expected",
            label_ref="Run2 JSON expected",
            formats=formats,
            logy=logy,
            lumi_fb=args.lumi,
            sqrts_tev=13.6,
            header_left="#bf{CMS} #it{Preliminary}",
            ratio_title="(this-run2)/run2",
        )

        if failed_xs or failed_w:
            sys.stderr.write(
                "[compare-run2] missing combine points:"
                f" xs_failed={failed_xs} wilson_failed={failed_w}\n"
            )

if __name__ == "__main__":
    main()
