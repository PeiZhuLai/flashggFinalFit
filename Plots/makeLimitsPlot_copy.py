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

# 質量點
# massPoints = [5, 15, 30]  # 若需修改質量點，直接編輯此列表
massPoints = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]  # 若需修改質量點，直接編輯此列表
# 新增: 將原列表抽出為預設, 允許以 --masses 覆寫
DEFAULT_MASS_POINTS = massPoints[:]

# 常數
# https://arxiv.org/pdf/2402.09955
# ggf-xs 51960
# mh = 125.38  # H 質量 (GeV)
mh = 125.18  # H 質量 (GeV) Zebing
mz = 91.1876  # Z 質量 (GeV)
#https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR#TotalWidthAnchor
# gamma_HToSM = 4.143e-3  # H 總寬度 (GeV)
gamma_HToSM = 3.2e-3  # H 總寬度 (GeV) Zebing
decoupling_energy_scale = 1000.0  # GeV (1 TeV)
ZToll_br = 0.06729

def lamda_formula(x,y):
    return (1-x-y)**2 - 4*x*y

def calc_Wilson_coupling(br, ma):
    """
    根據 Branching Ratio 計算 Wilson coupling
    br: Branching Ratio (e.g. 1e-4)
    return: coupling (float)
    """
    # 反解 Wilson_Zh
    # gamma_HToZa = (mh**3 / (16 * math.pi)) * Wilson_Zh**2 / decoupling_energy_scale**2 * lamda_formula((mz/mh)**2, (ma/mh)**2)**1.5
    # br = gamma_HToZa / (gamma_HToSM + gamma_HToZa)

    # 反解 Wilson_Zh_divided_lambda
    # br * (gamma_HToSM + gamma_HToZa) = gamma_HToZa
    # br * gamma_HToSM = gamma_HToZa * (1 - br)
    # br * gamma_HToSM / (1 - br) = gamma_HToZa

    gamma_HToZa = (br * gamma_HToSM) / (1. - br)
    # Wilson_Zh = ( gamma_HToZa / ((mh**3 / (16. * math.pi)) * lamda_formula((mz/mh)**2, (ma/mh)**2)**1.5) * decoupling_energy_scale**2 ) ** 0.5 
    ratio = mh**3 / (16. * math.pi)
    Wilson_Zh = ( ( gamma_HToZa * decoupling_energy_scale**2 )  / ( ratio * lamda_formula((mz/mh)**2, (ma/mh)**2)**1.5) ) ** 0.5 

    return Wilson_Zh

def calc_indirect_Wilson(ma):
    """
    根據 Branching Ratio 計算 Wilson coupling
    br: Branching Ratio (e.g. 1e-4)
    return: coupling (float)
    """

    br = 0.08
    # 反解 Wilson_Zh
    # gamma_HToZa = (mh**3 / (16 * math.pi)) * Wilson_Zh**2 / decoupling_energy_scale**2 * lamda_formula((mz/mh)**2, (ma/mh)**2)**1.5
    # br = gamma_HToZa / (gamma_HToSM + gamma_HToZa)

    # 反解 Wilson_Zh_divided_lambda
    # br * (gamma_HToSM + gamma_HToZa) = gamma_HToZa
    # br * gamma_HToSM = gamma_HToZa * (1 - br)
    # br * gamma_HToSM / (1 - br) = gamma_HToZa

    gamma_HToZa = (br * gamma_HToSM) / (1. - br)
    Wilson_Zh = ( gamma_HToZa / ((mh**3 / (16. * math.pi)) * lamda_formula((mz/mh)**2, (ma/mh)**2)**1.5) * decoupling_energy_scale**2 ) ** 0.5 

    return Wilson_Zh

def read_limits_from_file(fname: str) -> Tuple[bool, float, float, float, float, float, float, bool]:
    """
    讀取 AsymptoticLimits ROOT 檔案:
    回傳 (ok, q025, q16, q50, q84, q975, obs, hasObs)
    若無法取得中位數，ok = False
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

    # 改為直接使用 PyROOT 動態屬性，不用 SetBranchAddress
    has_quant = bool(t.GetBranch("quantileExpected"))

    n = t.GetEntries()
    eps = 1e-3  # 放寬容差
    entries = []  # 收集 (quant, limit, idx)
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
        # 無 quant 分支時後面用 fallback

    # 若 has_quant 但仍未取得 q50，嘗試 fallback (依 Combine 標準順序)
    if q50 < 0:
        if debug:
            print(f"[read_limits_from_file][DEBUG] Fallback triggered for {fname}")
            for qv, lv, idx in entries:
                print(f"  Entry {idx}: quant={qv} limit={lv}")
        if has_quant and len(entries) >= 5:
            # 預期順序: (obs?) q025 q16 q50 q84 q975
            # 判斷第一個是否觀察值：quant < 0
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
        elif not has_quant:
            # 已在原邏輯處理，這裡不需再做
            pass

    f.Close()

    if q50 < 0:
        if debug:
            print(f"[read_limits_from_file][DEBUG] Entries summary (no median found):")
            for qv, lv, idx in entries:
                print(f"  {idx}: quant={qv} limit={lv}")
        print(f"[read_limits_from_file] Median (q50) not found in file: {fname}", file=sys.stderr)
        return False, q025, q16, q50, q84, q975, obs, hasObs
    return True, q025, q16, q50, q84, q975, obs, hasObs


def BrazilianPlots(sample: int = 0,
                   isInt: bool = True,
                   year: int = 0,
                   APV: bool = False,
                   drawObs: bool = True,
                   setLimitsOnBR: bool = False,
                   setLimitsOnWilsonCoefficient: bool = False,  # 新增: 以 Wilson 係數呈現
                   masses: List[int] = None,
                   outdir: str = "output_plots",
                   assume_xs: float = 100.0,
                   ggF_xs: float = 52170.0,
                   lumi_fb: float = 61.89,
                   formats: List[str] = None,
                   logy: bool = True,
                   save_root: bool = False,
                   tag_suffix: str = "",
                   indirect_min: float = None,      # 新增: 紅線最小 ma
                   indirect_max: float = None,      # 新增: 紅線最大 ma
                   indirect_step: float = None):    # 新增: 紅線步長
    """
    生成巴西圖 (Brazilian plot)
    setLimitsOnBR: True 則轉成 BR
    setLimitsOnWilsonCoefficient: True 則將 (sigma*BR)/ggF_xs 轉換成 Wilson coupling
    若為 Wilson 模式, 並提供 indirect_* 參數, 會用平滑 (min,max,step) 生成間接限制紅線
    """
    def make_file_name(m):
        return f"/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results/higgsCombine{m}.AsymptoticLimits.mH125.38.root"

    g_exp = TGraph()
    g_exp_1s = TGraphAsymmErrors()
    g_exp_2s = TGraphAsymmErrors()
    g_obs = TGraph()
    # 新增: 間接限制紅線 (僅 Wilson 模式)
    if setLimitsOnWilsonCoefficient:
        g_indirect = TGraph()

    idx = 0
    failed = []

    if masses is None:
        masses = massPoints
    if formats is None:
        formats = ["png", "pdf"]
    os.makedirs(outdir, exist_ok=True)

    for m in masses:
        fname = make_file_name(m)
        ok, q025, q16, q50, q84, q975, obs, hasObs = read_limits_from_file(fname)
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

        # 修正: Wilson 誤差用「轉換後差值」而非直接把差值丟進公式
        if setLimitsOnWilsonCoefficient:
            mid_br   = q50 * assume_xs / ( ggF_xs * ZToll_br )
            up1_br   = q84 * assume_xs / ( ggF_xs * ZToll_br )
            dn1_br   = q16 * assume_xs / ( ggF_xs * ZToll_br )
            up2_br   = q975 * assume_xs / ( ggF_xs * ZToll_br )
            dn2_br   = q025 * assume_xs / ( ggF_xs * ZToll_br )
            mid_coup = calc_Wilson_coupling(mid_br, m)
            up1_coup = calc_Wilson_coupling(up1_br, m)
            dn1_coup = calc_Wilson_coupling(dn1_br, m)
            up2_coup = calc_Wilson_coupling(up2_br, m)
            dn2_coup = calc_Wilson_coupling(dn2_br, m)
            exp = mid_coup
            p1s = max(0.0, up1_coup - mid_coup)
            m1s = max(0.0, mid_coup - dn1_coup)
            p2s = max(0.0, up2_coup - mid_coup)
            m2s = max(0.0, mid_coup - dn2_coup)

        g_exp.SetPoint(idx, m, exp)
        g_exp_1s.SetPoint(idx, m, exp)
        g_exp_2s.SetPoint(idx, m, exp)
        g_exp_1s.SetPointError(idx, 0, 0, m1s, p1s)
        g_exp_2s.SetPointError(idx, 0, 0, m2s, p2s)

        if drawObs and hasObs:
            if not setLimitsOnBR:
                g_obs.SetPoint(idx, m, obs * assume_xs)
            else:
                g_obs.SetPoint(idx, m, obs * assume_xs / ggF_xs)

        print(f"[BrazilianPlots] m={m} ({'BR' if setLimitsOnBR else 'XS'}) median={exp:.3g}")
        idx += 1

    if idx == 0:
        print("[BrazilianPlots] No valid points. Abort.", file=sys.stderr)
        if failed:
            print("  Tried masses:", failed, file=sys.stderr)
        return

    if failed:
        print("[BrazilianPlots] Missing masses:", failed)

    # 先決定 tag 以便使用唯一畫布名稱
    tag = "BR" if setLimitsOnBR else "XS"
    if setLimitsOnWilsonCoefficient:
        tag = "Wilson"
    if tag_suffix:
        tag = f"{tag}_{tag_suffix}"

    c = TCanvas(f"cLimits_{tag}", "", 800, 600)  # 改成唯一名稱避免 ROOT 警告
    c.SetBottomMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.15)

    xmin = min(masses) - 1
    xmax = max(masses) + 1

    if setLimitsOnWilsonCoefficient:
        ytitle = "|C^{eff}_{ZH}| [#frac{\Lambda}{1 TeV}]"
        ymax = 10.
        ymin = 1e-3
    else:
        if not setLimitsOnBR:
            ytitle = "#sigma(pp #rightarrow H) #times B(#rightarrow Za #rightarrow 2l + 2#gamma) [fb]"
            ymax = 100
            ymin = 8e-2
        else:
            ytitle = "Br(H #rightarrow Za #rightarrow 2l + 2#gamma)"
            ymax = 3e-2
            ymin = 1e-6

    frame = TH1F("frame", f";m_{{a}} (GeV);{ytitle}", 100, xmin, xmax)
    frame.SetStats(0)
    frame.SetMaximum(ymax)
    frame.SetMinimum(ymin)
    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.GetYaxis().CenterTitle(True)

    # 樣式
    g_exp.SetMarkerStyle(24)
    g_exp.SetMarkerColor(1)
    g_exp.SetMarkerSize(0.8)
    g_exp.SetLineColor(1)
    g_exp.SetLineWidth(3)
    g_exp.SetLineStyle(2)

    g_exp_1s.SetFillColor(3)
    g_exp_1s.SetLineColor(3)
    g_exp_1s.SetFillStyle(1001)

    g_exp_2s.SetFillColor(5)
    g_exp_2s.SetLineColor(5)
    g_exp_2s.SetFillStyle(1001)

    if drawObs:
        g_obs.SetMarkerStyle(20)
        g_obs.SetLineColor(1)
        g_obs.SetLineWidth(3)

    if setLimitsOnWilsonCoefficient:
        # 紅線樣式
        g_indirect.SetLineColor(2)
        g_indirect.SetLineWidth(3)
        g_indirect.SetLineStyle(1)

    if setLimitsOnWilsonCoefficient:
        x1, y1, x2, y2 = 0.59, 0.62, 0.97, 0.86
    else:
        if not setLimitsOnBR:
            x1, y1, x2, y2 = 0.59, 0.68, 0.97, 0.87
        else:
            x1, y1, x2, y2 = 0.58, 0.65, 0.96, 0.87

    leg = TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.045)
    leg.AddEntry(g_exp, "Median expected", "l")
    leg.AddEntry(g_exp_1s, "68% expected", "f")
    leg.AddEntry(g_exp_2s, "95% expected", "f")
    if setLimitsOnWilsonCoefficient:
        leg.AddEntry(g_indirect, "Indirect search", "l")
    # if drawObs: leg.AddEntry(g_obs, "Observed", "lp")

    frame.Draw()
    g_exp_2s.Draw("3 same")
    g_exp_1s.Draw("3 same")
    g_exp.Draw("L same")
    # 僅在有點時才畫觀察值，避免 "illegal number of points (0)"
    if drawObs and g_obs.GetN() > 0:
        g_obs.Draw("LP same")
    elif drawObs:
        print("[BrazilianPlots] No observed points to draw (skipping).")
    if setLimitsOnWilsonCoefficient:
        # 平滑紅線取樣
        _imin = indirect_min if indirect_min is not None else float(min(masses))
        _imax = indirect_max if indirect_max is not None else float(max(masses))
        _istep = indirect_step if indirect_step is not None else 0.25
        if _istep <= 0 or _imax <= _imin:
            # fallback
            _imin, _imax, _istep = float(min(masses)), float(max(masses)), max(0.25, (max(masses)-min(masses))/200.)
        j = 0
        x = _imin
        while x <= _imax + 1e-9:
            g_indirect.SetPoint(j, x, calc_indirect_Wilson(x))
            j += 1
            x += _istep
        g_indirect.Draw("L same")
    leg.Draw()

    if _has_cms_lumi and CMS_lumi:
        cms = CMS_lumi()
        # 與原本呼叫一致，iPosX=0 代表左上；可依需求調整
        cms.set_lumi(c, lumi_fb, 0, text="Preliminary", drawLumi=True)
    else:
        c.cd()
        import ROOT
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.DrawLatex(0.16, 0.92, f"CMS (unofficial)  {lumi_fb} fb^{{-1}} (13 TeV)")

    if setLimitsOnWilsonCoefficient:
        c.cd()
        import ROOT
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextFont(42)
        latex.DrawLatex(0.19, 0.81, f"B(a #rightarrow #gamma#gamma) = 100%")

    if logy:
        c.SetLogy()
    gPad.SetTicks(1, 1)
    gPad.RedrawAxis()

    for ext in formats:
        c.SaveAs(os.path.join(outdir, f"Limits_{tag}.{ext}"))

    if save_root:
        from ROOT import TFile
        rpath = os.path.join(outdir, f"Limits_{tag}.root")
        rf = TFile(rpath, "RECREATE")
        g_exp.Write("g_exp")
        g_exp_1s.Write("g_exp_1s")
        g_exp_2s.Write("g_exp_2s")
        if drawObs:
            g_obs.Write("g_obs")
        if setLimitsOnWilsonCoefficient:
            g_indirect.Write("g_indirect")
        rf.Close()
        print(f"[BrazilianPlots] Saved ROOT graphs to {rpath}")

def LimitPlots(year: int = 2018, APV: bool = False, setLimitsOnWilsonCoefficient: bool = False, **kwargs):
    BrazilianPlots(sample=-1, isInt=True, year=year, APV=APV, drawObs=True,
                   setLimitsOnBR=False, setLimitsOnWilsonCoefficient=setLimitsOnWilsonCoefficient, **kwargs)
    BrazilianPlots(sample=-1, isInt=True, year=year, APV=APV, drawObs=True,
                   setLimitsOnBR=True, setLimitsOnWilsonCoefficient=setLimitsOnWilsonCoefficient, **kwargs)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Make Brazilian limit plots (自動輸出 XS / BR / Wilson 三種)")
    parser.add_argument("-y", "--year", type=int, default=2018, help="Year")
    parser.add_argument("--apv", action="store_true", help="APV flag")
    parser.add_argument("--no-obs", action="store_true", help="Do not draw observed")
    # 移除 --br 與 --wilson 參數
    parser.add_argument("--outdir", default="output_plots", help="Output directory")
    parser.add_argument("--assume-xs", type=float, default=100.0, help="Assumed signal cross section (fb)")
    parser.add_argument("--ggf-xs", type=float, default=52170.0, help="ggF total cross section (fb)")
    parser.add_argument("--lumi", type=float, default=61.89, help="Luminosity in fb^-1 for label")
    parser.add_argument("--formats", default="png,pdf", help="Output formats, e.g. png,pdf,root")
    parser.add_argument("--linear-y", action="store_true", help="Use linear y-axis (default log)")
    parser.add_argument("--save-root", action="store_true", help="Also save TGraphs to ROOT file")
    parser.add_argument("--tag-suffix", default="", help="Extra tag suffix for output filenames")
    parser.add_argument("--masses", default="", help="逗號分隔質量點 (例: 5,15,30) 留空使用內建")
    parser.add_argument("--only", default="", help="只畫哪些: xs,br,wilson (逗號分隔), 留空=全部")
    parser.add_argument("--indirect-min", type=float, default=5, help="Wilson 紅線最小 ma (預設=min(masses))")
    parser.add_argument("--indirect-max", type=float, default=30, help="Wilson 紅線最大 ma (預設=max(masses))")
    parser.add_argument("--indirect-step", type=float, default=0.25, help="Wilson 紅線步長 (預設=0.25)")
    args = parser.parse_args()

    # 新增: masses override
    if args.masses.strip():
        try:
            masses = [int(x) for x in args.masses.split(",") if x.strip()]
        except ValueError:
            print("[main] 解析 --masses 失敗, 使用預設", flush=True)
            masses = DEFAULT_MASS_POINTS
    else:
        masses = DEFAULT_MASS_POINTS

    # 覆寫全域 massPoints 供後續 fallback 使用
    global massPoints
    massPoints = masses

    # 新增: 選擇輸出類型
    only_req = [s.strip().lower() for s in args.only.split(",") if s.strip()]
    valid_types = ["xs", "br", "wilson"]
    if only_req:
        plot_types = [t for t in only_req if t in valid_types]
        if not plot_types:
            print("[main] --only 無有效項目, 使用全部 (xs, br, wilson)")
            plot_types = valid_types
    else:
        plot_types = valid_types

    fmts = [f.strip() for f in args.formats.split(",") if f.strip()]
    save_root = args.save_root or ("root" in fmts)
    if "root" in fmts:
        fmts = [f for f in fmts if f != "root"]

    common_kwargs = dict(
        masses=masses,
        outdir=args.outdir,
        assume_xs=args.assume_xs,
        ggF_xs=args.ggf_xs,
        lumi_fb=args.lumi,
        formats=fmts,
        logy=not args.linear_y,
        save_root=save_root,
        tag_suffix=args.tag_suffix,
        indirect_min=args.indirect_min,       # 新增
        indirect_max=args.indirect_max,       # 新增
        indirect_step=args.indirect_step      # 新增
    )

    for pt in plot_types:
        if pt == "xs":
            BrazilianPlots(sample=-1, isInt=True, year=args.year, APV=args.apv,
                           drawObs=not args.no_obs, setLimitsOnBR=False,
                           setLimitsOnWilsonCoefficient=False, **common_kwargs)
        elif pt == "br":
            BrazilianPlots(sample=-1, isInt=True, year=args.year, APV=args.apv,
                           drawObs=not args.no_obs, setLimitsOnBR=True,
                           setLimitsOnWilsonCoefficient=False, **common_kwargs)
        elif pt == "wilson":
            BrazilianPlots(sample=-1, isInt=True, year=args.year, APV=args.apv,
                           drawObs=not args.no_obs, setLimitsOnBR=True,
                           setLimitsOnWilsonCoefficient=True, **common_kwargs)

if __name__ == "__main__":
    main()
