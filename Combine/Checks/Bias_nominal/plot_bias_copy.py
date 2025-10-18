# plot_pull_overlay_paths.py
import os, glob, re, argparse
import json
import ROOT
ROOT.gROOT.SetBatch(True)

# 新增：對齊 makeLimitsPlot.py 的 cmsLumi 匯入與旗標
try:
    from cmsLumi import CMS_lumi
    _has_cms_lumi = True
except Exception:
    CMS_lumi = None
    _has_cms_lumi = False

# === 1) 在這裡直接寫你的絕對路徑 ===
BIAS_BASE = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal"
PATH_BIAS_FITS = os.path.join(BIAS_BASE, "BiasFits")
PATH_BIAS_TOYS = os.path.join(BIAS_BASE, "BiasToys")   # 若之後要用到 toys，可用這個路徑

# 參數：注入真值（Asimov 常用 1.0）
R_TRUE = 1.0
TITLE  = "Bias study pull comparison"
OUTTAG = "pull_overlay"

# mean_vs_r 的軸範圍（寫在前面、直接改數值）
MEAN_VS_R_XMIN = -0.2
MEAN_VS_R_XMAX =  0.2
# 若想自動，設為 None；否則給數值（例如 0.6, 1.4）
MEAN_VS_R_YMIN = 0.1
MEAN_VS_R_YMAX = 1.9

# 新增：第一張 pull 疊圖的 X 軸範圍（宣告在前面即可調整；None 表示沿用預設 -5~5）
PULL_HIST_XMIN = -4.8   # 例如 -2.0
PULL_HIST_XMAX = 4.8   # 例如  2.0

def file_label(path):
    base = os.path.basename(path)
    m = re.search(r"biasStudy_(.+?)_fits\.root", base)
    return m.group(1) if m else base

def parse_truth_from_filename(path, default=R_TRUE):
    """
    嘗試從檔名解析 truth r（支援如 r1p0、r1.0、truth1p0、mu1p0 等），失敗則回傳 default。
    """
    base = os.path.basename(path)
    patterns = [
        r"(?:r(?:True|Inject(?:ed)?)?|truth|mu)[:=_-]?([0-9]+(?:p[0-9]+)?|\d+\.\d+)",
        r"(?:r)[:=_-]?([0-9]+)"
    ]
    for pat in patterns:
        m = re.search(pat, base, flags=re.IGNORECASE)
        if m:
            s = m.group(1).replace('p', '.')
            try:
                return float(s)
            except Exception:
                pass
    return default

# 新增：將 token（例如 "5", "1p0", "1.0", "r1p0"）解析為真實 mu 值
def parse_mu_from_token(token, default=R_TRUE):
    s = str(token)
    s = s.replace('p', '.')
    s = re.sub(r"^[a-zA-Z_]+", "", s)  # 去除可能的前綴字母
    try:
        # 純數字視為以 0.1 為單位（例如 "5" -> 0.5，與現有 "5_gaussfit.json" 命名相容）
        if re.fullmatch(r"\d+", s):
            return float(s) / 10.0
        return float(s)
    except Exception:
        return default

# 新增：讀入 gaussfit JSON，取出每個模型的 mean
def load_gaussfit_json(json_path):
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
        # 期望結構：{ "Bern1": {"mean": ..., "sigma": ...}, ... }
        means = {k: v.get("mean", None) for k, v in data.items() if isinstance(v, dict) and "mean" in v}
        # 過濾 None
        return {k: m for k, m in means.items() if m is not None}
    except Exception as e:
        print(f"[WARN] cannot load JSON: {json_path} ({e})")
        return {}

def _resolve_val_err(tree):
    """
    依據可用欄位決定用來算 pull 的 (value, error) 欄位名稱。
    優先順序（偏好 r，因為我們研究的是 r 的 bias）：
      1) r + limitErr
      2) r + (rErr 或 rErrP/rErrM)
      3) r（無逐事件誤差，後續以全體 RMS 做 fallback）
      4) limit + limitErr
      5) limit（僅最後備用）
    回傳: (val_name, err_name or None, has_quantile)
    """
    leaves = [l.GetName() for l in tree.GetListOfLeaves()]
    has_q = ("quantileExpected" in leaves)
    # 優先用 r
    if "r" in leaves and "limitErr" in leaves:
        return ("r", "limitErr", has_q)
    if "r" in leaves and ("rErr" in leaves or ("rErrP" in leaves and "rErrM" in leaves)):
        return ("r", None, has_q)
    if "r" in leaves:
        return ("r", None, has_q)
    # 最後才考慮 limit
    if "limit" in leaves and "limitErr" in leaves:
        return ("limit", "limitErr", has_q)
    if "limit" in leaves:
        return ("limit", None, has_q)
    return (None, None, has_q)

def _list_leaves(tree):
    try:
        if not tree:
            return []
        return [l.GetName() for l in tree.GetListOfLeaves()] or []
    except Exception:
        return []

def get_pulls(fname, r_true=1.0):
    import math
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        print(f"[WARN] cannot open {fname}")
        return []
    t = f.Get("limit")
    if not t:
        print(f"[WARN] no TTree 'limit' in {fname}")
        f.Close()
        return []

    leaves_list = _list_leaves(t)
    vname, ename, has_q = _resolve_val_err(t)
    if not vname:
        print(f"[WARN] no suitable (value,error) branches in {fname}. Leaves={leaves_list}")
        f.Close()
        return []

    def compute_pulls(apply_quantile_filter=True):
        pulls = []
        raw_vals = []
        used_r_vals = []
        n_total = int(t.GetEntriesFast())
        n_kept = 0
        n_qrej = 0
        n_poserr = 0
        has_r_branch = ("r" in leaves_list)
        for i in range(n_total):
            t.GetEntry(i)
            if apply_quantile_filter and has_q:
                try:
                    if getattr(t, "quantileExpected") >= 0:
                        n_qrej += 1
                        continue
                except Exception:
                    pass
            val = getattr(t, vname, float("nan"))
            raw_vals.append(val)
            if ename:
                err = getattr(t, ename, 0.0)
            else:
                err = 0.0
                if hasattr(t, "rErr"):
                    err = getattr(t, "rErr", 0.0)
                else:
                    rp = getattr(t, "rErrP", 0.0) if hasattr(t, "rErrP") else 0.0
                    rm = abs(getattr(t, "rErrM", 0.0)) if hasattr(t, "rErrM") else 0.0
                    if (rp > 0) or (rm > 0):
                        err = 0.5 * (rp + rm)
            if err and err > 0 and math.isfinite(err) and math.isfinite(val):
                pulls.append((val - r_true) / err)
                n_kept += 1
                n_poserr += 1
                if has_r_branch:
                    rv = getattr(t, "r", float("nan"))
                    if math.isfinite(rv):
                        used_r_vals.append(rv)
        return pulls, raw_vals, used_r_vals, n_total, n_kept, n_qrej, n_poserr

    # 第一輪：若有 quantileExpected，採用 quantile 過濾；否則不過濾
    apply_q = bool(has_q)
    pulls, raw_vals, r_used, n_total, n_kept, n_qrej, n_poserr = compute_pulls(apply_quantile_filter=apply_q)

    # 若全被 quantile 過濾掉，第二輪不套 quantile 過濾
    used_fallback = False
    mode = "with-quantile" if apply_q else "no-quantile"
    if has_q and n_kept == 0:
        pulls, raw_vals, r_used, n_total, n_kept, _, n_poserr = compute_pulls(apply_quantile_filter=False)
        mode = "fallback(no-quantile)"

    # 若仍無逐事件誤差可用（或全為 0），對 r 使用全體 RMS 做標準化
    if n_kept == 0 and vname == "r":
        # 計算 r 相對 r_true 的樣本標準差
        vals = [v for v in raw_vals if math.isfinite(v)]
        if len(vals) > 1:
            mean_diff_sq = sum((v - r_true) * (v - r_true) for v in vals) / float(len(vals))
            std = math.sqrt(mean_diff_sq) if mean_diff_sq > 0 else 0.0
            if std > 0:
                pulls = [ (v - r_true) / std for v in vals ]
                n_kept = len(pulls)
                used_fallback = True
                mode = "global-rms"
                r_used = list(vals)  # 同步記錄對應的 r 值
        # 若 std==0 或無有效值則保持空陣列

    f.Close()
    tag = f"{mode}"
    if used_fallback:
        tag += " (no per-entry err)"
    print(f"[INFO] {os.path.basename(fname)}: entries={n_total}, pulls_kept={n_kept}, mode={tag}")
    if n_kept == 0:
        print(f"[HINT] Available leaves: {leaves_list}")
    return pulls, r_used

def make_hist(values, name, bins=48, xmin=-5, xmax=5):
    # 依據前段常數覆蓋 X 範圍（若為 None 則保留函式預設）
    xmin_final = PULL_HIST_XMIN if PULL_HIST_XMIN is not None else xmin
    xmax_final = PULL_HIST_XMAX if PULL_HIST_XMAX is not None else xmax

    h = ROOT.TH1F(name, name, bins, xmin_final, xmax_final)
    for v in values: h.Fill(v)
    if h.Integral() > 0:
        # 歸一化（比較形狀）
        h.Scale(1.0 / h.Integral("width"))
    h.GetXaxis().SetTitle(r"(#mu - #mu_{True}) / #sigma_{#mu}")
    h.GetYaxis().SetTitle(f"A.U. / {h.GetBinWidth(1):.2f}")
    # 對齊 makeLimitsPlot.py 的字型大小與版面
    h.GetXaxis().SetTitleSize(0.055)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetYaxis().CenterTitle(True)
    h.GetXaxis().CenterTitle(True)
    return h

# 新增：與 makeLimitsPlot.py 一致的 CMS/Lumi 繪製輔助函式
def draw_cms_lumi(canvas, lumi_fb):
    """
    在指定的 canvas 上繪製 'CMS Preliminary' 與亮度文字。
    若有 cmsLumi，使用 cmsLumi；否則以 TLatex fallback。
    """
    if _has_cms_lumi and CMS_lumi:
        cms = CMS_lumi()
        cms.set_lumi(canvas, float(lumi_fb), 0, text="Preliminary", drawLumi=True)
    else:
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.DrawLatex(0.16, 0.92, f"CMS (unofficial)  {float(lumi_fb)} fb^{{-1}} (13 TeV)")

def parse_args():
    p = argparse.ArgumentParser(description="Plot pull overlay")
    p.add_argument("--mA", type=str, default=None, help="Prefix for output plots")
    # 新增：亮度參數（預設與 makeLimitsPlot.py 相同）
    p.add_argument("--lumi", type=float, default=61.89, help="Luminosity in fb^-1 for label")
    # 移除先前加入的 --mean-xmin/--mean-xmax
    return p.parse_args()

def main():
    args = parse_args()

    # 新增：在解析參數後才組 JSON 路徑，並嘗試讀取
    json_means = {}
    if args.mA:
        json_path = os.path.join(BIAS_BASE, "BiasPlots", f"{args.mA}_gaussfit.json")
        json_means = load_gaussfit_json(json_path)
        if json_means:
            print(f"[INFO] loaded JSON means from {json_path}")

    if not os.path.isdir(PATH_BIAS_FITS):
        raise SystemExit(f"[ERR] Not found: {PATH_BIAS_FITS}")

    # 抓 BiasFits 下的所有結果
    files = sorted(glob.glob(os.path.join(PATH_BIAS_FITS, "biasStudy_*_fits.root")))
    # 新增：排除含有 'split' 的檔案
    files = [f for f in files if "split" not in os.path.basename(f)]
    if not files:
        raise SystemExit(f"[ERR] No files in {PATH_BIAS_FITS}/biasStudy_*_fits.root")

    c = ROOT.TCanvas("c","c",800,600)
    # 與 makeLimitsPlot.py 對齊的邊界
    c.SetTopMargin(0.09)
    c.SetBottomMargin(0.14)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.13)

    leg = ROOT.TLegend(0.52, 0.62, 0.80, 0.87)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.SetTextSize(0.043)
    ROOT.gStyle.SetOptStat(0)

    colors = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kMagenta+1,
              ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kViolet, ROOT.kTeal+2]
    styles = [1,2,3,4,5,6,7,8]

    hists = []
    scatter_points = []  # [(mu_pull, r_truth, label, color)]
    for i, fpath in enumerate(files):
        pulls, _ = get_pulls(fpath, R_TRUE)
        if not pulls:
            print(f"[INFO] skip empty: {fpath}")
            continue
        h = make_hist(pulls, f"h_{i}")
        color = colors[i % len(colors)]
        h.SetLineColor(color)
        h.SetLineWidth(3)
        h.SetLineStyle(styles[i % len(styles)])
        label = file_label(fpath)
        hists.append((h, label))

    if not hists:
        raise SystemExit("[ERR] nothing to plot.")

    ymax = max(h.GetMaximum() for h,_ in hists) * 1.7
    for idx, (h, label) in enumerate(hists):
        h.SetMaximum(ymax)
        if idx == 0:
            h.SetTitle("")
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")

        fit = ROOT.TF1(f"g_{idx}", "gaus", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        fit.SetLineColor(h.GetLineColor()); fit.SetLineStyle(1)
        h.Fit(fit, "RQ")
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)
        leg.AddEntry(h, f"{label}  (#mu={mean:+.2f}, #sigma={sigma:.2f})", "l")

        # 收集散點：x= pull mean，y= 該檔案的 truth r（從檔名解析，找不到則用 R_TRUE）
        truth_y = parse_truth_from_filename(files[idx], R_TRUE)
        scatter_points.append((mean, truth_y, label, h.GetLineColor()))

    leg.Draw()
    z = ROOT.TLine(0, 0, 0, ymax); z.SetLineStyle(7); z.SetLineColor(ROOT.kGray+2); z.Draw()

    outtag = f"{args.mA}_{OUTTAG}" if args.mA else OUTTAG
    os.makedirs("plots_bias", exist_ok=True)

    # 與 makeLimitsPlot.py 一致的 ticks 與 redraw
    ROOT.gPad.SetTicks(1, 1)
    ROOT.gPad.RedrawAxis()

    # 新增：在第一張圖加上 CMS/Lumi
    draw_cms_lumi(c, args.lumi)

    c.SaveAs(f"plots_bias/{outtag}.png")
    c.SaveAs(f"plots_bias/{outtag}.pdf")

    # ===== 第二張圖：pull mean vs r truth =====
    # 改為優先使用 JSON 的 mean；若 JSON 讀不到，才使用上面 histogram 擬合的結果
    points = []
    if json_means:
        mu_truth = parse_mu_from_token(args.mA, R_TRUE)
        # 將 JSON 中各模型的 mean 映射為點 (mean, mu_truth)
        for j, (lab, mean) in enumerate(json_means.items()):
            col = colors[j % len(colors)]
            points.append((mean, mu_truth, lab, col))
    else:
        # 回退：使用先前的 histogram 擬合結果（保持舊邏輯）
        points = scatter_points

    if points:
        c2 = ROOT.TCanvas("c2","c2",800,600)
        # 與 makeLimitsPlot.py 對齊的邊界
        c2.SetTopMargin(0.09)
        c2.SetBottomMargin(0.14)
        c2.SetRightMargin(0.05)
        c2.SetLeftMargin(0.13)

        mg = ROOT.TMultiGraph()
        leg2 = ROOT.TLegend(0.65, 0.61, 0.93, 0.86)
        leg2.SetBorderSize(0); leg2.SetFillStyle(0)
        leg2.SetTextSize(0.045)

        for j, (xmu, yr, lab, col) in enumerate(points):
            g = ROOT.TGraph(1)
            g.SetPoint(0, xmu, yr)
            g.SetMarkerColor(col)
            g.SetMarkerStyle(20 + (j % 10))
            g.SetMarkerSize(2.0)
            mg.Add(g, "P")
            leg2.AddEntry(g, lab, "p")

        mg.SetTitle("")
        mg.Draw("A P")

        # 直接用常數套用 X/Y 軸範圍並重畫
        mg.GetXaxis().SetTitle("Pull Mean")
        mg.GetYaxis().SetTitle("True #mu")
        mg.GetXaxis().SetLimits(MEAN_VS_R_XMIN, MEAN_VS_R_XMAX)
        if (MEAN_VS_R_YMIN is not None) or (MEAN_VS_R_YMAX is not None):
            # 用 SetMinimum/SetMaximum 控制 Y 軸
            if MEAN_VS_R_YMIN is not None: mg.SetMinimum(MEAN_VS_R_YMIN)
            if MEAN_VS_R_YMAX is not None: mg.SetMaximum(MEAN_VS_R_YMAX)
        mg.Draw("A P")
        ROOT.gPad.Update()

        # 對齊 makeLimitsPlot.py 的字型大小與版面
        mg.GetXaxis().SetTitleSize(0.055)
        mg.GetXaxis().SetLabelSize(0.05)
        mg.GetXaxis().SetTitleOffset(1.1)
        mg.GetYaxis().SetTitleSize(0.05)
        mg.GetYaxis().SetLabelSize(0.05)
        mg.GetYaxis().SetTitleOffset(1.3)
        mg.GetYaxis().CenterTitle(True)
        mg.GetXaxis().CenterTitle(True)

        # 垂直線：x= -0.14, 0, +0.14（用目前 Y 範圍）
        y_min = mg.GetYaxis().GetXmin()
        y_max = mg.GetYaxis().GetXmax()
        vlines = []
        for xv in (-0.14, 0.0, 0.14):
            ln = ROOT.TLine(xv, y_min, xv, y_max)
            if abs(xv) < 1e-9:
                ln.SetLineStyle(7)
                ln.SetLineColor(ROOT.kGray+2)
            else:
                ln.SetLineStyle(1)
                ln.SetLineWidth(3)
                ln.SetLineColor(ROOT.kBlue+2)
            ln.Draw()
            vlines.append(ln)

        # 水平線：y = R_TRUE（用目前 X 範圍）
        x_min = mg.GetXaxis().GetXmin()
        x_max = mg.GetXaxis().GetXmax()
        hline = ROOT.TLine(x_min, R_TRUE, x_max, R_TRUE)
        hline.SetLineStyle(7); hline.SetLineColor(ROOT.kGray+2); hline.Draw()

        leg2.Draw()

        ROOT.gPad.SetTicks(1, 1)
        ROOT.gPad.RedrawAxis()

        # 新增：在第二張圖加上 CMS/Lumi
        draw_cms_lumi(c2, args.lumi)

        c2.SaveAs(f"plots_bias/{outtag}_mean_vs_r.png")
        c2.SaveAs(f"plots_bias/{outtag}_mean_vs_r.pdf")

    print(f"[OK] saved: {outtag}.png / {outtag}.pdf")

if __name__ == "__main__":
    main()
