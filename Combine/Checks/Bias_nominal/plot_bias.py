# plot_pull_overlay_paths.py (safe-ownership version)
import os, glob, re, argparse
import json
import ROOT
import math
ROOT.gROOT.SetBatch(True)

try:
    from cmsLumi import CMS_lumi
    _has_cms_lumi = True
except Exception:
    CMS_lumi = None
    _has_cms_lumi = False

# ---------- 小工具：统一托管 ROOT 对象，避免双重析构 ----------
def keep(obj, canvas=None, bucket_name="_keep"):
    try:
        ROOT.SetOwnership(obj, False)  # 交给 ROOT，避免 Python 也析构
    except Exception:
        pass
    if canvas is not None:
        if not hasattr(canvas, bucket_name):
            setattr(canvas, bucket_name, [])
        getattr(canvas, bucket_name).append(obj)
    return obj

# === 1) 在這裡直接寫你的絕對路徑 ===
BIAS_BASE = os.environ.get(
    "BIAS_BASE_OVERRIDE",
    "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal",
)
PATH_BIAS_FITS = os.path.join(BIAS_BASE, "BiasFits")
PATH_BIAS_TOYS = os.path.join(BIAS_BASE, "BiasToys")

R_TRUE = 1.0
TITLE  = "Bias study pull comparison"
OUTTAG = "pull_overlay_distribution"
OUTTAG_MEAN = "r_vs_bias"

MEAN_VS_R_XMIN = -0.5
MEAN_VS_R_XMAX =  0.5
MEAN_VS_R_YMIN = None
MEAN_VS_R_YMAX = None

PULL_HIST_NBINS = 80
PULL_HIST_XMIN = -4.0
PULL_HIST_XMAX =  4.0

def file_label(path):
    base = os.path.basename(path)
    m = re.search(r"biasStudy_(.+?)_fits\.root", base)
    return m.group(1) if m else base

def parse_truth_from_filename(path, default=R_TRUE):
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

def parse_mu_from_token(token, default=R_TRUE):
    s = str(token).replace('p','.')
    s = re.sub(r"^[a-zA-Z_]+", "", s)
    try:
        if re.fullmatch(r"\d+", s):
            return float(s)/10.0
        return float(s)
    except Exception:
        return default

def load_gaussfit_json(json_path):
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
        exp_val = data.get("exp", None)
        means = {}
        mean_errs = {}

        def _extract(container):
            for k, v in container.items():
                if isinstance(v, dict):
                    if "mean" in v:
                        try:
                            means[k] = float(v["mean"])
                        except Exception:
                            means[k] = v["mean"]
                    # 容錯支援多種 key 命名
                    err_val = v.get("mean_err", v.get("meanErr", v.get("meanerror", None)))
                    if err_val is not None:
                        try:
                            mean_errs[k] = float(err_val)
                        except Exception:
                            mean_errs[k] = err_val

        if isinstance(data.get("fit_results"), dict):
            _extract(data["fit_results"])
        else:
            # 兼容平鋪結構
            _extract({k: v for k, v in data.items() if isinstance(v, dict)})

        return means, mean_errs, exp_val
    except Exception as e:
        print(f"[WARN] cannot load JSON: {json_path} ({e})")
        return {}, {}, None

def _list_leaves(tree):
    try:
        if not tree:
            return []
        return [l.GetName() for l in tree.GetListOfLeaves()] or []
    except Exception:
        return []

def _quantile_matches(value, target, tol=1e-3):
    try:
        return abs(float(value) - float(target)) < tol
    except Exception:
        return False

def get_pulls(fname, r_true=1.0):
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        print(f"[WARN] cannot open {fname}")
        return [], []
    t = f.Get("limit")
    if not t:
        print(f"[WARN] no TTree 'limit' in {fname}")
        f.Close()
        return [], []

    leaves_list = _list_leaves(t)
    if "r" not in leaves_list or "quantileExpected" not in leaves_list:
        print(f"[WARN] missing required branches in {fname}. Leaves={leaves_list}")
        f.Close()
        return [], []

    pulls = []
    used_r_vals = []
    n_total = int(t.GetEntriesFast())
    n_triplets = n_total // 3
    n_bad_quantile = 0
    n_bad_unc = 0

    for itoy in range(n_triplets):
        i_bf = 3 * itoy
        i_lo = i_bf + 1
        i_hi = i_bf + 2

        t.GetEntry(i_bf)
        q_bf = getattr(t, "quantileExpected", float("nan"))
        bf = getattr(t, "r", float("nan"))

        t.GetEntry(i_lo)
        q_lo = getattr(t, "quantileExpected", float("nan"))
        lo = getattr(t, "r", float("nan"))

        t.GetEntry(i_hi)
        q_hi = getattr(t, "quantileExpected", float("nan"))
        hi = getattr(t, "r", float("nan"))

        if (not _quantile_matches(q_bf, -1.0) or
            not _quantile_matches(q_lo, -0.32) or
            not _quantile_matches(q_hi, 0.32)):
            n_bad_quantile += 1
            continue

        unc = 0.5 * (hi - lo)
        if unc <= 0.0 or not math.isfinite(unc) or not math.isfinite(bf):
            n_bad_unc += 1
            continue

        pulls.append((bf - r_true) / unc)
        used_r_vals.append(bf)

    f.Close()
    print(
        f"[INFO] {os.path.basename(fname)}: "
        f"entries={n_total}, triplets={n_triplets}, pulls_kept={len(pulls)}, "
        f"bad_quantile={n_bad_quantile}, bad_unc={n_bad_unc}"
    )
    if not pulls:
        print(f"[HINT] Available leaves: {leaves_list}")
    return pulls, used_r_vals

def make_hist(values, name, bins=48, xmin=-5, xmax=5):
    h = ROOT.TH1F(name, name, PULL_HIST_NBINS, PULL_HIST_XMIN, PULL_HIST_XMAX)
    h.Sumw2(False)
    for v in values:
        h.Fill(v)
    h.GetXaxis().SetTitle("Pull")
    h.GetYaxis().SetTitle("Entries")
    h.GetXaxis().SetTitleSize(0.055); h.GetXaxis().SetLabelSize(0.05); h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitleSize(0.05);  h.GetYaxis().SetLabelSize(0.05);  h.GetYaxis().SetTitleOffset(1.3)
    h.GetYaxis().CenterTitle(True);   h.GetXaxis().CenterTitle(True)
    return h

def draw_cms_lumi(canvas, lumi_fb):
    if _has_cms_lumi and CMS_lumi:
        cms = keep(CMS_lumi(), canvas)
        cms.set_lumi(canvas, float(lumi_fb), 0, text="Preliminary", drawLumi=True)
    else:
        latex = keep(ROOT.TLatex(), canvas)
        latex.SetNDC(); latex.SetTextSize(0.045)
        latex.DrawLatex(0.16, 0.92, f"CMS (unofficial)  {float(lumi_fb)} fb^{{-1}} (13 TeV)")

def parse_args():
    p = argparse.ArgumentParser(description="Plot pull overlay")
    p.add_argument("--mA", type=str, default=None, help="Prefix for output plots")
    p.add_argument("--lumi", type=float, default=172.13, help="Luminosity in fb^-1 for label")
    return p.parse_args()

def main():
    args = parse_args()

    json_means, json_errs, json_exp = {}, {}, None
    if args.mA:
        json_path = os.path.join(BIAS_BASE, "BiasJson", f"{args.mA}_gaussfit.json")
        json_means, json_errs, json_exp = load_gaussfit_json(json_path)
        if json_means: print(f"[INFO] loaded JSON means from {json_path}")
        if json_exp is not None: print(f"[INFO] JSON exp = {json_exp}")

    if not os.path.isdir(PATH_BIAS_FITS):
        raise SystemExit(f"[ERR] Not found: {PATH_BIAS_FITS}")

    files = sorted(glob.glob(os.path.join(PATH_BIAS_FITS, "biasStudy_*_fits.root")))
    files = [f for f in files if "split" not in os.path.basename(f)]
    if not files:
        raise SystemExit(f"[ERR] No files in {PATH_BIAS_FITS}/biasStudy_*_fits.root")

    c = keep(ROOT.TCanvas("c","c",800,600))
    c.SetTopMargin(0.09); c.SetBottomMargin(0.14)
    c.SetRightMargin(0.05); c.SetLeftMargin(0.13)

    leg = keep(ROOT.TLegend(0.52, 0.62, 0.80, 0.87), c)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.043)
    ROOT.gStyle.SetOptStat(0)

    colors = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kMagenta+1,
              ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kViolet, ROOT.kTeal+2]
    styles = [1,2,3,4,5,6,7,8]

    hists = []
    scatter_points = []
    fits = []

    for i, fpath in enumerate(files):
        pulls, _ = get_pulls(fpath, R_TRUE)
        if not pulls:
            print(f"[INFO] skip empty: {fpath}")
            continue
        h = keep(make_hist(pulls, f"h_{i}"), c)
        color = colors[i % len(colors)]
        h.SetLineColor(color); h.SetLineWidth(3); h.SetLineStyle(styles[i % len(styles)])
        label = file_label(fpath)
        hists.append((h, label))

    if not hists:
        raise SystemExit("[ERR] nothing to plot.")

    ymax = max(h.GetMaximum() for h,_ in hists) * 1.35
    for idx, (h, label) in enumerate(hists):
        h.SetMaximum(ymax)
        h.SetTitle("" if idx==0 else h.GetTitle())
        h.Draw("HIST" if idx==0 else "HIST SAME")

        fit = keep(ROOT.TF1(f"g_{idx}", "gaus", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()), c)
        fit.SetLineColor(h.GetLineColor()); fit.SetLineStyle(1)
        h.Fit(fit, "RQ")
        fits.append(fit)
        mean = fit.GetParameter(1); sigma = fit.GetParameter(2)
        leg.AddEntry(h, f"{label}  (#mu={mean:+.2f}, #sigma={sigma:.2f})", "l")

        truth_y = parse_truth_from_filename(files[idx], R_TRUE)
        scatter_points.append((mean, truth_y, label, h.GetLineColor()))

    leg.Draw()
    z = keep(ROOT.TLine(0, 0, 0, ymax), c); z.SetLineStyle(7); z.SetLineColor(ROOT.kGray+2); z.Draw()

    outtag = f"{args.mA}_{OUTTAG}" if args.mA else OUTTAG
    os.makedirs("plots_bias", exist_ok=True)

    ROOT.gPad.SetTicks(1, 1); ROOT.gPad.RedrawAxis()
    draw_cms_lumi(c, args.lumi)
    c.SaveAs(f"plots_bias/{outtag}.png"); c.SaveAs(f"plots_bias/{outtag}.pdf")

    # ===== 第二張圖：pull mean vs r truth =====
    points = []
    if json_means:
        y_val = json_exp if json_exp is not None else parse_mu_from_token(args.mA, R_TRUE)
        for j, (lab, mean) in enumerate(json_means.items()):
            col = colors[j % len(colors)]
            xerr = float(json_errs.get(lab, 0.0)) if lab in json_errs else 0.0
            points.append((mean, xerr, y_val, lab, col))
    else:
        # 從前一張圖的擬合均值，無 JSON 時誤差設為 0
        points = [(xmu, 0.0, yr, lab, col) for (xmu, yr, lab, col) in scatter_points]

    if points:
        c2 = keep(ROOT.TCanvas("c2","c2",800,600))
        c2.SetTopMargin(0.09); c2.SetBottomMargin(0.14)
        c2.SetRightMargin(0.05); c2.SetLeftMargin(0.13)

        mg = keep(ROOT.TMultiGraph(), c2)
        leg2 = keep(ROOT.TLegend(0.73, 0.58, 0.98, 0.87), c2)
        leg2.SetBorderSize(0); leg2.SetFillStyle(0); leg2.SetTextSize(0.045)

        graphs = []
        for j, (xmu, xerr, yr, lab, col) in enumerate(points):
            g = keep(ROOT.TGraphErrors(1), c2)
            g.SetPoint(0, xmu, yr)
            g.SetPointError(0, xerr, 0.0)  # x 向誤差棒；y 無誤差
            g.SetMarkerColor(col)
            g.SetMarkerStyle(20 + ((j + 1) % 10))
            g.SetMarkerSize(2.0)
            g.SetLineColor(col); g.SetLineWidth(3)
            mg.Add(g, "P")
            graphs.append(g)
            leg2.AddEntry(g, lab, "p")

        mg.SetTitle("")
        mg.GetXaxis().SetTitle("(#mu - #tilde{#mu}) / #sigma_{#mu}")
        mg.GetYaxis().SetTitle("#tilde{#mu}")
        mg.GetXaxis().SetLimits(MEAN_VS_R_XMIN, MEAN_VS_R_XMAX)

        if json_means:
            y_center = json_exp if json_exp is not None else parse_mu_from_token(args.mA, R_TRUE)
        else:
            y_vals = [yr for _, _, yr, _, _ in points]
            y_center = (sum(y_vals) / len(y_vals)) if y_vals else R_TRUE

        y_pad = 1.2 * (abs(y_center) if abs(y_center) > 1e-9 else 1.0)
        mg.SetMinimum(y_center - y_pad); mg.SetMaximum(y_center + y_pad)

        mg.Draw("A")
        ROOT.gPad.Modified(); ROOT.gPad.Update()

        mg.GetXaxis().SetTitleSize(0.055)
        mg.GetXaxis().SetLabelSize(0.05)
        mg.GetXaxis().SetTitleOffset(1.1)
        mg.GetYaxis().SetTitleSize(0.05)
        mg.GetYaxis().SetLabelSize(0.05)
        mg.GetYaxis().SetTitleOffset(1.3)
        mg.GetYaxis().CenterTitle(True)
        mg.GetXaxis().CenterTitle(True)

        x_min = ROOT.gPad.GetUxmin(); x_max = ROOT.gPad.GetUxmax()
        y_min = ROOT.gPad.GetUymin(); y_max = ROOT.gPad.GetUymax()

        vlines = []
        for xv in (-0.2, 0.0, 0.2):
            ln = keep(ROOT.TLine(xv, y_min, xv, y_max), c2)
            if abs(xv) < 1e-9:
                ln.SetLineStyle(7); ln.SetLineColor(ROOT.kGray+2)
            else:
                ln.SetLineStyle(1); ln.SetLineWidth(3); ln.SetLineColor(ROOT.kBlue+2)
            ln.Draw()
            vlines.append(ln)

        y_line = json_exp if json_exp is not None else y_center
        hline = keep(ROOT.TLine(x_min, y_line, x_max, y_line), c2)
        hline.SetLineStyle(7); hline.SetLineColor(ROOT.kGray+2); hline.Draw()

        for g in graphs:
            g.Draw("P SAME")

        leg2.Draw()
        ROOT.gPad.SetTicks(1, 1); ROOT.gPad.RedrawAxis()
        draw_cms_lumi(c2, args.lumi)
        ROOT.gPad.Modified(); ROOT.gPad.Update()

        lat = ROOT.TLatex()
        lat.SetTextFont(42)
        lat.SetTextAlign(13)
        lat.SetNDC()
        lat.SetTextSize(0.050)
        lat.DrawLatex(0.165, 0.86, f"m_{{a}} = {args.mA} GeV")

        outtag = f"{args.mA}_{OUTTAG_MEAN}" if args.mA else OUTTAG_MEAN
        c2.SaveAs(f"plots_bias/{outtag}.png")
        c2.SaveAs(f"plots_bias/{outtag}.pdf")


if __name__ == "__main__":
    main()
