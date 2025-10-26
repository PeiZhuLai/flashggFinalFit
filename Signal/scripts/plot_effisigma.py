channels = ["ele", "mu"]
years = ["2022preEE"]
mAs = [5, 15, 30]

import os
import json
from array import array
import numpy as np

try:
    import ROOT
except ImportError as e:
    raise SystemExit("PyROOT (ROOT) 未安裝或無法匯入") from e

# 與 plotSignalModel 一致的偏移量，用於 Y 軸 TitleOffset
OFFSET = 0.01
# 樣式設定：可自行調整
MARKER_SIZE = 1.3
LINE_WIDTH = 3

def _extract_effsigma(data, year=None):
    # 1) 先用年份當 key（例如 "2022preEE"）
    if isinstance(data, dict) and year and year in data and isinstance(data[year], (int, float)):
        return float(data[year])

    # 2) 後備：嘗試多種常見 key；也支援 value/val 包在子 dict
    candidate_keys = [
        "effSigma", "effective_sigma", "effsigma",
        "sigma_eff", "effSigma_GeV", "eff_sigma"
    ]
    if isinstance(data, dict):
        for k in candidate_keys:
            if k in data:
                v = data[k]
                if isinstance(v, (int, float)):
                    return float(v)
                if isinstance(v, dict):
                    for subk in ("value", "val", "mean", "mu"):
                        if subk in v and isinstance(v[subk], (int, float)):
                            return float(v[subk])
        # 若只有一個數值欄位也接受
        numeric_vals = [float(v) for v in data.values() if isinstance(v, (int, float))]
        if len(numeric_vals) == 1:
            return numeric_vals[0]
    if isinstance(data, list) and len(data) == 1 and isinstance(data[0], (int, float)):
        return float(data[0])
    return None

def read_effsigma(json_path, year=None):
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
        val = _extract_effsigma(data, year=year)
        if val is None:
            print(f"[WARN] JSON 內找不到有效的 effective sigma（嘗試以年份 '{year}' 作為 key）：{json_path}")
        return val
    except FileNotFoundError:
        print(f"[WARN] 找不到檔案: {json_path}")
        return None
    except Exception as e:
        print(f"[WARN] 讀取/解析失敗: {json_path} -> {e}")
        return None

def ensure_dir(d):
    if not os.path.isdir(d):
        os.makedirs(d, exist_ok=True)

def _quadratic_curve(xs, ys, n=300):
    # 優先使用 SciPy；若不可用則採用局部二次擬合作為後備
    try:
        from scipy.interpolate import make_interp_spline
        x = np.asarray(xs, dtype=float)
        y = np.asarray(ys, dtype=float)
        order = np.argsort(x)
        x = x[order]; y = y[order]
        if len(x) >= 3:
            spline = make_interp_spline(x, y, k=2)
            x_new = np.linspace(x.min(), x.max(), n)
            y_new = spline(x_new)
            return x_new, y_new
        # 少於3點時線性插值
        x_new = np.linspace(x.min(), x.max(), n)
        y_new = np.interp(x_new, x, y)
        return x_new, y_new
    except Exception:
        x = np.asarray(xs, dtype=float)
        y = np.asarray(ys, dtype=float)
        order = np.argsort(x)
        x = x[order]; y = y[order]
        if len(x) < 3:
            x_new = np.linspace(x.min(), x.max(), n)
            y_new = np.interp(x_new, x, y)
            return x_new, y_new
        x_new_list, y_new_list = [], []
        seg_pts = max(8, n // max(1, len(x) - 2))
        for i in range(len(x) - 2):
            xi = x[i:i+3]; yi = y[i:i+3]
            coeffs = np.polyfit(xi, yi, deg=2)
            left = x[i] if i == 0 else 0.5 * (x[i] + x[i+1])
            right = x[i+2] if i == len(x) - 3 else 0.5 * (x[i+1] + x[i+2])
            xs_seg = np.linspace(left, right, seg_pts)
            ys_seg = np.polyval(coeffs, xs_seg)
            if x_new_list and xs_seg[0] <= x_new_list[-1][-1]:
                mask = xs_seg > x_new_list[-1][-1]
                xs_seg = xs_seg[mask]; ys_seg = ys_seg[mask]
            x_new_list.append(xs_seg); y_new_list.append(ys_seg)
        x_new = np.concatenate(x_new_list) if x_new_list else x
        y_new = np.concatenate(y_new_list) if y_new_list else y
        return x_new, y_new

def draw_channel(channel, years, mAs, base_dir):
    ROOT.gStyle.SetOptStat(0)

    # 蒐集各 year 的圖
    graphs = []
    global_ymin, global_ymax = None, None

    colors = [ROOT.TColor.GetColor("#276FBF"), ROOT.TColor.GetColor("#183059"), ROOT.TColor.GetColor("#FC7A1E"), ROOT.TColor.GetColor("#33673B"), ROOT.TColor.GetColor("#34E4EA"), ROOT.TColor.GetColor("#F564A9")]
    markers = [20, 21, 22, 23, 24, 25, 26]

    for yi, year in enumerate(years):
        xs, ys = [], []
        for m in mAs:
            json_path = os.path.join(
                base_dir, "Signal", f"outdir_{channel}", "signalFit", "Plots",
                f"effSigma_{m}_{year}_{channel}.json"
            )
            val = read_effsigma(json_path, year=year)
            if val is None:
                continue
            xs.append(float(m))
            ys.append(float(val))

        if not xs:
            print(f"[INFO] 無可用點 (channel={channel}, year={year})，略過")
            continue

        # 依 mA 排序
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        xs = [xs[i] for i in order]
        ys = [ys[i] for i in order]

        # 平滑二次曲線
        sx, sy = _quadratic_curve(xs, ys, n=300)

        from array import array as carray
        # 平滑連線 graph
        gxl = carray('d', list(map(float, sx)))
        gyl = carray('d', list(map(float, sy)))
        gr_line = ROOT.TGraph(len(sx), gxl, gyl)
        gr_line.SetLineColor(colors[yi % len(colors)])
        gr_line.SetMarkerColor(colors[yi % len(colors)])
        gr_line.SetMarkerStyle(markers[yi % len(markers)])  # 供 legend 顯示
        gr_line.SetLineWidth(LINE_WIDTH)

        # 原始點 graph（只畫點）
        gx = carray('d', list(map(float, xs)))
        gy = carray('d', list(map(float, ys)))
        gr_pts = ROOT.TGraph(len(xs), gx, gy)
        gr_pts.SetLineColor(0)
        gr_pts.SetMarkerColor(colors[yi % len(colors)])
        gr_pts.SetMarkerStyle(markers[yi % len(markers)])
        gr_pts.SetMarkerSize(MARKER_SIZE)

        # 更新全域 y 範圍
        ymin, ymax = min(ys), max(ys)
        global_ymin = ymin if global_ymin is None else min(global_ymin, ymin)
        global_ymax = ymax if global_ymax is None else max(global_ymax, ymax)

        graphs.append((gr_line, gr_pts, year))

    if not graphs:
        print(f"[INFO] channel={channel} 無資料可畫")
        return

    # 與 plotSignalModel 一致的畫布與邊界
    c = ROOT.TCanvas(f"c_{channel}", "", 800, 600)
    c.SetMargin(0.12 + OFFSET, 0.035, 0.14, 0.09)
    c.SetTickx()
    c.SetTicky()

    # 使用第一個 year 的平滑連線做軸
    first_line, first_pts, first_year = graphs[0]
    first_line.SetTitle("")
    first_line.GetXaxis().SetTitle("m_{a} [GeV]")
    first_line.GetYaxis().SetTitle("Signal Resolution [GeV]")
    first_line.Draw("AL")

    # 套用與 plotSignalModel 一致的座標軸樣式
    ax, ay = first_line.GetXaxis(), first_line.GetYaxis()
    for axis in (ax, ay):
        axis.CenterTitle(True)
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(0.055)
        axis.SetLabelSize(0.05)
    ax.SetTitleOffset(1.15)
    ay.SetTitleOffset(1.1 + 10*OFFSET)
    ax.SetLabelOffset(0.009)

    # 固定座標範圍
    first_line.GetXaxis().SetLimits(0.0, 33.0)  # x-axis [0, 33]
    first_line.SetMinimum(0.1)                  # y-axis min = 0
    first_line.SetMaximum(4.2)                  # y-axis max = 5.2

    # 畫第一個年的原始點
    first_pts.Draw("P SAME")

    leg = ROOT.TLegend(0.65, 0.70, 0.93, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    # 新增：依 channel 設定 legend 標題
    channel_label = "Electron" if channel == "ele" else "Muon"
    leg.SetHeader(f"{channel_label} Channel", "L")
    leg.AddEntry(first_line, first_year, "lp")

    # 疊上其餘年：先線後點
    for gr_line, gr_pts, year in graphs[1:]:
        gr_line.Draw("L SAME")
        gr_pts.Draw("P SAME")
        leg.AddEntry(gr_line, year, "lp")

    leg.Draw()

    # 加上 CMS Preliminary 與 13.6 TeV（風格與 plottingTools.plotSignalModel 一致）
    lat0 = ROOT.TLatex()
    lat0.SetTextFont(42)
    lat0.SetTextAlign(11)
    lat0.SetNDC()
    lat0.SetTextSize(0.05)
    lat0.DrawLatex(0.12 + OFFSET, 0.92, "#bf{CMS} #it{Preliminary}")
    lat0.DrawLatex(0.82 + OFFSET, 0.92, "13.6 TeV")

    c.Update()

    outdir = os.path.join(base_dir, "Signal", f"outdir_{channel}", "signalFit", "Plots")
    ensure_dir(outdir)
    out_base = os.path.join(outdir, f"effSigmaVmA_{channel}")
    c.SaveAs(out_base + ".png")
    c.SaveAs(out_base + ".pdf")
    print(f"[OK] 儲存圖檔：{out_base}.png / .pdf")

def main():
    base_dir = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"

    for ch in channels:
        draw_channel(ch, years, mAs, base_dir)

if __name__ == "__main__":
    main()