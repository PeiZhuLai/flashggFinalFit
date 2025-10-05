import ROOT, os, sys, tempfile, shutil
import argparse

def flatten_one(path, wname="CMS_hza_workspace", inplace=False, output=None, force=False):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        print(f"[FAIL] 無法開啟: {path}")
        return False
    obj = f.Get(wname)
    if not obj:
        print(f"[SKIP] {path}: 不存在鍵 {wname}")
        f.Close()
        return False

    # 情況 1: 已經是扁平 RooWorkspace
    if obj.InheritsFrom("RooWorkspace"):
        print(f"[OK] {path}: 已是扁平結構，略過")
        f.Close()
        return True

    # 情況 2: 是目錄 → 嘗試取內層同名 RooWorkspace
    if obj.InheritsFrom("TDirectory"):
        f.cd(wname)
        inner = ROOT.gDirectory.Get(wname)
        if not inner or not inner.InheritsFrom("RooWorkspace"):
            print(f"[WARN] {path}: 目錄內找不到同名 RooWorkspace，略過")
            f.Close()
            return False
        ws = inner
    else:
        print(f"[WARN] {path}: 非 RooWorkspace / TDirectory 類型 {obj.ClassName()}，略過")
        f.Close()
        return False

    # 決定輸出檔名
    if output and inplace:
        print(f"[ERR] 不可同時指定 --output 與 --inplace")
        f.Close()
        return False
    if inplace:
        out_path = path
    else:
        if output:
            out_path = output
        else:
            base, ext = os.path.splitext(path)
            out_path = base + "_flat" + ext
        if os.path.exists(out_path) and not force:
            print(f"[SKIP] {out_path} 已存在 (加 --force 覆寫)")
            f.Close()
            return False

    # 寫出：直接含 RooWorkspace
    tmp_dir = None
    if inplace:
        # 建立暫存檔避免壞檔
        tmp_fd, tmp_name = tempfile.mkstemp(prefix="ws_flat_", suffix=".root",
                                            dir=os.path.dirname(path) or ".")
        os.close(tmp_fd)
        out_real = tmp_name
    else:
        out_real = out_path

    fout = ROOT.TFile(out_real, "RECREATE")
    ws.Write(wname)
    fout.Close()
    f.Close()

    if inplace:
        shutil.move(out_real, path)
        print(f"[FIXED][INPLACE] {path}")
    else:
        print(f"[CREATED] {out_path}")

    return True

def main():
    ap = argparse.ArgumentParser(description="扁平化雙層 RooWorkspace ROOT 檔案")
    ap.add_argument("files", nargs="+", help="輸入 ROOT 檔 (可多個)")
    ap.add_argument("-w","--wname", default="CMS_hza_workspace", help="Workspace 名稱")
    ap.add_argument("--inplace", action="store_true", help="原地覆寫 (謹慎使用)")
    ap.add_argument("-o","--output", help="單一輸入時指定輸出檔名")
    ap.add_argument("--force", action="store_true", help="允許覆寫已存在輸出檔")
    args = ap.parse_args()

    if args.output and len(args.files) != 1:
        print("[ERR] --output 只能搭配單一輸入檔")
        sys.exit(1)

    ok = 0
    for fp in args.files:
        if flatten_one(fp, wname=args.wname, inplace=args.inplace,
                       output=args.output, force=args.force):
            ok += 1
    print(f"[SUMMARY] 成功 {ok} / {len(args.files)}")

if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    main()
