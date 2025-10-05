import sys, ROOT

def inspect(path, wname="CMS_hza_workspace"):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"[FAIL] 無法開啟: {path}")
        return
    obj = f.Get(wname)
    if not obj:
        print(f"[MISS] {path}: 不存在鍵 '{wname}'")
        return
    cls = obj.ClassName()
    if obj.InheritsFrom("RooWorkspace"):
        print(f"[OK] {path}: 直接取得 RooWorkspace ({cls})")
        return
    if obj.InheritsFrom("TDirectory"):
        f.cd(wname)
        inner = ROOT.gDirectory.Get(wname)
        if inner and inner.InheritsFrom("RooWorkspace"):
            print(f"[NESTED] {path}: 外層 {cls}, 內層同名 RooWorkspace → 建議重寫檔案或修讀取端")
        else:
            print(f"[DIR] {path}: 是目錄 {cls}，內部未找到同名 RooWorkspace")
    else:
        print(f"[UNEXPECTED] {path}: 類型 {cls} 非 RooWorkspace/TDirectory")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python check_workspace_structure.py file1.root [file2.root ...]")
        sys.exit(0)
    for p in sys.argv[1:]:
        inspect(p)
