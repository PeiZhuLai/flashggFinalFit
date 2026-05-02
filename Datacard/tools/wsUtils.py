import ROOT

def _inherits(obj, class_name):
    if obj is None:
        return False
    try:
        return bool(obj.InheritsFrom(class_name))
    except Exception:
        return False

def _obj_name(obj):
    try:
        return obj.GetPath()
    except Exception:
        try:
            return obj.GetName()
        except Exception:
            return str(obj)

def _scan_dir_for_ws(tdir, ws_name):
    """在給定 TDirectoryFile 內尋找:
       1) 同名 RooWorkspace
       2) 否則第一個 RooWorkspace
    """
    if tdir is None:
        return None
    cand = tdir.Get(ws_name)
    if _inherits(cand, "RooWorkspace"):
        return cand

    nested = tdir.Get(f"{ws_name}/{ws_name}")
    if _inherits(nested, "RooWorkspace"):
        return nested

    for k in tdir.GetListOfKeys():
        try:
            obj = k.ReadObj()
        except Exception:
            obj = tdir.Get(k.GetName())
        if _inherits(obj, "RooWorkspace"):
            return obj
        if _inherits(obj, "TDirectory"):
            ws = _scan_dir_for_ws(obj, ws_name)
            if ws:
                return ws
    return None

def fetchWorkspace(file_or_obj, ws_name):
    """
    傳入:
      - TFile
      - TDirectoryFile
      - RooWorkspace
    取得 RooWorkspace 物件 (含處理 TDirectoryFile 巢狀結構)。
    回傳 (ws, logs)
      ws  : RooWorkspace or None
      logs: list of debug / error 字串
    """
    logs = []
    # 直接是 RooWorkspace
    if _inherits(file_or_obj, "RooWorkspace"):
        return file_or_obj, logs

    # ROOT.TFile 繼承自 TDirectoryFile，所以一定要先處理 TFile。
    if _inherits(file_or_obj, "TFile"):
        raw = file_or_obj.Get(ws_name)
        # 1) 直接取到
        if _inherits(raw, "RooWorkspace"):
            return raw, logs
        # 2) 取到的是目錄
        if _inherits(raw, "TDirectory"):
            ws = _scan_dir_for_ws(raw, ws_name)
            if ws:
                return ws, logs
        # 3) 明確處理常見的 dir/ws 結構
        raw = file_or_obj.Get(f"{ws_name}/{ws_name}")
        if _inherits(raw, "RooWorkspace"):
            return raw, logs
        # 4) 第一層鍵與其子目錄掃描
        top_keys = []
        for k in file_or_obj.GetListOfKeys():
            top_keys.append(k.GetName())
            try:
                obj = k.ReadObj()
            except Exception:
                obj = file_or_obj.Get(k.GetName())
            if _inherits(obj, "RooWorkspace"):
                return obj, logs
            if _inherits(obj, "TDirectory"):
                ws = _scan_dir_for_ws(obj, ws_name)
                if ws:
                    return ws, logs
        logs.append(f"[FETCH] 無法在檔案內找到 RooWorkspace='{ws_name}'；頂層鍵={top_keys}")
        return None, logs

    # 若是一般 TDirectoryFile
    if _inherits(file_or_obj, "TDirectory"):
        ws = _scan_dir_for_ws(file_or_obj, ws_name)
        if ws:
            return ws, logs
        logs.append(f"[FETCH] 於目錄 {_obj_name(file_or_obj)} 未找到 RooWorkspace='{ws_name}'")
        return None, logs

    logs.append(f"[FETCH] 不支援的物件型別: {type(file_or_obj)}")
    return None, logs
