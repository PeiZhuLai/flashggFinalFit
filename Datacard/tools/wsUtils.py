import ROOT

def _scan_dir_for_ws(tdir, ws_name):
    """在給定 TDirectoryFile 內尋找:
       1) 同名 RooWorkspace
       2) 否則第一個 RooWorkspace
    """
    if tdir is None:
        return None
    cand = tdir.Get(ws_name)
    if isinstance(cand, ROOT.RooWorkspace):
        return cand
    for k in tdir.GetListOfKeys():
        obj = tdir.Get(k.GetName())
        if isinstance(obj, ROOT.RooWorkspace):
            return obj
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
    if isinstance(file_or_obj, ROOT.RooWorkspace):
        return file_or_obj, logs

    # 若是 TDirectoryFile
    if isinstance(file_or_obj, ROOT.TDirectoryFile):
        ws = _scan_dir_for_ws(file_or_obj, ws_name)
        if ws:
            return ws, logs
        logs.append(f"[FETCH] 於目錄 {file_or_obj.GetName()} 未找到 RooWorkspace='{ws_name}'")
        return None, logs

    # 若是 TFile
    if isinstance(file_or_obj, ROOT.TFile):
        raw = file_or_obj.Get(ws_name)
        # 1) 直接取到
        if isinstance(raw, ROOT.RooWorkspace):
            return raw, logs
        # 2) 取到的是目錄
        if isinstance(raw, ROOT.TDirectoryFile):
            ws = _scan_dir_for_ws(raw, ws_name)
            if ws:
                return ws, logs
        # 3) 第一層鍵與其子目錄掃描
        top_keys = []
        for k in file_or_obj.GetListOfKeys():
            top_keys.append(k.GetName())
            obj = file_or_obj.Get(k.GetName())
            if isinstance(obj, ROOT.RooWorkspace):
                return obj, logs
            if isinstance(obj, ROOT.TDirectoryFile):
                ws = _scan_dir_for_ws(obj, ws_name)
                if ws:
                    return ws, logs
        logs.append(f"[FETCH] 無法在檔案內找到 RooWorkspace='{ws_name}'；頂層鍵={top_keys}")
        return None, logs

    logs.append(f"[FETCH] 不支援的物件型別: {type(file_or_obj)}")
    return None, logs
