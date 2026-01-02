import ROOT

def _find_rooworkspace_in_dir(tdir, preferred_name=None, logs=None):
  """在 TDirectory 內尋找 RooWorkspace（可指定偏好名稱）。"""
  logs = logs if logs is not None else []
  if tdir is None:
    return None

  # 1) 若 preferred_name 存在，先直取
  if preferred_name:
    try:
      obj = tdir.Get(preferred_name)
      if obj and obj.InheritsFrom("RooWorkspace"):
        logs.append(f"[FETCH] 在目錄 {tdir.GetPath()} 內找到 RooWorkspace='{preferred_name}' (direct)")
        return obj
      if obj:
        logs.append(f"[FETCH] 目錄 {tdir.GetPath()} 內 '{preferred_name}' class={obj.ClassName()} (不是 RooWorkspace)")
    except Exception as e:
      logs.append(f"[FETCH] 讀取目錄內 preferred_name 失敗: {e}")

  # 2) 掃描該目錄 keys，找第一個 RooWorkspace
  try:
    keys = tdir.GetListOfKeys()
    if not keys:
      logs.append(f"[FETCH] 目錄 {tdir.GetPath()} 無任何 keys")
      return None
    for k in keys:
      name = k.GetName()
      try:
        obj = k.ReadObj()
      except Exception:
        obj = None
      if obj and obj.InheritsFrom("RooWorkspace"):
        logs.append(f"[FETCH] 在目錄 {tdir.GetPath()} 內找到 RooWorkspace='{name}' (scan)")
        return obj
  except Exception as e:
    logs.append(f"[FETCH] 掃描目錄 keys 失敗: {e}")

  logs.append(f"[FETCH] 目錄 {tdir.GetPath()} 內未找到任何 RooWorkspace")
  return None


def fetchWorkspace(tfile, wsName):
  """
  Robust 取得 RooWorkspace:
    - 直接以 wsName 取 RooWorkspace
    - 若 wsName 指到的是 TDirectoryFile，進去找 RooWorkspace
    - 若找不到 wsName，嘗試在檔案根目錄掃描 (可選擇保守：只掃第一層)
  回傳 (ws, logs)
  """
  logs = []
  if tfile is None:
    return None, ["[FETCH] tfile is None"]

  # 0) 先嘗試直接 Get(wsName)
  obj = None
  try:
    obj = tfile.Get(wsName)
  except Exception as e:
    logs.append(f"[FETCH] tfile.Get('{wsName}') 失敗: {e}")

  if obj:
    logs.append(f"[FETCH] tfile.Get('{wsName}') -> class={obj.ClassName()}")
    if obj.InheritsFrom("RooWorkspace"):
      return obj, logs

    # 1) 若是目錄：進去找 RooWorkspace
    if obj.InheritsFrom("TDirectory"):
      ws = _find_rooworkspace_in_dir(obj, preferred_name=wsName, logs=logs)
      if ws:
        return ws, logs

  else:
    logs.append(f"[FETCH] tfile.Get('{wsName}') 回傳 None")

  # 2) 最後：在根目錄第一層掃描，找 RooWorkspace
  try:
    keys = tfile.GetListOfKeys()
    if keys:
      for k in keys:
        name = k.GetName()
        try:
          o = k.ReadObj()
        except Exception:
          o = None
        if not o:
          continue
        # 若第一層本來就是 RooWorkspace
        if o.InheritsFrom("RooWorkspace"):
          logs.append(f"[FETCH] 在檔案根目錄找到 RooWorkspace='{name}' (fallback scan)")
          return o, logs
        # 若第一層是目錄，且其名字剛好等於 wsName（你現在的情況）
        if o.InheritsFrom("TDirectory") and name == wsName:
          ws = _find_rooworkspace_in_dir(o, preferred_name=None, logs=logs)
          if ws:
            return ws, logs
  except Exception as e:
    logs.append(f"[FETCH] fallback scan 失敗: {e}")

  logs.append(f"[FETCH] 於目錄 {tfile.GetName()} 未找到 RooWorkspace='{wsName}'")
  return None, logs
