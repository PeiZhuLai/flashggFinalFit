import pickle
import pandas as pd

# 讀取 pkl 檔案
with open("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_ele/calcPhotonSyst/pkl/5_2022preEE.pkl", "rb") as f:   # 注意要用 "rb" 二进制模式
    data = pickle.load(f)

# 調整 pandas 顯示選項
pd.set_option("display.max_rows", None)       # 顯示所有列
pd.set_option("display.max_columns", None)    # 顯示所有欄
pd.set_option("display.width", None)          # 自動換行
pd.set_option("display.max_colwidth", None)   # 不截斷文字

# 確保為 DataFrame
if not isinstance(data, pd.DataFrame):
    try:
        data = pd.DataFrame(data)
    except Exception:
        raise TypeError("無法將載入的物件轉成 DataFrame，請檢查 pkl 內容。")

def print_formatted(df):
    if len(df) == 1:
        s = df.iloc[0]
        for c in df.columns:
            val = s[c]
            if isinstance(val, (int, float)):
                print(f"{c}: {val:.6f}")
            else:
                print(f"{c}: {val}")
    else:
        for idx, row in df.iterrows():
            print(f"row {idx}:")
            for c in df.columns:
                val = row[c]
                if isinstance(val, (int, float)):
                    print(f"  {c}: {val:.6f}")
                else:
                    print(f"  {c}: {val}")

# 改為列印全部欄位
print_formatted(data)
