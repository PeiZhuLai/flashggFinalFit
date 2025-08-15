import pandas as pd

# 直接读取
df = pd.read_pickle("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/yields_leptons/5_datacard_2022preEE_leptons.pkl")

# 看前几行
print(df.head())

# 看列名
print(df.columns)

# 如果只是想看某几列
print(df[['proc', 'cat', 'nominal_yield']])

pd.set_option('display.max_rows', None)   # 显示所有行
pd.set_option('display.max_columns', None) # 显示所有列
pd.set_option('display.width', None)      # 不自动换行
pd.set_option('display.max_colwidth', None)  # 列内容不截断

print(df[['year', 'type', 'procOriginal', 'proc', 'proc_s0', 'cat', 'inputWSFile',
       'nominalDataName', 'modelWSFile', 'model', 'rate', 'nominal_yield',
       'sumw2', 'nominal_yield_COWCorr']])
