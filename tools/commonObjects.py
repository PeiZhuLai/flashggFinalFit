import os

# Paths and directory
cmsswbase__ = os.environ['CMSSW_BASE']
cwd__ = os.environ['CMSSW_BASE']+"/src/flashggFinalFit"
swd__ = "%s/Signal"%cwd__
bwd__ = "%s/Background"%cwd__
dwd__ = "%s/Datacard"%cwd__
fwd__ = "%s/Combine"%cwd__
pwd__ = "%s/Plots"%cwd__
twd__ = "%s/Trees2WS"%cwd__

# Centre of mass energy string
sqrts__ = "13p6TeV"

# Luminosity map in fb^-1: for using UL 2018
# lumiMap = {
#     '2016':36.33, 
#     '2017':41.48, 
#     '2018':59.83, 
#     'combined':137.65, 
#     'merged':137.65,
#     '2022preEE':8.00,
#     '2022postEE':26.70
# }

lumiMap = {
    '16':16.81,
    '16APV':19.52, 
    '17':41.48, 
    '18':59.83, 
    'combined':137.65, 
    'merged':137.65,
    '2022preEE':62.50,
    '2022postEE':26.70,
    '2023preBPix':17.79,
    '2023postBPix':9.45
    # Set 2022preEE to 62.5 for lack of signal MC samples 
    # '2022preEE':7.98,
}

# If using ReReco samples then switch to lumiMap below (missing data in 2018 EGamma data set)
#lumiMap = {'2016':36.33, '2017':41.48, '2018':59.35, 'combined':137.17, 'merged':137.17}
lumiScaleFactor = 1000. # Converting from pb to fb

# Constants
BR_W_lnu = 3.*10.86*0.01
BR_Z_ll = 3*3.3658*0.01
BR_Z_nunu = 20.00*0.01
BR_Z_qq = 69.91*0.01
BR_W_qq = 67.41*0.01

# Production modes and decay channel: for extract XS from combine
# productionModes = ['ggH','qqH','ttH','tHq','tHW','ggZH','WH','ZH','bbH']
productionModes = ['ggH']
decayMode = 'hgg'

# List of years
# years_to_process = ['2016','2017','2018','2022preEE','2022postEE', '2023preBPix', '2023postBPix']
years_to_process = ['2022preEE']

# flashgg input WS objects
# inputWSName__ = "tagsDumper/cms_hza_13TeV"
# inputWSName__ = "tagsDumper/CMS_hza_workspace"
inputWSName__ = "CMS_hza_workspace/CMS_hza_workspace"
inputNuisanceExtMap = {'scales':'','scalesCorr':'','smears':''}
# Signal output WS objects
outputWSName__ = "wsig"
outputWSObjectTitle__ = "hggpdfsmrel"
outputWSNuisanceTitle__ = "CMS_hgg_nuisance"
outputNuisanceExtMap = {'scales':'%sscale'%sqrts__,'scalesCorr':'%sscaleCorr'%sqrts__,'smears':'%ssmear'%sqrts__,'scalesGlobal':'%sscale'%sqrts__}
# outputNuisanceExtMap = {'scales':'','scalesCorr':'','smears':'','scalesGlobal':''}
# Bkg output WS objects
bkgWSName__ = "multipdf"
