# Input config file for running trees2ws

trees2wsCfg = {

  # Name of RooDirectory storing input tree
  'inputTreeDir':'DiphotonTree',

  # Variables to be added to dataframe: use wildcard * for common strings
  'mainVars':["CMS_hza_mass","weight","dZ","weight_*"], # Vars to add to nominal RooDatasets
  'dataVars':["CMS_hza_mass","weight"], # Vars for data workspace (trees2ws_data.py script)
  'stxsVar':'stage1p2bin', # Var for STXS splitting: if using option doSTXSSplitting
  'systematicsVars':["CMS_hza_mass","weight"], # Variables to add to sytematic RooDataHists
  'theoryWeightContainers':{}, # Theory weights to add to nominal + NOTAG RooDatasets, value corresponds to number of weights (0-N)

  # List of shape systematics: use string YEAR for year-dependent systematics
  'systematics':['FNUF', 'Material', 
  'ElectronScale', 'ElectronSmear', 
  'MuonPtScale', 'MuonPtSmear', 
  'PhotonScale', 'PhotonSmear'],

  # Analysis categories: python list of cats or use 'auto' to extract from input tree
  'cats':'auto'

}
