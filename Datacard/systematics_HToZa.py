# Python file to store systematics: for STXS analysis

# Comment out all nuisances that you do not want to include

# THEORY SYSTEMATICS:

# For type:constant
#  1) specify same value for all processes
#  2) define process map json in ./theory_uncertainties (add process names where necessary!)

# For type:factory
# Tier system: adds different uncertainties to dataframe
#   1) shape: absolute yield of process kept constant, shape effects i.e. calc migrations across cats
#   2) ishape: as (1) but absolute yield for proc x cat is allowed to vary
#   3) norm: absolute yield of production mode (s0) kept constant but migrations across sub-processes e.g. STXS bins.Same value in each category.
#   4) inorm: as (3) but absolute yield of production mode (s0) can vary
#   5) inc: variations in production mode (s0), same value for each subprocess in each category
# Relations: shape = ishape/inorm
#            norm  = inorm/inc
# Specify as list in dict: e.g. 'tiers'=['inc','inorm','norm','ishape','shape']

theory_systematics = [
                # 若要改成 factory：示例 (請確認 workspace 內有對應權重或 hist)
                # {'name':'QCDscale_ggH','title':'QCDscale_ggH','type':'factory','prior':'lnN','correlateAcrossYears':1,'tiers':['inc','inorm','norm','ishape','shape']},
                # 'debug_constant': True
                {'name':'QCDscale_ggH','title':'QCDscale_ggH','type':'constant','prior':'lnN','correlateAcrossYears':1,'value':'theory_uncertainties/thu_ggh.json', 'json_key_field':'procOriginal'},
                {'name':'pdf_Higgs_ggH','title':'pdf_Higgs_ggH','type':'constant','prior':'lnN','correlateAcrossYears':1,'value':'theory_uncertainties/thu_ggh.json', 'json_key_field':'procOriginal'},
                {'name':'alphaS_ggH','title':'alphaS_ggH','type':'constant','prior':'lnN','correlateAcrossYears':1,'value':'theory_uncertainties/thu_ggh.json', 'json_key_field':'procOriginal'},
              ]
# PDF weight
# for i in range(1,60): theory_systematics.append( {'name':'pdfWeight_%g'%i, 'title':'CMS_hgg_pdfWeight_%g'%i, 'type':'factory','prior':'lnN','correlateAcrossYears':1,'tiers':['shape']} )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# EXPERIMENTAL SYSTEMATICS
# correlateAcrossYears = 0 : no correlation
# correlateAcrossYears = 1 : fully correlated
# correlateAcrossYears = -1 : partially correlated

experimental_systematics = [
                # Updated luminosity partial-correlation scheme: 13/5/21 (recommended simplified nuisances)
                # Ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiRecommendationsRun3#2024
                {'name':'lumi_13p6TeV_Uncorrelated','title':'lumi_13p6TeV_Uncorrelated','type':'constant','prior':'lnN','correlateAcrossYears':0,'value':{'2022preEE':'1.014','2022postEE':'1.014','2023preBPix':'1.013','2023postBPix':'1.013','2024':'1.016'}},
                {'name':'lumi_13p6TeV_Correlated_2223','title':'lumi_13p6TeV_Correlated_2223','type':'constant','prior':'lnN','correlateAcrossYears':1,'value':{'2022preEE':'1.0101','2022postEE':'1.0101','2023preBPix':'1.0101','2023postBPix':'1.0101'}},
                {'name':'lumi_13p6TeV_Correlated_2324','title':'lumi_13p6TeV_Correlated_2324','type':'constant','prior':'lnN','correlateAcrossYears':1,'value':{'2023preBPix':'1.0141','2023postBPix':'1.0141','2024':'1.0141'}},
                {'name':'lumi_13p6TeV_Correlated_222324','title':'lumi_13p6TeV_Correlated_222324','type':'constant','prior':'lnN','correlateAcrossYears':1,'value':{'2022preEE':'1.0120','2022postEE':'1.0120','2023preBPix':'1.0120','2023postBPix':'1.0120','2024':'1.0120'}},

                # Trigger
                # 命名規則: workspace 變數為 weight_<name>_{central,Up,Down} → 此處只需填 <name>
                # 對應: weight_hlt_sf_central / weight_hlt_sf_Up / weight_hlt_sf_Down
                {'name':'hlt_sf','title':'CMS_hza_trigger','type':'factory','prior':'lnN','correlateAcrossYears':1},

                # Pileup
                {'name':'pu_reweight_sf','title':'CMS_hza_pileup','type':'factory','prior':'lnN','correlateAcrossYears':1},

                # Photon
                {'name':'photon_id_sf_SelectedPhoton','title':'CMS_hza_photon_id','type':'factory','prior':'lnN','correlateAcrossYears':1},

                # MVA training variable reweighting.
                # Filled in root_MVAcut signal trees as:
                #   weight_mva_reweight_central / weight_mva_reweight_Up / weight_mva_reweight_Down
                {'name':'mva_reweight','title':'CMS_hza_mva_reweight','type':'factory','prior':'lnN','correlateAcrossYears':1},

                # Electrons
                {'name':'electron_reco_sf_SelectedElectron','title':'CMS_hza_electron_reco','type':'factory','prior':'lnN','correlateAcrossYears':1},
                {'name':'electron_wplid_sf_SelectedElectron','title':'CMS_hza_electron_id','type':'factory','prior':'lnN','correlateAcrossYears':1},
                {'name':'electron_wplid_sf_nomatch_SelectedGenNoRecoElectron','title':'CMS_hza_electron_id_nomatch','type':'factory','prior':'lnN','correlateAcrossYears':1},

                # Muons
                {'name':'muon_reco_sf_SelectedMuon','title':'CMS_hza_muon_reco','type':'factory','prior':'lnN','correlateAcrossYears':1},
                {'name':'muon_looseid_sf_SelectedMuon','title':'CMS_hza_muon_id','type':'factory','prior':'lnN','correlateAcrossYears':1},
                {'name':'muon_looseid_sf_nomatch_SelectedGenNoRecoMuon','title':'CMS_hza_muon_id_nomatch','type':'factory','prior':'lnN','correlateAcrossYears':1},
              ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Shape nuisances: effect encoded in signal model
# mode = (other,scalesGlobal,scales,scalesCorr,smears): match the definition in the signal models

signal_shape_systematics = [
                {'name':'FNUF','title':'FNUF','type':'signal_shape','mode':'scalesCorr','mean':'0.0','sigma':'1.0'},
                {'name':'Material','title':'Material','type':'signal_shape','mode':'scalesCorr','mean':'0.0','sigma':'1.0'},

                {'name':'ElectronScale','title':'ElectronScale','type':'signal_shape','mode':'scales','mean':'0.0','sigma':'1.0'},
                {'name':'ElectronSmear','title':'ElectronSmear','type':'signal_shape','mode':'smears','mean':'0.0','sigma':'1.0'},
                {'name':'MuonScale','title':'MuonScale','type':'signal_shape','mode':'scales','mean':'0.0','sigma':'1.0'},
                {'name':'MuonSmear','title':'MuonSmear','type':'signal_shape','mode':'smears','mean':'0.0','sigma':'1.0'},
                {'name':'PhotonScale','title':'PhotonScale','type':'signal_shape','mode':'scales','mean':'0.0','sigma':'1.0'},
                {'name':'PhotonSmear','title':'PhotonSmear','type':'signal_shape','mode':'smears','mean':'0.0','sigma':'1.0'},
              ]
