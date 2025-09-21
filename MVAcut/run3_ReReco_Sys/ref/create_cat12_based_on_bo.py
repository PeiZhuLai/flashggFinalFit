import awkward as ak
import sys
debug=False
import numpy as np
import os
import vector
vector.register_awkward()
import mplhep as hep
hep.style.use(hep.style.CMS)
import re
import pandas as pd
import glob
from scipy.optimize import minimize
import argparse
import random
from tqdm import tqdm
import matplotlib.pyplot as plt
import ast
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler
# 自定义函数，将参数字符串转换为列表
def str_to_list(arg):
    return ast.literal_eval(arg)
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--is_HH', type=str, default="False", help='is HH')
parser.add_argument('--inputFHFiles',type=str_to_list, help='inputFHFiles List')
parser.add_argument('--inputBKGFiles',type=str_to_list, help='input pp and dd Files List')
parser.add_argument('--data',type=str,default="/eos/user/s/shsong/HiggsDNA/UL18data/merged_nominal.parquet", help='data file')
parser.add_argument('--year',type=str,default="2017", help='year')
parser.add_argument('--local', default=False, action='store_true', help="run locally or on condor")
parser.add_argument('--model',type=str,default="/eos/user/z/zhenxuan/PNN_wwgg/boosted_FHSL/data/simple_DNN_train_boosted_100_epoch_all_lr000001/model.pth", help='model file') 
parser.add_argument('--scalar',type=str,default="/eos/user/z/zhenxuan/PNN_wwgg/boosted_FHSL/data/simple_DNN_train_boosted_100_epoch_all_lr000001/scaler_params.json", help='scalar file') 



args = parser.parse_args()
local = args.local

datapath = args.data
model_path = args.model
scalar_path = args.scalar
year = args.year
print(year)
if local:
    sys.path.append("/eos/user/s/shsong/pkgs_condor/parquet_to_root-0.3.0")
    sys.path.append("/eos/user/s/shsong/pkgs_condor/bayesian-optimization-1.4.3")
else:
    sys.path.append("./PBDT_HH_FHSL_combine_"+year+"/pkgs_condor/parquet_to_root-0.3.0")
    sys.path.append("./PBDT_HH_FHSL_combine_"+year+"/pkgs_condor/bayesian-optimization-1.4.3")
if debug:
    print(sys.path)
from parquet_to_root import parquet_to_root
from bayes_opt import BayesianOptimization
from bayes_opt.util import UtilityFunction

class MultiClassDNN_model(nn.Module):
    def __init__(self, input_size, output_size):
        super(MultiClassDNN_model, self).__init__()
        self.fc1 = nn.Sequential(
            nn.Linear(input_size, 256),
            nn.BatchNorm1d(256),
            nn.SiLU(),
            nn.Dropout(0.5)
        )
        self.fc2 = nn.Sequential(
            nn.Linear(256, 528),
            nn.BatchNorm1d(528),
            nn.SiLU(),
            nn.Dropout(0.5)
        )
        self.fc3 = nn.Sequential(
            nn.Linear(528, 528),
            nn.BatchNorm1d(528),
            nn.SiLU(),
            nn.Dropout(0.5)
        )
        self.fc4 = nn.Sequential(
            nn.Linear(528, 256),
            nn.BatchNorm1d(256),
            nn.SiLU(),
            nn.Dropout(0.5)
        )
        self.fc5 = nn.Sequential(
            nn.Linear(256, 64),
            nn.BatchNorm1d(64),
            nn.SiLU(),
            nn.Dropout(0.5)
        )
        self.fc6 = nn.Linear(64, output_size)

    def forward(self, x):
        out = self.fc1(x)
        out = self.fc2(out)
        out = self.fc3(out)
        out = self.fc4(out)
        out = self.fc5(out)
        out = self.fc6(out)
        return out
def kinematic_reweight( events_data, events_bkg_pp, events_bkg_dd, weight_data, weight_bkg_pp, weight_bkg_dd, var_name_list, bins_list):
    for i in range(len(var_name_list)):
        data_array = np.array(events_data[var_name_list[i]])
        bkg_array_pp = np.array(events_bkg_pp[var_name_list[i]])
        bkg_array_dd = np.array(events_bkg_dd[var_name_list[i]])
        weight_bkg_dd = np.array(weight_bkg_dd)
        weight_bkg_pp = np.array(weight_bkg_pp)
        print("data_array", data_array)
        hist_data, bins = np.histogram(data_array, bins=100, range=(bins_list[i][0], bins_list[i][1]), weights=np.array(weight_data))
        hist_dd, bins = np.histogram(bkg_array_dd, bins=100, range=(bins_list[i][0], bins_list[i][1]), weights=np.array(weight_bkg_dd))
        hist_pp, bins = np.histogram(bkg_array_pp, bins=100, range=(bins_list[i][0], bins_list[i][1]), weights=np.array(weight_bkg_pp))
        for i in range(len(bins)-1):
            reweight_factor = (hist_data[i]-hist_dd[i])/hist_pp[i]
            weight_bkg_pp[(bkg_array_pp>bins[i]) & (bkg_array_pp<bins[i+1])] *= reweight_factor
    return weight_bkg_pp

def target_function(bo1, bo2,s, b, bkg):
    if debug:
        print("bo1: ", bo1)
        print("bo2: ", bo2)
    events_sig1 = s[(s['PNN_score'] >= 1-bo1) & (s['PNN_score'] <=1)]
    events_datasideband_1 = b[(b['PNN_score'] >= 1-bo1) & (b['PNN_score'] <=1)]
    events_bkg1 = bkg[(bkg['PNN_score'] >= 1-bo1) & (bkg['PNN_score'] <=1)]
    events_sig2 = s[(s['PNN_score'] >= 1-bo2-bo1) & (s['PNN_score'] <1-bo1)]
    events_datasideband_2 = b[(b['PNN_score'] >= 1-bo2-bo1) & (b['PNN_score'] <1-bo1)]
    events_bkg2 = bkg[(bkg['PNN_score'] >= 1-bo2-bo1) & (bkg['PNN_score'] <1-bo1)]
    # force the number of data events in the mgg 100-115 to be at least 2
    if len(b[(b['PNN_score'] >= 1-bo1) & (b['Diphoton_mass'] < 115)]) < 2:
        significance=0
        return 0
    if len(events_datasideband_1[((events_datasideband_1.Diphoton_mass>135)|(events_datasideband_1.Diphoton_mass<115))]) < 5:
        significance=0
        return 0
    elif len(events_datasideband_2[((events_datasideband_2.Diphoton_mass>135)|(events_datasideband_2.Diphoton_mass<115))]) < 5:
        significance=0
        return 0
    if debug:
        print("len(events_datasideband_1): ", len(events_datasideband_1))
        print("len(events_datasideband_2): ", len(events_datasideband_2))
    # 1
    # check signal mass distribution
    hist_sig1, bin_sig1 = np.histogram(np.array(events_sig1['Diphoton_mass']), bins=100, range=(110,150), weights=np.array(events_sig1['weight_central']))
    if debug:
        print("hist_sig1: ", hist_sig1)
        print("bin_sig1: ", bin_sig1)
    # get signal mass FWHM
    # fwhm1 = get_fwhm(hist_sig1, bin_sig1)
    # get significance s/sqrt(s+b)
    s1 = np.sum(events_sig1['weight_central'][(events_sig1['Diphoton_mass'] > 115) & (events_sig1['Diphoton_mass'] < 135)])
    bkg1_weight = (events_bkg1['weight_central'])[(events_bkg1['Diphoton_mass'] > 115) & (events_bkg1['Diphoton_mass'] < 135)]
    b1 = np.sum(bkg1_weight) 
    d1 = np.sum(events_datasideband_1['weight_central'][(events_datasideband_1['Diphoton_mass'] > 135) | (events_datasideband_1['Diphoton_mass'] < 115)])
    if d1>10:
        return 0
    # significance1 = s1 / np.sqrt(b1)
    significance1 = s1 / np.sqrt(d1)
    # 2
    # check signal mass distribution
    hist_sig2, bin_sig2 = np.histogram(np.array(events_sig2['Diphoton_mass']), bins=100, range=(110,150), weights=np.array(events_sig2['weight_central']))
    # check signal mass FWHM
    # fwhm2 = get_fwhm(hist_sig2, bin_sig2)
    # get significance s/sqrt(s+b)
    s2 = np.sum(events_sig2['weight_central'][(events_sig2['Diphoton_mass'] > 115) & (events_sig2['Diphoton_mass'] < 135)])
    bkg2_weight = (events_bkg2['weight_central'])[(events_bkg2['Diphoton_mass'] > 115) & (events_bkg2['Diphoton_mass'] < 135)]
    b2 = np.sum(bkg2_weight) 
    d2 = np.sum(events_datasideband_2['weight_central'][(events_datasideband_2['Diphoton_mass'] > 135) | (events_datasideband_2['Diphoton_mass'] < 115)])
    # significance2 = s2 / np.sqrt(b2)
    significance2 = s2 / np.sqrt(d2)

    significance = np.sqrt(significance1**2 + significance2**2 )
    # fwhm = np.sqrt(fwhm1**2 + fwhm2**2)
    return significance

def get_jet_variables(events):
    # get fatjet pt
    events['fatjet_1_pt'] = events['fatjet_1_pt']
    events['fatjet_2_pt'] = events['fatjet_2_pt']
    events['fatjet_1_eta'] = events['fatjet_1_eta']
    events['fatjet_2_eta'] = events['fatjet_2_eta']
    events['fatjet_1_phi'] = events['fatjet_1_phi']
    events['fatjet_2_phi'] = events['fatjet_2_phi']
    events['fatjet_1_msoftdrop'] = events['fatjet_1_msoftdrop']
    events['fatjet_2_msoftdrop'] = events['fatjet_2_msoftdrop']
    # get fatjet and diphoton dR
    fatjet_1_4D = vector.obj(pt=events.fatjet_1_pt, eta=events.fatjet_1_eta, phi=events.fatjet_1_phi, mass=events.fatjet_1_msoftdrop)
    fatjet_2_4D = vector.obj(pt=events.fatjet_2_pt, eta=events.fatjet_2_eta, phi=events.fatjet_2_phi, mass=events.fatjet_2_msoftdrop)
    fatjet_3_4D = vector.obj(pt=events.fatjet_3_pt, eta=events.fatjet_3_eta, phi=events.fatjet_3_phi, mass=events.fatjet_3_msoftdrop)
    diphoton_4D = vector.obj(pt=events.Diphoton_pt, eta=events.Diphoton_eta, phi=events.Diphoton_phi, mass=events.Diphoton_mass)
    fatjet_1_diphoton_dR = fatjet_1_4D.deltaR(diphoton_4D)
    fatjet_2_diphoton_dR = fatjet_2_4D.deltaR(diphoton_4D)
    events['fatjet_1_diphoton_dR'] = np.where(events.fatjet_1_pt>0, fatjet_1_diphoton_dR, -999)
    events['fatjet_2_diphoton_dR'] = np.where(events.fatjet_2_pt>0, fatjet_2_diphoton_dR, -999)
    leadphoton=vector.obj(pt=events.LeadPhoton_pt,eta=events.LeadPhoton_eta,phi=events.LeadPhoton_phi,mass=events.LeadPhoton_mass)
    subleadphoton=vector.obj(pt=events.SubleadPhoton_pt,eta=events.SubleadPhoton_eta,phi=events.SubleadPhoton_phi,mass=events.SubleadPhoton_mass)
    # dR with photon and fatjet
    fatjet_1_leadphoton_dR = fatjet_1_4D.deltaR(leadphoton)
    events['fatjet_1_leadphoton_dR'] = np.where(events.fatjet_1_pt>0, fatjet_1_leadphoton_dR, -999)
    fatjet_1_subleadphoton_dR = fatjet_1_4D.deltaR(subleadphoton)
    events['fatjet_1_subleadphoton_dR'] = np.where(events.fatjet_1_pt>0, fatjet_1_subleadphoton_dR, -999)
    fatjet_2_leadphoton_dR = fatjet_2_4D.deltaR(leadphoton)
    events['fatjet_2_leadphoton_dR'] = np.where(events.fatjet_2_pt>0, fatjet_2_leadphoton_dR, -999)
    fatjet_2_subleadphoton_dR = fatjet_2_4D.deltaR(subleadphoton)
    events['fatjet_2_subleadphoton_dR'] = np.where(events.fatjet_2_pt>0, fatjet_2_subleadphoton_dR, -999)
    # get fatjet1 and fatjet2 dR
    fatjet_1_2_dR = fatjet_1_4D.deltaR(fatjet_2_4D)
    events['fatjet_1_2_dR'] = np.where(np.logical_and(events.fatjet_1_pt>0, events.fatjet_2_pt>0), fatjet_1_2_dR, -999)
    # get the maximum fatjets mass with the combination of 2 fatjets in 3 fatjets
    fatjet_12_4D = fatjet_1_4D+fatjet_2_4D
    fatjet_13_4D = fatjet_1_4D+fatjet_3_4D
    fatjet_23_4D = fatjet_2_4D+fatjet_3_4D
    fatjet_12_msoftdrop = np.where(((fatjet_1_4D.pt >0) & (fatjet_2_4D.pt >0)), fatjet_12_4D.mass, -999)
    fatjet_13_msoftdrop = np.where(((fatjet_1_4D.pt >0) & (fatjet_3_4D.pt >0)), fatjet_13_4D.mass, -999)
    fatjet_23_msoftdrop = np.where(((fatjet_2_4D.pt >0) & (fatjet_3_4D.pt >0)), fatjet_23_4D.mass, -999)
    max_fatjets_mass = np.maximum(fatjet_12_msoftdrop, np.maximum(fatjet_13_msoftdrop, fatjet_23_msoftdrop))
    events['max_fatjets_mass'] = max_fatjets_mass

    
    

    
    # get max WvsQCD score with three fatjets
    events['fatjet_1_WvsQCDMD'] = events['fatjet_1_WvsQCDMD']
    events['fatjet_2_WvsQCDMD'] = events['fatjet_2_WvsQCDMD']
    # get max H4qvsQCD score with three fatjets
    events['fatjet_1_Hqqqq_vsQCDTop'] = events['fatjet_1_Hqqqq_vsQCDTop']
    events['fatjet_2_Hqqqq_vsQCDTop'] = events['fatjet_2_Hqqqq_vsQCDTop']
    # get XbbvsQCD
    events['fatjet_1_XbbvsQCDMD'] = (events['fatjet_1_particleNetMD_Xbb']) / (events['fatjet_1_particleNetMD_Xbb'] + events['fatjet_1_particleNetMD_QCD'])
    events['fatjet_2_XbbvsQCDMD'] = (events['fatjet_2_particleNetMD_Xbb']) / (events['fatjet_2_particleNetMD_Xbb'] + events['fatjet_2_particleNetMD_QCD'])
    # add number of good AK4 jets
    events['nGoodAK4jets'] = events['nGoodAK4jets']
    # add number of good AK8 jets
    events['nGoodAK8jets'] = events['nGoodAK8jets']
    # get 4 ak4 jets 4D info
    events['jet_1_pt'] = events['jet_1_pt']
    events['jet_2_pt'] = events['jet_2_pt']
    events['jet_3_pt'] = events['jet_3_pt']
    events['jet_1_eta'] = events['jet_1_eta']
    events['jet_2_eta'] = events['jet_2_eta']
    events['jet_3_eta'] = events['jet_3_eta']
    events['jet_1_phi'] = events['jet_1_phi']
    events['jet_2_phi'] = events['jet_2_phi']
    events['jet_3_phi'] = events['jet_3_phi']
    events['jet_1_mass'] = events['jet_1_mass']
    events['jet_2_mass'] = events['jet_2_mass']
    events['jet_3_mass'] = events['jet_3_mass']
    

    return events


def get_sig_events_forApply(filename,mx,my):
    events = ak.from_parquet(filename)
    category_cut = ((events["category"]==1) | (events["category"]==2))
    photonID_cut = (events["Diphoton_minID_modified"]>-0.7)
    events = events[ category_cut & photonID_cut]
    events = get_jet_variables(events)
    events['mx'] = np.ones(len(events))*int(mx)
    events['my'] = np.ones(len(events))*int(my)
    return events
print('start to get input features')
input_features = ['Diphoton_pt','Diphoton_eta','Diphoton_phi','LeadPhoton_pt','LeadPhoton_eta','LeadPhoton_phi', 'SubleadPhoton_pt','SubleadPhoton_eta','SubleadPhoton_phi','Diphoton_dR','fatjet_1_pt','fatjet_2_pt','fatjet_1_eta','fatjet_2_eta','fatjet_1_phi','fatjet_2_phi','fatjet_1_diphoton_dR','fatjet_2_diphoton_dR','fatjet_1_leadphoton_dR','fatjet_1_subleadphoton_dR','fatjet_2_leadphoton_dR','fatjet_2_subleadphoton_dR','fatjet_1_2_dR','max_fatjets_mass','nGoodAK4jets','nGoodAK8jets','jet_1_pt','jet_2_pt','jet_3_pt','jet_1_eta','jet_2_eta','jet_3_eta','jet_1_phi','jet_2_phi','jet_3_phi','jet_1_mass','jet_2_mass','jet_3_mass','nGoodisoleptons','nGoodnonisoleptons','PuppiMET_pt','PuppiMET_sumEt','mx']
other_vars  = ['weight_central','Diphoton_mass']
FHsignal_path_list = []
SLsignal_path_list = []
ZZggsignal_path_list = []
TTggsignal_path_list = []
BBGGsignal_path_list = []
vbf_path_list = []
vh_path_list = []
tth_path_list = []
ggh_path_list = []
signal_output_name=[]
bbgg_output_name=[]
zzgg_output_name=[]
ttgg_output_name=[]
vbf_output_name=[]
vh_output_name=[]
tth_output_name=[]
ggh_output_name=[]
list_of_files = args.inputFHFiles #merged_nominal.parquet should be the first one
for FHfile in list_of_files:
    FHsignal_path_list.append(FHfile)
    SLfile=(FHfile.replace("2G4Q","2G2Q1L1Nu")).replace("HHFH","HHSL")
    bbggfile=(FHfile.replace("2G2WTo2G4Q","2B2G")).replace("HHFH","HHbbgg")
    zzggfile=(FHfile.replace("2G2W","2G2Z")).replace("HHFH","HHZZgg")
    ttggfile=(FHfile.replace("2G2WTo2G4Q","2G2Tau")).replace("HHFH","HHttgg")
    if "2018" in FHfile or "2017" in FHfile:
        path_year = year
    else:
        path_year = "2016UL_"+year.split("2016")[1]+"VFP"
    vbffile="/eos/user/s/shsong/HiggsDNA/SingleHiggs"+year.split("20")[1]+"/VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8_"+path_year+"/"+FHfile.split("/")[-1]
    vhfile="/eos/user/s/shsong/HiggsDNA/SingleHiggs"+year.split("20")[1]+"/VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_"+path_year+"/"+FHfile.split("/")[-1]
    tthfile="/eos/user/s/shsong/HiggsDNA/SingleHiggs"+year.split("20")[1]+"/ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_"+path_year+"/"+FHfile.split("/")[-1]
    gghfile="/eos/user/s/shsong/HiggsDNA/SingleHiggs"+year.split("20")[1]+"/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+path_year+"/"+FHfile.split("/")[-1]
    SLsignal_path_list.append(SLfile)
    ZZggsignal_path_list.append(zzggfile)
    TTggsignal_path_list.append(ttggfile)
    BBGGsignal_path_list.append(bbggfile)
    vbf_path_list.append(vbffile)
    vh_path_list.append(vhfile)
    tth_path_list.append(tthfile)
    ggh_path_list.append(gghfile)
    dir_name="CombineFHSL_MX" + FHfile.split("M-")[1].split("_")[0] + "_MH125_cat12_"+(FHfile.split("/")[-1]).split(".")[0]
    # dir_name = "CombineFHSL_MX1100_MH125_cat12_merged_FJER_down"
    signal_output_name.append(dir_name)
    bbgg_output_name.append(dir_name.replace("CombineFHSL","BBGG"))
    zzgg_output_name.append(dir_name.replace("CombineFHSL","ZZGG"))
    ttgg_output_name.append(dir_name.replace("CombineFHSL","TTGG"))
    vbf_output_name.append(dir_name.replace("CombineFHSL","VBF"))
    vh_output_name.append(dir_name.replace("CombineFHSL","VH"))
    tth_output_name.append(dir_name.replace("CombineFHSL","TTH"))
    ggh_output_name.append(dir_name.replace("CombineFHSL","GGH"))    
signal_samples = {'FHpath':FHsignal_path_list, 'sig_output_name':signal_output_name, 'SLpath':SLsignal_path_list, 'ZZggpath':ZZggsignal_path_list, 'TTggpath':TTggsignal_path_list, 'BBGGpath':BBGGsignal_path_list, 'VBFpath':vbf_path_list,'VHpath':vh_path_list,'ttHpath':tth_path_list,'ggHpath':ggh_path_list,'bbgg_output_name':bbgg_output_name, 'zzgg_output_name':zzgg_output_name,'ttgg_output_name':ttgg_output_name, 'vbf_output_name':vbf_output_name, 'vh_output_name':vh_output_name, 'tth_output_name':tth_output_name, 'ggh_output_name':ggh_output_name}

bkgfiles = args.inputBKGFiles

#-------------------------------------Get the boundary-------------------------------------#
#firstly use the merged_nominal.parquet to get the boundary
mx = signal_samples['sig_output_name'][0].split('_')[1].split('X')[1]
my = signal_samples['sig_output_name'][0].split('_')[2].split('H')[1]
events_sigFH = get_sig_events_forApply(signal_samples['FHpath'][0],mx,my)
events_sigFH['signal'] = 2*np.ones(len(events_sigFH))
events_sigSL = get_sig_events_forApply(signal_samples['SLpath'][0],mx,my)
events_sigSL['signal'] = np.ones(len(events_sigSL))

events_bbgg = get_sig_events_forApply(signal_samples['BBGGpath'][0],mx,my)
events_zzgg = get_sig_events_forApply(signal_samples['ZZggpath'][0],mx,my)
events_ttgg = get_sig_events_forApply(signal_samples['TTggpath'][0],mx,my)

events_data = get_sig_events_forApply(datapath,mx,my)
events_pp_cat1= get_sig_events_forApply(bkgfiles[0],mx,my)
events_pp_cat2= get_sig_events_forApply(bkgfiles[1],mx,my)
events_dd_cat1= get_sig_events_forApply(bkgfiles[2],mx,my)
events_dd_cat2= get_sig_events_forApply(bkgfiles[3],mx,my)
events_vbf= get_sig_events_forApply(signal_samples['VBFpath'][0],mx,my)
events_vh= get_sig_events_forApply(signal_samples['VHpath'][0],mx,my)
events_tth= get_sig_events_forApply(signal_samples['ttHpath'][0],mx,my)
events_ggh= get_sig_events_forApply(signal_samples['ggHpath'][0],mx,my)
import json
model = MultiClassDNN_model(len(input_features), 5)
state_dict = torch.load(model_path, map_location=torch.device('cpu'))
new_state_dict = {key.replace("module.", ""): value for key, value in state_dict.items()}
model.load_state_dict(new_state_dict)
model.eval()
with open(scalar_path, 'r') as f:
    loaded_params = json.load(f)
    loaded_scaler = StandardScaler()
    loaded_scaler.mean_ = np.array(loaded_params['mean'])
    loaded_scaler.scale_ = np.array(loaded_params['scale'])
print('successfully load the model and scalar')
def model_predict(event, model, loaded_scaler, input_features):
    df = ak.to_pandas(event[input_features + ['weight_central','Diphoton_mass']])
    df = df.replace(-999, 0)
    X_test = loaded_scaler.transform(df[input_features])
    X_test = torch.tensor(X_test).float()
    proba = model(X_test)
    proba = F.softmax(proba, dim=1)
    proba = proba.detach().numpy()
    pnn_score = (proba[:,4] + proba[:,2] + proba[:,3]) / (proba[:,0] + proba[:,1] + proba[:,2] + proba[:,3] + proba[:,4])

    del proba
    del df
    del X_test
    return pnn_score
# for sig, data, bkgmc, bbgg and zzgg
# evaluate the model


events_pp = ak.concatenate([events_pp_cat1,events_pp_cat2])
del events_pp_cat1
del events_pp_cat2
events_dd = ak.concatenate([events_dd_cat1,events_dd_cat2])
del events_dd_cat1
del events_dd_cat2
MC_pp_new_weight = kinematic_reweight(events_data=events_data, events_bkg_pp=events_pp, events_bkg_dd=events_dd, weight_data=events_data['weight_central'], weight_bkg_pp=events_pp['weight_central'], weight_bkg_dd=events_dd['weight_central'], var_name_list = ['Diphoton_pt'], bins_list = [[0, 300]])
events_pp['weight_central'] = MC_pp_new_weight

PNN_score = model_predict(events_pp, model, loaded_scaler, input_features)
events_pp['PNN_score'] = PNN_score

PNN_score = model_predict(events_dd, model, loaded_scaler, input_features)
events_dd['PNN_score'] = PNN_score
event_bkgmc=ak.concatenate([events_pp,events_dd])
del events_pp
del events_dd
print("get the PNN score for bkgmc")
PNN_score = model_predict(events_sigFH, model, loaded_scaler, input_features)
events_sigFH['PNN_score'] = PNN_score
PNN_score = model_predict(events_sigSL, model, loaded_scaler, input_features)
events_sigSL['PNN_score'] = PNN_score
print("get the PNN score for signal")
events_sig = ak.concatenate([events_sigFH,events_sigSL])
del events_sigFH
del events_sigSL

PNN_score = model_predict(events_bbgg, model, loaded_scaler, input_features)
events_bbgg['PNN_score'] = PNN_score
print("get the PNN score for bbgg")
del PNN_score

PNN_score = model_predict(events_data, model, loaded_scaler, input_features)
events_data['PNN_score'] = PNN_score
print("get the PNN score for data")

PNN_score = model_predict(events_zzgg, model, loaded_scaler, input_features)
events_zzgg['PNN_score'] = PNN_score

PNN_score = model_predict(events_ttgg, model, loaded_scaler, input_features)
events_ttgg['PNN_score'] = PNN_score
print('get all the PBDT score for merged nominal.parquet')
PNN_score = model_predict(events_ggh, model, loaded_scaler, input_features)
events_ggh['PNN_score'] = PNN_score
PNN_score = model_predict(events_vbf, model, loaded_scaler, input_features)
events_vbf['PNN_score'] = PNN_score
PNN_score = model_predict(events_vh, model, loaded_scaler, input_features)
events_vh['PNN_score'] = PNN_score
PNN_score = model_predict(events_tth, model, loaded_scaler, input_features)
events_tth['PNN_score'] = PNN_score
import mplhep as hep
import matplotlib.pyplot as plt
Xmass=FHfile.split("M-")[1].split("_")[0]
directory_path = "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/MX"+Xmass+"_MH125"
if os.path.exists(directory_path):
    rmcommand="rm -rf "+directory_path
    os.system(rmcommand)
    print(f"The directory {directory_path} exists.")
    os.mkdir(directory_path)
else:
    os.mkdir(directory_path)

plt.figure(figsize=(8, 6))
hep.style.use("CMS")
plt.hist(events_sig['PNN_score'],weights = 30*events_sig.weight_central ,range=(0,1),bins=20, histtype='step', label='30*signal',color='#FF0000')
plt.hist(events_bbgg['PNN_score'],weights = 30*events_bbgg.weight_central ,range=(0,1),bins=20, histtype='step', label='30*bbgg',color='#FFA500')
plt.hist(event_bkgmc['PNN_score'][(event_bkgmc.Diphoton_mass>135)|(event_bkgmc.Diphoton_mass<115)],weights = event_bkgmc.weight_central[(event_bkgmc.Diphoton_mass>135)|(event_bkgmc.Diphoton_mass<115)] ,range=(0,1),bins=20, histtype='stepfilled', label='pp+dd',color='#00BFFF')
hist, bins = np.histogram(events_data[(events_data.Diphoton_mass>135)|(events_data.Diphoton_mass<115)]['PNN_score'],bins=20)
non_zero_bins = hist > 0
bin_centers = 0.5 * (bins[:-1] + bins[1:])
plt.scatter(bin_centers[non_zero_bins], hist[non_zero_bins], marker='o', color='black', label='data')
plt.xlabel('PNN Score')
plt.ylabel('Events')
plt.yscale('log')
#y axis unit small
plt.legend()
plt.savefig("./PBDT_HH_FHSL_combine_"+year+"/flashgginput/MX"+Xmass+"_MH125/dnnscore.png", dpi=140)
print("save the PNN score distribution plot")
print("./PBDT_HH_FHSL_combine_"+year+"/flashgginput/MX"+Xmass+"_MH125/dnnscore.png")
plt.close()

json_name="/eos/cms/store/group/phys_b2g/shsong/flashggws/cat12/"+year+"/MX"+Xmass+"_MH125/boundaries.json"
with open(json_name, 'r') as f:
    boundaries = json.load(f)
cut1=boundaries[Xmass]['cut1']
cut2= boundaries[Xmass]['cut2']
boundary1=[]
boundary2=[]
boundary1.append(cut1)
boundary2.append(cut2)
print('high purity cut:',cut1)
print('low purity cut:',cut2)
FHhighpurityevents = events_sig[(events_sig.PNN_score>cut1)&(events_sig.signal==2)]
SLhighpurityevents = events_sig[(events_sig.PNN_score>cut1)&(events_sig.signal==1)]
FHlowpurityevents = events_sig[((events_sig.PNN_score>cut2)&(events_sig.PNN_score<=cut1)&(events_sig.signal==2))|(((events_sig.category==3)&(events_sig.PNN_score>cut2)&(events_sig.PNN_score<=cut1)&(events_sig.signal==2)))]
SLlowpurityevents = events_sig[((events_sig.PNN_score>cut2)&(events_sig.PNN_score<=cut1)&(events_sig.signal==1))|(((events_sig.category==3)&(events_sig.PNN_score>cut2)&(events_sig.PNN_score<=cut1)&(events_sig.signal==1)))]
highpurity_sigeff = ak.sum(events_sig[events_sig.PNN_score>cut1].weight_central)/ak.sum(events_sig.weight_central)
lowpurity_sigeff = ak.sum(events_sig[((events_sig.PNN_score>cut2)&(events_sig.PNN_score<=cut1))].weight_central)/ak.sum(events_sig.weight_central)

FHhighpurity_sigeff = ak.sum(FHhighpurityevents.weight_central)/ak.sum(events_sig[events_sig.signal==2].weight_central)
FHlowpurity_sigeff = ak.sum(FHlowpurityevents.weight_central)/ak.sum(events_sig[events_sig.signal==2].weight_central)
SLhighpurity_sigeff = ak.sum(SLhighpurityevents.weight_central)/ak.sum(events_sig[events_sig.signal==1].weight_central)
SLlowpurity_sigeff = ak.sum(SLlowpurityevents.weight_central)/ak.sum(events_sig[events_sig.signal==1].weight_central)
vbfhighpurity_sigeff = ak.sum(events_vbf[(events_vbf.PNN_score>cut1)].weight_central)/ak.sum(events_vbf.weight_central)
vbfhighpurity_sigeff = ak.sum(events_vbf[(events_vbf.PNN_score>cut1)].weight_central)/ak.sum(events_vbf.weight_central)
vhhighpurity_sigeff = ak.sum(events_vh[(events_vh.PNN_score>cut1)].weight_central)/ak.sum(events_vh.weight_central)
vhlowpurity_sigeff = ak.sum(events_vh[((events_vh.PNN_score>cut2)&(events_vh.PNN_score<=cut1))].weight_central)/ak.sum(events_vh.weight_central)
tthhighpurity_sigeff = ak.sum(events_tth[(events_tth.PNN_score>cut1)].weight_central)/ak.sum(events_tth.weight_central)
tthlowpurity_sigeff = ak.sum(events_tth[((events_tth.PNN_score>cut2)&(events_tth.PNN_score<=cut1))].weight_central)/ak.sum(events_tth.weight_central)
gghhighpurity_sigeff = ak.sum(events_ggh[(events_ggh.PNN_score>cut1)].weight_central)/ak.sum(events_ggh.weight_central)
gghlowpurity_sigeff = ak.sum(events_ggh[((events_ggh.PNN_score>cut2)&(events_ggh.PNN_score<=cut1))].weight_central)/ak.sum(events_ggh.weight_central)

bbggFHhighpurity_sigeff = ak.sum(events_bbgg[(events_bbgg.PNN_score>cut1)&(events_bbgg.category==2)].weight_central)/ak.sum(events_bbgg.weight_central)
bbggFHlowpurity_sigeff = ak.sum(events_bbgg[((events_bbgg.PNN_score>cut2)&(events_bbgg.PNN_score<=cut1)&(events_bbgg.category==2))].weight_central)/ak.sum(events_bbgg.weight_central)
zzggFHhighpurity_sigeff = ak.sum(events_zzgg[(events_zzgg.PNN_score>cut1)&(events_zzgg.category==2)].weight_central)/ak.sum(events_zzgg.weight_central)
zzggFHlowpurity_sigeff = ak.sum(events_zzgg[((events_zzgg.PNN_score>cut2)&(events_zzgg.PNN_score<=cut1)&(events_zzgg.category==2))].weight_central)/ak.sum(events_zzgg.weight_central)
ttggFHhighpurity_sigeff = ak.sum(events_ttgg[(events_ttgg.PNN_score>cut1)&(events_ttgg.category==2)].weight_central)/ak.sum(events_ttgg.weight_central)
ttggFHlowpurity_sigeff = ak.sum(events_ttgg[((events_ttgg.PNN_score>cut2)&(events_ttgg.PNN_score<=cut1)&(events_ttgg.category==2))].weight_central)/ak.sum(events_ttgg.weight_central)
bbggSLhighpurity_sigeff = ak.sum(events_bbgg[(events_bbgg.PNN_score>cut1)&(events_bbgg.category==1)].weight_central)/ak.sum(events_bbgg.weight_central)
bbggSLlowpurity_sigeff = ak.sum(events_bbgg[((events_bbgg.PNN_score>cut2)&(events_bbgg.PNN_score<=cut1)&(events_bbgg.category==1))].weight_central)/ak.sum(events_bbgg.weight_central)
zzggSLhighpurity_sigeff = ak.sum(events_zzgg[(events_zzgg.PNN_score>cut1)&(events_zzgg.category==1)].weight_central)/ak.sum(events_zzgg.weight_central)
zzggSLlowpurity_sigeff = ak.sum(events_zzgg[((events_zzgg.PNN_score>cut2)&(events_zzgg.PNN_score<=cut1)&(events_zzgg.category==1))].weight_central)/ak.sum(events_zzgg.weight_central)
ttggSLhighpurity_sigeff = ak.sum(events_ttgg[(events_ttgg.PNN_score>cut1)&(events_ttgg.category==1)].weight_central)/ak.sum(events_ttgg.weight_central)
ttggSLlowpurity_sigeff = ak.sum(events_ttgg[((events_ttgg.PNN_score>cut2)&(events_ttgg.PNN_score<=cut1)&(events_ttgg.category==1))].weight_central)/ak.sum(events_ttgg.weight_central)
highpurity_sidebandnum=len(events_data[(((events_data.Diphoton_mass>135)|(events_data.Diphoton_mass<115))&(events_data.PNN_score > cut1) & (events_data.PNN_score <= 1))])
lowpurity_sidebandnum=len(events_data[(((events_data.Diphoton_mass>135)|(events_data.Diphoton_mass<115))&(events_data.PNN_score > cut2) & (events_data.PNN_score <= cut1))])
def add_HtaggerSF(event,mass):
    mass_list=[500,550,600,650,700,750,800,850,900,1000,1100,1200,1250,1300,1400,1500,1600,1700,1750,1800,1900,2000,2200,2400,2500,2600,2800,3000]
    # mass_list=[1000]
    index=mass_list.index(mass)
    HvsQCD_SF=[0.993610,1.083977,0.983107,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,0.993610,1.083977,0.983107,1.134969,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,1.008370,0.962904,0.969104,1.002570,0.941431,1.025998,1.025998,1.004931,0.967810]
    HvsQCD_SFup=np.array([0.993610,1.083977,0.983107,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,0.993610,1.083977,0.983107,1.134969,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,1.008370,0.962904,0.969104,1.002570,0.941431,1.025998,1.025998,1.004931,0.967810])*1.22
    HvsQCD_SFdown=np.array([0.993610,1.083977,0.983107,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,0.993610,1.083977,0.983107,1.134969,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,1.008370,0.962904,0.969104,1.002570,0.941431,1.025998,1.025998,1.004931,0.967810])*0.75
    weight_PTransformer_up = ak.ones_like(event.category)
    weight_PTransformer_down = ak.ones_like(event.category)
    weight_PTransformer_central = ak.ones_like(event.category)
    Htagger_SF=event.category==2
    weight_PTransformer_central = ak.where(Htagger_SF, ak.ones_like(weight_PTransformer_central)*HvsQCD_SF[index], weight_PTransformer_central)
    weight_PTransformer_up = ak.where(Htagger_SF, ak.ones_like(weight_PTransformer_up)*HvsQCD_SFup[index], weight_PTransformer_up)
    weight_PTransformer_down = ak.where(Htagger_SF, ak.ones_like(weight_PTransformer_down)*HvsQCD_SFdown[index], weight_PTransformer_down)
    event['weight_PTransformer_up']=weight_PTransformer_up
    event['weight_PTransformer_down']=weight_PTransformer_down
    event['weight_PTransformer_central']=weight_PTransformer_central
    return event
events_sig = add_HtaggerSF(events_sig, int(mx))
events_bbgg = add_HtaggerSF(events_bbgg, int(mx))
events_zzgg = add_HtaggerSF(events_zzgg, int(mx))
events_ttgg = add_HtaggerSF(events_ttgg, int(mx))
events_vbf = add_HtaggerSF(events_vbf, int(mx))
events_vh = add_HtaggerSF(events_vh, int(mx))
events_ggh = add_HtaggerSF(events_ggh, int(mx))
events_tth = add_HtaggerSF(events_tth, int(mx))

def add_sf_branches(events):
    events["CMS_hgg_mass"]=events["Diphoton_mass"]
    events["weight"]=events["weight_central"]
    events["dZ"]=np.ones(len(events['CMS_hgg_mass']))
    events["muon_highptreco_sf_Down01sigma"]=events["weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_down"]/events["weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_central"]
    events["muon_highptreco_sf_Up01sigma"]=events["weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_up"]/events["weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_central"]
    events["muon_highptid_sf_Down01sigma"]=events["weight_nonisomuon_highptid_sf_SelectedMuon_noiso_down"]/events["weight_nonisomuon_highptid_sf_SelectedMuon_noiso_central"]
    events["muon_highptid_sf_Up01sigma"]=events["weight_nonisomuon_highptid_sf_SelectedMuon_noiso_up"]/events["weight_nonisomuon_highptid_sf_SelectedMuon_noiso_central"]
    events["L1_prefiring_sf_Down01sigma"]=events["weight_L1_prefiring_sf_down"]/events["weight_L1_prefiring_sf_central"]
    events["L1_prefiring_sf_Up01sigma"]=events["weight_L1_prefiring_sf_up"]/events["weight_L1_prefiring_sf_central"]
    events["puWeight_Up01sigma"]=events["weight_pu_reweight_sf_up"]/events["weight_pu_reweight_sf_central"]
    events["puWeight_Down01sigma"]=events["weight_pu_reweight_sf_down"]/events["weight_pu_reweight_sf_central"]
    events["jet_pu_id_sf_Up01sigma"]=events["weight_jet_puid_sf_SelectedJet_up"]/events["weight_jet_puid_sf_SelectedJet_central"]
    events["jet_pu_id_sf_Down01sigma"]=events["weight_jet_puid_sf_SelectedJet_down"]/events["weight_jet_puid_sf_SelectedJet_central"]
    events["electron_veto_sf_Diphoton_Photon_Up01sigma"]=events["weight_electron_veto_sf_Diphoton_Photon_up"]/events["weight_electron_veto_sf_Diphoton_Photon_central"]
    events["electron_veto_sf_Diphoton_Photon_Down01sigma"]=events["weight_electron_veto_sf_Diphoton_Photon_down"]/events["weight_electron_veto_sf_Diphoton_Photon_central"]
    events["isoelectron_id_sf_SelectedElectron_iso_Up01sigma"]=events["weight_isoelectron_id_sf_SelectedElectron_iso_up"]/events["weight_isoelectron_id_sf_SelectedElectron_iso_central"]
    events["isoelectron_id_sf_SelectedElectron_iso_Down01sigma"]=events["weight_isoelectron_id_sf_SelectedElectron_iso_down"]/events["weight_isoelectron_id_sf_SelectedElectron_iso_central"]
    events["isoelectron_id_sf_SelectedElectron_noiso_Up01sigma"]= events["weight_isoelectron_id_sf_SelectedElectron_noiso_up"]/events["weight_isoelectron_id_sf_SelectedElectron_noiso_central"]
    events["isoelectron_id_sf_SelectedElectron_noiso_Down01sigma"]= events["weight_isoelectron_id_sf_SelectedElectron_noiso_down"]/events["weight_isoelectron_id_sf_SelectedElectron_noiso_central"] 
    events["isomuon_id_sf_SelectedMuon_iso_Up01sigma"]=events["weight_isomuon_id_sf_SelectedMuon_iso_up"]/events["weight_isomuon_id_sf_SelectedMuon_iso_central"]
    events["isomuon_id_sf_SelectedMuon_iso_Down01sigma"]=events["weight_isomuon_id_sf_SelectedMuon_iso_down"]/events["weight_isomuon_id_sf_SelectedMuon_iso_central"]        
    events["isomuon_iso_sf_SelectedMuon_iso_Up01sigma"]=events["weight_isomuon_iso_sf_SelectedMuon_iso_up"]/events["weight_isomuon_iso_sf_SelectedMuon_iso_central"]
    events["isomuon_iso_sf_SelectedMuon_iso_Down01sigma"]=events["weight_isomuon_iso_sf_SelectedMuon_iso_down"]/events["weight_isomuon_iso_sf_SelectedMuon_iso_central"]        
    events["nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma"]=events["weight_nonisoelectron_id_sf_SelectedElectron_noiso_up"]/events["weight_nonisoelectron_id_sf_SelectedElectron_noiso_central"]
    events["nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma"]=events["weight_nonisoelectron_id_sf_SelectedElectron_noiso_down"]/events["weight_nonisoelectron_id_sf_SelectedElectron_noiso_central"]
    events["photon_id_sf_Diphoton_Photon_Up01sigma"]=events["weight_photon_id_sf_Diphoton_Photon_up"]/events["weight_photon_id_sf_Diphoton_Photon_central"]
    events["photon_id_sf_Diphoton_Photon_Down01sigma"]=events["weight_photon_id_sf_Diphoton_Photon_down"]/events["weight_photon_id_sf_Diphoton_Photon_central"]
    events["photon_presel_sf_Diphoton_Photon_Up01sigma"]=events["weight_photon_presel_sf_Diphoton_Photon_up"]/events["weight_photon_presel_sf_Diphoton_Photon_central"]
    events["photon_presel_sf_Diphoton_Photon_Down01sigma"]=events["weight_photon_presel_sf_Diphoton_Photon_down"]/events["weight_photon_presel_sf_Diphoton_Photon_central"]        
    events["trigger_sf_Up01sigma"]=events["weight_trigger_sf_up"]/events["weight_trigger_sf_central"]
    events["trigger_sf_Down01sigma"]=events["weight_trigger_sf_down"]/events["weight_trigger_sf_central"]
    if "weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_central" in events.fields:
        events["PNetWvsQCDW1_sf_Up01sigma"]=events["weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_up"]/events["weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_central"]
        events["PNetWvsQCDW1_sf_Down01sigma"]=events["weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_down"]/events["weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_central"]
        events["PNetWvsQCDW1_mistagging_sf_Up01sigma"]=events["weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_up"]/events["weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_central"]
        events["PNetWvsQCDW1_mistagging_sf_Down01sigma"]=events["weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_down"]/events["weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_central"]
    else:
        events["PNetWvsQCDW1_sf_Up01sigma"]=ak.ones_like(events.weight_central)
        events["PNetWvsQCDW1_sf_Down01sigma"]=ak.ones_like(events.weight_central)
        events["PNetWvsQCDW1_mistagging_sf_Up01sigma"]=ak.ones_like(events.weight_central)
        events["PNetWvsQCDW1_mistagging_sf_Down01sigma"]=ak.ones_like(events.weight_central)
    events["isoelectron_reco_sf_Up01sigma"]=events["weight_electron_reco_sf_SelectedElectron_iso_up"]/events["weight_electron_reco_sf_SelectedElectron_iso_central"]
    events["isoelectron_reco_sf_Down01sigma"]=events["weight_electron_reco_sf_SelectedElectron_iso_down"]/events["weight_electron_reco_sf_SelectedElectron_iso_central"]
    events["nonisoelectron_reco_sf_Up01sigma"]=events["weight_electron_reco_sf_SelectedElectron_noiso_up"]/events["weight_electron_reco_sf_SelectedElectron_noiso_central"]
    events["nonisoelectron_reco_sf_Down01sigma"]=events["weight_electron_reco_sf_SelectedElectron_noiso_down"]/events["weight_electron_reco_sf_SelectedElectron_noiso_central"]
    # events["btag_reshape_jes_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_jes"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_jes_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_jes"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_lf_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_lf"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_lf_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_lf"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_hfstats1_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_hfstats1"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_hfstats1_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_hfstats1"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_hfstats2_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_hfstats2"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_hfstats2_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_hfstats2"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_cferr1_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_cferr1"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_cferr1_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_cferr1"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_cferr2_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_cferr2"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_cferr2_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_cferr2"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_hf_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_hf"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_hf_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_hf"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_lfstats1_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_lfstats1"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_lfstats1_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_lfstats1"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_lfstats2_sf_Up01sigma"]=events["weight_btag_reshape_sf_SelectedJet_up_lfstats2"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    # events["btag_reshape_lfstats2_sf_Down01sigma"]=events["weight_btag_reshape_sf_SelectedJet_down_lfstats2"]/events["weight_btag_reshape_sf_SelectedJet_central"]
    events["PTransformerHtagger_sf_Up01sigma"]=events["weight_PTransformer_up"]/events["weight_PTransformer_central"]
    events["PTransformerHtagger_sf_Down01sigma"]=events["weight_PTransformer_down"]/events["weight_PTransformer_central"]
    if "weight_PNbb_veto_sf_SelectedFatJet_up" in events.fields:
        events["PNXbb_sf_Up01sigma"]=events["weight_PNbb_veto_sf_SelectedFatJet_up"]/events["weight_PNbb_veto_sf_SelectedFatJet_central"]
        events["PNXbb_sf_Down01sigma"]=events["weight_PNbb_veto_sf_SelectedFatJet_down"]/events["weight_PNbb_veto_sf_SelectedFatJet_central"]
        # events=events[['dZ','PNN_score','category','CMS_hgg_mass','weight','muon_highptreco_sf_Down01sigma','muon_highptreco_sf_Up01sigma','muon_highptid_sf_Down01sigma','muon_highptid_sf_Up01sigma','L1_prefiring_sf_Down01sigma','L1_prefiring_sf_Up01sigma','puWeight_Up01sigma','puWeight_Down01sigma','jet_pu_id_sf_Up01sigma','jet_pu_id_sf_Down01sigma','electron_veto_sf_Diphoton_Photon_Up01sigma','electron_veto_sf_Diphoton_Photon_Down01sigma','isoelectron_id_sf_SelectedElectron_iso_Up01sigma','isoelectron_id_sf_SelectedElectron_iso_Down01sigma','isoelectron_id_sf_SelectedElectron_noiso_Up01sigma','isoelectron_id_sf_SelectedElectron_noiso_Down01sigma','isomuon_id_sf_SelectedMuon_iso_Up01sigma','isomuon_id_sf_SelectedMuon_iso_Down01sigma','isomuon_iso_sf_SelectedMuon_iso_Up01sigma','isomuon_iso_sf_SelectedMuon_iso_Down01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma','photon_id_sf_Diphoton_Photon_Up01sigma','photon_id_sf_Diphoton_Photon_Down01sigma','photon_presel_sf_Diphoton_Photon_Up01sigma','photon_presel_sf_Diphoton_Photon_Down01sigma','trigger_sf_Up01sigma','trigger_sf_Down01sigma','PNetWvsQCDW1_sf_Up01sigma','PNetWvsQCDW1_sf_Down01sigma','PNetWvsQCDW1_mistagging_sf_Up01sigma','PNetWvsQCDW1_mistagging_sf_Down01sigma','isoelectron_reco_sf_Up01sigma','isoelectron_reco_sf_Down01sigma','nonisoelectron_reco_sf_Up01sigma','nonisoelectron_reco_sf_Down01sigma','btag_reshape_jes_sf_Up01sigma','btag_reshape_jes_sf_Down01sigma','btag_reshape_lf_sf_Up01sigma','btag_reshape_lf_sf_Down01sigma','btag_reshape_hfstats1_sf_Up01sigma','btag_reshape_hfstats1_sf_Down01sigma','btag_reshape_hfstats2_sf_Up01sigma','btag_reshape_hfstats2_sf_Down01sigma','btag_reshape_cferr1_sf_Up01sigma','btag_reshape_cferr1_sf_Down01sigma','btag_reshape_cferr2_sf_Up01sigma','btag_reshape_cferr2_sf_Down01sigma','btag_reshape_hf_sf_Up01sigma','btag_reshape_hf_sf_Down01sigma','btag_reshape_lfstats1_sf_Up01sigma','btag_reshape_lfstats1_sf_Down01sigma','btag_reshape_lfstats2_sf_Up01sigma','btag_reshape_lfstats2_sf_Down01sigma','PNXbb_sf_Up01sigma','PNXbb_sf_Down01sigma',"PTransformerHtagger_sf_Up01sigma","PTransformerHtagger_sf_Down01sigma"]]
        events=events[['dZ','PNN_score','category','CMS_hgg_mass','weight','muon_highptreco_sf_Down01sigma','muon_highptreco_sf_Up01sigma','muon_highptid_sf_Down01sigma','muon_highptid_sf_Up01sigma','L1_prefiring_sf_Down01sigma','L1_prefiring_sf_Up01sigma','puWeight_Up01sigma','puWeight_Down01sigma','jet_pu_id_sf_Up01sigma','jet_pu_id_sf_Down01sigma','electron_veto_sf_Diphoton_Photon_Up01sigma','electron_veto_sf_Diphoton_Photon_Down01sigma','isoelectron_id_sf_SelectedElectron_iso_Up01sigma','isoelectron_id_sf_SelectedElectron_iso_Down01sigma','isoelectron_id_sf_SelectedElectron_noiso_Up01sigma','isoelectron_id_sf_SelectedElectron_noiso_Down01sigma','isomuon_id_sf_SelectedMuon_iso_Up01sigma','isomuon_id_sf_SelectedMuon_iso_Down01sigma','isomuon_iso_sf_SelectedMuon_iso_Up01sigma','isomuon_iso_sf_SelectedMuon_iso_Down01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma','photon_id_sf_Diphoton_Photon_Up01sigma','photon_id_sf_Diphoton_Photon_Down01sigma','photon_presel_sf_Diphoton_Photon_Up01sigma','photon_presel_sf_Diphoton_Photon_Down01sigma','trigger_sf_Up01sigma','trigger_sf_Down01sigma','PNetWvsQCDW1_sf_Up01sigma','PNetWvsQCDW1_sf_Down01sigma','PNetWvsQCDW1_mistagging_sf_Up01sigma','PNetWvsQCDW1_mistagging_sf_Down01sigma','isoelectron_reco_sf_Up01sigma','isoelectron_reco_sf_Down01sigma','nonisoelectron_reco_sf_Up01sigma','nonisoelectron_reco_sf_Down01sigma','PNXbb_sf_Up01sigma','PNXbb_sf_Down01sigma',"PTransformerHtagger_sf_Up01sigma","PTransformerHtagger_sf_Down01sigma"]]
   
    else:
        # events=events[['dZ','PNN_score','category','CMS_hgg_mass','weight','muon_highptreco_sf_Down01sigma','muon_highptreco_sf_Up01sigma','muon_highptid_sf_Down01sigma','muon_highptid_sf_Up01sigma','L1_prefiring_sf_Down01sigma','L1_prefiring_sf_Up01sigma','puWeight_Up01sigma','puWeight_Down01sigma','jet_pu_id_sf_Up01sigma','jet_pu_id_sf_Down01sigma','electron_veto_sf_Diphoton_Photon_Up01sigma','electron_veto_sf_Diphoton_Photon_Down01sigma','isoelectron_id_sf_SelectedElectron_iso_Up01sigma','isoelectron_id_sf_SelectedElectron_iso_Down01sigma','isoelectron_id_sf_SelectedElectron_noiso_Up01sigma','isoelectron_id_sf_SelectedElectron_noiso_Down01sigma','isomuon_id_sf_SelectedMuon_iso_Up01sigma','isomuon_id_sf_SelectedMuon_iso_Down01sigma','isomuon_iso_sf_SelectedMuon_iso_Up01sigma','isomuon_iso_sf_SelectedMuon_iso_Down01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma','photon_id_sf_Diphoton_Photon_Up01sigma','photon_id_sf_Diphoton_Photon_Down01sigma','photon_presel_sf_Diphoton_Photon_Up01sigma','photon_presel_sf_Diphoton_Photon_Down01sigma','trigger_sf_Up01sigma','trigger_sf_Down01sigma','PNetWvsQCDW1_sf_Up01sigma','PNetWvsQCDW1_sf_Down01sigma','PNetWvsQCDW1_mistagging_sf_Up01sigma','PNetWvsQCDW1_mistagging_sf_Down01sigma','isoelectron_reco_sf_Up01sigma','isoelectron_reco_sf_Down01sigma','nonisoelectron_reco_sf_Up01sigma','nonisoelectron_reco_sf_Down01sigma','btag_reshape_jes_sf_Up01sigma','btag_reshape_jes_sf_Down01sigma','btag_reshape_lf_sf_Up01sigma','btag_reshape_lf_sf_Down01sigma','btag_reshape_hfstats1_sf_Up01sigma','btag_reshape_hfstats1_sf_Down01sigma','btag_reshape_hfstats2_sf_Up01sigma','btag_reshape_hfstats2_sf_Down01sigma','btag_reshape_cferr1_sf_Up01sigma','btag_reshape_cferr1_sf_Down01sigma','btag_reshape_cferr2_sf_Up01sigma','btag_reshape_cferr2_sf_Down01sigma','btag_reshape_hf_sf_Up01sigma','btag_reshape_hf_sf_Down01sigma','btag_reshape_lfstats1_sf_Up01sigma','btag_reshape_lfstats1_sf_Down01sigma','btag_reshape_lfstats2_sf_Up01sigma','btag_reshape_lfstats2_sf_Down01sigma',"PTransformerHtagger_sf_Up01sigma","PTransformerHtagger_sf_Down01sigma"]]
        events=events[['dZ','PNN_score','category','CMS_hgg_mass','weight','muon_highptreco_sf_Down01sigma','muon_highptreco_sf_Up01sigma','muon_highptid_sf_Down01sigma','muon_highptid_sf_Up01sigma','L1_prefiring_sf_Down01sigma','L1_prefiring_sf_Up01sigma','puWeight_Up01sigma','puWeight_Down01sigma','jet_pu_id_sf_Up01sigma','jet_pu_id_sf_Down01sigma','electron_veto_sf_Diphoton_Photon_Up01sigma','electron_veto_sf_Diphoton_Photon_Down01sigma','isoelectron_id_sf_SelectedElectron_iso_Up01sigma','isoelectron_id_sf_SelectedElectron_iso_Down01sigma','isoelectron_id_sf_SelectedElectron_noiso_Up01sigma','isoelectron_id_sf_SelectedElectron_noiso_Down01sigma','isomuon_id_sf_SelectedMuon_iso_Up01sigma','isomuon_id_sf_SelectedMuon_iso_Down01sigma','isomuon_iso_sf_SelectedMuon_iso_Up01sigma','isomuon_iso_sf_SelectedMuon_iso_Down01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma','photon_id_sf_Diphoton_Photon_Up01sigma','photon_id_sf_Diphoton_Photon_Down01sigma','photon_presel_sf_Diphoton_Photon_Up01sigma','photon_presel_sf_Diphoton_Photon_Down01sigma','trigger_sf_Up01sigma','trigger_sf_Down01sigma','PNetWvsQCDW1_sf_Up01sigma','PNetWvsQCDW1_sf_Down01sigma','PNetWvsQCDW1_mistagging_sf_Up01sigma','PNetWvsQCDW1_mistagging_sf_Down01sigma','isoelectron_reco_sf_Up01sigma','isoelectron_reco_sf_Down01sigma','nonisoelectron_reco_sf_Up01sigma','nonisoelectron_reco_sf_Down01sigma',"PTransformerHtagger_sf_Up01sigma","PTransformerHtagger_sf_Down01sigma"]]

    return events
events_sig = add_sf_branches(events_sig)
events_bbgg = add_sf_branches(events_bbgg)
events_zzgg = add_sf_branches(events_zzgg)
events_ttgg = add_sf_branches(events_ttgg)
events_vbf = add_sf_branches(events_vbf)
events_ggh = add_sf_branches(events_ggh)
events_vh = add_sf_branches(events_vh)
events_tth = add_sf_branches(events_tth)


events_data["CMS_hgg_mass"]=events_data["Diphoton_mass"]
events_data['weight']=events_data.weight_central
events_data=events_data[['CMS_hgg_mass','weight','PNN_score']]
events_data_highpurity = events_data[(events_data['PNN_score'] > cut1) & (events_data['PNN_score'] <= 1)]
events_data_lowpurity = events_data[(events_data['PNN_score'] > cut2) & (events_data['PNN_score'] <= cut1)]
massname="MX"+(signal_samples['sig_output_name'][0].split("MX"))[1].split("_cat12")[0]


dataA_rootname="Data_"+year+"_cat12highpurity_"+massname+".root"
dataA_treename="Data_13TeV_cat12highpurity"
dataB_rootname="Data_"+year+"_cat12lowpurity_"+massname+".root"
dataB_treename="Data_13TeV_cat12lowpurity"

data_Acat_output_path="./PBDT_HH_FHSL_combine_"+year+"/"+dataA_rootname.replace(".root",".parquet")
data_Bcat_output_path="./PBDT_HH_FHSL_combine_"+year+"/"+dataB_rootname.replace(".root",".parquet")
ak.to_parquet(events_data_highpurity, data_Acat_output_path)
ak.to_parquet(events_data_lowpurity, data_Bcat_output_path)
del events_data
events_sig_highpurity = events_sig[(events_sig['PNN_score'] > cut1) & (events_sig['PNN_score'] <= 1)]
events_sig_lowpurity = events_sig[(events_sig['PNN_score'] > cut2) & (events_sig['PNN_score'] <= cut1)]
ak.to_parquet(events_sig_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['sig_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_sig_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['sig_output_name'][0]+"_lowpurity.parquet")
events_bbgg_highpurity = events_bbgg[(events_bbgg['PNN_score'] > cut1) & (events_bbgg['PNN_score'] <= 1)]
events_bbgg_lowpurity = events_bbgg[(events_bbgg['PNN_score'] > cut2) & (events_bbgg['PNN_score'] <= cut1)]
ak.to_parquet(events_bbgg_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['bbgg_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_bbgg_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['bbgg_output_name'][0]+"_lowpurity.parquet")
events_zzgg_highpurity = events_zzgg[(events_zzgg['PNN_score'] > cut1) & (events_zzgg['PNN_score'] <= 1)]
events_zzgg_lowpurity = events_zzgg[(events_zzgg['PNN_score'] > cut2) & (events_zzgg['PNN_score'] <= cut1)]
ak.to_parquet(events_zzgg_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['zzgg_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_zzgg_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['zzgg_output_name'][0]+"_lowpurity.parquet")
events_ttgg_highpurity = events_ttgg[(events_ttgg['PNN_score'] > cut1) & (events_ttgg['PNN_score'] <= 1)]
events_ttgg_lowpurity = events_ttgg[(events_ttgg['PNN_score'] > cut2) & (events_ttgg['PNN_score'] <= cut1)]
ak.to_parquet(events_ttgg_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['ttgg_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_ttgg_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['ttgg_output_name'][0]+"_lowpurity.parquet")
events_vbf_highpurity = events_vbf[(events_vbf['PNN_score'] > cut1) & (events_vbf['PNN_score'] <= 1)]
events_vbf_lowpurity = events_vbf[(events_vbf['PNN_score'] > cut2) & (events_vbf['PNN_score'] <= cut1)]
ak.to_parquet(events_vbf_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['vbf_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_vbf_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['vbf_output_name'][0]+"_lowpurity.parquet")
events_vh_highpurity = events_vh[(events_vh['PNN_score'] > cut1) & (events_vh['PNN_score'] <= 1)]
events_vh_lowpurity = events_vh[(events_vh['PNN_score'] > cut2) & (events_vh['PNN_score'] <= cut1)]
ak.to_parquet(events_vh_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['vh_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_vh_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['vh_output_name'][0]+"_lowpurity.parquet")
events_tth_highpurity = events_tth[(events_tth['PNN_score'] > cut1) & (events_tth['PNN_score'] <= 1)]
events_tth_lowpurity = events_tth[(events_tth['PNN_score'] > cut2) & (events_tth['PNN_score'] <= cut1)]
ak.to_parquet(events_tth_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['tth_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_tth_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['tth_output_name'][0]+"_lowpurity.parquet")
events_ggh_highpurity = events_ggh[(events_ggh['PNN_score'] > cut1) & (events_ggh['PNN_score'] <= 1)]
events_ggh_lowpurity = events_ggh[(events_ggh['PNN_score'] > cut2) & (events_ggh['PNN_score'] <= cut1)]
ak.to_parquet(events_ggh_highpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['ggh_output_name'][0]+"_highpurity.parquet")
ak.to_parquet(events_ggh_lowpurity,"./PBDT_HH_FHSL_combine_"+year+"/"+signal_samples['ggh_output_name'][0]+"_lowpurity.parquet")


bbgg_highpurity=ak.sum(events_bbgg_highpurity['weight'])
bbgg_lowpurity=ak.sum(events_bbgg_lowpurity['weight'])
zzgg_highpurity=ak.sum(events_zzgg_highpurity['weight'])
zzgg_lowpurity=ak.sum(events_zzgg_lowpurity['weight'])
ttgg_highpurity=ak.sum(events_ttgg_highpurity['weight'])
ttgg_lowpurity=ak.sum(events_ttgg_lowpurity['weight'])
FHSL_highpurity=ak.sum(events_sig_highpurity['weight'])
FHSL_lowpurity=ak.sum(events_sig_lowpurity['weight'])
vbf_highpurity=ak.sum(events_vbf_highpurity['weight'])
vbf_lowpurity=ak.sum(events_vbf_lowpurity['weight'])
vh_highpurity=ak.sum(events_vh_highpurity['weight'])
vh_lowpurity=ak.sum(events_vh_lowpurity['weight'])
tth_highpurity=ak.sum(events_tth_highpurity['weight'])
tth_lowpurity=ak.sum(events_tth_lowpurity['weight'])
ggh_highpurity=ak.sum(events_ggh_highpurity['weight'])
ggh_lowpurity=ak.sum(events_ggh_lowpurity['weight'])
del events_sig
del events_bbgg
del events_zzgg
del events_ttgg



print('start to process all the samples')
def process_sig_samples(FHfile,SLfile,bbggfile,zzggfile,ttggfile,vbffile,vhfile,tthfile,gghfile,sigoutput,bbggoutput,zzggoutput,ttggoutput,vbfoutput,vhoutput,tthoutput,gghoutput,input_features,cut1,cut2):
    def get_jet_variables(events):
        events['fatjet_1_pt'] = events['fatjet_1_pt']
        events['fatjet_2_pt'] = events['fatjet_2_pt']
        events['fatjet_1_eta'] = events['fatjet_1_eta']
        events['fatjet_2_eta'] = events['fatjet_2_eta']
        events['fatjet_1_phi'] = events['fatjet_1_phi']
        events['fatjet_2_phi'] = events['fatjet_2_phi']
        events['fatjet_1_msoftdrop'] = events['fatjet_1_msoftdrop']
        events['fatjet_2_msoftdrop'] = events['fatjet_2_msoftdrop']
        # get fatjet and diphoton dR
        fatjet_1_4D = vector.obj(pt=events.fatjet_1_pt, eta=events.fatjet_1_eta, phi=events.fatjet_1_phi, mass=events.fatjet_1_msoftdrop)
        fatjet_2_4D = vector.obj(pt=events.fatjet_2_pt, eta=events.fatjet_2_eta, phi=events.fatjet_2_phi, mass=events.fatjet_2_msoftdrop)
        fatjet_3_4D = vector.obj(pt=events.fatjet_3_pt, eta=events.fatjet_3_eta, phi=events.fatjet_3_phi, mass=events.fatjet_3_msoftdrop)
        diphoton_4D = vector.obj(pt=events.Diphoton_pt, eta=events.Diphoton_eta, phi=events.Diphoton_phi, mass=events.Diphoton_mass)
        fatjet_1_diphoton_dR = fatjet_1_4D.deltaR(diphoton_4D)
        fatjet_2_diphoton_dR = fatjet_2_4D.deltaR(diphoton_4D)
        events['fatjet_1_diphoton_dR'] = np.where(events.fatjet_1_pt>0, fatjet_1_diphoton_dR, -999)
        events['fatjet_2_diphoton_dR'] = np.where(events.fatjet_2_pt>0, fatjet_2_diphoton_dR, -999)
        leadphoton=vector.obj(pt=events.LeadPhoton_pt,eta=events.LeadPhoton_eta,phi=events.LeadPhoton_phi,mass=events.LeadPhoton_mass)
        subleadphoton=vector.obj(pt=events.SubleadPhoton_pt,eta=events.SubleadPhoton_eta,phi=events.SubleadPhoton_phi,mass=events.SubleadPhoton_mass)
        # dR with photon and fatjet
        fatjet_1_leadphoton_dR = fatjet_1_4D.deltaR(leadphoton)
        events['fatjet_1_leadphoton_dR'] = np.where(events.fatjet_1_pt>0, fatjet_1_leadphoton_dR, -999)
        fatjet_1_subleadphoton_dR = fatjet_1_4D.deltaR(subleadphoton)
        events['fatjet_1_subleadphoton_dR'] = np.where(events.fatjet_1_pt>0, fatjet_1_subleadphoton_dR, -999)
        fatjet_2_leadphoton_dR = fatjet_2_4D.deltaR(leadphoton)
        events['fatjet_2_leadphoton_dR'] = np.where(events.fatjet_2_pt>0, fatjet_2_leadphoton_dR, -999)
        fatjet_2_subleadphoton_dR = fatjet_2_4D.deltaR(subleadphoton)
        events['fatjet_2_subleadphoton_dR'] = np.where(events.fatjet_2_pt>0, fatjet_2_subleadphoton_dR, -999)
        # get fatjet1 and fatjet2 dR
        fatjet_1_2_dR = fatjet_1_4D.deltaR(fatjet_2_4D)
        events['fatjet_1_2_dR'] = np.where(np.logical_and(events.fatjet_1_pt>0, events.fatjet_2_pt>0), fatjet_1_2_dR, -999)
        # get the maximum fatjets mass with the combination of 2 fatjets in 3 fatjets
        fatjet_12_4D = fatjet_1_4D+fatjet_2_4D
        fatjet_13_4D = fatjet_1_4D+fatjet_3_4D
        fatjet_23_4D = fatjet_2_4D+fatjet_3_4D
        fatjet_12_msoftdrop = np.where(((fatjet_1_4D.pt >0) & (fatjet_2_4D.pt >0)), fatjet_12_4D.mass, -999)
        fatjet_13_msoftdrop = np.where(((fatjet_1_4D.pt >0) & (fatjet_3_4D.pt >0)), fatjet_13_4D.mass, -999)
        fatjet_23_msoftdrop = np.where(((fatjet_2_4D.pt >0) & (fatjet_3_4D.pt >0)), fatjet_23_4D.mass, -999)
        max_fatjets_mass = np.maximum(fatjet_12_msoftdrop, np.maximum(fatjet_13_msoftdrop, fatjet_23_msoftdrop))
        events['max_fatjets_mass'] = max_fatjets_mass
        # get max WvsQCD score with three fatjets
        events['fatjet_1_WvsQCDMD'] = events['fatjet_1_WvsQCDMD']
        events['fatjet_2_WvsQCDMD'] = events['fatjet_2_WvsQCDMD']
        # get max H4qvsQCD score with three fatjets
        events['fatjet_1_Hqqqq_vsQCDTop'] = events['fatjet_1_Hqqqq_vsQCDTop']
        events['fatjet_2_Hqqqq_vsQCDTop'] = events['fatjet_2_Hqqqq_vsQCDTop']
        # get XbbvsQCD
        events['fatjet_1_XbbvsQCDMD'] = (events['fatjet_1_particleNetMD_Xbb']) / (events['fatjet_1_particleNetMD_Xbb'] + events['fatjet_1_particleNetMD_QCD'])
        events['fatjet_2_XbbvsQCDMD'] = (events['fatjet_2_particleNetMD_Xbb']) / (events['fatjet_2_particleNetMD_Xbb'] + events['fatjet_2_particleNetMD_QCD'])
        Hscore = ak.concatenate([ak.unflatten(events['fatjet_1_Hqqqq_vsQCDTop'], counts=1), ak.unflatten(events['fatjet_2_Hqqqq_vsQCDTop'], counts=1)], axis=1)
        events['max_fatjet_Hqqqq_vsQCDTop'] = Hscore[ak.argsort(Hscore, axis=-1, ascending=False)][:,0]
        events['nGoodAK4jets'] = events['nGoodAK4jets']
        # add number of good AK8 jets
        events['nGoodAK8jets'] = events['nGoodAK8jets']
        # get 4 ak4 jets 4D info
        events['jet_1_pt'] = events['jet_1_pt']
        events['jet_2_pt'] = events['jet_2_pt']
        events['jet_3_pt'] = events['jet_3_pt']
        events['jet_1_eta'] = events['jet_1_eta']
        events['jet_2_eta'] = events['jet_2_eta']
        events['jet_3_eta'] = events['jet_3_eta']
        events['jet_1_phi'] = events['jet_1_phi']
        events['jet_2_phi'] = events['jet_2_phi']
        events['jet_3_phi'] = events['jet_3_phi']
        events['jet_1_mass'] = events['jet_1_mass']
        events['jet_2_mass'] = events['jet_2_mass']
        events['jet_3_mass'] = events['jet_3_mass']
        return events    
    def get_sig_events_forApply(filename,mx,my):
        events = ak.from_parquet(filename)
        category_cut = ((events["category"]==1) | (events["category"]==2))
        photonID_cut = (events["Diphoton_minID_modified"]>-0.7)
        events = events[ category_cut & photonID_cut]
        events = get_jet_variables(events)
        events['mx'] = np.ones(len(events))*int(mx)
        events['my'] = np.ones(len(events))*int(my)
        return events
    def add_shape_uncertainty_br(events):
        events["CMS_hgg_mass"] = events["Diphoton_mass"]
        events["weight"] = events["weight_central"]
        events["dZ"] = np.ones(len(events['CMS_hgg_mass']))
        events = events[['CMS_hgg_mass', 'weight', 'dZ', 'PNN_score','category']]
        return events
    def model_predict(event, model, loaded_scaler, input_features):
        df = ak.to_pandas(event[input_features + ['weight_central','Diphoton_mass']])
        df = df.replace(-999, 0)
        X_test = loaded_scaler.transform(df[input_features])
        X_test = torch.tensor(X_test).float()
        proba = model(X_test)
        proba = F.softmax(proba, dim=1)
        proba = proba.detach().numpy()
        pnn_score = (proba[:,4] + proba[:,2] + proba[:,3]) / (proba[:,0] + proba[:,1] + proba[:,2] + proba[:,3] + proba[:,4])

        del proba
        del df
        del X_test
        return pnn_score
    mx = sigoutput.split('_')[1].split('X')[1]
    my = sigoutput.split('_')[2].split('H')[1]
    events_sigFH = get_sig_events_forApply(FHfile, mx, my)
    events_sigSL = get_sig_events_forApply(SLfile, mx, my)
    events_bbgg = get_sig_events_forApply(bbggfile, mx, my)
    events_zzgg = get_sig_events_forApply(zzggfile, mx, my)
    events_ttgg = get_sig_events_forApply(ttggfile, mx, my)
    events_vbf = get_sig_events_forApply(vbffile, mx, my)
    events_vh = get_sig_events_forApply(vhfile, mx, my)
    events_tth = get_sig_events_forApply(tthfile, mx, my)
    events_ggh = get_sig_events_forApply(gghfile, mx, my)
    PNN_score = model_predict(events_vbf, model, loaded_scaler, input_features)
    events_vbf['PNN_score'] = PNN_score
    PNN_score = model_predict(events_vh, model, loaded_scaler, input_features)
    events_vh['PNN_score'] = PNN_score
    PNN_score = model_predict(events_tth, model, loaded_scaler, input_features)
    events_tth['PNN_score'] = PNN_score
    PNN_score = model_predict(events_ggh, model, loaded_scaler, input_features)
    events_ggh['PNN_score'] = PNN_score
    PNN_score = model_predict(events_sigFH, model, loaded_scaler, input_features)
    events_sigFH['PNN_score'] = PNN_score
    del PNN_score
    PNN_score = model_predict(events_sigSL, model, loaded_scaler, input_features)
    events_sigSL['PNN_score'] = PNN_score
    del PNN_score
    events_sig = ak.concatenate([events_sigFH, events_sigSL])
    del events_sigFH
    del events_sigSL
    PNN_score= model_predict(events_bbgg, model, loaded_scaler, input_features)
    events_bbgg['PNN_score'] = PNN_score
    PNN_score= model_predict(events_zzgg, model, loaded_scaler, input_features)
    events_zzgg['PNN_score'] = PNN_score
    PNN_score= model_predict(events_ttgg, model, loaded_scaler, input_features)
    events_ttgg['PNN_score'] = PNN_score
    events_sig=add_shape_uncertainty_br(events_sig)
    events_bbgg=add_shape_uncertainty_br(events_bbgg)
    events_zzgg=add_shape_uncertainty_br(events_zzgg)
    events_ttgg=add_shape_uncertainty_br(events_ttgg)
    events_vbf=add_shape_uncertainty_br(events_vbf)
    events_vh=add_shape_uncertainty_br(events_vh)
    events_tth=add_shape_uncertainty_br(events_tth)
    events_ggh=add_shape_uncertainty_br(events_ggh)
    events_sig_highpurity = events_sig[(events_sig['PNN_score'] > cut1) & (events_sig['PNN_score'] <= 1)]
    events_sig_lowpurity = events_sig[(events_sig['PNN_score'] > cut2) & (events_sig['PNN_score'] <= cut1)]
    ak.to_parquet(events_sig_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + sigoutput + "_highpurity.parquet")
    ak.to_parquet(events_sig_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + sigoutput + "_lowpurity.parquet")
    events_bbgg_highpurity = events_bbgg[(events_bbgg['PNN_score'] > cut1) & (events_bbgg['PNN_score'] <= 1)]
    events_bbgg_lowpurity = events_bbgg[(events_bbgg['PNN_score'] > cut2) & (events_bbgg['PNN_score'] <= cut1)]
    ak.to_parquet(events_bbgg_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + bbggoutput + "_highpurity.parquet")
    ak.to_parquet(events_bbgg_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + bbggoutput + "_lowpurity.parquet")
    events_zzgg_highpurity = events_zzgg[(events_zzgg['PNN_score'] > cut1) & (events_zzgg['PNN_score'] <= 1)]
    events_zzgg_lowpurity = events_zzgg[(events_zzgg['PNN_score'] > cut2) & (events_zzgg['PNN_score'] <= cut1)]
    ak.to_parquet(events_zzgg_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + zzggoutput + "_highpurity.parquet")
    ak.to_parquet(events_zzgg_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + zzggoutput + "_lowpurity.parquet")
    events_ttgg_highpurity = events_ttgg[(events_ttgg['PNN_score'] > cut1) & (events_ttgg['PNN_score'] <= 1)]
    events_ttgg_lowpurity = events_ttgg[(events_ttgg['PNN_score'] > cut2) & (events_ttgg['PNN_score'] <= cut1)]
    ak.to_parquet(events_ttgg_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + ttggoutput + "_highpurity.parquet")
    ak.to_parquet(events_ttgg_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + ttggoutput + "_lowpurity.parquet")
    events_vbf_highpurity = events_vbf[(events_vbf['PNN_score'] > cut1) & (events_vbf['PNN_score'] <= 1)]
    events_vbf_lowpurity = events_vbf[(events_vbf['PNN_score'] > cut2) & (events_vbf['PNN_score'] <= cut1)]
    ak.to_parquet(events_vbf_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + vbfoutput + "_highpurity.parquet")
    ak.to_parquet(events_vbf_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + vbfoutput + "_lowpurity.parquet")
    events_vh_highpurity = events_vh[(events_vh['PNN_score'] > cut1) & (events_vh['PNN_score'] <= 1)]
    events_vh_lowpurity = events_vh[(events_vh['PNN_score'] > cut2) & (events_vh['PNN_score'] <= cut1)]
    ak.to_parquet(events_vh_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + vhoutput + "_highpurity.parquet")
    ak.to_parquet(events_vh_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + vhoutput + "_lowpurity.parquet")
    events_tth_highpurity = events_tth[(events_tth['PNN_score'] > cut1) & (events_tth['PNN_score'] <= 1)]
    events_tth_lowpurity = events_tth[(events_tth['PNN_score'] > cut2) & (events_tth['PNN_score'] <= cut1)]
    ak.to_parquet(events_tth_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + tthoutput + "_highpurity.parquet")
    ak.to_parquet(events_tth_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + tthoutput + "_lowpurity.parquet")
    events_ggh_highpurity = events_ggh[(events_ggh['PNN_score'] > cut1) & (events_ggh['PNN_score'] <= 1)]
    events_ggh_lowpurity = events_ggh[(events_ggh['PNN_score'] > cut2) & (events_ggh['PNN_score'] <= cut1)]
    ak.to_parquet(events_ggh_highpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + gghoutput + "_highpurity.parquet")
    ak.to_parquet(events_ggh_lowpurity, "./PBDT_HH_FHSL_combine_"+year+"/" + gghoutput + "_lowpurity.parquet")
    highpuritylen=len(events_sig_highpurity)

    return highpuritylen
for i in range(1, len(signal_samples['FHpath'])):
    result=process_sig_samples(signal_samples['FHpath'][i], signal_samples['SLpath'][i], signal_samples['BBGGpath'][i], signal_samples['ZZggpath'][i],signal_samples['TTggpath'][i], signal_samples['VBFpath'][i],signal_samples['VHpath'][i],signal_samples['ttHpath'][i],signal_samples['ggHpath'][i],signal_samples['sig_output_name'][i], signal_samples['bbgg_output_name'][i], signal_samples['zzgg_output_name'][i],signal_samples['ttgg_output_name'][i],signal_samples['vbf_output_name'][i],signal_samples['vh_output_name'][i],signal_samples['tth_output_name'][i],signal_samples['ggh_output_name'][i], input_features, cut1, cut2)
def process_highpurity_sigfile(file):
    signal = file.split("/")[-1].split("_highpurity")[0].split("merge")[0]
    signaltype = "wwgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhwwgg_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhwwgg_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhwwgg_125_13TeV_cat12highpurity"
    rootname="./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"+file.split("/")[-1].replace("parquet", "root")
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_highpurity_VBFfile(file):
    signal = file.split("/")[-1].split("_highpurity")[0].split("merge")[0]
    signaltype = "VBF"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhVBF_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhVBF_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhVBF_125_13TeV_cat12highpurity"
    rootname="./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"+file.split("/")[-1].replace("parquet", "root")
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")
    return tree_name
def process_highpurity_TTHfile(file):
    signal = file.split("/")[-1].split("_highpurity")[0].split("merge")[0]
    signaltype = "TTH"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhTTH_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhTTH_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhTTH_125_13TeV_cat12highpurity"
    rootname="./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"+file.split("/")[-1].replace("parquet", "root")
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")
    return tree_name
def process_highpurity_VHfile(file):
    signal = file.split("/")[-1].split("_highpurity")[0].split("merge")[0]
    signaltype = "VH"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhVH_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhVH_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhVH_125_13TeV_cat12highpurity"
    rootname="./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"+file.split("/")[-1].replace("parquet", "root")
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")
    return tree_name
def process_highpurity_GGHfile(file):
    signal = file.split("/")[-1].split("_highpurity")[0].split("merge")[0]
    signaltype = "GGH"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhGGH_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhGGH_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhGGH_125_13TeV_cat12highpurity"
    rootname="./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"+file.split("/")[-1].replace("parquet", "root")
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")
    return tree_name
def process_highpurity_zzggfile(file):
    signaltype = "zzgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhzzgg_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhzzgg_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhzzgg_125_13TeV_cat12highpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_highpurity_ttggfile(file):
    signaltype = "ttgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhttgg_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhttgg_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhttgg_125_13TeV_cat12highpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_highpurity_bbggfile(file):
    signaltype = "bbgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhbbgg_125_13TeV_cat12highpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhbbgg_125_13TeV_cat12highpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhbbgg_125_13TeV_cat12highpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_lowpurity_sigfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "wwgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhwwgg_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhwwgg_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhwwgg_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_lowpurity_VBFfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "VBF"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhVBF_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhVBF_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhVBF_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_lowpurity_VHfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "VH"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhVH_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhVH_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhVH_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_lowpurity_TTHfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "TTH"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhTTH_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhTTH_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhTTH_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_lowpurity_GGHfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "GGH"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]

    if "down" in file.split("/")[-1]:
        tree_name = "gghhGGH_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhGGH_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhGGH_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")

    return tree_name
def process_lowpurity_zzggfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "zzgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]
    if "down" in file.split("/")[-1]:
        tree_name = "gghhzzgg_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhzzgg_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhzzgg_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")
    return tree_name
def process_lowpurity_ttggfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "ttgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]
    if "down" in file.split("/")[-1]:
        tree_name = "gghhttgg_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhttgg_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhttgg_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")
    return tree_name
def process_lowpurity_bbggfile(file):
    signal = file.split("/")[-1].split("_lowpurity")[0].split("merge")[0]
    signaltype = "bbgg"
    name=file.split("/")[-1].split('merged_')[-1]
    if "down" in name:
        sys=name.split("_down")[0]
    elif "up" in name:
        sys=name.split("_up")[0]
    if "down" in file.split("/")[-1]:
        tree_name = "gghhbbgg_125_13TeV_cat12lowpurity_" + sys + "Down01sigma"
    elif "up" in file.split("/")[-1]:
        tree_name = "gghhbbgg_125_13TeV_cat12lowpurity_" + sys + "Up01sigma"
    elif "nominal" in file.split("/")[-1]:
        tree_name = "gghhbbgg_125_13TeV_cat12lowpurity"
    parquet_to_root(
        file,
        "./PBDT_HH_FHSL_combine_"+year+"/flashgginput/"
        + file.split("/")[-1].replace("parquet", "root"),
        treename=tree_name,
    )
    print("done")
    return tree_name
import multiprocessing

if __name__ == "__main__":
    print("starting to get root")
    Xmass=FHfile.split("M-")[1].split("_")[0]
    highpurity_sigfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/CombineFHSL_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    lowpurity_sigfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/CombineFHSL_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    lowpurity_bbggfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/BBGG_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    highpurity_bbggfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/BBGG_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    highpurity_zzggfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/ZZGG_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    lowpurity_zzggfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/ZZGG_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    highpurity_ttggfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/TTGG_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    lowpurity_ttggfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/TTGG_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    highpurity_VBFfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/VBF_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    lowpurity_VBFfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/VBF_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    highpurity_TTHfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/TTH_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    lowpurity_TTHfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/TTH_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    highpurity_VHfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/VH_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    lowpurity_VHfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/VH_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    highpurity_GGHfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/GGH_MX"+Xmass+"_MH125_cat12_m*_highpurity.parquet")
    lowpurity_GGHfiles = glob.glob("./PBDT_HH_FHSL_combine_"+year+"/GGH_MX"+Xmass+"_MH125_cat12_m*_lowpurity.parquet")
    #multiprocessing to convert parquet to root, and hadd the root files for highpurity and lowpurity WWggsiganl
    pool = multiprocessing.Pool(processes=10)
    results = []
    #highpurity signal
    for file in tqdm(highpurity_sigfiles):
        result = pool.apply_async(process_highpurity_sigfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_wwgg_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/CombineFHSL_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/CombineFHSL_MX'+Xmass+'*highpurity.root'
    os.system(command)
    pool = multiprocessing.Pool(processes=10)
    results = []
    #lowpurity signal
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_sigfiles):
        result = pool.apply_async(process_lowpurity_sigfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_wwgg_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/CombineFHSL_MX' + Xmass + '*lowpurity.root'
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/CombineFHSL_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    parquet_to_root(data_Acat_output_path,"./PBDT_HH_FHSL_combine_"+year+"/"+"flashgginput/MX"+Xmass+"_MH125/"+dataA_rootname,treename=dataA_treename,verbose=False)
    parquet_to_root(data_Bcat_output_path,"./PBDT_HH_FHSL_combine_"+year+"/"+"flashgginput/MX"+Xmass+"_MH125/"+dataB_rootname,treename=dataB_treename,verbose=False)
    #multiprocessing to convert parquet to root, and hadd the root files for highpurity and lowpurity bbggsignal
    #highpurity bbgg
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(highpurity_bbggfiles):
        result = pool.apply_async(process_highpurity_bbggfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_bbgg_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/BBGG_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/BBGG_MX'+Xmass+'*highpurity.root'
    os.system(command)
    #lowpurity bbgg
    #lowpurity signal
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_bbggfiles):
        result = pool.apply_async(process_lowpurity_bbggfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_bbgg_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/BBGG_MX' + Xmass + '*lowpurity.root'
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/BBGG_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    #multiprocessing to convert parquet to root, and hadd the root files for highpurity and lowpurity zzggsignal
    #highpurity zzgg
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(highpurity_zzggfiles):
        result = pool.apply_async(process_highpurity_zzggfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_zzgg_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/ZZGG_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/ZZGG_MX'+Xmass+'*highpurity.root'
    os.system(command)
    #lowpurity zzgg
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_zzggfiles):
        result = pool.apply_async(process_lowpurity_zzggfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_zzgg_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/ZZGG_MX' + Xmass + '*lowpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/ZZGG_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(highpurity_ttggfiles):
        result = pool.apply_async(process_highpurity_ttggfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_ttgg_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTGG_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTGG_MX'+Xmass+'*highpurity.root'
    os.system(command)
    #lowpurity ttgg
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_ttggfiles):
        result = pool.apply_async(process_lowpurity_ttggfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_ttgg_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTGG_MX' + Xmass + '*lowpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTGG_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(highpurity_VBFfiles):
        result = pool.apply_async(process_highpurity_VBFfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_VBF_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VBF_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VBF_MX'+Xmass+'*highpurity.root'
    os.system(command)
    #lowpurity VBF
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_VBFfiles):
        result = pool.apply_async(process_lowpurity_VBFfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_VBF_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VBF_MX' + Xmass + '*lowpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VBF_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(highpurity_VHfiles):
        result = pool.apply_async(process_highpurity_VHfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_VH_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VH_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VH_MX'+Xmass+'*highpurity.root'
    os.system(command)
    #lowpurity VH
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_VHfiles):
        result = pool.apply_async(process_lowpurity_VHfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_VH_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VH_MX' + Xmass + '*lowpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/VH_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(highpurity_GGHfiles):
        result = pool.apply_async(process_highpurity_GGHfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_GGH_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/GGH_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/GGH_MX'+Xmass+'*highpurity.root'
    os.system(command)
    #lowpurity GGH
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_GGHfiles):
        result = pool.apply_async(process_lowpurity_GGHfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_GGH_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/GGH_MX' + Xmass + '*lowpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/GGH_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(highpurity_TTHfiles):
        result = pool.apply_async(process_highpurity_TTHfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_TTH_cat12highpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTH_MX' + Xmass + '*highpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTH_MX'+Xmass+'*highpurity.root'
    os.system(command)
    #lowpurity TTH
    pool = multiprocessing.Pool(processes=10)
    results = []
    for file in tqdm(lowpurity_TTHfiles):
        result = pool.apply_async(process_lowpurity_TTHfile, args=(file,))
        results.append(result)
    pool.close()
    pool.join()
    tree_names = [result.get() for result in tqdm(results)]
    command = 'hadd ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/MX'+Xmass+'_MH125/MX'+Xmass+'_MH125_'+year+'_TTH_cat12lowpurity.root ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTH_MX' + Xmass + '*lowpurity.root'
    print(command)
    os.system(command)
    command = 'rm ./PBDT_HH_FHSL_combine_'+year+'/flashgginput/TTH_MX'+Xmass+'*lowpurity.root'
    os.system(command)
    
    
    