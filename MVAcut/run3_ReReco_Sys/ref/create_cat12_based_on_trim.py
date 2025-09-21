import awkward as ak
import sys
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
import json
import multiprocessing

# ------------------------------------------------------------
# Utils
# ------------------------------------------------------------

def str_to_list(arg):
    return ast.literal_eval(arg)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--is_HH', type=str, default="False", help='is HH')
parser.add_argument('--inputFHFiles', type=str_to_list, help='inputFHFiles List')
parser.add_argument('--inputBKGFiles', type=str_to_list, help='input pp and dd Files List')
parser.add_argument('--data', type=str, default="/eos/user/s/shsong/HiggsDNA/UL18data/merged_nominal.parquet", help='data file')
parser.add_argument('--year', type=str, default="2017", help='year')
parser.add_argument('--local', default=False, action='store_true', help="run locally or on condor")
parser.add_argument('--model', type=str, default="/eos/user/z/zhenxuan/PNN_wwgg/boosted_FHSL/data/simple_DNN_train_boosted_100_epoch_all_lr000001/model.pth", help='model file')
parser.add_argument('--scalar', type=str, default="/eos/user/z/zhenxuan/PNN_wwgg/boosted_FHSL/data/simple_DNN_train_boosted_100_epoch_all_lr000001/scaler_params.json", help='scalar file')

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
    sys.path.append(f"./PBDT_HH_FHSL_combine_{year}/pkgs_condor/parquet_to_root-0.3.0")
    sys.path.append(f"./PBDT_HH_FHSL_combine_{year}/pkgs_condor/bayesian-optimization-1.4.3")

from parquet_to_root import parquet_to_root
from bayes_opt import BayesianOptimization
from bayes_opt.util import UtilityFunction

# ------------------------------------------------------------
# Model
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Refactor helpers
# ------------------------------------------------------------

def build_jet_vars(events):
    fat1 = vector.obj(pt=events.fatjet_1_pt, eta=events.fatjet_1_eta, phi=events.fatjet_1_phi, mass=events.fatjet_1_msoftdrop)
    fat2 = vector.obj(pt=events.fatjet_2_pt, eta=events.fatjet_2_eta, phi=events.fatjet_2_phi, mass=events.fatjet_2_msoftdrop)
    fat3 = vector.obj(pt=events.fatjet_3_pt, eta=events.fatjet_3_eta, phi=events.fatjet_3_phi, mass=events.fatjet_3_msoftdrop)
    diph = vector.obj(pt=events.Diphoton_pt, eta=events.Diphoton_eta, phi=events.Diphoton_phi, mass=events.Diphoton_mass)
    lead = vector.obj(pt=events.LeadPhoton_pt, eta=events.LeadPhoton_eta, phi=events.LeadPhoton_phi, mass=events.LeadPhoton_mass)
    subl = vector.obj(pt=events.SubleadPhoton_pt, eta=events.SubleadPhoton_eta, phi=events.SubleadPhoton_phi, mass=events.SubleadPhoton_mass)

    events['fatjet_1_diphoton_dR'] = np.where(events.fatjet_1_pt>0, fat1.deltaR(diph), -999)
    events['fatjet_2_diphoton_dR'] = np.where(events.fatjet_2_pt>0, fat2.deltaR(diph), -999)
    events['fatjet_1_leadphoton_dR'] = np.where(events.fatjet_1_pt>0, fat1.deltaR(lead), -999)
    events['fatjet_1_subleadphoton_dR'] = np.where(events.fatjet_1_pt>0, fat1.deltaR(subl), -999)
    events['fatjet_2_leadphoton_dR'] = np.where(events.fatjet_2_pt>0, fat2.deltaR(lead), -999)
    events['fatjet_2_subleadphoton_dR'] = np.where(events.fatjet_2_pt>0, fat2.deltaR(subl), -999)
    events['fatjet_1_2_dR'] = np.where((events.fatjet_1_pt>0)&(events.fatjet_2_pt>0), fat1.deltaR(fat2), -999)

    m12 = (fat1+fat2).mass
    m13 = (fat1+fat3).mass
    m23 = (fat2+fat3).mass
    events['max_fatjets_mass'] = np.maximum(m12, np.maximum(m13, m23))

    eps = 1e-9
    events['fatjet_1_XbbvsQCDMD'] = events['fatjet_1_particleNetMD_Xbb'] / (events['fatjet_1_particleNetMD_Xbb'] + events['fatjet_1_particleNetMD_QCD'] + eps)
    events['fatjet_2_XbbvsQCDMD'] = events['fatjet_2_particleNetMD_Xbb'] / (events['fatjet_2_particleNetMD_Xbb'] + events['fatjet_2_particleNetMD_QCD'] + eps)
    return events

def load_events_parquet(filename, mx, my):
    events = ak.from_parquet(filename)
    sel = ((events["category"]==1) | (events["category"]==2)) & (events["Diphoton_minID_modified"]>-0.7)
    events = events[sel]
    events = build_jet_vars(events)
    events['mx'] = np.ones(len(events))*int(mx)
    events['my'] = np.ones(len(events))*int(my)
    return events

def predict_pnn_score(events, model, scaler, input_features):
    df = ak.to_pandas(events[input_features + ['weight_central','Diphoton_mass']])
    df = df.replace(-999, 0)
    X_test = torch.tensor(scaler.transform(df[input_features])).float()
    with torch.no_grad():
        proba = torch.softmax(model(X_test), dim=1).cpu().numpy()
    pnn_score = (proba[:,4] + proba[:,2] + proba[:,3]) / (proba[:,0] + proba[:,1] + proba[:,2] + proba[:,3] + proba[:,4] + 1e-12)
    return pnn_score

def infer_sys_tag(parquet_path: str):
    name = os.path.basename(parquet_path)
    if "merged_" in name:
        part = name.split("merged_")[-1]
        if "_down" in part:
            return part.split("_down")[0], "Down01sigma"
        if "_up" in part:
            return part.split("_up")[0], "Up01sigma"
    return "", ""

TREE_PREFIX = {
    "wwgg": "gghhwwgg_125_13TeV",
    "bbgg": "gghhbbgg_125_13TeV",
    "zzgg": "gghhzzgg_125_13TeV",
    "ttgg": "gghhttgg_125_13TeV",
    "VBF" : "gghhVBF_125_13TeV",
    "VH"  : "gghhVH_125_13TeV",
    "TTH" : "gghhTTH_125_13TeV",
    "GGH" : "gghhGGH_125_13TeV",
}

def make_tree_name(sample_key: str, purity: str, parquet_path: str):
    base = TREE_PREFIX[sample_key] + f"_cat12{purity}"
    sys_name, sys_dir = infer_sys_tag(parquet_path)
    if sys_name:
        return f"{base}_{sys_name}{sys_dir}"
    return base

def convert_parquet_to_root(parquet_path: str, out_dir: str, sample_key: str, purity: str, year: str):
    os.makedirs(out_dir, exist_ok=True)
    tree_name = make_tree_name(sample_key, purity, parquet_path)
    out_root  = os.path.join(out_dir, os.path.basename(parquet_path).replace(".parquet", ".root"))
    parquet_to_root(parquet_path, out_root, treename=tree_name)
    return out_root, tree_name

# ------------------------------------------------------------
# Original helper with fixes: kinematic reweight
# ------------------------------------------------------------

def kinematic_reweight(events_data, events_bkg_pp, events_bkg_dd, weight_data, weight_bkg_pp, weight_bkg_dd, var_name_list, bins_list):
    w_pp = np.array(weight_bkg_pp)  # copy
    for idx in range(len(var_name_list)):
        vname = var_name_list[idx]
        data_vals = np.array(events_data[vname])
        dd_vals   = np.array(events_bkg_dd[vname])
        pp_vals   = np.array(events_bkg_pp[vname])
        hist_data, bins = np.histogram(data_vals, bins=100, range=(bins_list[idx][0], bins_list[idx][1]), weights=np.array(weight_data))
        hist_dd, _      = np.histogram(dd_vals,   bins=100, range=(bins_list[idx][0], bins_list[idx][1]), weights=np.array(weight_bkg_dd))
        hist_pp, _      = np.histogram(pp_vals,   bins=100, range=(bins_list[idx][0], bins_list[idx][1]), weights=np.array(weight_bkg_pp))
        for j in range(len(bins)-1):
            denom = hist_pp[j]
            if denom <= 0:
                continue
            rw = (hist_data[j]-hist_dd[j]) / denom
            sel = (pp_vals>bins[j]) & (pp_vals<=bins[j+1])
            w_pp[sel] *= rw
    return w_pp

# ------------------------------------------------------------
# Original SF helpers (kept for compatibility)
# ------------------------------------------------------------

def add_HtaggerSF(event, mass):
    mass_list=[500,550,600,650,700,750,800,850,900,1000,1100,1200,1250,1300,1400,1500,1600,1700,1750,1800,1900,2000,2200,2400,2500,2600,2800,3000]
    index=mass_list.index(mass)
    HvsQCD_SF=[0.993610,1.083977,0.983107,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,0.993610,1.083977,0.983107,1.134969,1.134969,1.097200,0.973850,0.970535,1.039317,1.008370,1.008370,0.962904,0.969104,1.002570,0.941431,1.025998,1.025998,1.004931,0.967810]
    HvsQCD_SFup=np.array(HvsQCD_SF)*1.22
    HvsQCD_SFdown=np.array(HvsQCD_SF)*0.75
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

# NOTE: This is the original long add_sf_branches from user script (kept as-is with tiny safety guards)

def add_sf_branches(events):
    events["CMS_hgg_mass"]=events["Diphoton_mass"]
    events["weight"]=events["weight_central"]
    events["dZ"]=np.ones(len(events['CMS_hgg_mass']))
    # ratios with guards
    def ratio(up, cen):
        return events[up] / ak.where(events[cen]==0, 1, events[cen])
    events["muon_highptreco_sf_Down01sigma"]=ratio("weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_down","weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_central")
    events["muon_highptreco_sf_Up01sigma"]=ratio("weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_up","weight_nonisomuon_highptreco_sf_SelectedMuon_noiso_central")
    events["muon_highptid_sf_Down01sigma"]=ratio("weight_nonisomuon_highptid_sf_SelectedMuon_noiso_down","weight_nonisomuon_highptid_sf_SelectedMuon_noiso_central")
    events["muon_highptid_sf_Up01sigma"]=ratio("weight_nonisomuon_highptid_sf_SelectedMuon_noiso_up","weight_nonisomuon_highptid_sf_SelectedMuon_noiso_central")
    events["L1_prefiring_sf_Down01sigma"]=ratio("weight_L1_prefiring_sf_down","weight_L1_prefiring_sf_central")
    events["L1_prefiring_sf_Up01sigma"]=ratio("weight_L1_prefiring_sf_up","weight_L1_prefiring_sf_central")
    events["puWeight_Up01sigma"]=ratio("weight_pu_reweight_sf_up","weight_pu_reweight_sf_central")
    events["puWeight_Down01sigma"]=ratio("weight_pu_reweight_sf_down","weight_pu_reweight_sf_central")
    events["jet_pu_id_sf_Up01sigma"]=ratio("weight_jet_puid_sf_SelectedJet_up","weight_jet_puid_sf_SelectedJet_central")
    events["jet_pu_id_sf_Down01sigma"]=ratio("weight_jet_puid_sf_SelectedJet_down","weight_jet_puid_sf_SelectedJet_central")
    events["electron_veto_sf_Diphoton_Photon_Up01sigma"]=ratio("weight_electron_veto_sf_Diphoton_Photon_up","weight_electron_veto_sf_Diphoton_Photon_central")
    events["electron_veto_sf_Diphoton_Photon_Down01sigma"]=ratio("weight_electron_veto_sf_Diphoton_Photon_down","weight_electron_veto_sf_Diphoton_Photon_central")
    events["isoelectron_id_sf_SelectedElectron_iso_Up01sigma"]=ratio("weight_isoelectron_id_sf_SelectedElectron_iso_up","weight_isoelectron_id_sf_SelectedElectron_iso_central")
    events["isoelectron_id_sf_SelectedElectron_iso_Down01sigma"]=ratio("weight_isoelectron_id_sf_SelectedElectron_iso_down","weight_isoelectron_id_sf_SelectedElectron_iso_central")
    events["isoelectron_id_sf_SelectedElectron_noiso_Up01sigma"]= ratio("weight_isoelectron_id_sf_SelectedElectron_noiso_up","weight_isoelectron_id_sf_SelectedElectron_noiso_central")
    events["isoelectron_id_sf_SelectedElectron_noiso_Down01sigma"]= ratio("weight_isoelectron_id_sf_SelectedElectron_noiso_down","weight_isoelectron_id_sf_SelectedElectron_noiso_central")
    events["isomuon_id_sf_SelectedMuon_iso_Up01sigma"]=ratio("weight_isomuon_id_sf_SelectedMuon_iso_up","weight_isomuon_id_sf_SelectedMuon_iso_central")
    events["isomuon_id_sf_SelectedMuon_iso_Down01sigma"]=ratio("weight_isomuon_id_sf_SelectedMuon_iso_down","weight_isomuon_id_sf_SelectedMuon_iso_central")
    events["isomuon_iso_sf_SelectedMuon_iso_Up01sigma"]=ratio("weight_isomuon_iso_sf_SelectedMuon_iso_up","weight_isomuon_iso_sf_SelectedMuon_iso_central")
    events["isomuon_iso_sf_SelectedMuon_iso_Down01sigma"]=ratio("weight_isomuon_iso_sf_SelectedMuon_iso_down","weight_isomuon_iso_sf_SelectedMuon_iso_central")
    events["nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma"]=ratio("weight_nonisoelectron_id_sf_SelectedElectron_noiso_up","weight_nonisoelectron_id_sf_SelectedElectron_noiso_central")
    events["nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma"]=ratio("weight_nonisoelectron_id_sf_SelectedElectron_noiso_down","weight_nonisoelectron_id_sf_SelectedElectron_noiso_central")
    events["photon_id_sf_Diphoton_Photon_Up01sigma"]=ratio("weight_photon_id_sf_Diphoton_Photon_up","weight_photon_id_sf_Diphoton_Photon_central")
    events["photon_id_sf_Diphoton_Photon_Down01sigma"]=ratio("weight_photon_id_sf_Diphoton_Photon_down","weight_photon_id_sf_Diphoton_Photon_central")
    events["photon_presel_sf_Diphoton_Photon_Up01sigma"]=ratio("weight_photon_presel_sf_Diphoton_Photon_up","weight_photon_presel_sf_Diphoton_Photon_central")
    events["photon_presel_sf_Diphoton_Photon_Down01sigma"]=ratio("weight_photon_presel_sf_Diphoton_Photon_down","weight_photon_presel_sf_Diphoton_Photon_central")
    events["trigger_sf_Up01sigma"]=ratio("weight_trigger_sf_up","weight_trigger_sf_central")
    events["trigger_sf_Down01sigma"]=ratio("weight_trigger_sf_down","weight_trigger_sf_central")
    if "weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_central" in events.fields:
        events["PNetWvsQCDW1_sf_Up01sigma"]=ratio("weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_up","weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_central")
        events["PNetWvsQCDW1_sf_Down01sigma"]=ratio("weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_down","weight_PNet_WvsQCD_MD_sf_GenmatchendFatJet_1W_central")
        events["PNetWvsQCDW1_mistagging_sf_Up01sigma"]=ratio("weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_up","weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_central")
        events["PNetWvsQCDW1_mistagging_sf_Down01sigma"]=ratio("weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_down","weight_PNet_WvsQCD_mistagging_sf_UnmatchendFatJet_1W_central")
    else:
        one = ak.ones_like(events.weight_central)
        events["PNetWvsQCDW1_sf_Up01sigma"]=one
        events["PNetWvsQCDW1_sf_Down01sigma"]=one
        events["PNetWvsQCDW1_mistagging_sf_Up01sigma"]=one
        events["PNetWvsQCDW1_mistagging_sf_Down01sigma"]=one
    events["isoelectron_reco_sf_Up01sigma"]=ratio("weight_electron_reco_sf_SelectedElectron_iso_up","weight_electron_reco_sf_SelectedElectron_iso_central")
    events["isoelectron_reco_sf_Down01sigma"]=ratio("weight_electron_reco_sf_SelectedElectron_iso_down","weight_electron_reco_sf_SelectedElectron_iso_central")
    events["nonisoelectron_reco_sf_Up01sigma"]=ratio("weight_electron_reco_sf_SelectedElectron_noiso_up","weight_electron_reco_sf_SelectedElectron_noiso_central")
    events["nonisoelectron_reco_sf_Down01sigma"]=ratio("weight_electron_reco_sf_SelectedElectron_noiso_down","weight_electron_reco_sf_SelectedElectron_noiso_central")
    events["PTransformerHtagger_sf_Up01sigma"]=ratio("weight_PTransformer_up","weight_PTransformer_central")
    events["PTransformerHtagger_sf_Down01sigma"]=ratio("weight_PTransformer_down","weight_PTransformer_central")
    if "weight_PNbb_veto_sf_SelectedFatJet_up" in events.fields:
        events["PNXbb_sf_Up01sigma"]=ratio("weight_PNbb_veto_sf_SelectedFatJet_up","weight_PNbb_veto_sf_SelectedFatJet_central")
        events["PNXbb_sf_Down01sigma"]=ratio("weight_PNbb_veto_sf_SelectedFatJet_down","weight_PNbb_veto_sf_SelectedFatJet_central")
        events=events[['dZ','PNN_score','category','CMS_hgg_mass','weight','muon_highptreco_sf_Down01sigma','muon_highptreco_sf_Up01sigma','muon_highptid_sf_Down01sigma','muon_highptid_sf_Up01sigma','L1_prefiring_sf_Down01sigma','L1_prefiring_sf_Up01sigma','puWeight_Up01sigma','puWeight_Down01sigma','jet_pu_id_sf_Up01sigma','jet_pu_id_sf_Down01sigma','electron_veto_sf_Diphoton_Photon_Up01sigma','electron_veto_sf_Diphoton_Photon_Down01sigma','isoelectron_id_sf_SelectedElectron_iso_Up01sigma','isoelectron_id_sf_SelectedElectron_iso_Down01sigma','isoelectron_id_sf_SelectedElectron_noiso_Up01sigma','isoelectron_id_sf_SelectedElectron_noiso_Down01sigma','isomuon_id_sf_SelectedMuon_iso_Up01sigma','isomuon_id_sf_SelectedMuon_iso_Down01sigma','isomuon_iso_sf_SelectedMuon_iso_Up01sigma','isomuon_iso_sf_SelectedMuon_iso_Down01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma','photon_id_sf_Diphoton_Photon_Up01sigma','photon_id_sf_Diphoton_Photon_Down01sigma','photon_presel_sf_Diphoton_Photon_Up01sigma','photon_presel_sf_Diphoton_Photon_Down01sigma','trigger_sf_Up01sigma','trigger_sf_Down01sigma','PNetWvsQCDW1_sf_Up01sigma','PNetWvsQCDW1_sf_Down01sigma','PNetWvsQCDW1_mistagging_sf_Up01sigma','PNetWvsQCDW1_mistagging_sf_Down01sigma','isoelectron_reco_sf_Up01sigma','isoelectron_reco_sf_Down01sigma','nonisoelectron_reco_sf_Up01sigma','nonisoelectron_reco_sf_Down01sigma','PNXbb_sf_Up01sigma','PNXbb_sf_Down01sigma','PTransformerHtagger_sf_Up01sigma','PTransformerHtagger_sf_Down01sigma']]
    else:
        events=events[['dZ','PNN_score','category','CMS_hgg_mass','weight','muon_highptreco_sf_Down01sigma','muon_highptreco_sf_Up01sigma','muon_highptid_sf_Down01sigma','muon_highptid_sf_Up01sigma','L1_prefiring_sf_Down01sigma','L1_prefiring_sf_Up01sigma','puWeight_Up01sigma','puWeight_Down01sigma','jet_pu_id_sf_Up01sigma','jet_pu_id_sf_Down01sigma','electron_veto_sf_Diphoton_Photon_Up01sigma','electron_veto_sf_Diphoton_Photon_Down01sigma','isoelectron_id_sf_SelectedElectron_iso_Up01sigma','isoelectron_id_sf_SelectedElectron_iso_Down01sigma','isoelectron_id_sf_SelectedElectron_noiso_Up01sigma','isoelectron_id_sf_SelectedElectron_noiso_Down01sigma','isomuon_id_sf_SelectedMuon_iso_Up01sigma','isomuon_id_sf_SelectedMuon_iso_Down01sigma','isomuon_iso_sf_SelectedMuon_iso_Up01sigma','isomuon_iso_sf_SelectedMuon_iso_Down01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Up01sigma','nonisoelectron_id_sf_SelectedElectron_noiso_Down01sigma','photon_id_sf_Diphoton_Photon_Up01sigma','photon_id_sf_Diphoton_Photon_Down01sigma','photon_presel_sf_Diphoton_Photon_Up01sigma','photon_presel_sf_Diphoton_Photon_Down01sigma','trigger_sf_Up01sigma','trigger_sf_Down01sigma','PNetWvsQCDW1_sf_Up01sigma','PNetWvsQCDW1_sf_Down01sigma','PNetWvsQCDW1_mistagging_sf_Up01sigma','PNetWvsQCDW1_mistagging_sf_Down01sigma','isoelectron_reco_sf_Up01sigma','isoelectron_reco_sf_Down01sigma','nonisoelectron_reco_sf_Up01sigma','nonisoelectron_reco_sf_Down01sigma','PTransformerHtagger_sf_Up01sigma','PTransformerHtagger_sf_Down01sigma']]
    return events

# ------------------------------------------------------------
# Slice & dump (refactor of process_sig_samples)
# ------------------------------------------------------------

def slice_and_dump_all(samples_one_mass, model, scaler, input_features, cut1, cut2, year):
    mx = samples_one_mass['sig_output_name'].split('_')[1].split('X')[1]
    my = samples_one_mass['sig_output_name'].split('_')[2].split('H')[1]

    # load
    events_sig = ak.concatenate([
        load_events_parquet(samples_one_mass['FHpath'], mx, my),
        load_events_parquet(samples_one_mass['SLpath'], mx, my)
    ])
    ev_dict = {
        'bbgg': load_events_parquet(samples_one_mass['BBGGpath'], mx, my),
        'zzgg': load_events_parquet(samples_one_mass['ZZggpath'], mx, my),
        'ttgg': load_events_parquet(samples_one_mass['TTggpath'], mx, my),
        'VBF' : load_events_parquet(samples_one_mass['VBFpath'], mx, my),
        'VH'  : load_events_parquet(samples_one_mass['VHpath'],  mx, my),
        'TTH' : load_events_parquet(samples_one_mass['ttHpath'], mx, my),
        'GGH' : load_events_parquet(samples_one_mass['ggHpath'], mx, my)
    }

    # PNN
    events_sig['PNN_score'] = predict_pnn_score(events_sig, model, scaler, input_features)
    for k in ev_dict:
        ev_dict[k]['PNN_score'] = predict_pnn_score(ev_dict[k], model, scaler, input_features)

    # attach Htagger SF for mass-dependent ones then shape SF branches (保持原有输出列)
    mass_int = int(mx)
    events_sig = add_HtaggerSF(events_sig, mass_int)
    for k in ev_dict:
        ev_dict[k] = add_HtaggerSF(ev_dict[k], mass_int)

    events_sig = add_sf_branches(events_sig)
    for k in ev_dict:
        ev_dict[k] = add_sf_branches(ev_dict[k])

    outdir = f"./PBDT_HH_FHSL_combine_{year}/"
    def dump_pair(base_name, events):
        hp = events[(events['PNN_score'] > cut1) & (events['PNN_score'] <= 1)]
        lp = events[(events['PNN_score'] > cut2) & (events['PNN_score'] <= cut1)]
        ak.to_parquet(hp, os.path.join(outdir, base_name + "_highpurity.parquet"))
        ak.to_parquet(lp, os.path.join(outdir, base_name + "_lowpurity.parquet"))
        return len(hp)

    n_hp = dump_pair(samples_one_mass['sig_output_name'], events_sig)
    dump_pair(samples_one_mass['bbgg_output_name'], ev_dict['bbgg'])
    dump_pair(samples_one_mass['zzgg_output_name'], ev_dict['zzgg'])
    dump_pair(samples_one_mass['ttgg_output_name'], ev_dict['ttgg'])
    dump_pair(samples_one_mass['vbf_output_name'],  ev_dict['VBF'])
    dump_pair(samples_one_mass['vh_output_name'],   ev_dict['VH'])
    dump_pair(samples_one_mass['tth_output_name'],  ev_dict['TTH'])
    dump_pair(samples_one_mass['ggh_output_name'],  ev_dict['GGH'])

    return n_hp

# ------------------------------------------------------------
# Batch convert and hadd (includes DataA/B conversion)
# ------------------------------------------------------------

def batch_convert_and_hadd(year, Xmass, purity, sample_key, glob_prefix):
    patt = f"./PBDT_HH_FHSL_combine_{year}/{glob_prefix}_MX{Xmass}_MH125_cat12_m*_{purity}.parquet"
    files = glob.glob(patt)
    out_dir = f"./PBDT_HH_FHSL_combine_{year}/flashgginput/"
    dst_dir = f"./PBDT_HH_FHSL_combine_{year}/flashgginput/MX{Xmass}_MH125"
    os.makedirs(dst_dir, exist_ok=True)

    pool, results = multiprocessing.Pool(processes=10), []
    for fpath in files:
        results.append(pool.apply_async(convert_parquet_to_root, args=(fpath, out_dir, sample_key, purity, year)))
    pool.close(); pool.join()
    _ = [r.get() for r in results]

    hadd_out = os.path.join(dst_dir, f"MX{Xmass}_MH125_{year}_{sample_key}_cat12{purity}.root")
    os.system(f"hadd {hadd_out} {out_dir}/{glob_prefix}_MX{Xmass}*{purity}.root")
    os.system(f"rm   {out_dir}/{glob_prefix}_MX{Xmass}*{purity}.root")


def run_all_hadd(year, Xmass, dataA_parquet, dataB_parquet, dataA_rootname, dataB_rootname, dataA_treename, dataB_treename):
    # signal (wwgg)
    for purity in ("highpurity","lowpurity"):
        batch_convert_and_hadd(year, Xmass, purity, "wwgg", "CombineFHSL")

    # Data A/B parquet -> ROOT (放在 wwgg 后，其他之前都可；保持原输出位置与树名)
    dst_dir = f"./PBDT_HH_FHSL_combine_{year}/flashgginput/MX{Xmass}_MH125"
    parquet_to_root(dataA_parquet, os.path.join(dst_dir, dataA_rootname), treename=dataA_treename, verbose=False)
    parquet_to_root(dataB_parquet, os.path.join(dst_dir, dataB_rootname), treename=dataB_treename, verbose=False)

    # others
    for key, prefix in [
        ("bbgg","BBGG"), ("zzgg","ZZGG"), ("ttgg","TTGG"),
        ("VBF","VBF"), ("VH","VH"), ("TTH","TTH"), ("GGH","GGH")
    ]:
        for purity in ("highpurity","lowpurity"):
            batch_convert_and_hadd(year, Xmass, purity, key, prefix)

# ------------------------------------------------------------
# Start main pipeline (keeps original behavior)
# ------------------------------------------------------------

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

list_of_files = args.inputFHFiles # merged_nominal.parquet should be the first one
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
    vbffile=f"/eos/user/s/shsong/HiggsDNA/SingleHiggs{year.split('20')[1]}/VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8_{path_year}/"+FHfile.split("/")[-1]
    vhfile=f"/eos/user/s/shsong/HiggsDNA/SingleHiggs{year.split('20')[1]}/VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_{path_year}/"+FHfile.split("/")[-1]
    tthfile=f"/eos/user/s/shsong/HiggsDNA/SingleHiggs{year.split('20')[1]}/ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_{path_year}/"+FHfile.split("/")[-1]
    gghfile=f"/eos/user/s/shsong/HiggsDNA/SingleHiggs{year.split('20')[1]}/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_{path_year}/"+FHfile.split("/")[-1]
    SLsignal_path_list.append(SLfile)
    ZZggsignal_path_list.append(zzggfile)
    TTggsignal_path_list.append(ttggfile)
    BBGGsignal_path_list.append(bbggfile)
    vbf_path_list.append(vbffile)
    vh_path_list.append(vhfile)
    tth_path_list.append(tthfile)
    ggh_path_list.append(gghfile)
    dir_name="CombineFHSL_MX" + FHfile.split("M-")[1].split("_")[0] + "_MH125_cat12_"+(FHfile.split("/")[-1]).split(".")[0]
    signal_output_name.append(dir_name)
    bbgg_output_name.append(dir_name.replace("CombineFHSL","BBGG"))
    zzgg_output_name.append(dir_name.replace("CombineFHSL","ZZGG"))
    ttgg_output_name.append(dir_name.replace("CombineFHSL","TTGG"))
    vbf_output_name.append(dir_name.replace("CombineFHSL","VBF"))
    vh_output_name.append(dir_name.replace("CombineFHSL","VH"))
    tth_output_name.append(dir_name.replace("CombineFHSL","TTH"))
    ggh_output_name.append(dir_name.replace("CombineFHSL","GGH"))

signal_samples = {
    'FHpath':FHsignal_path_list, 'sig_output_name':signal_output_name,
    'SLpath':SLsignal_path_list, 'ZZggpath':ZZggsignal_path_list, 'TTggpath':TTggsignal_path_list,
    'BBGGpath':BBGGsignal_path_list, 'VBFpath':vbf_path_list,'VHpath':vh_path_list,'ttHpath':tth_path_list,'ggHpath':ggh_path_list,
    'bbgg_output_name':bbgg_output_name, 'zzgg_output_name':zzgg_output_name,'ttgg_output_name':ttgg_output_name,
    'vbf_output_name':vbf_output_name, 'vh_output_name':vh_output_name, 'tth_output_name':tth_output_name, 'ggh_output_name':ggh_output_name
}

bkgfiles = args.inputBKGFiles

# ------------------------------------- Get boundary using the first mass point ------------------------------------- #
mx = signal_samples['sig_output_name'][0].split('_')[1].split('X')[1]
my = signal_samples['sig_output_name'][0].split('_')[2].split('H')[1]

# load events for first point
get_sig_events_forApply = load_events_parquet  # alias for readability

events_sigFH = get_sig_events_forApply(signal_samples['FHpath'][0], mx, my)
events_sigFH['signal'] = 2*np.ones(len(events_sigFH))
events_sigSL = get_sig_events_forApply(signal_samples['SLpath'][0], mx, my)
events_sigSL['signal'] = np.ones(len(events_sigSL))

events_bbgg = get_sig_events_forApply(signal_samples['BBGGpath'][0], mx, my)
events_zzgg = get_sig_events_forApply(signal_samples['ZZggpath'][0], mx, my)
events_ttgg = get_sig_events_forApply(signal_samples['TTggpath'][0], mx, my)

events_data = get_sig_events_forApply(datapath, mx, my)
events_pp_cat1= get_sig_events_forApply(bkgfiles[0], mx, my)
events_pp_cat2= get_sig_events_forApply(bkgfiles[1], mx, my)
events_dd_cat1= get_sig_events_forApply(bkgfiles[2], mx, my)
events_dd_cat2= get_sig_events_forApply(bkgfiles[3], mx, my)
events_vbf= get_sig_events_forApply(signal_samples['VBFpath'][0], mx, my)
events_vh= get_sig_events_forApply(signal_samples['VHpath'][0], mx, my)
events_tth= get_sig_events_forApply(signal_samples['ttHpath'][0], mx, my)
events_ggh= get_sig_events_forApply(signal_samples['ggHpath'][0], mx, my)

# load model & scaler
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

# reweight MC pp using data - dd sideband
print('Reweighting pp background...')
events_pp = ak.concatenate([events_pp_cat1,events_pp_cat2]); del events_pp_cat1, events_pp_cat2

events_dd = ak.concatenate([events_dd_cat1,events_dd_cat2]); del events_dd_cat1, events_dd_cat2

MC_pp_new_weight = kinematic_reweight(
    events_data=events_data,
    events_bkg_pp=events_pp,
    events_bkg_dd=events_dd,
    weight_data=events_data['weight_central'],
    weight_bkg_pp=events_pp['weight_central'],
    weight_bkg_dd=events_dd['weight_central'],
    var_name_list = ['Diphoton_pt'],
    bins_list     = [[0, 300]]
)

events_pp['weight_central'] = MC_pp_new_weight

# predict PNN for bkg MC
print("get the PNN score for bkgmc")
events_pp['PNN_score'] = predict_pnn_score(events_pp, model, loaded_scaler, input_features)
events_dd['PNN_score'] = predict_pnn_score(events_dd, model, loaded_scaler, input_features)
event_bkgmc = ak.concatenate([events_pp, events_dd]); del events_pp, events_dd

# predict PNN for signals and others (first mass)
print("get the PNN score for signal & others (first mass)")
events_sigFH['PNN_score'] = predict_pnn_score(events_sigFH, model, loaded_scaler, input_features)
events_sigSL['PNN_score'] = predict_pnn_score(events_sigSL, model, loaded_scaler, input_features)

events_bbgg['PNN_score'] = predict_pnn_score(events_bbgg, model, loaded_scaler, input_features)

events_data['PNN_score']  = predict_pnn_score(events_data, model, loaded_scaler, input_features)

events_zzgg['PNN_score']  = predict_pnn_score(events_zzgg, model, loaded_scaler, input_features)

events_ttgg['PNN_score']  = predict_pnn_score(events_ttgg, model, loaded_scaler, input_features)

events_ggh['PNN_score']   = predict_pnn_score(events_ggh,  model, loaded_scaler, input_features)
events_vbf['PNN_score']   = predict_pnn_score(events_vbf,  model, loaded_scaler, input_features)
events_vh['PNN_score']    = predict_pnn_score(events_vh,   model, loaded_scaler, input_features)
events_tth['PNN_score']   = predict_pnn_score(events_tth,  model, loaded_scaler, input_features)

# plotting dir
Xmass = signal_samples['FHpath'][0].split("M-")[1].split("_")[0]
directory_path = f"./PBDT_HH_FHSL_combine_{year}/flashgginput/MX{Xmass}_MH125"
if os.path.exists(directory_path):
    os.system(f"rm -rf {directory_path}")
    print(f"The directory {directory_path} exists.")
    os.mkdir(directory_path)
else:
    os.mkdir(directory_path)

# plot PNN distribution (keep same style & filenames)
plt.figure(figsize=(8, 6))
hep.style.use("CMS")
plt.hist(ak.to_numpy(ak.flatten( (ak.Array(events_sigFH['PNN_score']).to_numpy()) )),
         weights = ak.to_numpy(ak.flatten( (ak.Array(events_sigFH['weight_central']).to_numpy()) ))*30,
         range=(0,1), bins=20, histtype='step', label='30*signal',)
plt.hist(ak.to_numpy(ak.flatten( (ak.Array(events_bbgg['PNN_score']).to_numpy()) )),
         weights = ak.to_numpy(ak.flatten( (ak.Array(events_bbgg['weight_central']).to_numpy()) ))*30,
         range=(0,1), bins=20, histtype='step', label='30*bbgg',)
mask_side = ((event_bkgmc.Diphoton_mass>135)|(event_bkgmc.Diphoton_mass<115))
plt.hist(ak.to_numpy(ak.flatten( (ak.Array(event_bkgmc['PNN_score'][mask_side]).to_numpy()) )),
         weights = ak.to_numpy(ak.flatten( (ak.Array(event_bkgmc['weight_central'][mask_side]).to_numpy()) )),
         range=(0,1), bins=20, histtype='stepfilled', label='pp+dd')

hist, bins = np.histogram(ak.to_numpy(ak.flatten( (ak.Array(events_data['PNN_score'][ (events_data.Diphoton_mass>135)|(events_data.Diphoton_mass<115) ]).to_numpy()) )), bins=20)
non_zero_bins = hist > 0
bin_centers = 0.5 * (bins[:-1] + bins[1:])
plt.scatter(bin_centers[non_zero_bins], hist[non_zero_bins], marker='o', color='black', label='data')
plt.xlabel('PNN Score')
plt.ylabel('Events')
plt.yscale('log')
plt.legend()
plt.savefig(f"./PBDT_HH_FHSL_combine_{year}/flashgginput/MX{Xmass}_MH125/dnnscore.png", dpi=140)
print("save the PNN score distribution plot")
print(f"./PBDT_HH_FHSL_combine_{year}/flashgginput/MX{Xmass}_MH125/dnnscore.png")
plt.close()

# boundaries
json_name=f"/eos/cms/store/group/phys_b2g/shsong/flashggws/cat12/{year}/MX{Xmass}_MH125/boundaries.json"
with open(json_name, 'r') as f:
    boundaries = json.load(f)
cut1 = boundaries[Xmass]['cut1']
cut2 = boundaries[Xmass]['cut2']
print('high purity cut:',cut1)
print('low purity cut:',cut2)

# merge signal FH+SL for efficiencies (kept)
events_sig_first = ak.concatenate([events_sigFH, events_sigSL])

# Htagger SF for samples then SF branches (as in original for first mass)
for ev in [events_sig_first, events_bbgg, events_zzgg, events_ttgg, events_vbf, events_vh, events_ggh, events_tth]:
    add_HtaggerSF(ev, int(mx))

# Attach SF branches (original)
events_sig_first = add_sf_branches(events_sig_first)
events_bbgg = add_sf_branches(events_bbgg)
events_zzgg = add_sf_branches(events_zzgg)
events_ttgg = add_sf_branches(events_ttgg)
events_vbf = add_sf_branches(events_vbf)
events_ggh = add_sf_branches(events_ggh)
events_vh = add_sf_branches(events_vh)
events_tth = add_sf_branches(events_tth)

# Prepare data A/B parquet
massname="MX"+(signal_samples['sig_output_name'][0].split("MX"))[1].split("_cat12")[0]

dataA_rootname=f"Data_{year}_cat12highpurity_{massname}.root"
dataA_treename="Data_13TeV_cat12highpurity"
dataB_rootname=f"Data_{year}_cat12lowpurity_{massname}.root"
dataB_treename="Data_13TeV_cat12lowpurity"

# data events
events_data_simple = ak.Array({
    'CMS_hgg_mass': events_data['Diphoton_mass'],
    'weight': events_data['weight_central'],
    'PNN_score': events_data['PNN_score']
})

hp_mask = (events_data_simple['PNN_score'] > cut1) & (events_data_simple['PNN_score'] <= 1)
lp_mask = (events_data_simple['PNN_score'] > cut2) & (events_data_simple['PNN_score'] <= cut1)

data_Acat_output_path=f"./PBDT_HH_FHSL_combine_{year}/"+dataA_rootname.replace(".root",".parquet")
data_Bcat_output_path=f"./PBDT_HH_FHSL_combine_{year}/"+dataB_rootname.replace(".root",".parquet")

ak.to_parquet(events_data_simple[hp_mask], data_Acat_output_path)
ak.to_parquet(events_data_simple[lp_mask], data_Bcat_output_path)

del events_data_simple

# First mass: also dump per-sample parquets for this mass (consistent with old behavior)
first_mass_samples = {
    'FHpath': signal_samples['FHpath'][0],
    'SLpath': signal_samples['SLpath'][0],
    'BBGGpath': signal_samples['BBGGpath'][0],
    'ZZggpath': signal_samples['ZZggpath'][0],
    'TTggpath': signal_samples['TTggpath'][0],
    'VBFpath': signal_samples['VBFpath'][0],
    'VHpath': signal_samples['VHpath'][0],
    'ttHpath': signal_samples['ttHpath'][0],
    'ggHpath': signal_samples['ggHpath'][0],
    'sig_output_name': signal_samples['sig_output_name'][0],
    'bbgg_output_name': signal_samples['bbgg_output_name'][0],
    'zzgg_output_name': signal_samples['zzgg_output_name'][0],
    'ttgg_output_name': signal_samples['ttgg_output_name'][0],
    'vbf_output_name': signal_samples['vbf_output_name'][0],
    'vh_output_name': signal_samples['vh_output_name'][0],
    'tth_output_name': signal_samples['tth_output_name'][0],
    'ggh_output_name': signal_samples['ggh_output_name'][0],
}
_ = slice_and_dump_all(first_mass_samples, model, loaded_scaler, input_features, cut1, cut2, year)

# Remaining masses
for i in range(1, len(signal_samples['FHpath'])):
    samples_one_mass = {
        'FHpath': signal_samples['FHpath'][i],
        'SLpath': signal_samples['SLpath'][i],
        'BBGGpath': signal_samples['BBGGpath'][i],
        'ZZggpath': signal_samples['ZZggpath'][i],
        'TTggpath': signal_samples['TTggpath'][i],
        'VBFpath': signal_samples['VBFpath'][i],
        'VHpath': signal_samples['VHpath'][i],
        'ttHpath': signal_samples['ttHpath'][i],
        'ggHpath': signal_samples['ggHpath'][i],
        'sig_output_name': signal_samples['sig_output_name'][i],
        'bbgg_output_name': signal_samples['bbgg_output_name'][i],
        'zzgg_output_name': signal_samples['zzgg_output_name'][i],
        'ttgg_output_name': signal_samples['ttgg_output_name'][i],
        'vbf_output_name': signal_samples['vbf_output_name'][i],
        'vh_output_name': signal_samples['vh_output_name'][i],
        'tth_output_name': signal_samples['tth_output_name'][i],
        'ggh_output_name': signal_samples['ggh_output_name'][i],
    }
    slice_and_dump_all(samples_one_mass, model, loaded_scaler, input_features, cut1, cut2, year)

# Convert parquet to root and hadd (now also includes DataA/B conversion)
print("starting to get root & hadd")
Xmass = signal_samples['FHpath'][0].split("M-")[1].split("_")[0]
run_all_hadd(year, Xmass, data_Acat_output_path, data_Bcat_output_path, dataA_rootname, dataB_rootname, dataA_treename, dataB_treename)

print("All done.")
