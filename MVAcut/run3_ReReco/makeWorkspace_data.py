# =============================================================================
# DEPRECATED -- 2025-12 vintage, kept for reference only. Does not run as-is.
#
# Superseded by the live chain in 1_grand.sh:
#     MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py   (reads the pre-scored
#     root_P2Root/run3_bdt_scored_nominal trees)  ->  Trees2WS/run_tree2ws.sh
#     -> trees2ws.py / trees2ws_data.py
# The call to this script is already commented out in Trees2WS/run_tree2ws.sh
# and in MVAcut/run3_ReReco/runjob.sh.
#
# Three things are stale (checked 2026-08-23); fixing only one is not enough:
#   1. Input tree  /eos/home-p/pelai/HZa/Root_Dataset/run3/Data/run3.root
#      no longer exists. Present equivalent:
#      /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal/Data/run3.root
#   2. BDT_filename HZaMVA/model_Za_BDT_run3.pkl no longer exists. Copies now
#      live in HZaMVA/{using,scripts}/ (2026-07, 16 features) and
#      HZaMVA/{keep,train3mA}/ (2025-12, 15 features).
#   3. The MVA_list below builds 15 features, matching only the 2025-12 model.
#      The current model needs 16: photon/H pT are normalized (pho1Pt_oHm,
#      pho2Pt_oHm, H_pt_oHm) and pho_pt_asym is added -- see
#      HiggsZaAna/HZaMVA/scripts/hza_features.py. Feeding it 15 raises
#      "ValueError: Feature shape mismatch, expected: 16, got 15".
#
# So this script did NOT produce the current MVAcut workspaces. Do not treat it
# as the reference for how the per-mA selection is defined today.
# =============================================================================
#############################
from ROOT import *
from ROOT import RooFit
import gc

from xgboost import XGBClassifier
import pickle

#############################

import argparse
parser = argparse.ArgumentParser(description="Prepare flashggFinalFit workspace: signal Fit")
parser.add_argument("-m", "--mass", dest="mass", type=float, default=1.0, help="ALP mass")
parser.add_argument('--interp', dest='interp', action='store_true', default=False, help='interpolation?')
args = parser.parse_args()

mass = int(args.mass)

BDT_filename="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/model_Za_BDT_run3.pkl"
# BDT_filename="/publicfs/cms/user/laipeizhu/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"
if args.interp:
    mvaCuts = {11:0.99, 12:0.99, 13:0.99, 14:0.99, 16:0.99, 17:0.99, 18:0.99, 19:0.99, 21:0.99, 22:0.99, 23:0.985, 24:0.985, 26:0.985, 27:0.985, 28:0.98, 29:0.98}
else:
    # mvaCuts = {1:0.955, 2:0.98, 3:0.985, 4:0.98, 5:0.985, 6:0.99, 7:0.985, 8:0.99, 9:0.99, 10:0.99, 15:0.99, 20:0.99, 25:0.985, 30:0.98} # zebing
    # mvaCuts = {1:0.945, 2:0.975, 3:0.985, 4:0.985, 5:0.99, 6:0.99, 7:0.99, 8:0.995, 9:0.995, 10:0.995, 15:0.99, 20:0.99, 25:0.99, 30:0.99} #  run2 trial 20 times 420
    # mvaCuts = {1:0.88, 2:0.94, 3:0.965, 4:0.955, 5:0.97, 6:0.975, 7:0.975, 8:0.975, 9:0.97, 10:0.975, 15:0.97, 20:0.97, 25:0.97, 30:0.96} # run2 trail 20 times 300
    mvaCuts = {1:0.88, 2:0.94, 3:0.965, 4:0.955, 5:0.97, 6:0.975, 7:0.975, 8:0.975, 9:0.97, 10:0.975, 15:0.97, 20:0.97, 25:0.97, 30:0.96} # use run2 to select run3
    # mvaCuts = {5:0.975, 15:0.97, 30:0.96} # train from run3 3 mass points


mvaCut = mvaCuts[mass]
model = pickle.load(open(BDT_filename, 'rb'))

############################
myfile = TFile('/eos/home-p/pelai/HZa/Root_Dataset/run3/Data/run3.root')
# myfile = TFile('/publicfs/cms/user/laipeizhu/ALP/Analysis_output/UL/run2/ALP_data.root')
mychain = myfile.Get('inclusive')
entries = mychain.GetEntriesFast()

print (entries)


w = RooWorkspace("CMS_hza_workspace")

Sqrts = RooRealVar("Sqrts","Sqrts",13.6)
CMS_hza_mass = RooRealVar("CMS_hza_mass","CMS_hza_mass",125.,95.,180.)

getattr(w,'import')(Sqrts)
getattr(w,'import')(CMS_hza_mass)


data_mass_cat0 = RooDataSet("data_mass_cat0","data_mass_cat0",RooArgSet(CMS_hza_mass))

data_mass_cat0_clone = data_mass_cat0.Clone()

for jentry in range(entries):
    nb = mychain.GetEntry(jentry)
    # if not mychain.passChaHadIso: continue
    # if not mychain.passNeuHadIso: continue
    # if not mychain.passdR_gl: continue
    # if not mychain.passHOverE: continue
    if mychain.H_m<95. or mychain.H_m>180.: continue


    param = (mychain.ALP_m - float(mass))/mychain.H_m
    MVA_list = [mychain.pho1Pt, mychain.pho1R9, mychain.pho1IetaIeta55, mychain.pho1PIso_noCorr ,mychain.pho2Pt, mychain.pho2R9, mychain.pho2IetaIeta55, mychain.pho2PIso_noCorr, mychain.ALP_calculatedPhotonIso, mychain.var_dR_Za, mychain.var_dR_g1g2, mychain.var_dR_g1Z, mychain.var_PtaOverMh, mychain.H_pt, param]
    MVA_value = model.predict_proba([MVA_list])[:, 1]
    if MVA_value < mvaCut:continue

    CMS_hza_mass.setVal(mychain.H_m)
    data_mass_cat0_clone.add(RooArgSet(CMS_hza_mass))

# c1 = TCanvas("c1","With Weight")
# c1.cd()
# massDist = CMS_hza_mass.frame(RooFit.Title("CMS_hza_mass"))
# data_mass_cat0.plotOn(massDist)
# massDist.Draw()
# c1.SaveAs("ALP_data_bkg_Am{0}_workspace.png".format(mass))

getattr(w,'import')(data_mass_cat0_clone)
# getattr(w, 'import')(data_mass_cat0_clone, RooFit.RecycleConflictNodes()) # PZ

w.writeToFile("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run3_ReReco/output/data/ALP_data_bkg_Am{0}_workspace.root".format(mass))

myfile.Close()           # PZ, Close the ROOT file
# del myfile             # PZ, Delete the file object to clear ROOT's memory

gc.collect()

# data_mass_cat0_clone.Delete()  # PZ, Explicitly delete the dataset
# del data_mass_cat0_clone       # PZ, Delete the object to release memory
# del w

