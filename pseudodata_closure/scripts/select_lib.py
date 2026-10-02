#!/usr/bin/env python3
"""Shared selection recipe for the HZa pseudo-data sculpting closure test.

Single source of truth for:
  - the adopted per-mA BDT working-point cuts (all 30 mass points),
  - the DY-background sample combination per era,
  - the fit-window / mass-range conventions (massLow=95 for all mA),
  - the event-selection loader (H_mass, factor weight, MVA_Score_mA_M{n}).

Both the real-data validation (validate_cuts.py) and the pseudo-data workspace
builder (build_pseudodata_ws.py) import from here so the selection is identical.

Reads scored P2Root under run3_bdt_scored_nominal (tree 'inclusive'):
  mass branch  H_mass
  weight       factor   (== 1 for Data)
  BDT score    MVA_Score_mA_M{mA}
"""
import os
import numpy as np
import uproot

SCORED_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal"
# CLOSURE_FSRFIX=1 (2026-09-26): read the FSR-fix re-production instead, where 2022 DY+gamma is
# the inclusive PTG-10to100 sample (the 10to50/50to100 slices were retired). Unset = unchanged.
FSRFIX = os.environ.get("CLOSURE_FSRFIX", "0") == "1"
if FSRFIX:
    SCORED_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
TREE = "inclusive"
MASS_BR = "H_mass"
WEIGHT_BR = "factor"
SCORE_FMT = "MVA_Score_mA_M{}"

ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

# DY background combination per era (identical to make_sculpt_R_table_from_wp.py /
# table_interpolate_bkgYield_1.py).
_bkg_2022 = ["DYGto2LG_10to100"] if FSRFIX else ["DYGto2LG_10to50", "DYGto2LG_50to100"]
_bkg_2023 = ["DYGto2LG_10to100"]
_bkg_dyll = ["DYJetsToLL"]
_bkg_dyll_2024 = ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
BKG_SAMPLES_BY_YEAR = {
    "2022preEE":    _bkg_2022 + _bkg_dyll,
    "2022postEE":   _bkg_2022 + _bkg_dyll,
    "2023preBPix":  _bkg_2023 + _bkg_dyll,
    "2023postBPix": _bkg_2023 + _bkg_dyll,
    "2024":         _bkg_2023 + _bkg_dyll_2024,
}

# Adopted per-mA BDT cuts (all 30 points):
#   trained points  -> Plot/output/MVAcut_points_run3.json
#   interpolated    -> MVAcut/run3_ReReco/makeWorkspace_data.py  --interp dict
# NOTE: mA1 is validated empirically against the real roohist (R=1 WP vs 0.9628).
# mA2 loosened 0.982 -> 0.975 (2026-07-26): the old WP had drifted to R=1.25 (a ~2.5 sigma
# fake in the pseudo-data closure); 0.975 is the local-R-minimum on the loose side,
# restoring R~1.08 and dropping the closure fake to ~1.5 sigma, consistent with the
# analysis's R~1 working-point criterion. See project_hza_pseudodata_sculpting_closure.
CUTS_TRAINED = {1: 0.9628, 2: 0.975, 3: 0.988, 4: 0.99, 5: 0.99, 6: 0.99,
                7: 0.99, 8: 0.99, 9: 0.99, 10: 0.99, 15: 0.99, 20: 0.985,
                25: 0.985, 30: 0.98}
# NOTE: 18,19,21,22 were validated empirically to use 0.985 (not the 0.99 in the
# stale makeWorkspace_data.py --interp dict); confirmed by reproducing the real
# roohist_data yields (validate_cuts.py).
CUTS_INTERP = {11: 0.99, 12: 0.99, 13: 0.99, 14: 0.99, 16: 0.99, 17: 0.99,
               18: 0.985, 19: 0.985, 21: 0.985, 22: 0.985, 23: 0.985, 24: 0.985,
               26: 0.985, 27: 0.985, 28: 0.98, 29: 0.98}
ADOPTED_CUTS = {**CUTS_TRAINED, **CUTS_INTERP}

if FSRFIX:
    # 2026-09-27: the hard-coded tables above are the pre-FSR-fix working points. With the
    # FSR-fix re-production the anchors come from the current MVAcut_points_run3.json and each
    # interpolated point takes the cut of its nearest anchor -- the same rule as
    # MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py::find_nearest_cut. Verified against the
    # real data roohists with validate_cuts.py.
    import json as _json
    _WP_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
    CUTS_TRAINED = {int(e["mA"]): float(e["MVAcut"]) for e in _json.load(open(_WP_JSON))["results"]}
    CUTS_INTERP = {m: CUTS_TRAINED[min(CUTS_TRAINED, key=lambda k: abs(k - m))]
                   for m in range(1, 31) if m not in CUTS_TRAINED}
    ADOPTED_CUTS = {**CUTS_TRAINED, **CUTS_INTERP}


def mass_low(mA):
    """Low edge of the CMS_hza_mass range: 95 for all mA.

    (2026-07-30: mA1,2 restored from 105 to 95 now that the low-mass sculpting is
    controlled -- R~1 working points, a clean background-only pseudo-data closure, and
    the fine 1 GeV mass-shape check all show no turn-on fake peak, so the turn-on region
    no longer needs to be skipped.)
    """
    return 95.0


MASS_HIGH = 180.0


def load_combined(mA, samples_by_year=BKG_SAMPLES_BY_YEAR):
    """Load combined (mass, weight, score) over all (sample, era) for one mA."""
    br = SCORE_FMT.format(int(mA))
    m, w, s = [], [], []
    for era, samples in samples_by_year.items():
        for sample in samples:
            path = os.path.join(SCORED_BASE, sample, f"{era}.root")
            if not os.path.exists(path):
                raise FileNotFoundError(path)
            a = uproot.open(path)[TREE].arrays([MASS_BR, WEIGHT_BR, br], library="np")
            m.append(a[MASS_BR]); w.append(a[WEIGHT_BR]); s.append(a[br])
    return np.concatenate(m), np.concatenate(w), np.concatenate(s)


def select(mA, m, w, s, cut):
    """Window + BDT-cut mask for one mA."""
    lo = mass_low(mA)
    sel = (m > lo) & (m < MASS_HIGH) & np.isfinite(m) & np.isfinite(w) & np.isfinite(s)
    return sel & (s > cut)
