#!/usr/bin/env python3
"""Rebuild an fTest input workspace from a previously written multipdf workspace.

The MVAcut workspaces under /eos/home-p/pelai/HZa/root_MVAcut/data/ only still
exist for mA_M1 and mA_M2, so fTest cannot be rerun for the other mass points
from the original inputs.  The multipdf file that fTest itself wrote, however,
still carries the exact dataset that was fitted (``Data_13p6TeV``), so an
equivalent input workspace can be reconstructed from it.

The rebuilt file contains the same unbinned dataset, renamed to the name fTest
looks for, inside a workspace called ``CMS_hza_workspace``.  Refitting it
reproduces the original fit by construction -- this is a way to rerun fTest with
different PDF parameter ranges, not a way to recover the lost MVAcut step.

Input
  <fitres>/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root   (or a backup copy)

Output
  <outdir>/mA_M<mA>/ws/run3.root

Usage
  cmsenv && python3 scripts/rebuild_ftest_input_from_multipdf.py --masses 29,30 \
      --source-suffix _backup_20260802
"""

import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BKG_DIR = os.path.dirname(SCRIPT_DIR)
DEF_FITRES = os.path.join(BKG_DIR, "ALP_BkgModel_ReReco", "fit_results_run3")
DEF_OUTDIR = os.path.join(BKG_DIR, "ALP_BkgModel_ReReco", "ftest_inputs_rebuilt")

# fTest with --isData 0 looks for these, in this order.
TARGET_DATASET = "data_mass_cat0"
SOURCE_CANDIDATES = ["Data_13p6TeV", "Data_13TeV", "data_mass_cat0", "data_mass"]


def rebuild(mA, fitres, suffix, outdir, mass_lo, mass_hi):
    src = os.path.join(fitres, "%d%s" % (mA, suffix), "CMS-HGG_mva_13p6TeV_multipdf.root")
    if not os.path.isfile(src):
        print("[skip] mA=%d: no multipdf at %s" % (mA, src))
        return None

    fin = ROOT.TFile.Open(src)
    if not fin or fin.IsZombie():
        print("[skip] mA=%d: cannot open %s" % (mA, src))
        return None

    ws_in = None
    for key in fin.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("RooWorkspace"):
            ws_in = obj
            break
    if ws_in is None:
        print("[skip] mA=%d: no workspace in %s" % (mA, src))
        fin.Close()
        return None

    data = None
    for name in SOURCE_CANDIDATES:
        data = ws_in.data(name)
        if data:
            break
    if data is None:
        print("[skip] mA=%d: none of %s found" % (mA, SOURCE_CANDIDATES))
        fin.Close()
        return None
    if data.ClassName() != "RooDataSet":
        print("[skip] mA=%d: %s is %s, need an unbinned RooDataSet"
              % (mA, data.GetName(), data.ClassName()))
        fin.Close()
        return None

    mass = ws_in.var("CMS_hza_mass")
    if mass is None:
        print("[skip] mA=%d: no CMS_hza_mass in the workspace" % mA)
        fin.Close()
        return None
    mass.setRange(mass_lo, mass_hi)
    mass.setBins(int(round(mass_hi - mass_lo)))

    ws_out = ROOT.RooWorkspace("CMS_hza_workspace", "CMS_hza_workspace")
    copy = data.Clone(TARGET_DATASET)
    copy.SetName(TARGET_DATASET)
    getattr(ws_out, "import")(copy)

    lumi = ws_in.var("IntLumi")
    if lumi:
        getattr(ws_out, "import")(lumi, ROOT.RooFit.RecycleConflictNodes())

    dst_dir = os.path.join(outdir, "mA_M%d" % mA, "ws")
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, "run3.root")
    ws_out.writeToFile(dst)

    n_in = data.numEntries()
    check = ROOT.TFile.Open(dst)
    ws_chk = check.Get("CMS_hza_workspace")
    d_chk = ws_chk.data(TARGET_DATASET) if ws_chk else None
    n_out = d_chk.numEntries() if d_chk else -1
    check.Close()
    fin.Close()

    ok = (n_out == n_in)
    print("[%s] mA=%-2d  %s (%d entries) -> %s (%d entries)"
          % ("ok  " if ok else "FAIL", mA, data.GetName(), n_in, dst, n_out))
    return dst if ok else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masses", default="29,30")
    ap.add_argument("--fitres", default=DEF_FITRES)
    ap.add_argument("--outdir", default=DEF_OUTDIR)
    ap.add_argument("--source-suffix", default="",
                    help="suffix of the source directory, e.g. _backup_20260802")
    ap.add_argument("--mass-lo", type=float, default=95.0)
    ap.add_argument("--mass-hi", type=float, default=180.0)
    opts = ap.parse_args()

    masses = [int(x) for x in opts.masses.split(",") if x.strip()]
    print("[input ] %s/<mA>%s/CMS-HGG_mva_13p6TeV_multipdf.root" % (opts.fitres, opts.source_suffix))
    print("[output] %s/mA_M<mA>/ws/run3.root" % opts.outdir)

    failures = 0
    for mA in masses:
        if rebuild(mA, opts.fitres, opts.source_suffix, opts.outdir,
                   opts.mass_lo, opts.mass_hi) is None:
            failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
