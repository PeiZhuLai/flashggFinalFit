#!/usr/bin/env python3
"""
Sanity script for the ALP H->Zg bkg-fit pipeline.

For every mA in [1..30] it collects:
  - envelope members from
      Background/ALP_BkgModel_ReReco/fit_results_run3/<mA>/HZAmassInde_fTest/EnvelopeResults.txt
  - expected limit + bands from
      Combine/output_combine_results/higgsCombine<mA>.AsymptoticLimits.mH125.38.root

Usage:
  # save current state as baseline CSV
  python3 sanity_envelope_limits.py --save sanity_old.csv

  # after re-running fTest + combine, snapshot again and diff vs baseline
  python3 sanity_envelope_limits.py --save sanity_new.csv --compare sanity_old.csv
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

FT_DIR = Path(
    "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/"
    "flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3"
)
COMBINE_DIR = Path(
    "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/"
    "flashggFinalFit/Combine/output_combine_results"
)
MASSES = list(range(1, 31))

ENV_PAT = re.compile(
    r"pdf\s*:\s*(\S+).*gof\s*:\s*([0-9.eE+\-]+).*isTruth\s*:\s*(\d+)"
)


def parse_envelope(mA):
    p = FT_DIR / str(mA) / "HZAmassInde_fTest" / "EnvelopeResults.txt"
    if not p.exists():
        return []
    members = []
    for line in p.read_text().splitlines():
        m = ENV_PAT.search(line)
        if m:
            members.append((m.group(1), float(m.group(2)), int(m.group(3))))
    return members


def read_limits():
    macro = f"""{{
  gROOT->SetBatch(kTRUE);
  for (int m = 1; m <= 30; ++m) {{
    TString fname = Form("{COMBINE_DIR}/higgsCombine%d.AsymptoticLimits.mH125.38.root", m);
    if (gSystem->AccessPathName(fname)) {{ printf("MISS %d\\n", m); continue; }}
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) {{ printf("MISS %d\\n", m); continue; }}
    TTree *t = (TTree*)f->Get("limit");
    if (!t) {{ printf("MISS %d\\n", m); f->Close(); continue; }}
    double lim; float q;
    t->SetBranchAddress("limit", &lim);
    t->SetBranchAddress("quantileExpected", &q);
    int n = t->GetEntries();
    printf("M %d", m);
    for (int i = 0; i < n; ++i) {{
      t->GetEntry(i);
      printf(" %.4f:%.6f", (double)q, lim);
    }}
    printf("\\n");
    f->Close();
  }}
}}"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".C", delete=False) as tmp:
        tmp.write(macro)
        macro_path = tmp.name
    try:
        res = subprocess.run(
            ["root", "-l", "-b", "-q", macro_path],
            capture_output=True, text=True, timeout=180,
        )
    finally:
        os.unlink(macro_path)
    if res.returncode != 0:
        sys.stderr.write(res.stderr)
    out = {}
    QREF = {0.025: "-2", 0.16: "-1", 0.5: "med", 0.84: "+1", 0.975: "+2", -1.0: "obs"}
    for line in res.stdout.splitlines():
        if not line.startswith("M "):
            continue
        toks = line.split()
        mA = int(toks[1])
        d = {k: None for k in ("-2", "-1", "med", "+1", "+2", "obs")}
        for pair in toks[2:]:
            qs, ls = pair.split(":")
            q, lim = float(qs), float(ls)
            for qref, key in QREF.items():
                if abs(q - qref) < 0.02:
                    d[key] = lim
                    break
        out[mA] = d
    return out


def collect():
    limits = read_limits()
    rows = []
    for mA in MASSES:
        env = parse_envelope(mA)
        truth = [p for p, _, t in env if t == 1]
        l = limits.get(mA, {})
        rows.append({
            "mA": mA,
            "n_env": len(env),
            "best_pdf": env[0][0] if env else "",
            "truth_members": "+".join(truth) if truth else "",
            "neg2sigma": l.get("-2"),
            "expected": l.get("med"),
            "pos2sigma": l.get("+2"),
        })
    return rows


def fmt(x, w=9, p=4):
    if x is None or x == "" or x == "None":
        return f"{'-':>{w}}"
    try:
        return f"{float(x):>{w}.{p}f}"
    except (TypeError, ValueError):
        return f"{'-':>{w}}"


def print_table(rows, title=""):
    if title:
        print(f"\n=== {title} ===")
    hdr = (f"{'mA':>3} {'n':>2} {'best':>10} {'truth_members':<32} "
           f"{'-2sig':>9} {'exp':>9} {'+2sig':>9}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['mA']:>3} {r['n_env']:>2} {r['best_pdf']:>10} "
              f"{r['truth_members']:<32} "
              f"{fmt(r['neg2sigma'])} {fmt(r['expected'])} {fmt(r['pos2sigma'])}")


def write_csv(rows, path):
    keys = ["mA", "n_env", "best_pdf", "truth_members",
            "neg2sigma", "expected", "pos2sigma"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in keys})


def read_csv(path):
    with open(path) as f:
        return [{k: v for k, v in r.items()} for r in csv.DictReader(f)]


def diff(old_rows, new_rows):
    old = {int(r["mA"]): r for r in old_rows}
    new = {int(r["mA"]): r for r in new_rows}
    print("\n=== diff (new vs old) ===")
    hdr = (f"{'mA':>3} {'truth_old':<28} {'truth_new':<28} "
           f"{'exp_old':>9} {'exp_new':>9} {'d_exp%':>8} "
           f"{'-2s_old':>9} {'-2s_new':>9}")
    print(hdr)
    print("-" * len(hdr))
    for m in sorted(set(old) | set(new)):
        o = old.get(m, {})
        n = new.get(m, {})
        to = o.get("truth_members", "") or "-"
        tn = n.get("truth_members", "") or "-"
        try:
            ex_o = float(o["expected"])
            ex_n = float(n["expected"])
            dpc = (ex_n - ex_o) / ex_o * 100.0
            ex_o_s = f"{ex_o:9.4f}"
            ex_n_s = f"{ex_n:9.4f}"
            dpc_s = f"{dpc:+7.1f}%"
        except (KeyError, TypeError, ValueError):
            ex_o_s = f"{'-':>9}"
            ex_n_s = f"{'-':>9}"
            dpc_s = f"{'-':>8}"
        flag = "  <-- env CHANGE" if to != tn else ""
        try:
            if abs(float(dpc_s.strip().rstrip('%'))) >= 20.0:
                flag += "  jump>=20%"
        except ValueError:
            pass
        print(f"{m:>3} {to:<28} {tn:<28} {ex_o_s} {ex_n_s} {dpc_s} "
              f"{fmt(o.get('neg2sigma'))} {fmt(n.get('neg2sigma'))}{flag}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--save", help="write current snapshot to CSV")
    ap.add_argument("--compare", help="diff current snapshot vs this CSV")
    args = ap.parse_args()

    rows = collect()
    print_table(rows, "current state")
    if args.save:
        write_csv(rows, args.save)
        print(f"\nwrote snapshot -> {args.save}")
    if args.compare:
        if not os.path.exists(args.compare):
            sys.exit(f"baseline not found: {args.compare}")
        diff(read_csv(args.compare), rows)


if __name__ == "__main__":
    main()
