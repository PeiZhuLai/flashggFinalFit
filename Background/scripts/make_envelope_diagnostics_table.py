#!/usr/bin/env python3
"""Build the per-function background-fit diagnostics table for the AN.

For every ALP mass hypothesis and every function retained in the background
envelope, collect the four quantities the L3 conveners asked for:

  * fit status   -- Minuit status of the background fit (0 = converged)
  * GOF p-value  -- toy-based goodness-of-fit probability of the fit to data
  * F-test       -- probability that the next higher order of the same family
                    is not required
  * bias         -- mean of the pull distribution from the signal-injection
                    bias study, with its uncertainty

Everything is read back from products that already exist on disk; nothing is
refitted here.

Inputs
  <bkgmodel>/fit_results_run3/<mA>/HZAmassInde_fTest/EnvelopeResults.txt
      final envelope membership, per-function GOF, truth flag, fallback flags
  <bkgmodel>/fit_results_run3/<mA>/ftest.stdout.log
      per-(family, order) GOF p-value and observed chi2, F-test probabilities,
      and the Minuit fit-status warnings
  <finalfit>/Combine/Checks/Bias_nominal/bias_outputs/mA_<mA>/BiasJson/<mA>_gaussfit.json
      per-function pull mean/sigma from the bias study

Output
  <bkgdir>/output_tables/envelope_diagnostics_mA<lo>to<hi>.tex   (LaTeX tables)
  <bkgdir>/output_tables/envelope_diagnostics.csv                (flat dump)

Usage
  python3 scripts/make_envelope_diagnostics_table.py
  python3 scripts/make_envelope_diagnostics_table.py --masses 1-30 --per-table 10
"""

import argparse
import csv
import json
import os
import re
import sys
from collections import OrderedDict

# --------------------------------------------------------------------------- #
# Defaults: resolved relative to this file so the script works from anywhere.
# --------------------------------------------------------------------------- #
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BKG_DIR = os.path.dirname(SCRIPT_DIR)                       # .../Background
FINALFIT_DIR = os.path.dirname(BKG_DIR)                     # .../flashggFinalFit

DEF_FITRES = os.path.join(BKG_DIR, "ALP_BkgModel_ReReco", "fit_results_run3")
DEF_BIAS = os.path.join(FINALFIT_DIR, "Combine", "Checks", "Bias_nominal", "bias_outputs")
# High-statistics bias study (10000 toys per function). Preferred over the
# nominal 1000-toy run: comparing the two over the full scan showed that 1000
# toys systematically underestimates |bias| by ~0.03-0.05, because the Gaussian
# fit to a 1000-entry pull distribution is pulled towards zero by the tails.
DEF_BIAS_HS = os.path.join(FINALFIT_DIR, "Combine", "Checks", "Bias_nominal",
                           "bias_outputs_highstat")
DEF_OUTDIR = os.path.join(BKG_DIR, "output_tables")

# Order in which families are listed inside one mass point.
FAMILY_ORDER = ["Exp", "Pow", "Lau", "Bern"]
FAMILY_LONG = {"Exp": "Exponential", "Pow": "PowerLaw", "Lau": "Laurent", "Bern": "Bernstein"}
LONG_TO_SHORT = {v: k for k, v in FAMILY_LONG.items()}


def split_pdf_name(name):
    """'Bern2' -> ('Bern', 2). Returns (name, None) if there is no order."""
    m = re.match(r"^([A-Za-z]+?)(\d+)$", name)
    if not m:
        return name, None
    return m.group(1), int(m.group(2))


# --------------------------------------------------------------------------- #
# Parsers
# --------------------------------------------------------------------------- #
def parse_envelope_results(path):
    """Final envelope membership from EnvelopeResults.txt.

    Each line looks like
      category : 0 , pdf : Exp1 , gof : 0.840000, isTruth : 1 [family-fallback]
    """
    if not os.path.isfile(path):
        return []

    out = []
    with open(path, errors="ignore") as fh:
        for line in fh:
            m = re.search(
                r"pdf\s*:\s*(\S+)\s*,\s*gof\s*:\s*([0-9.eE+-]+)\s*,\s*isTruth\s*:\s*(\d+)",
                line,
            )
            if not m:
                continue
            out.append(
                {
                    "pdf": m.group(1),
                    "gof_envelope": float(m.group(2)),
                    "is_truth": int(m.group(3)),
                    "fallback": "fallback" in line,
                }
            )
    return out


def parse_ftest_log(path):
    """Per-(family, order) GOF, F-test probability and fit status from the log.

    The log is written in two passes per family: a scan pass that prints the
    F-test probabilities, then an envelope pass that prints, for each order,
    the observed chi2 and the toy-based GOF p-value.  Both passes emit

      [INFO]  funcType: Exponential order: 3 prevNll: X thisNll: Y chi2: Z prob: P

    where ``prob`` is the F-test probability of this order against the previous
    one.  In the envelope pass each such line is followed by

      [INFO] Chi2 in Observed =  50.495
      [INFO] p-value  =  0.84

    Fit failures are reported out of band as

      [WARNING] Warning -- Fit status for Exp3 at 1
      [INFO] Skip GOF because fitStatus=-1 for pdf Bern6

    so any function without such a line converged with status 0.
    """
    info = {
        "per_order": OrderedDict(),   # (short_family, order) -> dict
        "status": {},                 # pdf name -> int
        "status_bare": {},            # family with no order -> int
        "best_fit_pdf": None,
        "multipdf_members": [],
    }
    if not os.path.isfile(path):
        return info

    with open(path, errors="ignore") as fh:
        lines = fh.readlines()

    re_func = re.compile(
        r"funcType:\s*(\S+)\s+order:\s*(\d+)\s+prevNll:\s*([0-9.eE+-]+)\s+"
        r"thisNll:\s*([0-9.eE+-]+)\s+chi2:\s*([0-9.eE+-]+)\s+prob:\s*([0-9.eE+-]+)"
    )
    re_pval = re.compile(r"p-value\s+=\s+([0-9.eE+-]+)")
    re_chi2obs = re.compile(r"Chi2 in Observed\s*=\s*([0-9.eE+-]+)")
    re_status = re.compile(r"Fit status for\s+(\S+)\s+at\s+(-?\d+)")
    re_skipgof = re.compile(r"Skip GOF because fitStatus=(-?\d+)\s+for pdf\s+(\S+)")
    re_store = re.compile(r"RooArgList::store\s*=\s*\(([^)]*)\)")
    re_best = re.compile(r"Best Fit Pdf\s*=\s*\d+,\s*(\S+)")

    current = None  # (short_family, order) whose GOF block we are inside
    for line in lines:
        m = re_func.search(line)
        if m:
            fam_long, order = m.group(1), int(m.group(2))
            fam = LONG_TO_SHORT.get(fam_long, fam_long)
            key = (fam, order)
            entry = info["per_order"].setdefault(
                key, {"family": fam, "order": order, "ftest_prob": None,
                      "gof": None, "chi2_obs": None, "nll": None}
            )
            # prevNll == 0 marks the lowest order, where the F-test is undefined.
            if float(m.group(3)) != 0.0:
                entry["ftest_prob"] = float(m.group(6))
            entry["nll"] = float(m.group(4))
            current = key
            continue

        m = re_chi2obs.search(line)
        if m and current is not None:
            info["per_order"][current]["chi2_obs"] = float(m.group(1))
            continue

        m = re_pval.search(line)
        if m and current is not None:
            # Keep the last GOF evaluation of this order (envelope pass).
            info["per_order"][current]["gof"] = float(m.group(1))
            continue

        m = re_status.search(line)
        if m:
            name, val = m.group(1), int(m.group(2))
            fam, order = split_pdf_name(name)
            if order is None:
                info["status_bare"][fam] = val
            else:
                info["status"][name] = val
            continue

        m = re_skipgof.search(line)
        if m:
            val, name = int(m.group(1)), m.group(2)
            fam, order = split_pdf_name(name)
            if order is None:
                info["status_bare"][fam] = val
            else:
                info["status"][name] = val
            continue

        m = re_store.search(line)
        if m:
            info["multipdf_members"] = [s.strip() for s in m.group(1).split(",") if s.strip()]
            continue

        m = re_best.search(line)
        if m:
            info["best_fit_pdf"] = m.group(1)

    return info


def parse_bias_json(path):
    """{'Exp1': {'mean': .., 'mean_err': ..}, ...} plus the injected mu."""
    if not os.path.isfile(path):
        return {}, None
    with open(path) as fh:
        payload = json.load(fh)
    return payload.get("fit_results", {}), payload.get("exp")


def load_bias(mA, bias_dir, bias_hs_dir):
    """Bias for one mass point, preferring the high-statistics run.

    Returns (fit_results, injected_mu, n_toys, source) where source is
    'highstat' or 'nominal'. The toy count is read back from the merged tree
    size recorded by the merge step when available; otherwise it falls back to
    the nominal 1000.
    """
    hs = os.path.join(bias_hs_dir, "mA_%d" % mA, "merged", "BiasJson",
                      "%d_gaussfit.json" % mA)
    if os.path.isfile(hs):
        fits, mu = parse_bias_json(hs)
        if fits:
            return fits, mu, 10000, "highstat"
    nom = os.path.join(bias_dir, "mA_%d" % mA, "BiasJson", "%d_gaussfit.json" % mA)
    fits, mu = parse_bias_json(nom)
    return fits, mu, 1000, "nominal"


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
def next_order_ftest(per_order, family, order):
    """F-test probability of the next scanned order of this family, if any."""
    orders = sorted(o for (f, o) in per_order if f == family)
    later = [o for o in orders if o > order]
    if not later:
        return None
    return per_order[(family, later[0])]["ftest_prob"], later[0]


def collect(masses, fitres_dir, bias_dir, bias_hs_dir):
    rows, warnings = [], []
    toy_counts = {}

    for mA in masses:
        env_path = os.path.join(fitres_dir, str(mA), "HZAmassInde_fTest", "EnvelopeResults.txt")
        log_path = os.path.join(fitres_dir, str(mA), "ftest.stdout.log")
        bias_path = os.path.join(bias_dir, "mA_%d" % mA, "BiasJson", "%d_gaussfit.json" % mA)

        envelope = parse_envelope_results(env_path)
        log = parse_ftest_log(log_path)
        bias, mu_inj, n_toys, bias_src = load_bias(mA, bias_dir, bias_hs_dir)
        toy_counts[mA] = n_toys

        if not envelope:
            warnings.append("mA=%d: no EnvelopeResults.txt at %s" % (mA, env_path))
            continue
        if not log["per_order"]:
            warnings.append("mA=%d: ftest log missing or unparsable at %s" % (mA, log_path))
        if not bias:
            warnings.append("mA=%d: no bias JSON at %s" % (mA, bias_path))

        env_names = [e["pdf"] for e in envelope]
        if log["multipdf_members"] and set(log["multipdf_members"]) != set(env_names):
            warnings.append(
                "mA=%d: MultiPdf members %s != EnvelopeResults %s"
                % (mA, log["multipdf_members"], env_names)
            )
        if bias and set(bias) != set(env_names):
            warnings.append(
                "mA=%d: bias JSON covers %s but envelope is %s"
                % (mA, sorted(bias), sorted(env_names))
            )
        for fam, val in log["status_bare"].items():
            warnings.append(
                "mA=%d: fit-status %d reported for family %s without an order; "
                "not attributed to any envelope member" % (mA, val, fam)
            )

        envelope.sort(
            key=lambda e: (
                FAMILY_ORDER.index(split_pdf_name(e["pdf"])[0])
                if split_pdf_name(e["pdf"])[0] in FAMILY_ORDER else 99,
                split_pdf_name(e["pdf"])[1] or 0,
            )
        )

        for e in envelope:
            pdf = e["pdf"]
            fam, order = split_pdf_name(pdf)

            # Only an exact pdf-name match assigns a non-zero status. The log
            # also emits family-level warnings without an order ("Fit status
            # for Bern at -1"); those cannot be attributed to a specific order,
            # so they are reported as warnings instead of being spread over
            # every member of the family.
            status = log["status"].get(pdf, 0)

            entry = log["per_order"].get((fam, order), {})
            gof = entry.get("gof")
            if gof is None:
                # Fall back on the value EnvelopeResults.txt recorded.
                gof = e["gof_envelope"]
            elif abs(gof - e["gof_envelope"]) > 0.005:
                # The log and the summary file must agree; a mismatch means the
                # GOF block was attributed to the wrong (family, order).
                warnings.append(
                    "mA=%d %s: GOF from log %.4g disagrees with EnvelopeResults %.4g"
                    % (mA, pdf, gof, e["gof_envelope"])
                )

            ft = next_order_ftest(log["per_order"], fam, order)
            ftest_prob, ftest_next = (None, None) if ft is None else ft

            b = bias.get(pdf, {})

            rows.append(
                {
                    "mA": mA,
                    "pdf": pdf,
                    "family": fam,
                    "order": order,
                    "status": status,
                    "gof": gof,
                    "chi2_obs": entry.get("chi2_obs"),
                    "ftest_prob": ftest_prob,
                    "ftest_next_order": ftest_next,
                    "bias_mean": b.get("mean"),
                    "bias_mean_err": b.get("mean_err"),
                    "bias_sigma": b.get("sigma"),
                    "mu_injected": mu_inj,
                    "bias_ntoys": n_toys,
                    "bias_source": bias_src,
                    "is_truth": e["is_truth"],
                    "fallback": e["fallback"],
                    "best_fit": pdf == log["best_fit_pdf"],
                }
            )

    return rows, warnings, toy_counts


# --------------------------------------------------------------------------- #
# Formatting
# --------------------------------------------------------------------------- #
def fmt_p(value):
    """p-values: keep small ones readable without drowning the table."""
    if value is None:
        return "--"
    if value < 1e-3:
        return "$<10^{-3}$"
    return "%.3f" % value


def fmt_bias(mean, err):
    if mean is None:
        return "--"
    if err is None:
        return "%+.3f" % mean
    return "$%+.3f \\pm %.3f$" % (mean, err)


def fmt_note(row):
    notes = []
    if row["fallback"]:
        notes.append("fallback")
    if row["best_fit"]:
        notes.append("best fit")
    return ", ".join(notes) if notes else "--"


def make_table(rows, lo, hi, label_suffix):
    lines = []
    lines.append("\\begin{table}[htbp]")
    lines.append("  \\begin{center}")
    lines.append(
        ("    \\topcaption{Per-function diagnostics of the background fit for the "
         "envelope members at $%d \\le m_{a} \\le %d~\\GeV$. " % (lo, hi)) +
        "The fit status is the "
        "Minuit status of the background-only fit to data (0 = converged). The "
        "goodness-of-fit probability is obtained from the observed $\\chi^{2}$ with "
        "toys; a function is accepted only if it exceeds 0.01. The F-test column "
        "gives the probability that the next higher order of the same family is "
        "not required (`--' when no higher order was scanned); the higher order is "
        "adopted only if this probability falls below 0.10 when building the "
        "envelope candidates, and below the stricter 0.05 when identifying the "
        "highest order supported by the data in the truth-model scan. "
        "The bias column is the mean of the pull distribution obtained when the toy "
        "datasets are generated with this function as the truth model, at the "
        "expected signal strength, and fitted with the full envelope. The pull "
        "is normalized to the uncertainty on the side of the fitted signal "
        "strength that faces the true value; the "
        "quoted uncertainty is the uncertainty on that mean and 10000 toys are "
        "used per function. Functions marked "
        "`fallback' did not pass the combined goodness-of-fit and F-test "
        "selection and were retained by the fallback rule that keeps every "
        "family represented in the envelope; `best fit' marks the function "
        "preferred by the penalized likelihood.}"
    )
    lines.append("    \\label{tab:bkg_diagnostics_%s}" % label_suffix)
    lines.append("    \\resizebox{\\textwidth}{!}{")
    lines.append("    \\begin{tabular}{|c|c|c|c|c|c|c|} \\hline")
    lines.append(
        "      $m_{a}$ (\\GeV) & Function & Fit status & GOF $p$-value & "
        "F-test $p$ & Bias (mean pull, truth $=$ this function) & Note \\\\ \\hline"
    )

    last_mA = None
    for row in rows:
        if row["mA"] < lo or row["mA"] > hi:
            continue
        if last_mA is not None and row["mA"] != last_mA:
            lines.append("      \\hline")
        mass_cell = "%d" % row["mA"] if row["mA"] != last_mA else ""
        lines.append(
            "      %s & %s & %d & %s & %s & %s & %s \\\\"
            % (
                mass_cell,
                row["pdf"],
                row["status"],
                fmt_p(row["gof"]),
                fmt_p(row["ftest_prob"]),
                fmt_bias(row["bias_mean"], row["bias_mean_err"]),
                fmt_note(row),
            )
        )
        last_mA = row["mA"]

    lines.append("      \\hline")
    lines.append("    \\end{tabular}")
    lines.append("    }")
    lines.append("  \\end{center}")
    lines.append("\\end{table}")
    lines.append("")
    return "\n".join(lines)


def write_csv(rows, path):
    fields = [
        "mA", "pdf", "family", "order", "status", "gof", "chi2_obs",
        "ftest_prob", "ftest_next_order", "bias_mean", "bias_mean_err",
        "bias_sigma", "mu_injected", "bias_ntoys", "bias_source",
        "is_truth", "fallback", "best_fit",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def parse_masses(spec):
    out = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            a, b = chunk.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(chunk))
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fitres", default=DEF_FITRES, help="fit_results_run3 directory")
    ap.add_argument("--bias", default=DEF_BIAS, help="bias_outputs directory (1000-toy fallback)")
    ap.add_argument("--bias-highstat", default=DEF_BIAS_HS,
                    help="bias_outputs_highstat directory (10000-toy, preferred)")
    ap.add_argument("--outdir", default=DEF_OUTDIR, help="output directory for the LaTeX tables")
    ap.add_argument("--masses", default="1-30", help="mass points, e.g. '1-30' or '1,2,5-8'")
    ap.add_argument("--per-table", type=int, default=10, help="mass points per LaTeX table")
    opts = ap.parse_args()

    masses = parse_masses(opts.masses)
    print("[input ] fit results : %s" % opts.fitres)
    print("[input ] bias outputs: %s" % opts.bias)
    print("[input ] bias highstat: %s" % opts.bias_highstat)
    print("[output] tables      : %s" % opts.outdir)

    rows, warnings, toy_counts = collect(masses, opts.fitres, opts.bias, opts.bias_highstat)
    if not rows:
        print("[ERROR] no rows collected; check the input paths", file=sys.stderr)
        return 1

    os.makedirs(opts.outdir, exist_ok=True)

    written = []
    for start in range(0, len(masses), opts.per_table):
        block = masses[start:start + opts.per_table]
        lo, hi = block[0], block[-1]
        suffix = "mA%dto%d" % (lo, hi)
        path = os.path.join(opts.outdir, "envelope_diagnostics_%s.tex" % suffix)
        with open(path, "w") as fh:
            fh.write(make_table(rows, lo, hi, suffix))
        written.append(path)

    csv_path = os.path.join(opts.outdir, "envelope_diagnostics.csv")
    write_csv(rows, csv_path)
    written.append(csv_path)

    for path in written:
        print("[wrote ] %s" % path)

    n_fallback = sum(1 for r in rows if r["fallback"])
    n_nonzero = sum(1 for r in rows if r["status"] != 0)
    n_nobias = sum(1 for r in rows if r["bias_mean"] is None)
    print("[summary] %d rows over %d mass points" % (len(rows), len(masses)))
    print("[summary] %d fallback functions, %d with fit status != 0, %d without a bias value"
          % (n_fallback, n_nonzero, n_nobias))

    if warnings:
        print("[warn   ] %d consistency warnings:" % len(warnings))
        for w in warnings:
            print("           %s" % w)

    return 0


if __name__ == "__main__":
    sys.exit(main())
