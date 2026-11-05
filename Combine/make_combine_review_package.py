#!/usr/bin/env python3
"""Self-contained datacard package for the Combine review.

The production datacards point to their shape files with absolute AFS paths, which
CombineHarvester (ValidateDatacards.py) cannot parse and which a reviewer cannot read.
This copies, for every m_a, the cross-section card (no theory uncertainties on the
Higgs boson production cross section) and the Wilson card (with them), together with
every shape file they use, into one directory and rewrites the shapes lines to
relative paths. The production files are only read.

Output layout (<out>):
  mA<m>_xsec.txt, mA<m>_wilson.txt
  shapes/signal_ele/*.root, shapes/signal_mu/*.root, shapes/background/mA<m>/*.root
  MANIFEST.txt   source card -> packaged card, and adler32 of every copied file
"""
import argparse, os, shutil, zlib

HERE = os.path.dirname(os.path.abspath(__file__))
FF = os.path.dirname(HERE)
CARDS = os.path.join(FF, "Datacard", "output_Datacard_leptons")
CARD_NAME = {"xsec": "%d_pruned_datacard_leptons_xsec.txt", "wilson": "%d_pruned_datacard_leptons.txt"}


def adler32(path):
    v = 1
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            v = zlib.adler32(chunk, v)
    return "%08x" % (v & 0xffffffff)


def packaged_path(src):
    if "/Signal/outdir_ele/" in src:
        return os.path.join("shapes", "signal_ele", os.path.basename(src))
    if "/Signal/outdir_mu/" in src:
        return os.path.join("shapes", "signal_mu", os.path.basename(src))
    if "/fit_results_run3/" in src:
        m = src.split("/fit_results_run3/")[1].split("/")[0]
        return os.path.join("shapes", "background", "mA" + m, os.path.basename(src))
    raise SystemExit("unexpected shape file location: " + src)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--masses", default="1-30")
    args = ap.parse_args()
    lo, hi = (int(x) for x in args.masses.split("-"))
    if os.path.exists(args.out) and os.listdir(args.out):
        raise SystemExit("output directory is not empty: " + args.out)
    os.makedirs(args.out, exist_ok=True)

    copied, manifest = {}, []
    for m in range(lo, hi + 1):
        for kind, pattern in CARD_NAME.items():
            src_card = os.path.join(CARDS, pattern % m)
            out_lines = []
            for line in open(src_card):
                tok = line.split()
                if tok and tok[0] == "shapes":
                    src = tok[3]
                    rel = packaged_path(src)
                    if rel not in copied:
                        dst = os.path.join(args.out, rel)
                        os.makedirs(os.path.dirname(dst), exist_ok=True)
                        shutil.copy2(src, dst)
                        copied[rel] = (src, adler32(dst))
                        assert copied[rel][1] == adler32(src), src
                    line = line.replace(src, rel, 1)
                out_lines.append(line)
            name = "mA%d_%s.txt" % (m, kind)
            open(os.path.join(args.out, name), "w").writelines(out_lines)
            manifest.append("%s <- %s" % (name, src_card))
    with open(os.path.join(args.out, "MANIFEST.txt"), "w") as f:
        f.write("\n".join(manifest) + "\n\n")
        for rel in sorted(copied):
            f.write("%s %s <- %s\n" % (copied[rel][1], rel, copied[rel][0]))
    print("[package] %d cards, %d shape files -> %s" % (len(manifest), len(copied), args.out))


if __name__ == "__main__":
    main()
