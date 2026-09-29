"""Does the targeted fetch miss distal-junction reads? DJ reads in the fetch against the whole-file scan, per genome.
usage: python3 analysis/dj/fetch_capture.py [--samples analysis/dj/selection.tsv]   (from the repository root)"""
import argparse, csv, gzip, json, os
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--samples", help="TSV with a `sample` column (default: every genome)"); a = ap.parse_args()
samples = [r["sample"] for r in csv.DictReader(open(a.samples), delimiter="\t")] if a.samples else sorted(f.split(".")[0] for f in os.listdir("counts_scan"))
r = []
for s in samples:
    fs, ff = f"counts_scan/{s}.json.gz", f"counts_fetch/{s}.json.gz"
    if not (os.path.exists(fs) and os.path.exists(ff)):
        continue
    S = {c["name"]: c["reads"] for c in json.load(gzip.open(fs, "rt"))["classes"]}
    F = {c["name"]: c["reads"] for c in json.load(gzip.open(ff, "rt"))["classes"]}
    if S.get("DJ") and F.get("DJ"):
        r.append(F["DJ"] / S["DJ"])
r = np.array(r)
print(f"{len(r)} genomes; DJ reads fetch / scan: median {np.median(r):.4f}, mean {r.mean():.4f}, min {r.min():.4f}, max {r.max():.4f}")
