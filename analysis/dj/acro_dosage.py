"""Are the values between steps mosaic losses or gains of an acrocentric chromosome in the cell line? Per-chromosome dosage
from each genome's 800 control regions (the estimate's control_qc), against its distal-junction step.
usage: python3 analysis/dj/acro_dosage.py [--samples analysis/dj/selection.tsv] [--cache cache/scan]"""
import argparse, csv, gzip, json, os
import numpy as np
from ngsdose import resources
ap = argparse.ArgumentParser(); ap.add_argument("--samples", default="analysis/dj/selection.tsv"); ap.add_argument("--cache", default="cache/scan"); a = ap.parse_args()
from pathlib import Path
root = Path(resources.default_bundle()); B = resources.Bundle(root)
bed = [l.rstrip("\n").split("\t") for l in open(root / B.meta.get("controls_bed", "controls.bed")) if l.strip() and not l.startswith("#")]
chrom = np.array([b[0] for b in bed if b[3] == "control"])
acro = ["chr13", "chr14", "chr15", "chr21", "chr22"]
step = {r["sample"]: float(r["DJ.step"]) for r in csv.DictReader(open("docs/data/cohort.tsv"), delimiter="\t") if r.get("DJ.step") not in ("", "NA")}
rows = []
for r in csv.DictReader(open(a.samples), delimiter="\t"):
    s = r["sample"]; f = f"{a.cache}/{s}.estimate.json.gz"
    if s not in step or not os.path.exists(f):
        continue
    lr = np.array([np.nan if x is None else x for x in json.load(gzip.open(f, "rt"))["control_qc"]["region_log_ratio"]], float)
    rows.append((s, step[s], {c: 2 * float(np.exp(np.nanmedian(lr[chrom == c]))) for c in acro}))
X = np.array([[st] + [d[c] for c in acro] for _, st, d in rows]); off = X[:, 1:].sum(axis=1) - 10
print(f"{len(rows)} genomes; corr(DJ step, sum of acrocentric dosage offsets) = {np.corrcoef(X[:, 0], off)[0, 1]:+.3f}; per-chromosome dosage SD: "
      + ", ".join(f"{c} {X[:, 1 + i].std():.3f}" for i, c in enumerate(acro)))
print("genomes with an acrocentric off by more than 0.15 copies (chr22 has few control regions):")
for s, st, d in sorted(rows, key=lambda t: -max(abs(v - 2) for v in t[2].values())):
    o = {c: round(v - 2, 2) for c, v in d.items() if abs(v - 2) > 0.15}
    if o:
        print(f"  {s}: DJ step {st:+.2f}, offsets {o}")
