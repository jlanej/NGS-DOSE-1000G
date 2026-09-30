"""The offsets of the polymorphic intervals' efficiencies, anchored on assemblies.

In an interval that a common deletion removes from some junction copies, the cohort's median genome does not hold the
interval at its level, so the cohort's median gives the interval's windows the wrong efficiency. The offset is the
median, over the genomes whose HPRC assembly resolves the junction, of log(calibrated estimate of the interval /
the assembly's copies of it), the estimate taken on the genome's own scale against its assembly (its core level over
the assembly's core copies). Run with the class's `polymorphic` rules emptied, so that what is measured is the
cohort-median calibration itself.

usage (from the results repository's root): python3 analysis/dj/polymorphic_offsets.py [--cache cache/scan] [--assemblies meta/dj_hprc]
"""
import argparse, csv, os, sys
import numpy as np
from ngsdose import cohort, resources
from ngsdose.tables import load_result

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="cache/scan"); ap.add_argument("--assemblies", default="meta/dj_hprc"); ap.add_argument("--table", default="docs/data/dj_hprc.tsv")
ap.add_argument("--leave-out", nargs="*", default=["HG00673"], help="resolved genomes whose whole level disagrees with their assembly")
a = ap.parse_args()
B = resources.Bundle()
rules = B.calibration()
intervals = [tuple(p["interval"]) for p in rules["DJ"]["polymorphic"]] + [(0, 5000), (23000, 30000), (128000, 137000), (160000, 165000), (225000, 230000), (300000, 340000)]
rules["DJ"] = dict(rules["DJ"], polymorphic=[], segments=False)
files = sorted(f for f in os.listdir(a.cache) if f.endswith(".estimate.json.gz"))
prof = {}
cohort.cohort_table((load_result(f"{a.cache}/{f}") for f in files), B.anchors(), rules=rules, profiles=prof, n_control_pcs=0, log=lambda m: print(m, file=sys.stderr))
P = prof["DJ"]; idx = {s: i for i, s in enumerate(P["samples"])}; starts = np.array(P["start"]); cn = P["cn"].astype(float)
level = np.exp(np.nanmedian(np.log(cn[:, P["level"]]), axis=1))
H = {}
for r in csv.DictReader(open(f"{a.assemblies}/haplotypes.tsv"), delimiter="\t"):
    H.setdefault(r["sample"], []).append([float(r[f"b{k}"]) if r[f"b{k}"] != "" else np.nan for k in range(80)])
T = {r["sample"]: r for r in csv.DictReader(open(a.table), delimiter="\t")}
res = [s for s in H if len(H[s]) == 2 and T.get(s, {}).get("resolved") == "True" and s not in a.leave_out and s in idx]
print(f"{len(res)} genomes with a resolved assembly")
for a0, b0 in intervals:
    m = (starts >= a0) & (starts < b0) & np.isfinite(cn).any(axis=0)
    ks = [k for k in range(80) if a0 <= k * 5000 < b0]
    v = []
    for s in res:
        i = idx[s]; asm = np.nanmedian(np.sum(H[s], axis=0)[ks])
        if np.isfinite(asm) and asm > 0 and m.any():
            own = level[i] / float(T[s]["assembly_core"])
            v.append(np.log(np.nanmedian(cn[i, m]) / own / asm))
    v = np.array(v)
    print(f"{a0 // 1000:3d}-{b0 // 1000:3d} kb  {int(m.sum()):3d} windows  offset {np.median(v):+.3f}  (quartiles {np.percentile(v, 25):+.3f} .. {np.percentile(v, 75):+.3f}; n = {len(v)})")
