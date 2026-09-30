"""The polymorphic intervals by release batch and by generation.

In the trios the children read the distal intervals of the junction lower than the mean of their parents, which
Mendel does not allow. Most children were sequenced in the later release batch (698 genomes) and most parents in
the earlier (2,504), so the difference is either the batch's or the generation's. About a hundred parents sit in the
later batch and tell the two apart: this prints each interval's value (the median of its calibrated windows, less
the genome's level) by batch, for all genomes, for parents and for children.

usage (from the results repository's root): python3 analysis/dj/intervals_by_batch.py [--cache cache/scan]
"""
import argparse, csv, os, sys, warnings
import numpy as np
from ngsdose import cohort, resources
from ngsdose.tables import load_result

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="cache/scan"); ap.add_argument("--qc", default="meta/ngspca_sample_qc.tsv")
ap.add_argument("--ped", default="meta/20130606_g1k_3202_samples_ped_population.txt")
a = ap.parse_args()
warnings.simplefilter("ignore", RuntimeWarning)
B = resources.Bundle()
rules = B.calibration()
intervals = [tuple(p["interval"]) for p in rules["DJ"]["polymorphic"]] + [(40000, 60000), (240000, 262000), (300000, 340000)]
rules["DJ"] = dict(rules["DJ"], segments=False)
files = sorted(f for f in os.listdir(a.cache) if f.endswith(".estimate.json.gz"))
prof = {}
cohort.cohort_table((load_result(f"{a.cache}/{f}") for f in files), B.anchors(), rules=rules, profiles=prof, n_control_pcs=0, log=lambda m: print(m, file=sys.stderr))
P = prof["DJ"]; names = list(P["samples"]); idx = {s: i for i, s in enumerate(names)}; starts = np.array(P["start"]); cn = P["cn"].astype(float)
level = np.exp(np.nanmedian(np.log(cn[:, P["level"]]), axis=1))
has = np.isfinite(cn).any(axis=0)
qc = {r["SAMPLE_ID"]: r for r in csv.DictReader(open(a.qc), delimiter="\t")}
ped = {p["SampleID"]: p for p in csv.DictReader(open(a.ped), delimiter=" ")}
trio = {s: p for s, p in ped.items() if s in idx and p["FatherID"] in idx and p["MotherID"] in idx}
parents = {p["FatherID"] for p in trio.values()} | {p["MotherID"] for p in trio.values()}
batch = np.array([qc.get(s, {}).get("RELEASE_BATCH", "") for s in names])
is_par = np.array([s in parents for s in names]); is_kid = np.array([s in trio and s not in parents for s in names])
groups = [("all, batch 2504", batch == "2504"), ("all, batch 698", batch == "698"), ("parents, batch 2504", is_par & (batch == "2504")),
          ("parents, batch 698", is_par & (batch == "698")), ("children, batch 698", is_kid & (batch == "698"))]
print("interval\twindows\t" + "\t".join(f"{g} (n {int(m.sum())})" for g, m in groups) + "\tlater less earlier, all\tlater less earlier, parents")
for a0, b0 in intervals:
    m = (starts >= a0) & (starts < b0) & has
    d = np.nanmedian(cn[:, m], axis=1) - level
    cell = lambda sel, d=d: f"{d[sel].mean():+.3f} ± {d[sel].std(ddof=1) / sel.sum() ** 0.5:.3f}"
    diff = lambda x, y, d=d: f"{d[x].mean() - d[y].mean():+.3f} ± {(d[x].var(ddof=1) / x.sum() + d[y].var(ddof=1) / y.sum()) ** 0.5:.3f}"
    print(f"{a0 // 1000}-{b0 // 1000} kb\t{int(m.sum())}\t" + "\t".join(cell(sel) for _, sel in groups)
          + f"\t{diff(groups[1][1], groups[0][1])}\t{diff(groups[3][1], groups[2][1])}")
