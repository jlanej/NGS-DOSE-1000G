"""Is a genome's reading its own, or the cohort's? A model learned on the 2,504 genomes of the first release is applied,
unchanged, to the 698 added later (another sequencing batch), and set against what the whole cohort's own model reads.

usage (from the repository's root; about eight minutes): python3 analysis/karyotype/fixed_model.py
"""
import csv
import sys

import numpy as np

sys.path.insert(0, "analysis/karyotype")
from common import ARMS, GC, KAR, K, cohort, mad  # noqa: E402

C = cohort()
samples, names, Y, full = C["samples"], C["names"], C["Y"], C["readings"]
qc = {r["SAMPLE_ID"]: r for r in csv.DictReader(open("meta/ngspca_sample_qc.tsv"), delimiter="\t")}
batch = np.array([qc.get(s, {}).get("RELEASE_BATCH", "2504") for s in samples])
first, later = np.flatnonzero(batch == "2504"), np.flatnonzero(batch == "698")
_, m1, _ = K.cohort([(names, Y[i]) for i in first], ARMS, rules=KAR.get("rules"), gc=GC)
rb, _, _ = K.cohort([(names, Y[i]) for i in later], ARMS, model=m1, fit_own=False)
same = sum(1 for j, i in enumerate(later) if rb[j].karyotype() == full[i].karyotype())
print(f"model learned on the {len(first):,} genomes of the first release, applied to the {len(later)} added later: {same} of {len(later)} karyotypes are "
      "written as the whole cohort's model writes them")
for j, i in enumerate(later):
    if rb[j].karyotype() != full[i].karyotype():
        print(f"   {samples[i]}: {full[i].karyotype()} | {rb[j].karyotype()}")
d = {c: np.array([rb[j].chromosomes[c].copies - full[i].chromosomes[c].copies for j, i in enumerate(later)]) for c in K.AUTOSOMES}
se = {c: np.median([full[i].chromosomes[c].se for i in later]) for c in K.AUTOSOMES}
print(f"levels, fixed model less the cohort's own: robust SD {min(mad(d[c]) for c in d):.4f} to {max(mad(d[c]) for c in d):.4f} copies over the autosomes; "
      f"as a share of the level's standard error {min(mad(d[c]) / se[c] for c in d):.2f} to {max(mad(d[c]) / se[c] for c in d):.2f}")
