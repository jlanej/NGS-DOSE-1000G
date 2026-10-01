"""Three checks of the calls from inside the measurement: places in random order, lower depth, and the fetch.

  1. With the places of every chromosome in random order (a window's pieces kept together), a real stretch is scattered
     and cannot be found: what the chain still calls is chance. Clean genomes, with the windows and without.
  2. Counting noise as at a lower depth is added to clean genomes with the windows, read against the model learned at
     the cohort's depth: the genome's own noise factor must widen the errors, not make calls.
  3. The fetch-mode counts of the same files (the regions the bundle held before the windows), read against the
     scan's model, set against the scan's reading of the same regions.

usage (from the repository's root; about fifteen minutes): python3 analysis/karyotype/checks.py [--cache cache/scan] [--fetch cache/fetch]
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, "analysis/karyotype")
from common import K, _gather, aligned, cohort, mad, shuffled_places  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="cache/scan")
ap.add_argument("--fetch", default="cache/fetch")
ap.add_argument("--seed", type=int, default=5)
a = ap.parse_args()
C = cohort(a.cache)
model, samples, win = C["model"], C["samples"], C["with_windows"]
Ym, tab, clean = aligned(C)
rng = np.random.default_rng(a.seed)

for label, group, reps in (("with the windows", [i for i in clean if samples[i] in win], 12), ("without", [i for i in clean if samples[i] not in win][:1500], 1)):
    if not group:
        continue
    present = np.isfinite(Ym[group[0]])
    found, n = {"whole": 0, "stretch": 0}, 0
    for _ in range(reps):
        t = shuffled_places(tab, model, present, rng)
        for i in group:
            n += 1
            for e in K.read(model, model.names, Ym[i], tab=t).events:
                found["whole" if e.span == "whole" else "stretch"] += 1
    print(f"1. {label}: {len(group):,} clean genomes x {reps} random orders = {n:,} readings: {found['stretch']} stretches and {found['whole']} whole chromosomes called")

group = [i for i in clean if samples[i] in win]
if group:
    length = np.where(tab.start >= 0, tab.end - tab.start, 12000).astype(float)
    depth0 = 37.2
    print(f"2. {len(group)} clean genomes with the windows at lower depth (counting noise added), the model unchanged:")
    for depth in (30, 20, 15, 10, 6, 4):
        extra = np.where(tab.kind == "A", 1.0, 2.0) * 300.0 / length * (1 / depth - 1 / depth0)
        ev, g, z, se = 0, [], [], []
        for i in group:
            rd = K.read(model, model.names, Ym[i] + rng.normal(0, np.sqrt(np.maximum(extra, 0))))
            ev += len(rd.events)
            g.append(rd.noise)
            z += [rd.chromosomes[c].z for c in K.AUTOSOMES if rd.chromosomes[c].z is not None]
            se.append(np.median([rd.chromosomes[c].se for c in K.AUTOSOMES]))
        print(f"   {depth:>2}x: noise factor {np.median(g):.2f}; z of the autosomes, robust SD {mad(np.array(z)):.2f}; events {ev}; SE of a typical autosome {np.median(se):.4f} copies")

fetch = _gather(a.fetch)
by = {s: i for i, s in enumerate(samples)}
same, n = 0, 0
for s, v in fetch.items():
    if s not in by or v is None:
        continue
    rf = K.read(model, *v)
    i = by[s]
    pos = {q: k for k, q in enumerate(C["names"])}
    keep = [q for q in v[0] if q in pos]
    rs = K.read(model, keep, C["Y"][i][[pos[q] for q in keep]])            # the scan, on the regions the fetch holds
    n += 1
    same += rf is not None and rs is not None and rf.karyotype() == rs.karyotype()
print(f"3. read from the fetch-mode counts against the cohort's model: {same:,} of {n:,} karyotypes are written as from the scan's counts of the same regions")
