#!/usr/bin/env python3
"""Where a parent's partial copy is not called in the child, does the child's profile hold it?

The page counts a parent's partial copy as passed on when the child's call has the parent's breakpoint. A call that
missed the child's copy would count as "not passed on". This reads the child's calibrated profile itself: the mean of
the 20-kb core blocks beyond the parent's breakpoint less the mean of those before it (the block that holds the
breakpoint left out). A child that holds the copy steps as its parent does, one copy down or up; one that does not
reads 0. The pairs are the page's (report.report.partial_copies): the parent's call has one breakpoint in the core,
the other parent is at ten with none near it, the three calls are settled.
Reads docs/data/cohort.tsv, docs/data/dj_blocks.tsv and the pedigree; run from the repository root.
"""
import csv
import statistics as st
import sys
from collections import Counter

sys.path.insert(0, ".")
from report.report import DJ_LEFT_OUT, dj_breakpoints  # noqa: E402

PED = "meta/20130606_g1k_3202_samples_ped_population.txt"
TOL = 12
rows = {r["sample"]: r for r in csv.DictReader(open("docs/data/cohort.tsv"), delimiter="\t")}
B = {r["sample"]: r for r in csv.DictReader(open("docs/data/dj_blocks.tsv"), delimiter="\t")}
ped = {p["SampleID"]: p for p in csv.DictReader(open(PED), delimiter=" ")}
core = [k for k in range(0, 400, 20) if not any(k < b and k + 20 > a for a, b in DJ_LEFT_OUT)]
settled = lambda s: s in rows and rows[s]["DJ.call"] == "settled"
bp = {s: dj_breakpoints(rows[s], 10) for s in rows if settled(s)}


def step(s, at):
    """The profile's step across `at` kb: the core blocks beyond it less those before it."""
    before = [float(B[s][f"b{k}"]) for k in core if k + 20 <= at]
    beyond = [float(B[s][f"b{k}"]) for k in core if k >= at]
    return None if len(before) < 2 or len(beyond) < 2 else st.mean(beyond) - st.mean(before)


pairs, few = [], 0
for child, p in ped.items():
    f, m = p["FatherID"], p["MotherID"]
    if not (settled(child) and settled(f) and settled(m)):
        continue
    for who, par, other in (("father", f, m), ("mother", m, f)):
        if len(bp[par]) != 1 or int(rows[other]["DJ.copies"]) != 10 or any(abs(x - bp[par][0][0]) <= 2 * TOL for x, _ in bp[other]):
            continue
        at, s = bp[par][0]
        found = any(abs(x - at) <= TOL and t * s > 0 for x, t in bp[child])
        sp, sc = step(par, at), step(child, at)
        if sp is None or sc is None:
            few += 1
            continue
        pairs.append(dict(who=who, parent=par, child=child, at=at, step=s, found=found, parent_step=sp, child_step=sc * (1 if s > 0 else -1)))

print(f"{len(pairs)} pairs read ({few} more whose breakpoint leaves fewer than two core blocks on a side)")
print(f"parents' own step across their breakpoint: median {st.median(abs(r['parent_step']) for r in pairs):.2f} copies")
for found in (True, False):
    v = [r["child_step"] for r in pairs if r["found"] == found]
    print(f"children whose call {'has' if found else 'lacks'} the breakpoint: {len(v)}; their step, signed by the parent's: median {st.median(v):+.2f} "
          f"(quartiles {st.quantiles(v, n=4)[0]:+.2f} to {st.quantiles(v, n=4)[2]:+.2f})")
held = lambda r: r["child_step"] > 0.5
print(f"by the profile (a step of more than half a copy the parent's way): passed on in {sum(map(held, pairs))} of {len(pairs)}; by the calls in {sum(r['found'] for r in pairs)}")
print(f"  the call lacks it and the profile holds it: {sum(1 for r in pairs if not r['found'] and held(r))}; the call has it and the profile does not: {sum(1 for r in pairs if r['found'] and not held(r))}")
t = Counter((r["who"], held(r)) for r in pairs)
for who in ("father", "mother"):
    print(f"  {who}'s: {t[(who, True)]} of {t[(who, True)] + t[(who, False)]} by the profile")
for r in sorted(pairs, key=lambda r: r["child_step"]):
    if r["found"] != held(r):
        print(f"    {r['who']:6s} {r['parent']} -> {r['child']}  breakpoint at {r['at']:.0f} kb ({r['step']:+d}); parent's step {r['parent_step']:+.2f}, child's {r['child_step']:+.2f} (signed); "
              f"child's call: {rows[r['child']]['DJ.variants']}")
