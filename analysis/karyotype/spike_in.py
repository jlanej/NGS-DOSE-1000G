"""What the calls make of a chromosome or an arm gained or lost in a share of the cells: the change is added to the regions
of clean genomes counted with the karyotype windows, and read back against the cohort's model, from all their regions and
from those of a lighter control set. The noise is the genomes' own.

usage (from the repository's root; about fifteen minutes): python3 analysis/karyotype/spike_in.py [--set screen]
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, "analysis/karyotype")
from common import BUNDLE, K, aligned, cohort  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--set", default="screen", help="the lighter control set to read the same genomes from as well (`ngsdose control-sets`)")
ap.add_argument("--seed", type=int, default=7)
a = ap.parse_args()
C = cohort()
model, samples, win, full = C["model"], C["samples"], C["with_windows"], C["readings"]
Ym, tab, clean = aligned(C, [i for i, s in enumerate(samples) if s in win])
men = [i for i in clean if full[i].complement == "XY"]
women = [i for i in clean if full[i].complement == "XX"]
u = float(np.nanmedian(model.u))
shares = (0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.5, 1.0)
gain, loss = lambda f: np.log(1 + f / 2), lambda f: np.log(1 - f / 2)
sets = {"every region": np.ones(len(model.names), bool)}
f = BUNDLE.dir / f"controls.{a.set}.bed"
if f.exists():
    keep = {"%s:%s-%s" % tuple(line.split("\t")[:3]) for line in open(f) if line.strip() and not line.startswith("#")}
    sets[a.set] = np.array([n in keep for n in model.names])


def test(tag, genomes, where, shift, want, chrom, use, fr=shares):
    out = []
    for share in fr:
        hit, est = 0, []
        for i in genomes:
            y = Ym[i].copy()
            y[where] += shift(share)
            y[~use] = np.nan
            ev = [e for e in K.read(model, model.names, y).events if e.chrom == chrom and e.span in want]
            hit += bool(ev)
            est += [abs(e.delta) for e in ev[:1]]
        out.append(f"{share:.2f}: {100 * hit / len(genomes):3.0f}%" + (f" ({np.mean(est):.2f})" if est else "       "))
    print(f"   {tag:18s} regions {int((where & use & np.isfinite(model.sd)).sum()):4d}  " + "  ".join(out), flush=True)


for name, use in sets.items():
    print(f"== {name}: share of the cells: found in this share of {len(clean)} clean genomes (the share read back)")
    for c in ("chr1", "chr13", "chr18", "chr21", "chr19", "chr22"):
        test("+" + c[3:], clean, tab.chrom == c, gain, ("whole",), c, use)
    for c in ("chr7", "chr21"):
        test("-" + c[3:], clean, tab.chrom == c, loss, ("whole",), c, use)
    for c, arm in (("chr5", "p"), ("chr18", "p"), ("chr17", "p"), ("chr1", "q")):
        test(f"+{c[3:]}{arm}", clean, (tab.chrom == c) & (tab.arm == arm), gain, (arm, "pter" if arm == "p" else "qter"), c, use)
    X, Y = tab.kind == "X", tab.kind == "Y"
    test("X lost, a woman", women, X, lambda f: np.log((1 + (1 - f) * u) / (1 + u)), ("whole",), "chrX", use, shares[:-1])
    test("X gained, a man", men, X, lambda f: np.log(1 + f * u), ("whole",), "chrX", use, shares[:-1])
    test("Y lost, a man", men, Y, lambda f: np.log(max(1 - f, 1e-3)), ("whole",), "chrY", use, shares[:-1])
    test("Y gained, a man", men, Y, lambda f: np.log(1 + f), ("whole",), "chrY", use, shares[:-1])
