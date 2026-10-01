"""Chromosomes in copies on the page (report/karyotype.py, report/karyotype_page.py): the summary of a cohort whose
truth is known, the flags, and the section's text."""
import re
import sys
from pathlib import Path

import numpy as np
import pytest

K = pytest.importorskip("ngsdose.karyotype", reason="ngsdose before 0.3.0 does not read chromosomes")
from ngsdose import resources  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from report import karyotype as kmod  # noqa: E402
from report import karyotype_page  # noqa: E402
from report.report_page import Page  # noqa: E402

BUNDLE = resources.Bundle()
ARMS, LENGTHS = BUNDLE.karyotype()["arms"], BUNDLE.contig_lengths()
ACRO = kmod.ACROCENTRIC
U2 = 0.95


def layout(per_arm=10):
    names = []
    for c in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
        b, n = ARMS[c], LENGTHS[c]
        arms = [(b + 3_000_000, min(n, 27_000_000 if c == "chrY" else n) - 3_000_000)]
        if c not in ACRO and c != "chrY":
            arms.insert(0, (3_000_000, b - 3_000_000))
        for lo, hi in arms:
            names += [f"{c}:{s}-{s + 12_000}" for s in np.linspace(lo, hi, per_arm).astype(int)]
    return names


@pytest.fixture(scope="module")
def cohort():
    """120 genomes: plain men and women, and by construction a trisomy 21, a 45,X, an XXY, chromosome 12 gained in 30%
    of the cells, an X lost in a quarter of a woman's, the long arm of 5 gained in a child and in her father, and one
    clone that gained three chromosomes in a fifth of its cells."""
    rng = np.random.default_rng(3)
    names = layout()
    tab = K.table(names, ARMS)
    n = 120
    Y = rng.normal(0, 0.1, (n, 1)) + rng.normal(0, 0.05, len(names))[None, :] + rng.normal(0, 0.03, (n, len(names)))
    men = np.arange(n) % 2 == 0

    def sex(i, nx, ny):
        Y[i, tab.kind == "X"] += np.log(max(1 + (nx - 1) * U2, 1e-3) / 2)
        Y[i, tab.kind == "Y"] += np.log(ny / 2) if ny else np.log(1.5e-3)
    for i in range(n):
        if i not in (3, 4, 9):
            sex(i, 1 if men[i] else 2, 1 if men[i] else 0)
    sex(3, 1, 0)                                                       # 45,X (a woman by the pedigree)
    sex(4, 2, 1)                                                       # 47,XXY
    Y[9, tab.kind == "X"] += np.log((1 + 0.75 * U2) / 2)               # an X lost in a quarter of the cells
    Y[9, tab.kind == "Y"] += np.log(1.5e-3)
    Y[5, tab.chrom == "chr21"] += np.log(1.5)
    Y[6, tab.chrom == "chr12"] += np.log(1.15)
    q5 = (tab.chrom == "chr5") & (tab.arm == "q")
    Y[10, q5] += np.log(1.5)                                           # a child (S10) and her father (S12)
    Y[12, q5] += np.log(1.5)
    for c in ("chr9", "chr12", "chr15"):
        Y[20, tab.chrom == c] += np.log(1.1)
    samples = [f"S{i}" for i in range(n)]
    rd, model, info = K.cohort([(names, y) for y in Y], ARMS)
    chroms = K.table(model.names, model.arms).chromosomes()
    rows = [dict(sample=s, pop="POP", flagged_chromosomes="chr21" if i == 5 else "chr5" if i == 10 else None, **{"DJ.cn": 11.02 if i == 5 else 10.0, "DJ.call": "settled", "DJ.copies": 11 if i == 5 else 10},
                 **{"ngspca.chrX": r.x + rng.normal(0, 0.01)}, **K.columns(r, chroms)) for i, (s, r) in enumerate(zip(samples, rd))]
    ped = {s: dict(sex="M" if men[i] else "F", father="0", mother="0", pop="POP") for i, s in enumerate(samples)}
    ped["S3"]["sex"], ped["S10"]["father"], ped["S10"]["mother"] = "F", "S12", "S13"
    kout = dict(samples=samples, readings=rd, model=model, info=info)
    return rows, kout, ped


def test_the_summary_counts_what_the_genomes_hold(cohort):
    rows, kout, ped = cohort
    said = []
    k = kmod.run(rows, kout, ped, log=said.append)
    assert k["n"] == 120 and k["model"]["regions"] == 420 and k["model"]["x"] == 20 and k["model"]["y"] == 10 and abs(k["model"]["u"] - U2) < 0.02
    assert k["plain"] == 120 - 8 and k["status"]["fractional"] == 3 and "chromosomes: 120 genomes read" in said[0]
    sx = k["sex"]
    assert sx["by_complement"] == {"X": 1, "XX": 59, "XXY": 1, "XY": 59} and sx["mismatch"] == []
    assert {t["sample"]: t["karyotype"] for t in sx["not_plain"] if t["complement"] not in ("XX", "XY")} == {"S3": "45,X", "S4": "47,XXY"}
    part = sx["part"]["x_lost"]
    assert part["n"] == 1 and abs(part["q"][1] - 0.25) < 0.04 and sx["part"]["y_lost"]["n"] == 0 and sx["ngspca"]["r"] > 0.99
    au = k["autosomes"]
    by = {c["chrom"]: c for c in au["by_chromosome"]}
    assert by["chr21"]["whole"] == 1 and by["chr12"]["part"] == 2 and by["chr9"]["part"] == 1 and au["gained"] == 5 and au["lost"] == 0
    assert [(s_["sample"], s_["chrom"], s_["span"], s_["cells"]) for s_ in sorted(au["stretches"], key=lambda x: x["sample"])] == [("S10", "chr5", "q", "every cell"), ("S12", "chr5", "q", "every cell")]
    (clone,) = au["clones"]
    assert clone["sample"] == "S20" and len(clone["chromosomes"]) == 3 and abs(clone["share"] - 0.2) < 0.04
    # the table of chromosomes: every autosome once, chrX by its copies, chrY for the genomes that hold one
    labels = [c["label"] for c in k["chromosomes"]]
    assert labels[:22] == [str(i) for i in range(1, 23)] and labels[22:] == ["X, one copy", "X, two copies", "Y, one copy"]
    c21 = next(c for c in k["chromosomes"] if c["label"] == "21")
    assert c21["regions"] == 10 and c21["whole"] == 1 and 0.01 < c21["se"] < 0.04 and c21["seen_from"] == max(0.03, 5 * c21["se"])


def test_the_checks_use_the_arms_the_junction_and_the_family(cohort):
    rows, kout, ped = cohort
    ck = kmod.run(rows, kout, ped)["checks"]
    a = ck["arms"]
    assert a["n"] == 3 and a["same_side"] == 3 and a["r"] > 0.95                 # chr12 twice and chr9 (chr15 and chr21 have one arm)
    t = ck["tails"]["rows"][0]
    assert t["z"] == 3 and t["below"] <= 8 and abs(t["expected"] - 0.00135 * ck["tails"]["n"]) < 0.1
    assert ck["flags"] == {"n": 2, "whole": 1, "stretch": 1, "neither": 0, "stretches": ["S10 +5q"]}
    d, e = [d for d in ck["dj"] if d["whole_chromosome"]]                         # the largest excess first
    assert d["sample"] == "S5" and abs(d["extra"] - 1) < 0.06 and d["dj"] == 11.02 and e["sample"] == "S20" and abs(e["extra"] - 0.2) < 0.06
    rel = ck["relatives"]
    assert rel["pairs"] == 2 and [(s_["child"], s_["parent"], s_["role"]) for s_ in rel["shared"]] == [("S10", "S12", "father")]
    ev = kmod.event_rows(kout["samples"], kout["readings"])
    assert {e["sample"] for e in ev} == {"S5", "S6", "S9", "S10", "S12", "S20"} and next(e for e in ev if e["sample"] == "S5")["label"] == "+21"


def test_flags_say_what_is_not_plain():
    assert kmod.flags({"karyotype": "46,XX"}) == [] and kmod.flags({}) == [] and kmod.flags({"karyotype": "NA"}) == []
    assert kmod.flags({"karyotype": "47,XY,+21", "karyotype.status": "settled"}) == ["karyotype 47,XY,+21"]
    f = kmod.flags({"karyotype": "46,XX,-X[0.20]", "karyotype.status": "fractional", "karyotype.note": "X and Y are off"})
    assert f == ["karyotype 46,XX,-X[0.20] (part of the cells in brackets)", "X and Y are off"]
    assert kmod.flags({"karyotype": "46,XY", "karyotype.status": "uncertain"}) == ["chromosomes read too coarsely to settle"]


def test_the_section_says_it_in_words(cohort):
    rows, kout, ped = cohort
    data = dict(meta=dict(title="t", as_of="2026-10-01"), karyotype=kmod.run(rows, kout, ped))
    data["karyotype"]["fetch"] = dict(n=120, identical=120, differ=[])
    P = Page(data)
    karyotype_page.render(P, data, rows)
    html = "".join(P.parts)
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
    assert 'id="karyotype"' in html and "karyo_xy" in P.charts and P.charts["karyo_xy"]["x"] == "chrX.copies"
    for phrase in ("59 XY", "59 XX", "S3 45,X", "S4 47,XXY", "5 gained and 0 lost", "S20 reads", "120 of 120 genomes are written exactly as from the scan",
                   "S10 +5q and father S12 +5q", "What depth does not see", "S5: +21, junction 11.02"):
        assert phrase in text, phrase
    assert "None" not in text and "nan" not in text.lower().replace("nanopore", "")
    # nothing to say: nothing is written
    P2 = Page(dict(meta=data["meta"], karyotype=None))
    karyotype_page.render(P2, dict(karyotype=None), rows)
    assert P2.parts == [] and kmod.run(rows, {}, ped) is None and kmod.run(rows, dict(model=None), ped) is None
