"""The distal junction in profile (report/dj.py), the screen-to-tables script (pipeline/hprc_dj.py), and the
page's section, on synthetic profiles and on the pilot counts with a synthetic pair of assembly tables."""
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
from report import dj  # noqa: E402
import hprc_dj  # noqa: E402

PILOT = ROOT / "pilot"
UNIT, W = dj.UNIT, 250


def windows(level, rng, deletion=None, gain=None, noise=0.5):
    """A genome's usable windows: `level` copies everywhere, `deletion`/`gain` = (start, end, copies) on part of the unit."""
    st = np.arange(0, UNIT, W, dtype=float)
    cn = np.full(len(st), float(level))
    for iv, sign in ((deletion, -1), (gain, 1)):
        if iv:
            a, b, k = iv
            cn[(st >= a) & (st < b)] += sign * k
    cn = cn * np.exp(rng.normal(0, noise / level, len(st)))
    return np.column_stack([st, cn])


def test_segments_find_the_polymorphic_part_and_the_core_level_ignores_it():
    rng = np.random.default_rng(1)
    S, Wn = [], []
    for i in range(300):
        state = rng.choice([-1, 0, 1], p=[0.25, 0.55, 0.2])           # a 200-220 kb polymorphism in half the genomes
        S.append(f"S{i}")
        Wn.append(windows(10, rng, deletion=(200000, 220000, -state) if state else None))
    P = dj.profiles(S, Wn, None)
    assert P["block"].shape == (300, dj.NBLOCK) and P["sub"].shape == (300, dj.NSUB)
    seg = dj.segments(P)
    assert seg["hyper_intervals"] == [(200000, 220000)]
    assert 10 not in seg["core_blocks"] and len(seg["core_blocks"]) == dj.NBLOCK - 1
    st = dj.cohort_stats(P, seg, 10.0)
    assert abs(st["pin"] - 1) < 0.02 and st["core"]["within_0_2"] > 0.95
    assert st["breakpoints"][200] > 100 and st["breakpoints"][220] > 100 and st["breakpoints"][100] < 10
    states = st["states"][200]
    assert states.get(9, 0) > 40 and states.get(11, 0) > 30


def test_the_level_is_calibrated_by_the_efficiencies():
    rng = np.random.default_rng(2)
    st = np.arange(0, UNIT, W, dtype=float)
    a = np.where(st < 200000, 0.1, -0.1)                              # windows read 10% high in the first half, 10% low in the second
    Wn = [np.column_stack([st, 10 * np.exp(a) * np.exp(rng.normal(0, 0.03, len(st)))]) for _ in range(20)]
    eff = dict(start=st.tolist(), a=a.tolist())
    P = dj.profiles([f"S{i}" for i in range(20)], Wn, eff)
    assert np.allclose(P["level"], 10, atol=0.1)
    assert np.allclose(np.nanmean(P["block"], axis=0), 10, atol=0.15)
    raw = dj.profiles([f"S{i}" for i in range(20)], Wn, None)
    assert abs(np.nanmean(raw["block"][:, :10]) - 11.05) < 0.2       # uncalibrated, the first half reads high


def test_group_of():
    assert dj.group_of(-1.02) == "carrier" and dj.group_of(0.05) == "zero" and dj.group_of(0.55) == "between" and dj.group_of(None) == ""
    assert dj.group_of(-0.31) == "between" and dj.group_of(-0.29) == "zero" and dj.group_of(1.9) == "carrier"


def test_copies_from_a_paf_and_their_classes(tmp_path):
    """A complete copy, a partial copy in tandem with it on one contig, and a copy cut by a contig end."""
    def paf(ctg, clen, q0, t0, t1, ident=0.99):
        aln = t1 - t0
        return f"{ctg}\t{clen}\t{q0}\t{q0 + aln}\t+\tDJ\t{UNIT}\t{t0}\t{t1}\t{int(aln * ident)}\t{aln}\t60\n"
    lines = []
    for t0 in range(3000, 400000, 10000):                                 # the complete copy, pieces of 8 kb with 2-kb masked gaps
        lines.append(paf("ctgA", 5_000_000, 1_000_000 + t0, t0, min(t0 + 8000, UNIT)))
    for t0 in range(3000, 316000, 10000):                                 # the partial copy, 400 kb downstream on the same contig
        lines.append(paf("ctgA", 5_000_000, 1_600_000 + t0, t0, min(t0 + 8000, 316000)))
    for t0 in range(110000, 400000, 10000):                               # a copy whose first 110 kb fell off the contig's start
        lines.append(paf("ctgB", 300_000, t0 - 110000, t0, min(t0 + 8000, UNIT)))
    lines.append(paf("ctgC", 1_000_000, 500_000, 59900, 65900, 0.95))     # a repeat element hit: too short to be a copy
    p = tmp_path / "x.masked.paf"
    p.write_text("".join(lines))
    cps = hprc_dj.copies_from_paf(p, UNIT)
    complete = max(c["covered_core_bp"] for c in cps)
    kinds = sorted((c["contig"], hprc_dj.classify(c, UNIT, complete)) for c in cps)
    assert kinds == [("ctgA", "complete"), ("ctgA", "partial"), ("ctgB", "truncated at contig end")]
    part = next(c for c in cps if c["dj_end"] < 320000)
    assert part["dj_start"] == 3000 and part["dj_end"] == 316000 and abs(part["blocks"][0] - 0.8 * 17 / 20) < 0.1


def test_haplotype_row_medians_are_robust_to_single_kmer_losses():
    rng = np.random.default_rng(3)
    core = np.sort(rng.choice(UNIT - 31, 60000, replace=False))
    bgc = np.full(len(core), 5, dtype=np.int32)
    bgc[rng.random(len(core)) < 0.06] -= 1                                # 6% of k-mers lost to nucleotide differences
    bgc[(core >= 200000) & (core < 220000)] -= 1                          # one copy lacks 200-220 kb
    row = hprc_dj.haplotype_row("S", "hap1", bgc, core, UNIT)
    assert row["unit_median"] == 5 and row["b0"] == 5 and row["b40"] == 4 and row["b43"] == 4 and row["b44"] == 5
    assert 0.85 < row["frac_at_median"] < 0.96                        # 6% lost, and the 5% of k-mers in the deleted part


@pytest.mark.skipif(len(list((PILOT / "counts_nygc").glob("*.json.gz"))) < 6, reason="pilot counts not present")
def test_the_page_compares_the_pilot_with_synthetic_assemblies(tmp_path):
    """Two pilot genomes given assemblies: one with five complete copies per haplotype, one whose paternal
    haplotype has four, one of them lacking 200-220 kb; the page's section and tables follow."""
    asm = tmp_path / "dj_hprc"
    asm.mkdir()
    cols = ["sample", "haplotype", "n_core", "unit_median", "frac_at_median", "frac_above_median", "mean"] + [f"b{i}" for i in range(dj.NSUB)]
    rows = []
    for s, hap, n, hole in (("HG00733", "mat", 5, None), ("HG00733", "pat", 5, None), ("NA19240", "mat", 5, None), ("NA19240", "pat", 4, (40, 44))):
        sub = [n] * dj.NSUB
        if hole:
            for i in range(*hole):
                sub[i] = n - 1
        rows.append([s, hap, 169808, n, 0.8, 0.01, n - 0.2] + sub)
    with open(asm / "haplotypes.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows))
    ccols = ["sample", "haplotype", "contig", "contig_len", "q_start", "q_end", "dj_start", "dj_end", "covered_core_bp", "to_contig_start", "to_contig_end", "strands", "blocks", "class"]
    crow = [["HG00733", "mat", "c1", 3000000, 1000000, 1400000, 3000, 400000, 305000, 1000000, 1600000, "+", ",".join(["1.00"] * 20), "complete"]] * 5 \
        + [["HG00733", "pat", "c2", 3000000, 1000000, 1400000, 3000, 400000, 305000, 1000000, 1600000, "+", ",".join(["1.00"] * 20), "complete"]] * 5 \
        + [["NA19240", "mat", "c3", 3000000, 1000000, 1400000, 3000, 400000, 305000, 1000000, 1600000, "+", ",".join(["1.00"] * 20), "complete"]] * 5 \
        + [["NA19240", "pat", "c4", 3000000, 1000000, 1400000, 3000, 400000, 290000, 1000000, 1600000, "+", ",".join(["1.00"] * 20), "complete, internal gaps"]] * 4
    with open(asm / "copies.tsv", "w") as fh:
        fh.write("\t".join(ccols) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in crow))
    out = tmp_path / "report"
    ped = tmp_path / "ped.txt"
    ped.write_text("FamilyID SampleID FatherID MotherID Sex Population Superpopulation\nPR05 HG00733 HG00731 HG00732 2 PUR AMR\nPR05 HG00731 0 0 1 PUR AMR\nPR05 HG00732 0 0 2 PUR AMR\n"
                   "Y117 NA19240 NA19239 NA19238 2 YRI AFR\nY117 NA19239 0 0 1 YRI AFR\nY117 NA19238 0 0 2 YRI AFR\n")
    subprocess.run([sys.executable, "-m", "report", "--fetch", str(PILOT / "counts_nygc"), "-p", str(ped), "--pilot", str(PILOT), "-o", str(out), "-j", "2",
                    "--as-of", "2026-09-28", "--dj-assemblies", str(asm)], check=True, cwd=ROOT, capture_output=True)
    html = (out / "index.html").read_text()
    for must in ("The junction in profile", "Against the HPRC release-2 assemblies of 2 cohort members", "What the assemblies get wrong", "What follows for the method",
                 'id="chart-djasm"', "data/dj_hprc.tsv", "data/dj_blocks.tsv"):
        assert must in html, must
    visible = html.split('<script id="report-data"')[0]
    assert "NaN" not in visible and "None" not in visible.replace("None of", "")
    with open(out / "data" / "dj_hprc.tsv") as fh:
        T = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert set(T) == {"HG00733", "NA19240"}
    assert float(T["HG00733"]["assembly_core"]) == 10 and T["HG00733"]["resolved"] == "True" and T["HG00733"]["haplotypes"] == "mat=5 + pat=5"
    assert float(T["NA19240"]["assembly_core"]) == 9 and float(T["NA19240"]["assembly_mean"]) < 9
    with open(out / "data" / "dj_hprc_blocks.tsv") as fh:
        B = [r for r in csv.DictReader(fh, delimiter="\t") if r["sample"] == "NA19240"]
    assert len(B) == dj.NBLOCK and float(next(r for r in B if r["block_kb"] == "200")["assembly"]) == 8 and float(next(r for r in B if r["block_kb"] == "40")["assembly"]) == 9
    # the pilot genomes read near ten on the pinned scale, and the profile table has every genome
    with open(out / "data" / "dj_blocks.tsv") as fh:
        P = list(csv.DictReader(fh, delimiter="\t"))
    assert len(P) == 12 and all(8.5 < float(r["core_level_pinned"]) < 11.5 for r in P)
    if (out / "dj_assemblies.png").exists():
        assert (out / "dj_assemblies.png").stat().st_size > 10000 and 'src="dj_assemblies.png"' in html
