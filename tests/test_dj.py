"""The distal junction in profile (report/dj.py, on the cohort layer's profiles), the screen-to-tables script
(pipeline/hprc_dj.py), and the page's section, on synthetic profiles and on the pilot counts with a synthetic
pair of assembly tables."""
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
STARTS = np.arange(0, UNIT, W)


def library(rng, n=300, poly=(200000, 220000), exclude=((0, 30000), (190000, 232000)), calls=True):
    """What `cohort_table(..., profiles=...)` returns for a class: calibrated window profiles of a cohort at ten copies,
    a polymorphic interval that one copy lacks in a quarter of the genomes and holds twice in a fifth, and the calls."""
    from ngsdose import segments
    cn = np.full((n, len(STARTS)), 10.0)
    state = rng.choice([-1, 0, 1], n, p=[0.25, 0.55, 0.2])
    inside = (STARTS >= poly[0]) & (STARTS < poly[1])
    cn[:, inside] += state[:, None]
    cn = cn * np.exp(rng.normal(0, 0.05, cn.shape))
    level = np.ones(len(STARTS), bool)
    for a, b in exclude:
        level &= ~((STARTS + W > a) & (STARTS < b))
    names = [f"S{i}" for i in range(n)]
    cl = {}
    if calls:
        for i, s in enumerate(names[:40]):
            cl[s] = segments.describe(segments.segment(STARTS, cn[i], unit_length=UNIT), list(exclude), expected=10)
    return dict(samples=names, start=STARTS.tolist(), end=(STARTS + W).tolist(), cn=cn.astype(np.float32), level=level, calls=cl,
                scale=dict(rule="mode", expected=10.0, factor=1.02, mode_on_anchors=9.8, n_main=280, n=n), offsets=[]), state


def test_segments_find_the_polymorphic_part_and_the_core_is_the_librarys():
    rng = np.random.default_rng(1)
    lib, state = library(rng)
    P = dj.profiles(lib)
    assert P["block"].shape == (300, dj.NBLOCK) and P["sub"].shape == (300, dj.NSUB)
    assert np.allclose(P["level"], 10, atol=0.08)                                      # the level is the core's: the polymorphism does not move it
    seg = dj.segments(P, dict(level_exclude=[(0, 30000), (190000, 232000)]))
    assert seg["hyper_intervals"] == [(200000, 220000)]
    assert seg["excluded"] == [(0, 30000), (190000, 232000)] and dj.segments(P)["excluded"] == [] and len(dj.segments(P)["core_blocks"]) == dj.NBLOCK
    assert set(range(dj.NBLOCK)) - set(seg["core_blocks"]) == {0, 1, 9, 10, 11}
    st = dj.cohort_stats(P, seg)
    assert st["core"]["within_0_2"] > 0.99 and st["scale"]["rule"] == "mode"
    assert st["breakpoints"][200] > 100 and st["breakpoints"][220] > 100 and st["breakpoints"][100] < 10
    c = st["calls"]
    assert c["n"] == 40 and c["copies"] == {10: 40} and c["whole"] == {10: 40} and c["with_partial_copy"] == 0 and c["scale_uncertain"] == 0
    assert {b["kb"] for b in c["breakpoints"]} <= {195, 200, 215, 220}
    # the calls' states in the polymorphic interval are the genomes'
    got = np.array([P["calls"][s].state_at(210000) - 10 for s in P["samples"][:40]])
    assert (got == state[:40]).mean() > 0.95
    assert len(dj.call_table(P)) == sum(len(P["calls"][s].segments) for s in P["calls"])


def test_without_profiles_there_is_nothing():
    assert dj.profiles(None) is None and dj.run([], None, {}, {}, None) is None


def test_group_of():
    g = dj.group_of
    assert g({"DJ.copies": 9, "DJ.partial": "none", "DJ.call": "settled"}) == "carrier" and g({"DJ.copies": 10, "DJ.partial": "+1:0-316kb", "DJ.call": "settled"}) == "partial copy"
    assert g({"DJ.copies": 11, "DJ.partial": "none", "DJ.call": "uncertain"}) == "uncertain" and g({"DJ.copies": 10, "DJ.partial": "none", "DJ.call": "settled"}) == "ten"
    assert g({"DJ.step": -1.02}) == "carrier" and g({"DJ.step": 0.05}) == "zero" and g({"DJ.step": 0.55}) == "between" and g({"DJ.step": None}) == ""


def _row(sample, sex, cn, copies, variants="none", partial=None, call="settled"):
    return dict(sample=sample, sex_inferred=sex, **{"DJ.cn": cn, "DJ.copies": copies, "DJ.variants": variants, "DJ.partial": partial or "none", "DJ.scale_f": 1.0, "DJ.call": call})


def test_the_step_logic_by_calls():
    """The copies a settled call is described against give the steps; a copy that holds an end of the unit is tallied and its
    transmission tested apart; an uncertain call stays out; the calls are tested position by position in the trios."""
    from report.report import dj_mendel, dj_states, dj_steps, partial_class
    ped = {f"kid{i}": dict(father=f"dad{i}", mother=f"mum{i}") for i in range(1, 6)}
    rows = [_row("dad1", "M", 9.0, 9), _row("mum1", "F", 10.0, 10), _row("kid1", "F", 9.05, 9),                       # a father's loss, passed on
            _row("dad2", "M", 10.0, 10), _row("mum2", "F", 9.0, 9), _row("kid2", "M", 10.0, 10),                       # a mother's, not
            _row("dad3", "M", 10.0, 10), _row("mum3", "F", 10.0, 10, "-1:197-217kb"), _row("kid3", "M", 11.0, 11),     # a gain neither parent has
            _row("dad4", "M", 10.8, 10, "+1:16-318kb", "+1:16-318kb"), _row("mum4", "F", 10.0, 10), _row("kid4", "M", 10.75, 10, "+1:0-315kb", "+1:0-315kb"),
            _row("dad5", "M", 10.45, 11, call="uncertain"), _row("mum5", "F", 10.0, 10), _row("kid5", "F", 10.0, 10)] \
        + [_row(f"x{i}", "F", 10.0 + 0.01 * (i % 7 - 3), 10) for i in range(30)]
    d = dj_steps(rows, ped)
    assert d["basis"] == "calls" and d["n_called"] == 45 and d["n_settled"] == 44 and d["uncertain"] == ["dad5"]
    assert d["near"][-1] == 3 and d["near"][1] == 1 and d["near"][0] == 40 and d["between"] == 2
    assert d["whole"] == {9: 3, 10: 38, 11: 1} and d["plain"] == 38
    assert d["transmitted"] == 1 and d["not_transmitted"] == 1 and d["unclassified"] == 0 and d["de_novo"] == ["kid3"]
    assert d["by_parent"] == {"father, loss": dict(transmitted=1, not_transmitted=0), "mother, loss": dict(transmitted=0, not_transmitted=1)}
    assert {c["sample"] for c in d["carriers"]} == {"dad1", "kid1", "mum2", "kid3"}
    p = d["partial"]
    assert p["n_carriers"] == 2 and p["n_gain"] == 2 and p["n_loss"] == 0 and p["de_novo"] == []
    assert p["kinds"] == [dict(kind="copy of 0-316 kb", n=2, transmitted=1, not_transmitted=0)]
    assert p["transmitted"] == 1 and p["not_transmitted"] == 0 and p["several"] == 0 and p["by_parent"] == dict(father=dict(transmitted=1, not_transmitted=0))
    from report.report import dj_breakpoints
    assert dj_breakpoints(rows[9], 10) == [(320.0, -1)] and dj_breakpoints(rows[0], 10) == []                   # dad4's copy ends at 318 kb, between the blocks at 315 and 320 kb; dad1 holds nine throughout
    assert dj_breakpoints(rows[7], 10) == []                                                                       # mum3's deletion lies in a stretch the level leaves out
    # one copy, two descriptions: ten with a copy that lacks the last 84 kb, and ten with a copy of the first 316 kb
    a, b = _row("a", "F", 9.8, 10, "-1:316-400kb", "-1:316-400kb"), _row("b", "F", 10.8, 10, "+1:0-317kb", "+1:0-317kb")
    (qa, sa), (qb, sb) = dj_breakpoints(a, 10)[0], dj_breakpoints(b, 10)[0]
    assert abs(qa - qb) <= 5 and sa == sb == -1
    trio = [_row("dad", "M", 10.0, 10), a | dict(sample="mum"), b | dict(sample="kid")] + [_row(f"x{i}", "F", 10.0, 10) for i in range(5)]
    q = dj_steps(trio, dict(kid=dict(father="dad", mother="mum")))["partial"]
    assert q["transmitted"] == 1 and q["not_transmitted"] == 0 and q["by_parent"] == dict(mother=dict(transmitted=1, not_transmitted=0))
    assert partial_class(123, 400) == "copy of 122-400 kb" and partial_class(0, 263) == "copy of 0-262 kb" and partial_class(0, 124, -1) == "loss of 0-122 kb"
    # position by position: kid3's gain is in no parent; the polymorphic deletion of mum3 is not in her son, which is Mendelian
    st = dj_states(rows[9], 10)
    assert st[0] == 10 and st[10] == 11 and st[62] == 11 and st[64] == 10 and len(st) == 80
    m = d["mendel"]
    assert m["n_trios"] == 4 and m["clean"] == 3 and m["worst"][0]["child"] == "kid3"
    assert 0.7 < m["core"]["share"] < 0.8 and m["polymorphic"]["informative"] > m["polymorphic"]["passed"] > 0
    assert m["core"]["informative"] > 0 and abs(m["core"]["passed_share"] - 2 / 3) < 0.05                      # of dad1's, mum2's and dad4's, two passed on
    assert m["core"]["child_dev"] > m["core"]["child_explained"] > 0
    # without calls, by the level
    for r in rows:
        for k in ("DJ.copies", "DJ.partial", "DJ.scale_f", "DJ.variants", "DJ.call"):
            r.pop(k)
    d = dj_steps(rows, ped)
    assert d["basis"] == "level" and d["near"][-1] == 3 and d["near"][1] == 3 and "partial" not in d            # by the level, a partial copy reads as a step


def _families(rng, n, passed, noise=0.05, poly=(200000, 220000), exclude=((0, 30000), (190000, 232000)), calls=True):
    """Trios whose parents lack the polymorphic interval in none, one or two of their copies (at two loci); a child
    takes each deleted copy of a parent with probability `passed` (one half is Mendel's)."""
    from ngsdose import segments
    inside = (STARTS >= poly[0]) & (STARTS < poly[1])
    level = np.ones(len(STARTS), bool)
    for a, b in exclude:
        level &= ~((STARTS + W > a) & (STARTS < b))
    names, state, ped = [], [], {}
    for i in range(n):
        f, m = rng.choice([0, 1, 2], 2, p=[0.55, 0.35, 0.10])                   # copies that lack it
        kid = int(rng.binomial(f, passed) + rng.binomial(m, passed))
        names += [f"dad{i}", f"mum{i}", f"kid{i}"]
        state += [-f, -m, -kid]
        ped[f"kid{i}"] = dict(father=f"dad{i}", mother=f"mum{i}")
    cn = np.full((len(names), len(STARTS)), 10.0)
    cn[:, inside] += np.array(state)[:, None]
    cn = cn * np.exp(rng.normal(0, noise, cn.shape))
    calls = {s: segments.describe(segments.segment(STARTS, cn[i], unit_length=UNIT), list(exclude), expected=10) for i, s in enumerate(names)} if calls else {}
    lib = dict(samples=names, start=STARTS.tolist(), end=(STARTS + W).tolist(), cn=cn.astype(np.float32), level=level, calls=calls, scale={}, offsets=[])
    return lib, ped


def test_what_varies_is_passed_on_is_measured_without_a_call():
    """The child's value in a polymorphic interval against the mean of its parents': at Mendel's rate the slope is the
    value's reliability; a variant passed on less often has a slope below it; too few trios give nothing."""
    rules = dict(polymorphic=[dict(name="200-220 kb", interval=(200000, 220000))], level_exclude=[(0, 30000), (190000, 232000)])
    lib, ped = _families(np.random.default_rng(5), 600, 0.5, calls=False)
    (r,) = dj.inheritance(dj.profiles(lib), ped, rules, n_boot=300)
    assert r["n_trios"] == 600 and r["windows"] == 80 and r["reliability"] > 0.97 and r["noise_sd"] < 0.1
    assert abs(r["child_minus_midparent"]) < 3 * r["child_minus_midparent_se"] < 0.1
    assert r["slope_lo"] < r["reliability"] < r["slope_hi"] and abs(r["of_expected"] - 1) < 0.1
    v = r["by_value"]
    assert v["pairs"] > 150 and abs(v["passed"] / v["pairs"] - 0.5) < 0.1 and v["p"] > 0.01 and v["new"] == 0 and v["both_at_zero"] > 100
    assert r["by_calls"] is None                                              # no calls, no count by them
    lib, ped = _families(np.random.default_rng(7), 40, 0.5)
    (r,) = dj.inheritance(dj.profiles(lib), ped, rules, n_boot=100)
    assert r["by_calls"] == r["by_value"] and r["by_calls"]["pairs"] > 5      # a deletion of 20 kb at this noise is called in every genome
    lib, ped = _families(np.random.default_rng(6), 250, 0.25, calls=False)
    (r,) = dj.inheritance(dj.profiles(lib), ped, rules, n_boot=300)
    assert r["slope_hi"] < r["reliability"] and abs(r["of_expected"] - 0.5) < 0.15
    assert r["by_value"]["passed"] / r["by_value"]["pairs"] < 0.35 and r["by_value"]["p"] < 0.01 and r["by_calls"] is None
    assert dj.inheritance(dj.profiles(lib), dict(list(ped.items())[:10]), rules) is None
    assert dj.inheritance(dj.profiles(lib), ped, None) is None and dj.inheritance(dj.profiles(lib), None, rules) is None


def test_events_that_arose_in_the_lines_would_show_in_the_children():
    """Carrier parents pass their events to a quarter of the children; had the excess arisen in culture, the children of
    two parents at ten would hold new ones, and the count says how many."""
    from report.report import dj_steps
    ped, rows = {}, []
    for i in range(40):                                                    # a father at nine, a mother at ten; one child in four at nine
        ped[f"kid{i}"] = dict(father=f"dad{i}", mother=f"mum{i}")
        rows += [_row(f"dad{i}", "M", 9.0, 9), _row(f"mum{i}", "F", 10.0, 10), _row(f"kid{i}", "F", 9.0 if i % 4 == 0 else 10.0, 9 if i % 4 == 0 else 10)]
    for i in range(40, 140):                                               # both parents at ten; one child with a new gain
        ped[f"kid{i}"] = dict(father=f"dad{i}", mother=f"mum{i}")
        rows += [_row(f"dad{i}", "M", 10.0, 10), _row(f"mum{i}", "F", 10.0, 10), _row(f"kid{i}", "M", 11.0 if i == 40 else 10.0, 11 if i == 40 else 10)]
    d = dj_steps(rows, ped)
    ln = d["lines"]
    assert d["transmitted"] == 10 and d["not_transmitted"] == 30 and d["de_novo"] == ["kid40"]
    assert ln["trios"] == 140 and ln["parents"] == 280 and ln["parents_carrying"] == 40 and ln["both_plain"] == 100 and ln["new"] == ["kid40"]
    assert ln["transmitted"] == 10 and ln["pairs"] == 40 and ln["p_half"] < 0.01
    assert abs(ln["excess_share"] - 0.5) < 1e-9 and abs(ln["rate"] - 0.5 * 40 / 280) < 1e-9 and abs(ln["expected_new"] - 100 * 0.5 * 40 / 280) < 1e-9
    assert 0 < ln["p_new"] < 0.01
    n = d["not_passed"]
    assert n["parent"]["n"] == 30 and abs(n["parent"]["mean"]) < 0.05 and abs(n["child"]["mean"]) < 0.05 and n["parent"]["within_0_3"] == 30 == n["child"]["within_0_3"]


def test_copies_from_a_paf_and_their_classes(tmp_path):
    """A complete copy, a partial copy in tandem with it on one contig, and a copy cut by a contig end."""
    def paf(ctg, clen, q0, t0, t1, ident=0.99):
        aln = t1 - t0
        return f"{ctg}\t{clen}\t{q0}\t{q0 + aln}\t+\tDJ\t{UNIT}\t{t0}\t{t1}\t{int(aln * ident)}\t{aln}\t60\n"
    lines = []
    for t0 in range(3000, 400000, 10000):                                 # the complete copy, pieces of 8 kb with 2-kb masked gaps
        lines.append(paf("ctgA", 5_000_000, 1_000_000 + t0, t0, min(t0 + 8000, UNIT)))
    for t0 in range(3000, 316000, 10000):                                 # the partial copy, 200 kb downstream on the same contig
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
    haplotype has four, one of them lacking 200-220 kb; the page's section and tables follow. Twelve genomes
    are too few for a pin of the scale: the page says so, and the calls are made all the same."""
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
    for must in ("The junction in profile", "Whole numbers of copies along the unit", "Against the HPRC release-2 assemblies of 2 cohort members", "What the assemblies get wrong",
                 "What the method does now, and what remains", "too small to pin the scale", 'id="chart-djasm"', "data/dj_hprc.tsv", "data/dj_blocks.tsv", "data/dj_calls.tsv"):
        assert must in html, must
    visible = html.split('<script id="report-data"')[0]
    assert "NaN" not in visible and "None" not in visible.replace("None of", "")
    with open(out / "data" / "dj_hprc.tsv") as fh:
        T = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert set(T) == {"HG00733", "NA19240"}
    assert float(T["HG00733"]["assembly_core"]) == 10 and T["HG00733"]["resolved"] == "True" and T["HG00733"]["haplotypes"] == "mat=5 + pat=5"
    assert float(T["NA19240"]["assembly_core"]) == 9 and float(T["NA19240"]["assembly_mean"]) < 9
    assert T["HG00733"]["copies"] in ("9", "10", "11") and T["HG00733"]["scale_f"]
    with open(out / "data" / "dj_hprc_blocks.tsv") as fh:
        B = [r for r in csv.DictReader(fh, delimiter="\t") if r["sample"] == "NA19240"]
    assert len(B) == dj.NBLOCK and float(next(r for r in B if r["block_kb"] == "200")["assembly"]) == 8 and float(next(r for r in B if r["block_kb"] == "40")["assembly"]) == 9
    # the pilot genomes read near ten, and the tables have every genome
    with open(out / "data" / "dj_blocks.tsv") as fh:
        P = list(csv.DictReader(fh, delimiter="\t"))
    assert len(P) == 12 and all(8.5 < float(r["cn"]) < 11.5 for r in P) and all(r["copies"] in ("9", "10", "11") for r in P)
    with open(out / "data" / "dj_calls.tsv") as fh:
        C = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["sample"] for r in C} == {r["sample"] for r in P} and all(int(r["start"]) < int(r["end"]) for r in C)
    with open(out / "data" / "cohort.tsv") as fh:
        K = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert K["HG00733"]["DJ.copies"] == T["HG00733"]["copies"] and "DJ.cn_unit" in K["HG00733"] and "DJ.variants" in K["HG00733"] and K["HG00733"]["DJ.call"] in ("settled", "uncertain")
    if (out / "dj_assemblies.png").exists():
        assert (out / "dj_assemblies.png").stat().st_size > 10000 and 'src="dj_assemblies.png"' in html
