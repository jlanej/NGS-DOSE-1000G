"""`python -m report` on the committed pilot counts: the page and its tables exist, and the numbers in
them are the pilot's numbers."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "pilot"
pytestmark = pytest.mark.skipif(len(list((PILOT / "counts_nygc").glob("*.json.gz"))) < 6, reason="pilot counts not present")


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    out = tmp_path_factory.mktemp("report")
    ped = out / "ped.txt"
    # the four pilot trios, with sex and population, in the 1000 Genomes pedigree format
    lines = ["FamilyID SampleID FatherID MotherID Sex Population Superpopulation"]
    for fam, c, f, m, cs, pop, sp in (("1463", "NA12878", "NA12891", "NA12892", 2, "CEU", "EUR"), ("Y117", "NA19240", "NA19239", "NA19238", 2, "YRI", "AFR"),
                                       ("SH032", "HG00514", "HG00512", "HG00513", 2, "CHS", "EAS"), ("PR05", "HG00733", "HG00731", "HG00732", 2, "PUR", "AMR")):
        lines += [f"{fam} {c} {f} {m} {cs} {pop} {sp}", f"{fam} {f} 0 0 1 {pop} {sp}", f"{fam} {m} 0 0 2 {pop} {sp}"]
    ped.write_text("\n".join(lines) + "\n")
    subprocess.run([sys.executable, "-m", "report", "--fetch", str(PILOT / "counts_nygc"), "-p", str(ped), "--hall", str(PILOT / "hall2021_MOESM1.txt"), "--pilot", str(PILOT),
                    "--pcs", str(ROOT / "meta/ngspca/svd.pcs.txt"), "-o", str(out), "-j", "2", "--as-of", "2026-09-22"], check=True, cwd=ROOT, capture_output=True)
    return out


def test_outputs_exist_and_are_reproducible(report):
    for f in ("index.html", "report.json", "data/cohort.tsv", "data/fetch.tsv", "data/transmission.tsv", "data/flags.tsv", "data/efficiencies.json"):
        assert (report / f).stat().st_size > 0, f
    html = (report / "index.html").read_text()
    for must in ("Every number and figure on this page is recomputed", 'id="chart-auto"', 'id="chart-trio"', "Hall, Turner", "report-data", "<table"):
        assert must in html
    # every figure stands alone in a screenshot: numbered, captioned, with its source; the 18S ratio is named as the published method
    n_fig = html.count('<figure id="fig-')
    n_supp = html.count('<div class="title">Figure S')
    assert n_fig >= 5 and html.count("<figcaption>") == n_fig and html.count('<span class="src">') == n_fig
    assert all(f'<div class="title">Figure {i}. ' in html for i in range(1, n_fig - n_supp + 1))
    assert all(f'<div class="title">Figure S{i}. ' in html for i in range(1, n_supp + 1))
    assert "18S depth ratio (published" in html and "18S depth ratio (literature)" not in html
    visible = html.split('<script id="report-data"')[0]                # the prose and tables, not the code
    assert "NaN" not in visible and "None" not in visible.replace("None of", "") and "nan" not in visible.split("<style>")[1].split("</style>")[0]
    # the second run reuses the cached estimates and produces the same page
    first = html
    subprocess.run([sys.executable, "-m", "report", "--fetch", str(PILOT / "counts_nygc"), "-p", str(report / "ped.txt"), "--hall", str(PILOT / "hall2021_MOESM1.txt"), "--pilot", str(PILOT),
                    "--pcs", str(ROOT / "meta/ngspca/svd.pcs.txt"), "-o", str(report), "-j", "2", "--as-of", "2026-09-22"], check=True, cwd=ROOT, capture_output=True)
    assert (report / "index.html").read_text() == first


def test_the_numbers_are_the_pilots(report):
    d = json.loads((report / "report.json").read_text())
    assert d["meta"]["n"] == 12 and d["meta"]["primary_mode"] == "fetch" and d["meta"]["n_fetch"] == 12
    kt = d["known_truth"]
    assert abs(kt["auto"]["mean"] - 1.996) < 0.003 and abs(kt["chrX"]["M"]["mean"] - 0.995) < 0.01 and 9.6 < kt["DJ"]["mean"] < 9.8
    assert kt["chrY"]["M"]["n"] == 4 and kt["chrY"]["F"]["n"] == 8 and kt["chrY"]["F"]["max"] < 0.01 and kt["sex"]["mismatch"] == []
    assert d["trios"]["n_complete"] == 4 and {t["column"] for t in d["trios"]["table"]} >= {"rDNA45S.cn", "rDNA45S.18S.flat", "truth.auto", "chrM.copies"}
    h = d["hall"]
    assert h["n"] == 5 and h["flat"]["r"] > 0.97 and 1.05 < h["flat_ratio"] < 1.1 and 1.0 < h["dup_corrected_ratio"] < 1.05
    assert any("HG00732" == s and "chrX" in f for s, f in d["flags"])            # the culture that lost an X
    rep = d["replicates"]                                                          # the same twelve people on an older technology
    assert rep["n"] == 12 and rep["table"]["calibrated"]["icc"] > 0.95 and rep["table"]["flat"]["icc"] < 0.5 and rep["table"]["flat_centred"]["icc"] > 0.8
    assert len(d["samples"]) == 12 and all(s["sex_inferred"] in ("M", "F") for s in d["samples"])
    cols = (report / "data" / "cohort.tsv").read_text().splitlines()[0].split("\t")
    assert {"sample", "rDNA45S.cn", "truth.chrY", "chrM.copies", "flags", "sex_inferred"} <= set(cols)


def test_the_coverage_qc_comparison(report):
    """NGS-PCA's QC table for the same files: a synthetic one built from the pilot's own numbers with a
    constant duplicate-flag offset and noise must come back at r > 0.99 with the offset recovered."""
    import numpy as np

    from report import report as R
    d = json.loads((report / "report.json").read_text())
    rows = [dict(s) for s in d["samples"]]
    rng = np.random.default_rng(3)
    qc = report / "sample_qc.tsv"
    with open(qc, "w") as fh:
        fh.write("SAMPLE_ID\tMEAN_AUTOSOMAL_COV\tX_COV_RATIO\tY_COV_RATIO\tINFERRED_SEX\tMTDNA_CN\tRELEASE_BATCH\n")
        for r in rows:
            fh.write(f'{r["sample"]}\t{0.88 * r["depth"]:.3f}\t{r["truth.chrX"] / 2:.4f}\t{0.85 * r["truth.chrY"] / 2:.4f}\t{r["sex_inferred"]}\t{0.92 * r["chrM.copies"] * np.exp(rng.normal(0, 0.02)):.2f}\t2504\n')
        fh.write("NA00000\t30\t0.5\t0.4\tM\t500\t698\n")                               # a sample that was not counted
    q = R.ngspca_qc_comparison(rows, qc)
    assert q["n"] == 12 and q["n_qc"] == 13 and q["sex_agree"] == q["sex_n"] == 12
    assert q["mtdna"]["r"] > 0.99 and abs(q["mtdna"]["ratio"]["median"] - 0.92) < 0.02 and q["chrX"]["r"] > 0.999 and abs(q["depth"]["ratio"]["median"] - 0.88) < 0.01
    assert [m["sample"] for m in q["mosaic_X"]] == ["HG00732"] and rows[0]["ngspca.batch"] == "2504"


def test_a_man_with_two_x_chromosomes_is_set_apart():
    """A man whose reads show two X chromosomes and a Y (47,XXY) is not a failure of the X model: the men's
    figures are those of men with one X, and he is listed on his own."""
    from report.report import known_truth
    row = lambda s, sex, x, y: {"sample": s, "sex_inferred": sex, "truth.auto": 2.0, "truth.chrX": x, "truth.chrY": y, "DJ.cn": 10.0}
    rows = [row(f"M{i}", "M", 0.99 + 0.005 * (i % 3), 1.0) for i in range(20)] + [row(f"F{i}", "F", 1.95, 0.0) for i in range(20)]
    rows.append(row("XXY", "M", 1.97, 0.97))
    sx = known_truth(rows)["sex"]
    assert [d["sample"] for d in sx["men_extra_x"]] == ["XXY"]
    assert sx["men_one_x"]["n"] == 20 and sx["men_one_x"]["max"] < 1.01


def test_a_dj_step_is_classified_by_the_nearer_hypothesis():
    """A child nearer zero than the carrier parent's step did not receive it, even within half a copy of the step; a
    child within half a copy of neither is counted as unclassified, not dropped."""
    from report.report import dj_steps
    rows = [{"sample": f"B{i}", "DJ.cn": 10.0} for i in range(20)]
    ped = {}
    for p, pv, c, cv in (("P1", 10.75, "C1", 10.32), ("P2", 11.47, "C2", 10.59), ("P3", 11.0, "C3", 10.95)):
        rows += [{"sample": p, "DJ.cn": pv}, {"sample": c, "DJ.cn": cv}]
        ped[c] = dict(father=p, mother="0")
    d = dj_steps(rows, ped)
    assert (d["transmitted"], d["not_transmitted"], d["unclassified"], d["n_pairs"]) == (1, 1, 1, 3)
    assert [(u["parent"], u["child"]) for u in d["unclassified_pairs"]] == [("P2", "C2")]
    assert set(d["near"]) == {-2, -1, 0, 1, 2}


def _trios(n=24, seed=5):
    import numpy as np
    from ngsdose.trios import Trio
    rng = np.random.default_rng(seed)
    rows, trios, population = [], [], {}
    for i in range(n):
        f, m = rng.normal(300, 40, 2)
        c = (f + m) / 2 + rng.normal(0, 20)
        pop = ("AAA", "BBB")[i % 2]
        for s, v, sex in ((f"C{i}", c, "MF"[i % 2]), (f"F{i}", f, "M"), (f"M{i}", m, "F")):
            rows.append({"sample": s, "rDNA45S.cn": float(v), "chrEBV.copies": 0.0, "sex": sex})
            population[s] = pop
        trios.append(Trio(f"C{i}", f"F{i}", f"M{i}", pop))
    return rows, trios, population


def test_a_trio_metric_with_no_variance_among_the_parents_is_listed_not_tested():
    """chrEBV in blood-derived DNA is zero in everyone: no slope, no reliability, and no row among the tested metrics."""
    from report.report import transmission_by_sex, trio_analysis
    rows, trios, population = _trios()
    cols = [("rDNA45S.cn", "45S"), ("chrEBV.copies", "EBV episomes per cell")]
    t = trio_analysis(rows, trios, population, cols)
    assert [r["column"] for r in t["table"]] == ["rDNA45S.cn"] and t["scatter_column"] == "rDNA45S.cn"
    assert t["constant"] == [dict(column="chrEBV.copies", label="EBV episomes per cell", group="culture", n_trios=24, reason="no variance among the parents")]
    b = transmission_by_sex(rows, trios, population, cols)
    assert "chrEBV.copies" not in b["table"] and [c["column"] for c in b["constant"]] == ["chrEBV.copies"]


def test_the_trio_total_counts_the_pedigrees_own_rows():
    """One 1000 Genomes pedigree row names a father with no row of his own (HG02567): the release has 602 trios, although
    the population map, which gives that father his child's population, would complete 603."""
    from report.report import load_pedigree, trio_analysis
    ped = ROOT / "meta" / "20130606_g1k_3202_samples_ped_population.txt"
    if not ped.exists():
        pytest.skip("pedigree not present")
    info, trio_list, population = load_pedigree(ped)
    assert len(info) == 3202 and "HG02567" not in info
    assert trio_analysis([], trio_list, population, [], sequenced=set(info))["n_total"] == 602


def test_singular_values_in_any_float_notation(tmp_path):
    from report.report import load_singular_values
    f = tmp_path / "svd.singularvalues.txt"
    f.write_text("SINGULAR_VALUES\n338.7\n1.2e+03\n9.5E-4\n3.387055774857011E2\n\n")
    assert load_singular_values(f) == [338.7, 1200.0, 0.00095, 338.7055774857011]
    f.write_text("SINGULAR_VALUES\n338.7\n12,0\n")
    with pytest.raises(ValueError, match="line 3"):
        load_singular_values(f)


def _copy_counts(dst, names):
    import shutil
    dst.mkdir(parents=True, exist_ok=True)
    for src, name in names:
        shutil.copy2(PILOT / "counts_nygc" / f"{src}.json.gz", dst / name)
    return dst


def test_the_cache_follows_the_bundle_and_the_counts(tmp_path):
    """A cached estimate is reused only while the counts file (size, modification time), the ngsdose code and the bundle's
    files are those it was made with; a change to bundle.json's descriptive text alone keeps it."""
    import os
    import shutil

    from ngsdose import resources
    from ngsdose.tables import load_result

    from report import report as R
    bundle = tmp_path / "bundle"
    shutil.copytree(resources.default_bundle(), bundle, ignore=shutil.ignore_patterns("build_inputs"))
    d = _copy_counts(tmp_path / "counts", [("HG00512", "HG00512.json.gz")])
    paths, cache, logs = sorted(d.glob("*.json.gz")), tmp_path / "cache", []
    run = lambda: R.estimate_all(paths, cache, bundle, 1, logs.append)
    rows, bad = run()
    assert bad == [] and "0 cached, 1 new" in logs[-1] and rows[0]["_estimate"] == str(cache / "HG00512.estimate.json.gz")
    cn = rows[0]["rDNA45S.cn_single"]
    run()
    assert "1 cached, 0 new" in logs[-1]
    meta = json.loads((bundle / "bundle.json").read_text())
    meta["reference"] += " (reworded)"
    (bundle / "bundle.json").write_text(json.dumps(meta))
    run()
    assert "1 cached, 0 new" in logs[-1]
    anchors = json.loads((bundle / "anchors.json").read_text())
    anchors["rDNA45S"]["intervals"] = [[1000, 1250], [2000, 2250], [40000, 40500]]
    (bundle / "anchors.json").write_text(json.dumps(anchors))
    rows, _ = run()
    assert "cache invalidated for 1 (ngsdose or bundle changed)" in logs[-1] and rows[0]["rDNA45S.cn_single"] != cn
    assert load_result(cache / "HG00512.estimate.json.gz")["report_extras"]["fingerprint"] == R.fingerprint(bundle)["key"]
    st = os.stat(paths[0])
    os.utime(paths[0], ns=(st.st_atime_ns, st.st_mtime_ns - 10**9))              # an older copy of a recount: older, but not the same
    run()
    assert "cache invalidated for 1 (counts file changed)" in logs[-1]


def test_an_unreadable_counts_file_is_left_out_and_named(tmp_path):
    from report import report as R
    d = _copy_counts(tmp_path / "counts", [("HG00512", "HG00512.json.gz"), ("HG00513", "HG00513.json.gz")])
    (d / "HG00513.json.gz").write_bytes((d / "HG00513.json.gz").read_bytes()[:20000])        # a partial copy
    (d / "junk.json.gz").write_text("not gzip")
    logs = []
    rows, bad = R.estimate_all(sorted(d.glob("*.json.gz")), tmp_path / "cache", None, 1, logs.append)
    assert [r["sample"] for r in rows] == ["HG00512"]
    assert sorted((b["counts"], b["stage"]) for b in bad) == [(str(d / "HG00513.json.gz"), "read"), (str(d / "junk.json.gz"), "read")]
    assert any("HG00513.json.gz" in m and "EOFError" in m for m in logs) and "2 left out" in logs[-1]
    (d / "HG00512.json.gz").unlink()
    with pytest.raises(SystemExit, match="could be used"):
        R.estimate_all(sorted(d.glob("*.json.gz")), tmp_path / "cache", None, 1, logs.append)


def test_estimates_are_found_through_their_file_and_a_sample_counted_twice_stops_the_run(tmp_path):
    """A counts file need not be named after its sample (a biobank's <eid>_23193_0_0); an unreadable one is listed in
    meta.unreadable; two files of one mode naming the same sample stop the run, naming both."""
    names = [("HG00512", "HG00512.json.gz"), ("HG00513", "HG00513.json.gz"), ("HG00514", "HG00514.json.gz"), ("NA12878", "NA12878_23193_0_0.json.gz")]
    d = _copy_counts(tmp_path / "fetch", names)
    (d / "HG00731.json.gz").write_bytes((PILOT / "counts_nygc" / "HG00731.json.gz").read_bytes()[:20000])
    args = [sys.executable, "-m", "report", "--fetch", str(d), "-o", str(tmp_path / "out"), "--cache", str(tmp_path / "cache"), "-j", "1", "--as-of", "2026-09-22"]
    subprocess.run(args, check=True, cwd=ROOT, capture_output=True)
    m = json.loads((tmp_path / "out" / "report.json").read_text())["meta"]
    assert m["n"] == 4 and "NA12878" in {s["sample"] for s in json.loads((tmp_path / "out" / "report.json").read_text())["samples"]}
    assert [(u["mode"], u["stage"], u["counts"]) for u in m["unreadable"]] == [("fetch", "read", str(d / "HG00731.json.gz"))]
    _copy_counts(d, [("NA12878", "NA12878.rerun.json.gz")])
    p = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    assert p.returncode != 0 and "NA12878 is in two fetch counts files" in p.stderr and "NA12878_23193_0_0.json.gz" in p.stderr and "NA12878.rerun.json.gz" in p.stderr


def test_a_correlation_lies_inside_its_own_interval():
    """The Fisher interval is computed on r itself, not on r clipped; at |r| = 1 (to within rounding) there is none."""
    import math

    import numpy as np

    from report.report import corr
    x = np.arange(1.0, 101.0)
    d = corr(x, 2 * x)
    assert d["r"] > 1 - 1e-12 and math.isnan(d["r_lo"]) and math.isnan(d["r_hi"])
    d = corr(x, x + 1e-2 * np.sin(x))
    assert 0.999999 < d["r"] < 1 and d["r_lo"] <= d["r"] <= d["r_hi"]


def test_a_column_with_too_few_usable_values_is_left_unadjusted(tmp_path):
    """chrEBV at 0 in everyone (blood-derived DNA) has no positive value to regress on the PCs: its .adj and
    .adj_ngspca stay NA and it is listed as skipped, and the other columns are still adjusted."""
    import numpy as np

    from report.report import pc_analysis
    n, n_pc = 60, 3
    rows = [{"sample": f"S{i}", "chrEBV.copies": 0.0, "chrM.copies": (300.0 + i if i < 5 else 0.0), "truth.auto": 2 + 0.01 * np.sin(3 * i),
             "ctrlPC1": np.sin(i), "ctrlPC2": np.cos(i)} for i in range(n)]
    rng = np.random.default_rng(1)
    (tmp_path / "svd.pcs.txt").write_text("SAMPLE\t" + "\t".join(f"PC{j + 1}" for j in range(n_pc)) + "\n"
                                          + "".join(f"S{i}.by1000.\t" + "\t".join(f"{v:.6f}" for v in rng.normal(size=n_pc)) + "\n" for i in range(n)))
    (tmp_path / "svd.singularvalues.txt").write_text("\n".join(str(v) for v in [400.0, 200.0, 100.0] + list(np.linspace(20, 10, 57))) + "\n")
    (tmp_path / "svd.bins.txt").write_text("".join(f"chr1\t{j * 1000}\t{(j + 1) * 1000}\n" for j in range(5000)))
    logged = []
    out = pc_analysis(rows, {"mp": 2}, [], {}, str(tmp_path / "svd.pcs.txt"), logged.append)
    adj = out["adjusted"]
    assert {d["column"] for d in adj["skipped"]} == {"chrEBV.copies", "chrM.copies"}
    assert "0 usable values" in next(d["reason"] for d in adj["skipped"] if d["column"] == "chrEBV.copies")
    assert "truth.auto" in adj["columns"] and "chrEBV.copies" not in adj["columns"]
    assert all(r.get("chrEBV.copies.adj") is None for r in rows) and all(r["truth.auto.adj"] is not None for r in rows)
    assert any("chrEBV.copies not adjusted" in m for m in logged)
    ng = out["ngspca"]
    assert ng["n_pc"] == n_pc and {d["column"] for d in ng["skipped"]} == {"chrEBV.copies", "chrM.copies"}
    assert all(r.get("chrEBV.copies.adj_ngspca") is None for r in rows) and all(r["truth.auto.adj_ngspca"] is not None for r in rows)
