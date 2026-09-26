"""The page, the trio PDF and the evidence figure say what the numbers say: every qualitative sentence follows the comparison it
makes, a metric with nothing to test is listed as not tested instead of breaking the page, and every step, pair and path the
data hold is shown. Built once from a slice of the cohort's own counts (24 trios), then re-rendered from altered copies of it."""
import copy
import html as H
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCAN, FETCH = ROOT / "counts_scan", ROOT / "counts_fetch"
PED = ROOT / "meta" / "20130606_g1k_3202_samples_ped_population.txt"
N_TRIOS = 24
have_counts = len(list(SCAN.glob("*.json.gz"))) >= 200
needs_counts = pytest.mark.skipif(not have_counts, reason="the cohort's counts are not present")


def text(page_html: str) -> str:
    """The visible prose of a page, tags stripped."""
    vis = page_html.split('<script id="report-data"')[0].split("</style>")[1]
    return H.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", vis)))


def args_for(tmp, out):
    return ["--scan", str(tmp / "scan"), "--fetch", str(tmp / "fetch"), "-p", str(PED), "--hall", str(ROOT / "meta/hall2021_MOESM1.txt"),
            "--pilot", str(ROOT / "pilot"), "--pcs", str(ROOT / "meta/ngspca/svd.pcs.txt"), "--qc", str(ROOT / "meta/ngspca_sample_qc.tsv"),
            "--censat", str(ROOT / "hprc_censat"), "--ddpcr", str(ROOT / "assembly_rdna/tables/potapova_comparison.tsv"),
            "-o", str(out), "--cache", str(tmp / "cache"), "--as-of", "2026-09-24"]


@pytest.fixture(scope="module")
def subset(tmp_path_factory):
    """The first 24 complete trios, scanned and fetched; returns (tmp dir, data, rows) as `render` received them."""
    if not have_counts:
        pytest.skip("the cohort's counts are not present")
    tmp = tmp_path_factory.mktemp("page")
    have = {p.name.split(".")[0] for p in SCAN.glob("*.json.gz")} & {p.name.split(".")[0] for p in FETCH.glob("*.json.gz")}
    trios = []
    with open(PED) as fh:
        next(fh)
        for line in fh:
            fam, s, f, m = line.split()[:4]
            if f != "0" and m != "0" and {s, f, m} <= have:
                trios.append((s, f, m))
    chosen = sorted({x for t in trios[:N_TRIOS] for x in t})
    for mode, src in (("scan", SCAN), ("fetch", FETCH)):
        (tmp / mode).mkdir()
        for s in chosen:
            os.symlink(src / f"{s}.json.gz", tmp / mode / f"{s}.json.gz")
    # the estimates in a separate process (its worker pool), then the build again in this one from the cache, to keep what render gets
    subprocess.run([sys.executable, "-m", "report", *args_for(tmp, tmp / "first"), "-j", "4"], check=True, cwd=ROOT, capture_output=True)
    import report.report as R
    got = {}
    render = R.render

    def keep(data, rows):
        got.update(data=copy.deepcopy(data), rows=copy.deepcopy(rows))
        return render(data, rows)
    R.render = keep
    try:
        R.main(args_for(tmp, tmp / "out") + ["-j", "1"])
    finally:
        R.render = render
    return tmp, got["data"], got["rows"]


@needs_counts
def test_a_constant_trio_metric_is_listed_not_tested_and_the_page_is_built(subset, tmp_path):
    """chrEBV in blood-derived DNA is zero in everyone: the trio tables send it to `constant`, and the page, the PDF and the figure
    are built from them, naming the metric as not tested rather than showing it with empty statistics."""
    import report.report as R
    from ngsdose.tables import write_table
    from report.report_page import page
    _, data, rows = subset
    rows = copy.deepcopy(rows)
    for r in rows:
        r["chrEBV.copies"] = 0.0
    ped, trio_list, population = R.load_pedigree(PED)
    d = copy.deepcopy(data)
    tr = R.trio_analysis(rows, trio_list, population, R.TRIO_COLUMNS, set(ped))
    tr["by_sex"] = R.transmission_by_sex(rows, trio_list, population, R.TRIO_COLUMNS)
    tr["points"], tr["batches"] = R.trio_points(tr.pop("values"), rows, R.TRIO_COLUMNS), d["trios"].get("batches")
    d["trios"] = tr
    d = R.rnd(d, 5)
    assert [c["column"] for c in d["trios"]["constant"]] == ["chrEBV.copies"] and "chrEBV.copies" not in d["trios"]["by_sex"]["table"]
    t = text(page(d, rows))
    assert "Not tested (no slope): EBV episomes per cell (no variance among the parents)." in t
    assert "Not tested by sex: EBV episomes per cell" in t
    assert "NaN" not in t and "None" not in t.replace("None of", "")
    (tmp_path / "data").mkdir()
    write_table(d["trios"]["table"], tmp_path / "data" / "transmission.tsv")
    (tmp_path / "report.json").write_text(json.dumps(d))
    for mod, extra in (("report.trio_report", ["--data", str(tmp_path / "data"), "-o", str(tmp_path / "trio.pdf")]),
                       ("report.evidence_figure", ["--pilot", str(ROOT / "pilot"), "-o", str(tmp_path / "evidence.png")])):
        subprocess.run([sys.executable, "-m", mod, "--report", str(tmp_path / "report.json"), *extra], check=True, cwd=ROOT, capture_output=True)
    assert (tmp_path / "trio.pdf").stat().st_size > 0 and (tmp_path / "evidence.png").stat().st_size > 0


@needs_counts
def test_what_the_page_shows_is_what_the_data_hold(subset):
    from report.report_page import page
    _, data, rows = subset
    html = page(copy.deepcopy(data), rows)
    t = text(html)
    charts = json.loads(html.split('<script id="report-data" type="application/json">')[1].split("</script>")[0].replace("<\\/", "</"))["charts"]
    # the by-sex heatmap holds correlations, and its hover label says so
    assert all(c.startswith("child–parent correlation (Pearson r)") for c in charts["heat_sex"]["col_titles"])
    # the S-phase figure draws the women its caption and r describe
    sph = data["biology"]["DJ_vs_chrX_female"]
    assert "where" not in charts["sphase"] and len(charts["sphase"]["points"]) == sph["n"]
    assert all(1.85 <= p["x"] <= 2.15 for p in charts["sphase"]["points"])
    # figures name the repository that holds their data; the reproduce block is the one regenerate.sh runs
    assert "github.com/jlanej/NGS-DOSE-1000G" in t
    reg, repro = (ROOT / "regenerate.sh").read_text(), t.split("Data and reproducibility")[1]
    for flag in ("--scan", "--fetch", "--hall", "--pilot", "--pcs", "--censat", "--qc", "--ddpcr", "--cache", "meta/20130606_g1k_3202_samples_ped_population.txt"):
        assert flag in reg and flag in repro, flag
    assert "regenerate.sh" in t and "results repository's" not in t and "the 2,504 unrelated samples)" not in t
    # nothing claims more than the data: no unqualified 'minute' or '0.5 GB', no 'nothing is typed in'
    assert "One minute of the file" not in t and "nothing is typed in" not in t and "within what a dozen lines can resolve" not in t
    assert "0.5 GB of a 15-GB CRAM" not in t or "by the pipeline's estimate" in t
    # the path table: every class the counts carry, with the path the data give it
    paths = {r["class"]: r["path"] for r in data["fetch_paths"]["rows"]}
    assert paths["rDNA45S"] == "fetch-direct" and paths["HSat2"] == "scan-only"
    table = html.split('id="paths"')[1].split("</table>")[0]
    assert all(f"<td class=\"\">{c}</td>" in table for c in paths)
    assert "no sinks in the bundle" in t and "DRAGEN" in t
    # the fetch's sinks, from the sinks files the counts record
    f0 = data["meta"]["sinks"]["fetch"][0]
    assert f"{f0['intervals']:,} intervals" in t and f0["sha256"] in t
    # every whole-copy step of the distal junction that a genome reaches is listed, and the pairs add up
    dj = data["known_truth"]["DJ_steps"]
    if dj.get("n_pairs") and any(int(k) and v for k, v in dj["near"].items()):
        assert f"transmitted in {dj['transmitted']} of {dj['n_pairs']}" in t


@needs_counts
def test_the_culture_paragraph_counts_its_intervals(subset):
    """Two culture metrics with intervals wholly below zero, one far from it with a spousal correlation: no 'one such interval is chance'."""
    from report.report_page import page
    _, data, rows = subset
    d = copy.deepcopy(data)
    by = {t["column"]: t for t in d["trios"]["table"]}
    by["depth"].update(R=-0.4, R_lo=-0.6, R_hi=-0.2, spousal_r=0.3, spousal_lo=0.1, spousal_hi=0.45, perm_p=0.9)
    by["gc_rel_65"].update(R=-0.2, R_lo=-0.39, R_hi=-0.01, perm_p=0.7)
    t = text(page(d, rows))
    para = t.split("Culture and library.")[1].split("A printable assessment")[0]
    assert "one such interval" not in para and "that of sequencing depth" in para and "that of library GC bias" in para
    assert "standard errors below zero, which chance does not give" in para and "spousal correlation is 0.30" in para
    # with none below zero, the paragraph says every interval includes zero
    d = copy.deepcopy(data)
    for t_ in d["trios"]["table"]:
        if t_["group"] == "culture":
            t_.update(R=0.05, R_lo=-0.2, R_hi=0.3, perm_p=0.5)
    para = text(page(d, rows)).split("Culture and library.")[1]
    assert "every interval includes zero" in para and "shuffling children among families" not in para


@needs_counts
def test_qualitative_sentences_follow_their_numbers(subset):
    from report.report_page import page
    _, data, rows = subset
    d = copy.deepcopy(data)
    # a +2 step and an unclassified pair appear in the counts the page prints
    dj = d["known_truth"]["DJ_steps"]
    dj["near"] = {"-2": 0, "-1": 3, "0": 60, "1": 2, "2": 2}
    dj.update(transmitted=2, not_transmitted=1, unclassified=1, n_pairs=4,
              unclassified_pairs=[dict(parent="HGX1", parent_step=1.47, child="HGX2", child_step=0.59)])
    dj["carriers"] = dj.get("carriers") or [dict(sample="HGX1", pop="CEU", sex="F", step=1.47, relatives=[])]
    # mitochondrial content varies more than the 45S here
    d["biology"].update(cn45_cv=0.2, chrM_cv=0.25)
    t = text(page(d, rows))
    assert "3 genomes at −1, 60 at 0, 2 at +1 and 2 at +2" in t and "at −2" not in t.split("Within ±0.3 of a step:")[1].split(";")[0]
    assert "transmitted in 2 of 4 (1 unclassified)" in t and "HGX1 (+1.47) and HGX2 (+0.59) fits neither" in t
    assert "mitochondrial content more than it" in t
    # the satellite spread sentence is written only when the least varying classes do read the lowest R
    low_claim = "read the lowest R" in t
    sats = [x for x in data["trios"]["table"] if x["group"] == "satellites" and x["column"] not in ("TEL.mass_Mb", "HSat1B.mass_Mb")]
    assert low_claim or "R does not follow the spread between people" in t or len(sats) < 3


@needs_counts
def test_the_methods_and_assembly_text_describe_the_current_rules(subset):
    """The control-QC rule, the DJ k-mer rule, the two gap conventions and placeholder gaps, DRAGEN by version and the hedged S-phase
    and Hall readings, as the code and the data now have them."""
    from ngsdose.estimate import MIN_CHROM_REGIONS
    from report.report_page import page
    _, data, rows = subset
    d = copy.deepcopy(data)
    base = dict(ratio_median=1.0, sd_log=0.04, sd_log_robust=0.04, cv_assembly=0.04, pearson=0.6, spearman=0.6, n_far=0, n_far_assembly_short=0)
    d["satellites"]["hprc"] = dict(n_samples=10, stats={"aSatHOR": dict(base, n=10, n_gapped=1, n_unsized_gap=8), "HSat3": dict(base, n=10, n_gapped=0, n_unsized_gap=1)},
                                   rows=[dict(sample=f"S{i}", cls=c, assembly_Mb=140.0 + i, assembly_gapped_Mb=0.0, ngsdose_Mb=141.0 + i) for i in range(10) for c in ("aSatHOR", "HSat3")])
    d["biology"]["DJ_vs_chrX_female"] = dict(n=628, r=-0.06, r_lo=-0.14, r_hi=0.02)
    t = text(page(d, rows))
    assert f"fewer than {MIN_CHROM_REGIONS} regions" in t and "taken without the chromosomes already flagged" in t
    assert "at least five remaining regions" not in t
    assert "five times in CHM13, all inside the five distal junctions (once in chr21's)" in t and "exactly once in each" not in t
    # both gap conventions, placeholder gaps counted and shown, and no claim that gaps are immaterial
    assert "standalone GAP record" in t and "placeholder gaps" in t
    assert "8 of the 10 have a 100-bp placeholder gap" in t and "1 of the 10 has a 100-bp placeholder gap" in t
    assert "Placeholder gaps are common next to aSatHOR arrays (8 of 10 compared assemblies)" in t
    assert "immaterial" not in t and "without a marked gap" not in t
    assert "DRAGEN 3.7.6" in t and "are untested; in the one" not in t
    assert t.count("3.7 family that UK Biobank and All of Us are thought to use (3.7.8; not checked against their headers)") == 2
    assert "the version UK Biobank" not in t and "the family UK Biobank and All of Us use" not in t
    assert "weakens a shared S-phase explanation without ruling out a small one" in t and "does not explain them" not in t
    if (d.get("hall") or {}).get("dup_denominator_ratio") is not None and d["hall"].get("n", 0) >= 3:
        assert "not their stated method" in t and "leaves duplicate-flagged reads out" not in t


def test_the_evidence_figure_bins_every_step():
    import numpy as np
    from report.evidence_figure import step_bins, step_counts
    steps = np.array([-2.1, -1.0, 0.0, 0.1, 1.0, 1.94, 2.0])
    h, _ = np.histogram(steps, bins=step_bins(steps))
    assert h.sum() == len(steps)
    assert step_counts({"-2": 4, "-1": 32, "0": 1030, "1": 36, "2": 2, "3": 0}) == [(-2, 4), (-1, 32), (0, 1030), (1, 36), (2, 2)]


def _report_json():
    p = ROOT / "docs" / "report.json"
    if not p.exists():
        pytest.skip("docs/report.json is not present")
    return json.loads(p.read_text())


def test_the_pdf_summary_survives_families_it_cannot_judge():
    """Early in a run no satellite family has four assembled people, or none that differ by 10%: the summary says so."""
    from report.trio_report import summary_blocks
    d = _report_json()
    for change in (dict(n=3), dict(cv_assembly=0.05)):
        e = copy.deepcopy(d)
        for st in e["satellites"]["hprc"]["stats"].values():
            st.update(change)
        res = dict(summary_blocks(e, [], []))["Results"]
        assert "Assemblies:" in res


def test_the_pdf_counts_men_with_one_x_and_lists_the_others():
    from report.trio_report import truth_rows
    d = _report_json()
    sx = d["known_truth"]["sex"]
    sx["men_extra_x"] = [dict(sample="HGXXY", chrX=1.97, chrY=0.97)]
    sx["men_one_x"] = dict(n=9, mean=0.99, sd=0.007)
    rows = truth_rows(d)
    assert rows[1][0] == "chrX in men with one X (60 regions)" and rows[1][3] == "9" and rows[1][2] == "0.990 ± 0.007"
    assert rows[2][0].startswith("chrX / chrY, HGXXY") and rows[2][2] == "1.97 / 0.97"


def test_the_pdf_is_built_below_twenty_trios(tmp_path):
    """The pilot's four trios: no bootstrap, so no interval columns in transmission.tsv; the PDF is still written."""
    pilot = ROOT / "pilot"
    if len(list((pilot / "counts_nygc").glob("*.json.gz"))) < 6:
        pytest.skip("pilot counts not present")
    ped = tmp_path / "ped.txt"
    lines = ["FamilyID SampleID FatherID MotherID Sex Population Superpopulation"]
    for fam, c, f, m, pop, sp in (("1463", "NA12878", "NA12891", "NA12892", "CEU", "EUR"), ("Y117", "NA19240", "NA19239", "NA19238", "YRI", "AFR"),
                                  ("SH032", "HG00514", "HG00512", "HG00513", "CHS", "EAS"), ("PR05", "HG00733", "HG00731", "HG00732", "PUR", "AMR")):
        lines += [f"{fam} {c} {f} {m} 2 {pop} {sp}", f"{fam} {f} 0 0 1 {pop} {sp}", f"{fam} {m} 0 0 2 {pop} {sp}"]
    ped.write_text("\n".join(lines) + "\n")
    out = tmp_path / "out"
    subprocess.run([sys.executable, "-m", "report", "--fetch", str(pilot / "counts_nygc"), "-p", str(ped), "--pilot", str(pilot), "-o", str(out),
                    "-j", "2", "--as-of", "2026-09-22"], check=True, cwd=ROOT, capture_output=True)
    assert "R_lo" not in (out / "data" / "transmission.tsv").read_text().splitlines()[0].split("\t")
    subprocess.run([sys.executable, "-m", "report.trio_report", "--report", str(out / "report.json"), "--data", str(out / "data"), "-o", str(out / "trio.pdf")],
                   check=True, cwd=ROOT, capture_output=True)
    assert (out / "trio.pdf").stat().st_size > 0
