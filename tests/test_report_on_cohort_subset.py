"""The page on a slice of the cohort's own counts: enough complete trios for the bootstrap intervals, the
paired comparisons and the PC sweep, so that every section the full page has is built here - the fetch
check, transmission by sex, the batch count, Hall, the coverage QC, the assemblies, ddPCR - from real data,
in about a minute."""
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCAN, FETCH = ROOT / "counts_scan", ROOT / "counts_fetch"
N_TRIOS = 24
pytestmark = pytest.mark.skipif(len(list(SCAN.glob("*.json.gz"))) < 200, reason="the cohort's counts are not present")


@pytest.fixture(scope="module")
def page(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("subset")
    have = {p.name.split(".")[0] for p in SCAN.glob("*.json.gz")} & {p.name.split(".")[0] for p in FETCH.glob("*.json.gz")}
    ped = ROOT / "meta" / "20130606_g1k_3202_samples_ped_population.txt"
    trios = []
    with open(ped) as fh:
        next(fh)
        for line in fh:
            fam, s, f, m = line.split()[:4]
            if f != "0" and m != "0" and {s, f, m} <= have:
                trios.append((s, f, m))
    chosen = sorted({x for t in trios[:N_TRIOS] for x in t})
    assert len(chosen) >= 60, "fewer than sixty genomes in the first trios: the sweep needs sixty"
    for mode, src in (("scan", SCAN), ("fetch", FETCH)):
        d = tmp / mode
        d.mkdir()
        for s in chosen:
            os.symlink(src / f"{s}.json.gz", d / f"{s}.json.gz")
    args = ["--scan", str(tmp / "scan"), "--fetch", str(tmp / "fetch"), "-p", str(ped), "--hall", str(ROOT / "meta/hall2021_MOESM1.txt"),
            "--pilot", str(ROOT / "pilot"), "--pcs", str(ROOT / "meta/ngspca/svd.pcs.txt"), "--qc", str(ROOT / "meta/ngspca_sample_qc.tsv"),
            "--censat", str(ROOT / "hprc_censat"), "--ddpcr", str(ROOT / "assembly_rdna/tables/potapova_comparison.tsv"),
            "-o", str(tmp / "out"), "--cache", str(tmp / "cache"), "-j", "4", "--as-of", "2026-09-24"]
    subprocess.run([sys.executable, "-m", "report", *args], check=True, cwd=ROOT, capture_output=True)
    return tmp / "out", len(chosen)


def test_every_section_is_built_from_real_data(page):
    out, n = page
    d = json.loads((out / "report.json").read_text())
    assert d["meta"]["n"] == n and d["meta"]["n_both"] == n
    tr = d["trios"]
    assert tr["n_complete"] >= 20 and all("R_lo" in t for t in tr["table"]) and tr["compare"]
    R = {t["column"]: t for t in tr["table"]}
    assert R["rDNA45S.cn"]["R"] > 0.5 and R["truth.auto"]["R"] < 0.5 and R["rDNA45S.cn"].get("perm_p") is not None
    assert tr.get("by_sex") and tr.get("batches", {}).get("n") == tr["n_complete"]
    fc = d["fetch_check"]
    assert fc["n"] == n and fc["agreement"]["rDNA45S.cn"]["r"] > 0.99 and {t["column"] for t in fc["trios"]["table"]} >= {"rDNA45S.cn", "truth.auto"}
    assert "sweep" in d["pcs"] and d["hall"]["n"] >= 3 and d["ngspca_qc"]["n"] == n and d["replicates"]["n"] == 12
    assert d["ddpcr"]["n"] >= 9 and d["ddpcr"]["ngsdose"]["r"] > 0.9
    assert d["known_truth"]["auto"]["n"] == n and 1.98 < d["known_truth"]["auto"]["mean"] < 2.02
    html = (out / "index.html").read_text()
    visible = html.split('<script id="report-data"')[0]
    for must in ('id="chart-heat"', 'id="chart-ddpcr"', 'id="chart-trio"', "Hall, Turner", "python -m report"):
        assert must in html, must
    assert "chart failed" not in visible and "NaN" not in visible and "None" not in visible.replace("None of", "")
    for f in ("data/cohort.tsv", "data/transmission.tsv", "data/trios.tsv", "data/fetch_check.tsv", "data/ddpcr.tsv", "data/modes.tsv"):
        assert (out / f).stat().st_size > 0, f
    rows = list(csv.DictReader(open(out / "data/transmission.tsv"), delimiter="\t"))
    assert {r["column"] for r in rows} >= {"rDNA45S.cn", "rDNA5S.cn", "HSat3.mass_Mb", "chrM.copies", "insert_median"}
