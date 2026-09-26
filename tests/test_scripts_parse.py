"""Every shell script in the pipeline and the pilot parses; and the pipeline's guards hold, run with stub
sbatch / squeue / curl / sleep on PATH: a refused sbatch leaves no dispatch mark, the stager ends (rather than
waits for ever) when a failed job has kept the only CRAM the disk has room for, but neither it nor 01a takes a
job for ended while squeue does not answer, a cut download never takes its final name, a fetch that does not
match its scan is caught, counts are not removed when there is no CRAM to make them again from, and the default
FETCH_PANELS leaves out a panel the bundle has no sinks for."""
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PIPE = ROOT / "pipeline"
needs_bash = pytest.mark.skipif(shutil.which("bash") is None or shutil.which("gzip") is None, reason="needs bash and gzip")

sys.path.insert(0, str(PIPE))
import check_counts  # noqa: E402


@needs_bash
def test_shell_scripts_parse():
    scripts = sorted(list((ROOT / "pipeline").glob("*.sh")) + list((ROOT / "pilot").glob("*.sh")) + [ROOT / "regenerate.sh"])
    assert len(scripts) >= 10
    for sh in scripts:
        subprocess.run(["bash", "-n", str(sh)], check=True)


def _gz(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def test_check_counts(tmp_path):
    panel = tmp_path / "tel.tsv.gz"
    _gz(panel, "##ngs-dose-panel v1\n##k=31\n##class\tid=0\tname=TEL\tkind=compositional\tlength=1200\n#kmer\tclass\tpos\tstrand\nAAA\t0\t0\t+\n")
    old, new = tmp_path / "old.bed", tmp_path / "new.bed"
    old.write_text("chr21\t0\t10\tDJ\nchrUn\t0\t10\trDNA45S\n")
    new.write_text("track name=x\nchr21\t0\t10\tDJ\nchr1\t0\t21000\tTEL\n")
    assert check_counts.main(["fetch-panels", str(old), str(panel)]) == 1
    assert check_counts.main(["fetch-panels", str(new), str(panel)]) == 0

    regions = [{"name": "chr1:1-100", "role": "control", "obs": 7}, {"name": "chrM:1-100", "role": "dosage", "obs": 90}]
    scan, fetch, dj = tmp_path / "s.json.gz", tmp_path / "f.json.gz", tmp_path / "dj.bed"
    dj.write_text("chr21\t0\t1000\tDJ\nchrUn\t0\t900\tDJ\n")        # chrUn's bin ends at the contig's end, 900
    _gz(scan, json.dumps(_scan(regions)))
    _gz(fetch, json.dumps(_fetch(regions[::-1], dj, dj=100)))
    assert check_counts.main(["same-reads", str(scan), str(fetch)]) == 0
    assert check_counts.main(["same-reads", str(scan), str(fetch), str(dj)]) == 0
    assert check_counts.main(["same-reads", str(scan), str(fetch), str(new)]) == 0              # another sinks file: not checked
    _gz(fetch, json.dumps(_fetch(regions, dj, dj=60)))                    # chrUn's DJ reads are gone (an index without it)
    assert check_counts.main(["same-reads", str(scan), str(fetch)]) == 0
    assert check_counts.main(["same-reads", str(scan), str(fetch), str(dj)]) == 1
    _gz(fetch, json.dumps(_fetch(regions, dj, dj=121)))                   # more than the scan has
    assert check_counts.main(["same-reads", str(scan), str(fetch), str(dj)]) == 1
    co = _scan(regions)                                                   # a scan that also loaded HSat2 may have given it
    co["classes"].append({"name": "HSat2", "kind": "compositional", "reads": 5})   # a read the fetch gives to DJ: a note
    _gz(scan, json.dumps(co))
    assert check_counts.main(["same-reads", str(scan), str(fetch), str(dj)]) == 0
    _gz(fetch, json.dumps(_fetch(regions, dj, dj=60)))                    # but fewer than inside the sinks still fails
    assert check_counts.main(["same-reads", str(scan), str(fetch), str(dj)]) == 1
    _gz(scan, json.dumps(_scan(regions)))
    _gz(fetch, json.dumps(_fetch([regions[0], dict(regions[1], obs=0)], dj, dj=100)))   # a fetch through a cut index
    assert check_counts.main(["same-reads", str(scan), str(fetch)]) == 1
    fetch.write_bytes(fetch.read_bytes()[:20])
    assert check_counts.main(["same-reads", str(scan), str(fetch)]) == 2


def _scan(regions, build="fae112470d74"):
    """120 DJ reads: 60 in chr21's sink bin, 40 in chrUn's, 20 outside the sinks."""
    return {"engine_build": build, "ctrl_reads": 97, "regions": regions, "placement_bin": 1000, "placement_bin_compositional": 10000,
            "contigs": [{"name": "chr21", "len": 1000}, {"name": "chrUn", "len": 900}, {"name": "chr22", "len": 5000}],
            "classes": [{"name": "DJ", "kind": "positional", "reads": 120}],
            "placements": [{"class": "DJ", "contig": "chr21", "start": 0, "reads": 60}, {"class": "DJ", "contig": "chrUn", "start": 0, "reads": 40},
                           {"class": "DJ", "contig": "chr22", "start": 0, "reads": 20}]}


def _fetch(regions, sinks, dj, build="fae112470d74"):
    return {"engine_build": build, "ctrl_reads": 97, "regions": regions, "classes": [{"name": "DJ", "kind": "positional", "reads": dj}],
            "sinks_sha256": hashlib.sha256(sinks.read_bytes()).hexdigest()}


STUBS = {
    "sbatch": """#!/usr/bin/env bash
q=$(cat "$STUB/quota" 2>/dev/null || echo 99)
[ "$q" -le 0 ] && { echo "sbatch: error: QOSMaxSubmitJobPerUserLimit" >&2; exit 1; }
echo $((q - 1)) > "$STUB/quota"; n=$(( $(cat "$STUB/next" 2>/dev/null || echo 1000) + 1 )); echo $n > "$STUB/next"; echo $n
""",
    "squeue": """#!/usr/bin/env bash
[ -e "$STUB/sqfail" ] && { echo "slurm_load_jobs error: Socket timed out on send/recv operation" >&2; exit 1; }
cat "$STUB/alive" 2>/dev/null; exit 0
""",                                                                     # lists the jobs in $STUB/alive
    "sleep": "#!/usr/bin/env bash\nexit 0\n",
    "curl": """#!/usr/bin/env bash
[ "$1" = --help ] && exit 0
out=""; url=""; while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift ;; -*) ;; *) url=$1 ;; esac; shift; done
src="$STUB/src/$(basename "$url")"; [ -e "$src" ] || exit 22
if [ -e "$src.cut" ]; then head -c 10 "$src" > "$out"; exit 18; fi
cp "$src" "$out"
""",
}


def _harness(tmp_path):
    stubs, stub = tmp_path / "stubs", tmp_path / "stub"
    stubs.mkdir()
    (stub / "src").mkdir(parents=True)
    for name, text in STUBS.items():
        (stubs / name).write_text(text)
        (stubs / name).chmod(0o755)
    work = tmp_path / "work"
    (work / "crams").mkdir(parents=True)
    (work / "dispatched").mkdir()
    env = dict(os.environ, PATH=f"{stubs}:{os.environ['PATH']}", STUB=str(stub), WORK_DIR=str(work), CRAM_DIR=str(work / "crams"),
               FETCH_PANELS="", SIF="", NGSDOSE_BIN="/bin/false")
    return env, stub, work


def _cram(path):
    path.write_bytes(b"CRAM" + os.urandom(64))
    _gz(Path(f"{path}.crai"), "0\t1\t100\t0\t0\t0\n")
    return hashlib.md5(path.read_bytes()).hexdigest()


@needs_bash
def test_refused_sbatch_leaves_no_mark(tmp_path):
    env, stub, work = _harness(tmp_path)
    lines = []
    for s in ("S1", "S2", "S3"):
        lines.append(f"{s}\thttps://x/{s}.cram\t{_cram(work / 'crams' / f'{s}.cram')}\n")
    (work / "manifest.tsv").write_text("".join(lines))
    (stub / "quota").write_text("1")
    r = subprocess.run(["bash", str(PIPE / "01a_dispatch_staged.sh")], env=env, capture_output=True, text=True)
    assert r.returncode == 1 and "QOSMaxSubmitJobPerUserLimit" in r.stderr
    assert (work / "dispatched" / "S1").read_text().strip() == "1001" and not (work / "dispatched" / "S2").exists()
    (stub / "quota").write_text("9")
    r = subprocess.run(["bash", str(PIPE / "01a_dispatch_staged.sh")], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert all((work / "dispatched" / s).read_text().strip() for s in ("S1", "S2", "S3"))      # S1's job ended without counts: again


@needs_bash
def test_stager_ends_when_a_kept_cram_fills_the_disk(tmp_path):
    env, stub, work = _harness(tmp_path)
    md5 = _cram(work / "crams" / "FAILED1.cram")
    shutil.copy(work / "crams" / "FAILED1.cram", stub / "src" / "NEXT2.cram")
    (work / "dispatched" / "FAILED1").write_text("111\n")
    (work / "manifest.tsv").write_text(f"FAILED1\thttps://x/FAILED1.cram\t{md5}\nNEXT2\thttps://x/NEXT2.cram\t{md5}\n")
    r = subprocess.run(["bash", str(PIPE / "01_stage_and_dose.sh")], env=dict(env, MAX_LOCAL_CRAMS="1"), capture_output=True, text=True, timeout=60)
    assert r.returncode == 1 and "STOPPED" in r.stderr and "FAILED1" in r.stderr, r.stdout + r.stderr
    assert "[FAILED1] staged, job 1001" in r.stdout                         # submitted again before anything waited
    assert not (work / "crams" / "NEXT2.cram").exists()


@needs_bash
def test_cut_download_never_takes_its_name(tmp_path):
    env, stub, work = _harness(tmp_path)
    (stub / "src" / "a.cenSat.bed").write_text("track name=x\nc1\t0\t10\tHSat2\t0\t.\n")
    (stub / "src" / "a.cenSat.bed.cut").write_text("")
    (stub / "src" / "err.cenSat.bed").write_text("<?xml version='1.0'?><Error><Code>NoSuchKey</Code></Error>\n")
    script = f"""source {PIPE}/config.sh
download https://h/a.cenSat.bed {work}/a.cenSat.bed bed4 && echo CUT-ACCEPTED
download https://h/err.cenSat.bed {work}/err.cenSat.bed bed4 && echo ERROR-PAGE-ACCEPTED
rm {stub}/src/a.cenSat.bed.cut; download https://h/a.cenSat.bed {work}/a.cenSat.bed bed4 && echo WHOLE"""
    r = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True, timeout=60)
    assert r.stdout.split() == ["WHOLE"], r.stdout + r.stderr
    assert not (work / "err.cenSat.bed").exists() and (work / "a.cenSat.bed").read_text().endswith("HSat2\t0\t.\n")


@needs_bash
def test_no_job_is_taken_for_ended_while_squeue_does_not_answer(tmp_path):
    env, stub, work = _harness(tmp_path)
    md5 = _cram(work / "crams" / "S1.cram")
    (work / "dispatched" / "S1").write_text("111\n")
    (work / "manifest.tsv").write_text(f"S1\thttps://x/S1.cram\t{md5}\nS2\thttps://x/S2.cram\t{md5}\n")
    (stub / "sqfail").write_text("")
    r = subprocess.run(["bash", str(PIPE / "01a_dispatch_staged.sh")], env=env, capture_output=True, text=True)
    assert r.returncode == 0 and "dispatched 0 sample(s); 1 with a job left for the next run" in r.stdout, r.stdout + r.stderr
    assert (work / "dispatched" / "S1").read_text().strip() == "111"
    try:           # the disk is full and squeue cannot say whether S1's job is alive: the stager waits, it does not stop
        r = subprocess.run(["bash", str(PIPE / "01_stage_and_dose.sh")], env=dict(env, MAX_LOCAL_CRAMS="1"), capture_output=True, text=True, timeout=5)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode() + (e.stderr or b"").decode()
    assert "waiting for room" in out and "STOPPED" not in out and "staged, job" not in out, out


@needs_bash
def test_counts_are_kept_without_a_cram_to_redo_them(tmp_path):
    env, stub, work = _harness(tmp_path)
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    (bundle / "sinks.bed").write_text("chr21\t0\t1000\tDJ\nchrUn\t0\t900\tDJ\n")
    env = dict(env, BUNDLE=str(bundle), EXPECTED_ENGINE_BUILD="zzz")
    regions = [{"name": "chr1:1-100", "role": "control", "obs": 7}]
    scan, fetch = work / "counts_scan" / "S1.json.gz", work / "counts_fetch" / "S1.json.gz"
    scan.parent.mkdir()
    fetch.parent.mkdir()
    _gz(scan, json.dumps(_scan(regions)))
    _gz(fetch, json.dumps(_fetch(regions, bundle / "sinks.bed", dj=100)))
    run = lambda e: subprocess.run(["bash", str(PIPE / "01b_dose_sample.sh"), "S1", "1"], env=e, capture_output=True, text=True)  # noqa: E731
    r = run(env)                                                  # another build, and no CRAM: left in place
    assert r.returncode == 1 and "left in place" in r.stderr and scan.exists() and fetch.exists(), r.stderr
    r = run(dict(env, EXPECTED_ENGINE_BUILD="fae1124"))
    assert r.returncode == 0 and "already done" in r.stdout, r.stderr
    _gz(fetch, json.dumps(_fetch(regions, bundle / "sinks.bed", dj=60)))
    r = run(dict(env, EXPECTED_ENGINE_BUILD="fae1124"))          # a fetch that lost reads, and no CRAM: moved aside
    assert r.returncode == 1 and not fetch.exists() and Path(f"{fetch}.rejected").exists() and scan.exists(), r.stderr


@needs_bash
def test_the_default_fetch_panels_leave_out_a_panel_without_sinks(tmp_path):
    env, stub, work = _harness(tmp_path)
    src = tmp_path / "NGS-DOSE"
    (src / "resources" / "GRCh38").mkdir(parents=True)
    (src / "resources" / "experimental").mkdir()
    _gz(src / "resources" / "experimental" / "telomere.k31.panel.tsv.gz", "##k=31\n##class\tid=0\tname=TEL\tkind=compositional\n#kmer\n")
    env.pop("FETCH_PANELS")
    env["NGSDOSE_SRC"] = str(src)
    script = f"source {PIPE}/config.sh; check_fetch_panels; echo \"rc=$? [$FETCH_PANELS]\""
    for sinks, default, explicit in (("chr21\t0\t10\tDJ\n", "rc=0 []", "rc=1"), ("chr21\t0\t10\tDJ\nchr1\t0\t10\tTEL\n", "rc=0 [/", "rc=0 [/")):
        (src / "resources" / "GRCh38" / "sinks.bed").write_text(sinks)
        r = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
        assert r.stdout.startswith(default), r.stdout + r.stderr
        r = subprocess.run(["bash", "-c", script], env=dict(env, FETCH_PANELS=str(src / "resources" / "experimental" / "telomere.k31.panel.tsv.gz")),
                           capture_output=True, text=True)
        assert r.stdout.startswith(explicit), r.stdout + r.stderr
