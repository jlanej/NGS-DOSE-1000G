"""Every shell script in the pipeline and the pilot parses; and the pipeline's guards hold, run with stub
sbatch / squeue / curl / sleep on PATH: a refused sbatch leaves no dispatch mark, the stager ends (rather than
waits for ever) when a failed job has kept the only CRAM the disk has room for, but neither it nor 01a takes a
job for ended while squeue does not answer, a cut download never takes its final name, a fetch that does not
match its scan is caught, counts are not removed when there is no CRAM to make them again from, and the default
FETCH_PANELS leaves out a panel the bundle has no sinks for. With a stub engine that records its command lines:
the default scan and fetch are exactly those of the fae1124 run; a fetch plan (FETCH_PRESET ...) is built once,
reused, refused when its settings or files change, and replaces the bundle's sinks, panels, flags and controls (a
lighter controls file is compared with the scan on the regions it holds; --classes needs an engine that takes it; a
controls file the bundle does not name is refused); candidate panels join the scans only, and a candidate that shares
a k-mer or a class name is stopped; 06_learn_sinks.sh refuses fewer than 30 scans, writes the named classes' sinks,
statistics and menu, and marks the held-out statistics as such; estimate is told a fetch plan's sinks."""
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


def test_check_counts_with_a_lighter_controls_file(tmp_path):
    """A fetch plan's controls (controls.lite200): the fetch must have used them, and its regions are compared with the
    scan's counts of the same regions, ctrl_reads with the scan's reads in the fetch's control regions."""
    full = [{"name": "c1", "role": "control", "obs": 7}, {"name": "c2", "role": "control", "obs": 5}, {"name": "t1", "role": "test", "obs": 90}]
    full_fa, lite_fa, dj = tmp_path / "controls.fa.gz", tmp_path / "controls.lite200.fa.gz", tmp_path / "dj.bed"
    _gz(full_fa, ">c1\n>c2\n>t1 role=test\n")
    _gz(lite_fa, ">c1\n>t1 role=test\n")
    dj.write_text("chr21\t0\t1000\tDJ\nchrUn\t0\t900\tDJ\n")
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
    scan, fetch = tmp_path / "s.json.gz", tmp_path / "f.json.gz"
    _gz(scan, json.dumps(dict(_scan(full), ctrl_reads=12, controls="controls.fa.gz", controls_sha256=sha(full_fa))))

    def fetched(regions, ctrl_reads=7, fa=lite_fa, dj_reads=100):
        _gz(fetch, json.dumps(dict(_fetch(regions, dj, dj=dj_reads), ctrl_reads=ctrl_reads, controls=fa.name, controls_sha256=sha(fa))))

    same = lambda *extra: check_counts.main(["same-reads", str(scan), str(fetch), str(dj), *map(str, extra)])  # noqa: E731
    fetched([full[2], full[0]])
    assert same(lite_fa) == 0
    assert same() == 1                                                   # without the plan's controls: all regions compared
    assert same(full_fa) == 1                                            # not the controls the plan says
    fetched([full[0], dict(full[2], obs=89)])                             # a region that lost a read
    assert same(lite_fa) == 1
    fetched([full[0], full[2]], ctrl_reads=8)
    assert same(lite_fa) == 1
    fetched([full[0], full[2]], dj_reads=60)                              # class reads are still checked
    assert same(lite_fa) == 1
    fetched([full[0], full[2], {"name": "c9", "role": "control", "obs": 1}])   # not a subset of the scan's regions
    assert same(lite_fa) == 2
    fetched([full[0], dict(full[2], role="control")])
    assert same(lite_fa) == 2
    fetched([full[2]], ctrl_reads=0)                                      # no control region at all
    assert same(lite_fa) == 2
    fetched(full, ctrl_reads=12, fa=full_fa)                              # a plan with the scan's own controls: compared whole
    assert same(full_fa) == 0
    fetched([full[0], full[2]], fa=full_fa)
    assert same(full_fa) == 1


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


# --- the fetch menu and the candidates -------------------------------------------------------------------------

# The stub engine: `count --help` names --classes unless $STUB/noclasses exists; a count records its command line and
# writes a counts file whose regions are those of its controls FASTA (role from the header, 5 reads each)
ENGINE = r"""#!/usr/bin/env bash
[ "${2:-}" = --help ] && { [ -e "$STUB/noclasses" ] || echo "      --classes <CLASSES>"; exit 0; }
printf '%s\t' "$@" >> "$STUB/engine.log"; echo >> "$STUB/engine.log"
out=""; mode=""; ctrl=""; prev=""
for a in "$@"; do [ "$prev" = -o ] && out=$a; [ "$prev" = -m ] && mode=$a; [ "$prev" = -c ] && ctrl=$a; prev=$a; done
sha=$( (sha256sum "$ctrl" 2>/dev/null || shasum -a 256 "$ctrl") | cut -d' ' -f1)
IFS=$'\t' read -r regions ctrl_reads < <(gzip -dc "$ctrl" | awk '/^>/ {role = "control"; for (i = 2; i <= NF; i++) if ($i ~ /^role=/) role = substr($i, 6)
  printf "%s{\"name\": \"%s\", \"role\": \"%s\", \"obs\": 5}", (n++ ? ", " : ""), substr($1, 2), role; if (role == "control") c += 5}
  END {printf "\t%d\n", c}')
printf '{"engine_build": "fae112470d74", "mode": "%s", "controls": "%s", "controls_sha256": "%s", "ctrl_reads": %s, "regions": [%s], "classes": [], "placements": [], "sinks_sha256": "x"}' \
  "$mode" "$(basename "$ctrl")" "$sha" "$ctrl_reads" "$regions" | gzip -c > "$out"
"""
# ngsdose, as far as the scripts use it: `fetchplan` (absent while $STUB/nofetchplan exists) and `sinks`
NGSDOSE = r"""#!/usr/bin/env bash
echo "$*" >> "$STUB/ngsdose.log"
cmd=$1; shift
if [ "$cmd" = fetchplan ]; then
  [ -e "$STUB/nofetchplan" ] && { echo "ngsdose: error: argument cmd: invalid choice: 'fetchplan'" >&2; exit 2; }
  [ "$1" = --help ] && { echo "  --engine ENGINE"; exit 0; }
  prev=""; ctrl=""; for a in "$@"; do [ "$prev" = -o ] && out=$a; [ "$prev" = --controls ] && ctrl=$a; [ "$prev" = --menu ] && menu=$a; prev=$a; done
  [ -n "$ctrl" ] || ctrl="$(dirname "$menu")/GRCh38/controls.bed"
  printf '# plan\nchr1\t0\t21000\tTEL\n' > "$out.sinks.bed"; cat "$STUB/plan_panels" > "$out.panels.txt"; echo "${ctrl%.bed}.fa.gz" > "$out.controls.txt"
  if [ -e "$STUB/plan_flags" ]; then cat "$STUB/plan_flags" > "$out.count_flags.txt"; else printf -- '--unmapped\n' > "$out.count_flags.txt"; fi
  : > "$out.scan_panels.txt"; echo "option" > "$out.plan.tsv"; echo "option"
  exit 0
fi
if [ "$1" = --help ]; then
  case "$cmd" in sinks) echo "--stats FILE"; [ -e "$STUB/noheldout" ] || echo "--held-out" ;; estimate) echo "--fetch-sinks BED" ;; esac
  exit 0
fi
files=(); ev=""; o=""; st=""; prev=""
for a in "$@"; do
  case "$prev" in --evaluate) ev=$a ;; -o) o=$a ;; --stats) st=$a ;; --classes) ;; *) case "$a" in *.json.gz) files+=("$a") ;; esac ;; esac
  prev=$a
done
H="class\trank\tcontig\tstart\tend\tbp\tscans\tshare_median\tshare_p10\treads_median\tshare_per_read\tcum_capture_median\tcum_capture_p10"
if [ -n "$ev" ]; then
  printf "$H\nNEWC\t1\tchr1\t0\t20000\t20000\t9\t0.99\t0.98\t10\t0.1\t0.99\t0.98\n" > "$st"
  printf 'sample\tclass\tscan_reads\tcaptured\tfraction\tunmapped\n'
  for f in "${files[@]}"; do printf '%s\tNEWC\t100\t99\t0.99000\t1\n%s\tDJ\t100\t100\t1.00000\t0\n' "$(basename "$f")" "$(basename "$f")"; done
else
  printf 'chr1\t0\t20000\tNEWC\nchr21\t0\t1000\tDJ\n' > "$o"
  printf "# in-sample\n$H\nDJ\t1\tchr21\t0\t1000\t1000\t9\t1\t1\t10\t0.1\t1\t1\nNEWC\t1\tchr1\t0\t20000\t20000\t9\t0.99\t0.98\t10\t0.1\t0.99\t0.98\n" > "$st"
fi
"""
TIMEOUT = "#!/usr/bin/env bash\nshift; exec \"$@\"\n"
MENU = ("# a menu\n##preset\tcore_tel\tthe test\n"
        "class\tgroup\tkind\tstatus\ttier\tpresets\tpanel\tsinks\tmeasures\ttruth\tnotes\tstats\n"
        "controls\tcontrols\tcontrols\tshipped\tA\t-\t-\tGRCh38/controls.bed\tcontrols\t-\t-\t-\n"
        "DJ\tacro\tpositional\tshipped\tA\tcore_tel\tGRCh38/panel.k31.tsv.gz\tGRCh38/sinks.bed\tDJ\t-\t-\t-\n"
        "TEL\ttelomere\tcompositional\tshipped\tB\tcore_tel\texperimental/telomere.k31.panel.tsv.gz\tGRCh38/sinks.bed\tTEL\t-\t-\t-\n")


def _panel(path, classes, kmers):
    _gz(path, "##ngs-dose-panel v1\n##k=31\n" + "".join(f"##class\tid={i}\tname={n}\tkind=positional\tlength=100\tcircular=0\n"
                                                         for i, n in enumerate(classes))
        + "#kmer\tclass\tpos\tstrand\n" + "".join(f"{k}\t0\t0\t+\n" for k in kmers))


def _image(tmp_path, env, stub):
    """A checkout standing in for the image (SIF unset): the bundle without TEL sinks, as in fae1124, the two
    experimental panels, a fetch menu; and a stub engine, ngsdose and timeout."""
    src = tmp_path / "NGS-DOSE"
    b, x = src / "resources" / "GRCh38", src / "resources" / "experimental"
    b.mkdir(parents=True)
    x.mkdir()
    _panel(b / "panel.k31.tsv.gz", ["DJ"], ["A" * 31])
    _panel(x / "satellites.CHM13v2.k31.panel.tsv.gz", ["HSat2"], ["C" * 30 + "A"])
    _panel(x / "telomere.k31.panel.tsv.gz", ["TEL"], ["AACCCT" * 5 + "A"])
    (b / "sinks.bed").write_text("chr21\t0\t1000\tDJ\n")
    (b / "controls.bed").write_text("chr1\t100\t200\tc1\n")
    _gz(b / "controls.fa.gz", ">c1\nACGT\n")
    (src / "resources" / "fetch_menu.tsv").write_text(MENU)
    stubs = Path(env["PATH"].split(":")[0])
    for name, text in (("ngs-dose", ENGINE), ("ngsdose", NGSDOSE), ("timeout", TIMEOUT)):
        (stubs / name).write_text(text)
        (stubs / name).chmod(0o755)
    (stub / "plan_panels").write_text(f"{x / 'telomere.k31.panel.tsv.gz'}\n")
    env = dict(env, NGSDOSE_SRC=str(src), NGSDOSE_BIN=str(stubs / "ngs-dose"), KEEP_CRAMS="1", THREADS="8")
    for k in ("FETCH_PANELS", "EXTRA_PANELS", "BUNDLE", "SLURM_CPUS_PER_TASK", "SLURM_SUBMIT_DIR", "CANDIDATE_PANELS", "FETCH_CLASSES",
              "FETCH_PRESET", "FETCH_BUDGET_MB", "FETCH_CAPTURE", "FETCH_MENU", "FETCH_PLAN_ARGS", "FETCH_CRAI", "FETCH_PLAN_DIR"):
        env.pop(k, None)
    return env, src


def _calls(stub):
    log = stub / "engine.log"
    calls = [line.rstrip("\t").split("\t") for line in log.read_text().splitlines()] if log.exists() else []
    log.unlink(missing_ok=True)
    return calls


def _dose(env, work, sample):
    _cram(work / "crams" / f"{sample}.cram")
    return subprocess.run(["bash", str(PIPE / "01b_dose_sample.sh"), sample, "1"], env=env, capture_output=True, text=True, timeout=120)


@needs_bash
def test_default_scan_and_fetch_are_those_of_the_cohort_run(tmp_path):
    env, stub, work = _harness(tmp_path)
    env, src = _image(tmp_path, env, stub)
    b, x = src / "resources" / "GRCh38", src / "resources" / "experimental"
    r = _dose(env, work, "S1")
    assert r.returncode == 0, r.stdout + r.stderr
    cram, ref = work / "crams" / "S1.cram", work / "reference" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    common = ["count", "-i", str(cram), "--index", f"{cram}.crai", "-T", str(ref), "-c", f"{b}/controls.fa.gz", "-p", f"{b}/panel.k31.tsv.gz",
              "-s", "S1", "-@", "8"]
    assert _calls(stub) == [common + ["-m", "scan", "-p", str(x / "satellites.CHM13v2.k31.panel.tsv.gz"), "-p", str(x / "telomere.k31.panel.tsv.gz"),
                                      "-o", str(work / "counts_scan" / "S1.json.gz")],
                            common + ["-m", "fetch", "--sinks", f"{b}/sinks.bed", "-o", str(work / "counts_fetch" / "S1.json.gz")]]
    assert not (work / "fetchplan").exists() and not (stub / "ngsdose.log").exists()
    # 01_count.sh, fetch mode, from a local CRAM
    (work / "manifest.tsv").write_text(f"S9\t{cram}\n")
    r = subprocess.run(["bash", str(PIPE / "01_count.sh")], env=dict(env, MODE="fetch"), capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    out = work / "counts_fetch" / "S9.json.gz.tmp.gz"
    assert _calls(stub) == [["count", "-i", str(cram), "-T", str(ref), "-c", f"{b}/controls.fa.gz", "-p", f"{b}/panel.k31.tsv.gz", "-@", "8", "-s", "S9",
                             "-o", str(out), "-m", "fetch", "--sinks", f"{b}/sinks.bed", "--index", f"{cram}.crai"]]


@needs_bash
def test_a_fetch_plan_is_built_once_and_kept(tmp_path):
    env, stub, work = _harness(tmp_path)
    env, src = _image(tmp_path, env, stub)
    b, x = src / "resources" / "GRCh38", src / "resources" / "experimental"
    plan = work / "fetchplan"
    env = dict(env, FETCH_PRESET="core_tel")
    r = _dose(env, work, "S1")
    assert r.returncode == 0, r.stdout + r.stderr
    scan, fetch = _calls(stub)
    assert "-m" in scan and scan[scan.index("-m") + 1] == "scan" and f"{b}/panel.k31.tsv.gz" in scan         # the scan is the scan
    assert fetch[fetch.index("--sinks") + 1] == str(plan / "plan.sinks.bed") and f"{b}/panel.k31.tsv.gz" not in fetch
    assert fetch[fetch.index("-p") + 1] == str(x / "telomere.k31.panel.tsv.gz") and "--unmapped" in fetch
    assert "fetchplan --menu" in (stub / "ngsdose.log").read_text() and (plan / "settings.txt").exists()
    r = _dose(env, work, "S2")                                            # built once, then reused
    assert r.returncode == 0, r.stdout + r.stderr
    assert (stub / "ngsdose.log").read_text().count(" -o ") == 1 and len(_calls(stub)) == 2
    (work / "manifest.tsv").write_text(f"S3\t{work / 'crams' / 'S1.cram'}\n")
    r = subprocess.run(["bash", str(PIPE / "01_count.sh")], env=dict(env, MODE="fetch"), capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    (call,) = _calls(stub)
    assert call[call.index("--sinks") + 1] == str(plan / "plan.sinks.bed") and "--unmapped" in call and f"{b}/panel.k31.tsv.gz" not in call
    for change in ({"FETCH_CAPTURE": "0.99"}, {"FETCH_PRESET": "core_tel truths"}):
        r = _dose(dict(env, **change), work, "S4")
        assert r.returncode == 1 and "built from other settings" in r.stderr and not _calls(stub), r.stderr
    (b / "sinks.bed").write_text("chr21\t0\t1001\tDJ\n")                  # a file the menu names has changed
    r = _dose(env, work, "S4")
    assert r.returncode == 1 and "built from other settings" in r.stderr and "sinks.bed" in r.stderr, r.stderr
    r = _dose(dict(env, FETCH_PRESET=""), work, "S4")                    # a job that did not inherit the settings: said so
    assert r.returncode == 0 and "holds a fetch plan" in r.stderr, r.stderr
    assert _calls(stub)[1][-3:-2] == [str(b / "sinks.bed")]
    shutil.rmtree(plan)
    for bad, msg in (({"FETCH_PRESET": "", "FETCH_CAPTURE": "0.995"}, "choose its options"),
                     ({"FETCH_PANELS": str(x / "telomere.k31.panel.tsv.gz")}, "unset FETCH_PANELS")):
        r = _dose(dict(env, **bad), work, "S5")
        assert r.returncode == 1 and msg in r.stderr and not _calls(stub), r.stderr
    (stub / "nofetchplan").write_text("")                                 # an image whose ngsdose predates fetchplan
    r = _dose(env, work, "S5")
    assert r.returncode == 1 and "does not have" in r.stderr and not _calls(stub) and not plan.exists(), r.stderr


@needs_bash
def test_a_fetch_plan_uses_its_controls_and_count_flags(tmp_path):
    env, stub, work = _harness(tmp_path)
    env, src = _image(tmp_path, env, stub)
    b = src / "resources" / "GRCh38"
    _gz(b / "controls.fa.gz", ">c1\nACGT\n>c2\nACGT\n>t1 flank=1000 role=test\nACGT\n")
    (b / "controls.lite200.bed").write_text("chr1\t100\t200\tc1\n")
    _gz(b / "controls.lite200.fa.gz", ">c1\nACGT\n>t1 flank=1000 role=test\nACGT\n")
    (stub / "plan_flags").write_text("--classes=DJ,TEL\n")
    plan, full, lite = work / "fetchplan", str(b / "controls.fa.gz"), str(b / "controls.lite200.fa.gz")
    env = dict(env, FETCH_PRESET="core_tel", FETCH_PLAN_ARGS=f"--controls {b}/controls.lite200.bed")
    r = _dose(env, work, "S1")
    assert r.returncode == 0, r.stdout + r.stderr
    scan, fetch = _calls(stub)
    assert scan[scan.index("-c") + 1] == full and fetch[fetch.index("-c") + 1] == lite and fetch[-3] == "--classes=DJ,TEL"
    assert "the fetch counted 2 of the scan's 3 regions" in r.stderr and "controls controls.lite200.fa.gz" in r.stderr
    assert f"--engine {env['NGSDOSE_BIN']}" in (stub / "ngsdose.log").read_text()      # the engine the fetch runs
    assert (plan / "plan.controls.txt").read_text().strip() == lite
    (work / "manifest.tsv").write_text(f"S9\t{work / 'crams' / 'S1.cram'}\n")               # 01_count.sh fetches the same way
    r = subprocess.run(["bash", str(PIPE / "01_count.sh")], env=dict(env, MODE="fetch"), capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    (call,) = _calls(stub)
    assert call[call.index("-c") + 1] == lite and "--classes=DJ,TEL" in call and call[call.index("--sinks") + 1] == str(plan / "plan.sinks.bed")
    # a fetch made with the full controls while the plan says lite200: removed and made again, from the CRAM that is here
    for d in ("counts_scan", "counts_fetch"):
        shutil.copy(work / "counts_scan" / "S1.json.gz", work / d / "S2.json.gz")
    r = _dose(env, work, "S2")
    assert r.returncode == 0 and "not with" in r.stderr and "removed" in r.stderr, r.stdout + r.stderr
    (call,) = _calls(stub)
    assert call[call.index("-m") + 1] == "fetch" and call[call.index("-c") + 1] == lite
    # a lighter file that is not a subset of the scan's regions: the two cannot be compared, and nothing is removed
    _gz(b / "controls.lite200.fa.gz", ">c1\nACGT\n>c9\nACGT\n")
    r = _dose(env, work, "S3")
    assert r.returncode == 1 and "could not compare" in r.stderr and "c9" in r.stderr, r.stdout + r.stderr
    assert (work / "counts_fetch" / "S3.json.gz").exists() and (work / "crams" / "S3.cram").exists()
    _calls(stub)
    (stub / "noclasses").write_text("")                                   # an engine without --classes (fae1124)
    r = _dose(env, work, "S4")
    assert r.returncode == 1 and "does not take" in r.stderr and not _calls(stub), r.stderr
    (stub / "noclasses").unlink()
    shutil.rmtree(plan)
    other = tmp_path / "mine"
    other.mkdir()
    (other / "controls.x.bed").write_text("chr1\t100\t200\tc1\n")
    _gz(other / "controls.x.fa.gz", ">c1\nACGT\n")
    r = _dose(dict(env, FETCH_PLAN_ARGS=f"--controls {other}/controls.x.bed"), work, "S5")    # not a subset the bundle names
    assert r.returncode == 1 and "would refuse every counts file" in r.stderr and not _calls(stub), r.stderr
    shutil.rmtree(plan)
    (b / "controls.lite200.fa.gz").unlink()                               # the FASTA of the plan's controls is not there
    r = _dose(env, work, "S5")
    assert r.returncode == 1 and "not there" in r.stderr and not _calls(stub), r.stderr


@needs_bash
def test_estimate_is_told_the_fetch_plans_sinks(tmp_path):
    env, stub, work = _harness(tmp_path)
    env, src = _image(tmp_path, env, stub)
    script = f"source {PIPE}/config.sh; estimate_sinks_args; echo \"rc=$? [${{ESTIMATE_SINKS_ARGS[*]}}]\""
    run = lambda **e: subprocess.run(["bash", "-c", script], env=dict(env, **e), capture_output=True, text=True).stdout.strip()  # noqa: E731
    assert run() == "rc=0 []" and not (stub / "ngsdose.log").exists()   # no plan: nothing, and nothing asked of the image
    (work / "fetchplan").mkdir()
    (work / "fetchplan" / "plan.sinks.bed").write_text("chr1\t0\t10\tTEL\n")
    assert run() == f"rc=0 [--fetch-sinks {work}/fetchplan/plan.sinks.bed]"
    (tmp_path / "old.bed").write_text("chr1\t0\t10\tTEL\n")
    assert run(FETCH_SINKS_KNOWN=f"{tmp_path}/old.bed {work}/fetchplan/plan.sinks.bed") == \
        f"rc=0 [--fetch-sinks {tmp_path}/old.bed {work}/fetchplan/plan.sinks.bed]"
    assert run(FETCH_SINKS_KNOWN=f"{tmp_path}/none.bed").startswith("rc=1 []")


@needs_bash
def test_candidate_panels_are_scanned_only_and_checked(tmp_path):
    env, stub, work = _harness(tmp_path)
    env, src = _image(tmp_path, env, stub)
    cand = tmp_path / "candidates"
    cand.mkdir()
    _panel(cand / "new.panel.tsv.gz", ["NEWC"], ["G" * 31, "T" * 30 + "G"])
    (cand / "notes.tsv").write_text("not a panel\n")
    env = dict(env, CANDIDATE_PANELS=str(cand))
    script = f"source {PIPE}/config.sh; echo \"[$EXTRA_PANELS] [$APPTAINER_BINDS]\"; check_candidate_panels; echo rc=$?"
    r = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert r.stdout.split()[-1] == "rc=0" and f"{cand}/new.panel.tsv.gz]" in r.stdout and f",{cand}]" in r.stdout, r.stdout + r.stderr
    assert "NEWC" in r.stderr and (work / "candidate_panels.checked").exists()
    r = _dose(env, work, "S1")
    assert r.returncode == 0, r.stdout + r.stderr
    scan, fetch = _calls(stub)
    assert scan[-3:-2] == [str(cand / "new.panel.tsv.gz")] and str(cand / "new.panel.tsv.gz") not in fetch
    for classes, kmers in ((["NEWC"], ["AACCCT" * 5 + "A"]), (["TEL"], ["G" * 31])):             # a shared k-mer; a class name
        _panel(cand / "new.panel.tsv.gz", classes, kmers)
        r = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
        assert r.stdout.split()[-1] == "rc=1" and "telomere.k31.panel.tsv.gz" in r.stderr, r.stdout + r.stderr
        r = _dose(env, work, "S2")
        assert r.returncode == 1 and "clash" in r.stderr and not _calls(stub), r.stderr
    r = subprocess.run(["bash", "-c", script], env=dict(env, CANDIDATE_PANELS=str(tmp_path / "nowhere")), capture_output=True, text=True)
    assert r.stdout.split()[-1] == "rc=1" and "neither a panel file nor a directory" in r.stderr


def _scans(d, n, classes=("DJ", "NEWC")):
    d.mkdir(parents=True, exist_ok=True)
    for i in range(n):
        _gz(d / f"S{i:03d}.json.gz", json.dumps({"sample": f"S{i:03d}", "mode": "scan",
                                                  "classes": [{"name": c, "kind": "positional", "reads": 100} for c in classes]}))


@needs_bash
def test_learn_sinks_needs_30_scans_and_writes_the_menu(tmp_path):
    env, stub, work = _harness(tmp_path)
    env, src = _image(tmp_path, env, stub)
    cand = tmp_path / "candidates"
    cand.mkdir()
    _panel(cand / "new.panel.tsv.gz", ["NEWC"], ["G" * 31])
    env = dict(env, CANDIDATE_PANELS=str(cand))
    _scans(work / "counts_scan", 29)
    _scans(work / "counts_scan", 0)
    for i in range(29, 40):                                               # scans made before the candidate: not counted
        _gz(work / "counts_scan" / f"S{i:03d}.json.gz", json.dumps({"sample": f"S{i:03d}", "mode": "scan", "classes": [{"name": "DJ", "kind": "positional"}]}))
    learn = lambda: subprocess.run(["bash", str(PIPE / "06_learn_sinks.sh"), "cand1", "NEWC"], env=env, capture_output=True, text=True, timeout=120)  # noqa: E731
    r = learn()
    assert r.returncode == 1 and "29 of 40 scans" in r.stderr and "at least 30" in r.stderr and not (work / "sinks" / "cand1.bed").exists(), r.stderr
    _gz(work / "counts_scan" / "S039.json.gz", json.dumps({"sample": "S039", "mode": "scan", "classes": [{"name": "NEWC", "kind": "positional"}]}))
    r = subprocess.run(["bash", str(PIPE / "06_learn_sinks.sh"), "cand1", "NEWC"], env=dict(env, CANDIDATE_PANELS=""), capture_output=True, text=True)
    assert r.returncode == 0 and "is not written" in r.stderr and (work / "sinks" / "cand1.bed").exists(), r.stderr    # no panel, no menu
    assert not (work / "sinks" / "cand1.menu.tsv").exists()
    r = learn()
    assert r.returncode == 0, r.stdout + r.stderr
    out = work / "sinks"
    assert (out / "cand1.bed").read_text() == "chr1\t0\t20000\tNEWC\n"                     # DJ, which the learner adds, left out
    assert "DJ\t" not in (out / "cand1.train.stats.tsv").read_text() and "NEWC\t1" in (out / "cand1.stats.tsv").read_text()
    roles = [line.split("\t")[0] for line in (out / "cand1.scans.txt").read_text().splitlines() if not line.startswith("#")]
    assert roles.count("train") == 15 and roles.count("heldout") == 15
    assert "NEWC\tpositional\t1\t20000\t15\t0.99000\t0.99000\t0.99000\t0.01000" in r.stdout
    log = (stub / "ngsdose.log").read_text().splitlines()
    assert any("--evaluate" in line and "--held-out" in line for line in log)                  # held-out statistics say so
    assert not any("--evaluate" not in line and "--held-out" in line for line in log)
    menu = (out / "cand1.menu.tsv").read_text()
    assert "##preset\tcand1\t" in menu and f"{src}/resources/GRCh38/sinks.bed" in menu            # the image's paths made absolute
    row = next(line.split("\t") for line in menu.splitlines() if line.startswith("NEWC\t"))
    assert row[3].startswith("experimental:") and row[5] == "cand1" and row[6] == str(out / "cand1.new.panel.tsv.gz") and row[7] == str(out / "cand1.bed")
    assert (out / "cand1.new.panel.tsv.gz").read_bytes() == (cand / "new.panel.tsv.gz").read_bytes()
    try:
        from ngsdose import fetchplan
    except ImportError:
        return
    if hasattr(fetchplan, "read_menu"):                                   # the menu is one that `ngsdose fetchplan` reads
        m = fetchplan.read_menu(out / "cand1.menu.tsv")
        assert m.members("cand1") == ["NEWC"] and m.options["NEWC"].status == "experimental"
