#!/usr/bin/env bash
# Shared configuration for the 1000 Genomes 30x run. Override any variable by exporting it.
EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"      # this directory
REPO="$(cd "$EX_DIR/.." && pwd)"                               # this repository (meta/, pilot/, report/)
NGSDOSE_SRC="${NGSDOSE_SRC:-$REPO/../NGS-DOSE}"                # a checkout of the method, for runs without the container

WORK_DIR="${WORK_DIR:-/scratch/${USER}/ngs_dose_1000G}"
LOG_DIR="${LOG_DIR:-$WORK_DIR/logs}"
MANIFEST="${MANIFEST:-$WORK_DIR/manifest.tsv}"      # SAMPLE <TAB> CRAM (https:// URL or local path)

# Two ways to run. With SIF set, everything - engine, python package, resource bundle, aria2c -
# comes from the container image and nothing is installed on the host; these scripts come with this
# repository. Without it: a checkout of NGS-DOSE at $NGSDOSE_SRC with `cargo build --release` and
# `pip install .`.
SIF="${SIF:-}"
# The image 00_setup.sh pulls into $SIF when that file does not exist yet. :latest follows every push to
# NGS-DOSE's main; a cohort should keep one engine, so pin it: :sha-<first 7 characters of the commit>
# (published for every push to main), :X.Y.Z, or @sha256:<digest>. The cohort's counts so far were made
# with docker://ghcr.io/jlanej/ngs-dose:sha-fae1124.
NGSDOSE_IMAGE="${NGSDOSE_IMAGE:-docker://ghcr.io/jlanej/ngs-dose:latest}"
# The engine commit (or a prefix of it) every new counts file must record as engine_build, e.g. fae1124.
# A file from another build is removed and its job fails, with the CRAM kept. Empty: not checked.
EXPECTED_ENGINE_BUILD="${EXPECTED_ENGINE_BUILD:-}"
NGSDOSE_BIN="${NGSDOSE_BIN:-$NGSDOSE_SRC/target/release/ngs-dose}"
if [ -n "$SIF" ]; then IMAGE_ROOT=/opt/ngs-dose; else IMAGE_ROOT="$NGSDOSE_SRC"; fi
BUNDLE="${BUNDLE:-$IMAGE_ROOT/resources/GRCh38}"
# Further panels for whole-file scans (space-separated): the experimental satellite families and the
# telomeric repeat (NGS-DOSE's resources/experimental/README.md). A scan counts their reads wherever they
# were aligned and records where that was, which is what sinks are learned from. EXTRA_PANELS="" scans
# with the bundle's classes only.
EXTRA_PANELS="${EXTRA_PANELS-$IMAGE_ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz $IMAGE_ROOT/resources/experimental/telomere.k31.panel.tsv.gz}"
# Panels loaded in fetch mode as well. A fetch reads only the controls and the bundle's sinks, so every
# class of these panels needs intervals in $BUNDLE/sinks.bed. The bundle has sinks for the telomeric repeat
# (92% of its reads lie within 25 kb of a chromosome end), added in NGS-DOSE 45c9f9a (2026-09-22). It has
# none yet for the satellite families. In these CRAMs sinks learned from 30 scans hold >= 99.8% of each
# family's reads in every one of 200 other genomes (99.85% in two of three random draws; HSat1B >= 96.6%), but no such sinks ship, so here the
# satellites are measured by the scan alone.
# Left unset, FETCH_PANELS is the telomere panel when the image's bundle has its sinks, and nothing (with
# a note) when it has not, as in the fae1124 image that counted the cohort's first 1,748 genomes: the
# default never loads a class a fetch would miss. Set explicitly (exported, or FETCH_PANELS="" for the
# bundle's classes alone), a panel whose classes lack sinks stops every script that fetches.
if [ -z "${FETCH_PANELS+set}" ]; then FETCH_PANELS_DEFAULT=1; else FETCH_PANELS_DEFAULT=0; fi
FETCH_PANELS="${FETCH_PANELS-$IMAGE_ROOT/resources/experimental/telomere.k31.panel.tsv.gz}"

# New classes are scanned first (pipeline/README.md, "New classes, and what a fetch reads"). CANDIDATE_PANELS:
# k-mer panels on the host, space-separated, each a file or a directory whose *panel.tsv.gz files are taken. They
# join EXTRA_PANELS, so every scan counts their reads wherever they were aligned, and their directories join
# APPTAINER_BINDS. Before anything is counted, check_candidate_panels (every script that counts runs it) makes sure
# that no candidate shares a k-mer or a class name with the panels loaded beside it, or with another candidate: the
# engine drops a k-mer that two loaded panels claim from both, which would change the counts of classes the cohort
# already measures. A candidate is never fetched: it has no sinks until 06_learn_sinks.sh learns them from the scans
# that carry it and checks them on held-out ones.
CANDIDATE_PANELS="${CANDIDATE_PANELS:-}"

# What a fetch reads. Unset, it is what it has always been: the controls and $BUNDLE/sinks.bed, with the bundle panel
# and FETCH_PANELS. Set any of FETCH_CLASSES (options by name), FETCH_PRESET (presets of the menu), FETCH_BUDGET_MB
# (options in tier order until the plan would read more than this many MB of a CRAM) - and FETCH_CAPTURE (per
# class, keep the sink intervals that yield most per read until the held-out capture reaches this fraction) - and the
# fetch follows a plan that `ngsdose fetchplan` builds once, from FETCH_MENU, into FETCH_PLAN_DIR: its sinks BED, its
# panels, its flags (--classes=... where the panels define classes the plan leaves out and the engine takes it) and its
# controls FASTA (FETCH_PLAN_ARGS="--controls $BUNDLE/controls.lite200.bed" reads 200 of the 800 control regions and
# all 182 truth and dosage regions; the scan-vs-fetch check before a CRAM is removed then compares the regions both
# counts files hold). Every later run with the same settings and files reuses that plan; one that would build
# another stops, as a cohort's fetches should share one plan. The image must have `ngsdose fetchplan` (the fae1124
# image has not). `ngsdose fetchplan --menu $FETCH_MENU --list` shows the options and presets.
FETCH_CLASSES="${FETCH_CLASSES:-}"
FETCH_PRESET="${FETCH_PRESET:-}"
FETCH_BUDGET_MB="${FETCH_BUDGET_MB:-}"
FETCH_CAPTURE="${FETCH_CAPTURE:-}"
FETCH_MENU_GIVEN="${FETCH_MENU:+1}"                # the default menu is the image's own and is never bound
FETCH_MENU="${FETCH_MENU:-$IMAGE_ROOT/resources/fetch_menu.tsv}"   # or one written by 06_learn_sinks.sh
FETCH_PLAN_ARGS="${FETCH_PLAN_ARGS:-}"             # further fetchplan flags, e.g. "--fill" or "--capture-class TEL=0.99"
FETCH_PLAN_DIR="${FETCH_PLAN_DIR:-$WORK_DIR/fetchplan}"
# The CRAM indexes the plan is costed on (a budget needs them): FETCH_CRAI, or the indexes of the manifest's first
# FETCH_COST_SAMPLES samples (downloaded to $WORK_DIR/crai, as 01_count.sh does), with the contigs of $REF_FASTA.fai.
FETCH_CRAI="${FETCH_CRAI:-}"
FETCH_COST_SAMPLES="${FETCH_COST_SAMPLES:-5}"
# staged CRAMs for the full-accuracy path (01b_dose_sample.sh): $CRAM_DIR/<sample>.cram(.crai)
CRAM_DIR="${CRAM_DIR:-$WORK_DIR/crams}"

# CRAM decoding needs the reference: a local FASTA (recommended) or htslib's REF_PATH/REF_CACHE
REF_FASTA="${REF_FASTA:-$WORK_DIR/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa}"

# what the container has to see (bind sources must exist: 00_setup.sh creates them)
APPTAINER_BINDS="${APPTAINER_BINDS:-$WORK_DIR,$REPO,$(dirname "$REF_FASTA"),$CRAM_DIR}"
bind_dir() {  # dir: add it to APPTAINER_BINDS unless a bound directory holds it already, or it is the image's own
  local b src
  # a host directory bound over the image's /opt/ngs-dose would hide the image's files (the bundle among them)
  if [ -n "$SIF" ]; then case "$1/" in "$IMAGE_ROOT"/*) return 0 ;; esac; fi
  for b in ${APPTAINER_BINDS//,/ }; do
    src="${b%%:*}"; src="${src%/}"
    [ -n "$src" ] || continue
    case "$1/" in "$src"/*) return 0 ;; esac
  done
  APPTAINER_BINDS="${APPTAINER_BINDS:+$APPTAINER_BINDS,}$1"
}
# the candidate panels, as absolute paths (the same inside the container, through the bind); a problem is kept for
# check_candidate_panels to report, since this file is sourced, not run
SHIPPED_SCAN_PANELS="$EXTRA_PANELS"; CANDIDATE_PANEL_FILES=""; CANDIDATE_PANELS_ERROR=""
for _c in $CANDIDATE_PANELS; do
  if [ -d "$_c" ]; then
    _n=0
    for _f in "$_c"/*panel.tsv.gz; do [ -s "$_f" ] && CANDIDATE_PANEL_FILES="$CANDIDATE_PANEL_FILES $(cd "$(dirname "$_f")" && pwd)/$(basename "$_f")" && _n=$((_n + 1)); done
    [ "$_n" -gt 0 ] || CANDIDATE_PANELS_ERROR="$CANDIDATE_PANELS_ERROR; directory $_c has no *panel.tsv.gz file"
  elif [ -s "$_c" ]; then
    CANDIDATE_PANEL_FILES="$CANDIDATE_PANEL_FILES $(cd "$(dirname "$_c")" && pwd)/$(basename "$_c")"
  else
    CANDIDATE_PANELS_ERROR="$CANDIDATE_PANELS_ERROR; $_c is neither a panel file nor a directory on this host"
  fi
done
for _f in $CANDIDATE_PANEL_FILES; do
  case "$_f" in *,*|*:*) CANDIDATE_PANELS_ERROR="$CANDIDATE_PANELS_ERROR; $_f: a path with a comma or a colon cannot be bound" ;; esac
  case " $SHIPPED_SCAN_PANELS " in *" $_f "*) CANDIDATE_PANELS_ERROR="$CANDIDATE_PANELS_ERROR; $_f is in EXTRA_PANELS already" ;; esac
  bind_dir "$(dirname "$_f")"
done
CANDIDATE_PANEL_FILES="${CANDIDATE_PANEL_FILES# }"; CANDIDATE_PANELS_ERROR="${CANDIDATE_PANELS_ERROR#; }"
[ -z "$CANDIDATE_PANEL_FILES" ] || EXTRA_PANELS="${EXTRA_PANELS:+$EXTRA_PANELS }$CANDIDATE_PANEL_FILES"
unset _c _f _n
# a menu of the user's own, on the host, is bound too (the panels and sinks it names must be where the container sees them)
if [ -n "$SIF" ] && [ -n "$FETCH_MENU_GIVEN" ] && [ -f "$FETCH_MENU" ]; then bind_dir "$(cd "$(dirname "$FETCH_MENU")" && pwd)"; fi
# what the fetch reads; check_fetch_panels replaces these with the plan's when a plan is asked for
FETCH_SINKS_BED="$BUNDLE/sinks.bed"; FETCH_CONTROLS="$BUNDLE/controls.fa.gz"; FETCH_FLAGS=""; FETCH_PLAN=""
ngsdose_engine() { if [ -n "$SIF" ]; then apptainer exec --bind "$APPTAINER_BINDS" "$SIF" ngs-dose "$@"; else "$NGSDOSE_BIN" "$@"; fi; }
# any other command of the image (ngsdose, python3, aria2c, test): run where the bundle paths are valid
ngsdose_py() { if [ -n "$SIF" ]; then apptainer exec --bind "$APPTAINER_BINDS" "$SIF" "$@"; else "$@"; fi; }
have_file() { ngsdose_py test -s "$1"; }

# scan  = read every record: the full-accuracy mode. Placement-independent, measures the
#         (experimental) satellite classes too, records where every class read was aligned and
#         what else sits in those bins (1 kb for positional classes, 10 kb for compositional ones),
#         and is what sinks are learned from and checked against. ~5 min wall per 30x genome on
#         8 CPUs with every panel loaded (-@ 8: 8 classifier + 8 decode threads; median elapsed_sec
#         of the cohort's 1,748 scans 304 s, 10-90% 217-450 s; CPU time not recorded); ~15 GB per
#         CRAM, so stage the CRAMs (01b_dose_sample.sh) rather than stream them if the whole cohort
#         is to be scanned.
# fetch = controls + learned sinks through the index: ~0.5 GB and ~1 min per sample straight from
#         the public bucket (~5 s from a local file), no CRAM on disk: 45S, 5S, DJ, the telomeric
#         repeat (FETCH_PANELS), the known-truth regions, chrM and chrEBV. (The 1,748 fetch counts in
#         counts_fetch/ predate the TEL sinks and hold 45S, 5S and DJ only.) The biobank-scale mode;
#         01b_dose_sample.sh runs it on every staged CRAM as well, so the cohort says what it costs
#         in accuracy, sample by sample.
# Counts go to $WORK_DIR/counts_$MODE; MODE selects which of them 01_count.sh makes and
# 02_cohort.sh analyses.
MODE="${MODE:-scan}"
COUNTS_DIR="${COUNTS_DIR:-$WORK_DIR/counts_$MODE}"
EST_DIR="${EST_DIR:-$WORK_DIR/estimates_$MODE}"
THREADS="${THREADS:-8}"
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-10}"
SAMPLE_TIMEOUT="${SAMPLE_TIMEOUT:-5400}"            # seconds; a stalled HTTPS connection can hang
S3_HTTPS_BASE="${S3_HTTPS_BASE:-https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage}"

# NGS-PCA outputs for this cohort (coverage PCs); the copy in this repository is the default
NGSPCA_DIR="${NGSPCA_DIR:-$REPO/meta/ngspca}"
PEDIGREE="${PEDIGREE:-$WORK_DIR/20130606_g1k_3202_samples_ped_population.txt}"
PEDIGREE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

# --- shared by the scripts -----------------------------------------------------------------------------

# Downloads (00_setup.sh, 01_stage_and_dose.sh, 01_count.sh, 04_hprc_satellites.sh) go to <dest>.part,
# resumed after an interruption, are checked, and only then take their own name: a file under its own
# name is whole. aria2c comes from the image when SIF is set and from the host otherwise; curl (on the
# host) is used when that aria2c cannot be run. Call `download` in a condition (`download ... || ...`).
have_aria2() { [ -n "${_ARIA2:-}" ] || { ngsdose_py aria2c --version >/dev/null 2>&1 && _ARIA2=1 || _ARIA2=0; }; [ "$_ARIA2" = 1 ]; }
md5_of() { (md5sum "$1" 2>/dev/null || md5 -r "$1") | awk '{print $1}'; }
file_ok() {  # path [check]: an MD5, "gzip" (gzip -t), "bed4" (every data line has >= 4 tab-separated fields), or none (non-empty)
  [ -s "$1" ] || return 1
  case "${2:-}" in
    "") ;;
    gzip) gzip -t "$1" 2>/dev/null ;;
    bed4) awk -F'\t' '!/^(track|browser|#)/ && NF && NF < 4 {bad = 1; exit} END {exit bad}' "$1" ;;
    *) [ "$(md5_of "$1")" = "$2" ] ;;
  esac
}
crai_ok() { [ ! -e "$1.aria2" ] && file_ok "$1" gzip; }     # a CRAM index that is whole (htslib reads a cut one without a word)
download() {  # url dest [check, as for file_ok]
  local url="$1" dest="$2" check="${3:-}" part="$2.part" md5="" try rc aria retry=(--retry 5)
  case "$check" in ""|gzip|bed4) ;; *) md5="$check" ;; esac
  curl --help all 2>/dev/null | grep -q -- --retry-all-errors && retry+=(--retry-all-errors)
  for try in 1 2 3 4; do
    if have_aria2; then
      aria=1
      ngsdose_py aria2c --quiet=true --console-log-level=warn --summary-interval=0 -x "${DOWNLOAD_CONNECTIONS:-8}" -s "${DOWNLOAD_CONNECTIONS:-8}" -k 16M \
        --file-allocation=none --auto-file-renaming=false --allow-overwrite=true --continue=true --max-tries=5 --retry-wait=20 \
        ${md5:+--checksum=md5=$md5} -d "$(dirname "$part")" -o "$(basename "$part")" "$url"; rc=$?
    else
      aria=0
      curl -sSL --fail "${retry[@]}" -C - -o "$part" "$url"; rc=$?
    fi
    if [ "$aria" = 1 ] && [ -e "$part.aria2" ]; then
      :                                                   # unfinished: the next try resumes it
    elif [ "$rc" = 0 ] || [ "$aria" = 1 ] || [ -n "$md5" ] || [ "$check" = gzip ]; then
      # finished, or a file that can be checked on its own (curl fails to resume a file that is already whole)
      if { [ "$rc" = 0 ] && [ "$aria" = 1 ] && [ -n "$md5" ]; } || file_ok "$part" "$check"; then
        mv -f "$part" "$dest" && rm -f "$part.aria2" && return 0
      fi
      if [ "$rc" = 0 ] || [ "$aria" = 1 ]; then rm -f "$part"; fi      # whole but wrong: not worth resuming
    fi
    if [ "$aria" = 0 ] && { [ "$rc" = 22 ] || [ "$rc" = 33 ]; }; then rm -f "$part"; fi   # an HTTP error (e.g. a range past the end of a
                                                                        # whole file) or a server that cannot resume: start again from zero
    echo "  download attempt $try failed: $(basename "$dest")" >&2
    [ "$try" = 4 ] || sleep $(( try * 30 ))
  done
  return 1
}

# Per-sample jobs (01_stage_and_dose.sh, 01a_dispatch_staged.sh): $WORK_DIR/dispatched/<sample> holds the
# SLURM id of the sample's latest job, written only once sbatch has accepted it; <sample>.tries has a line
# per job. A sample whose job has ended without counts is submitted again, up to MAX_TRIES jobs.
MAX_TRIES="${MAX_TRIES:-3}"
# job_gone: 0 when the sample's latest job has left the queue (squeue answered and does not list it) or the
# sample has none; 1 while it is queued or running; 2 when squeue did not answer (a busy slurmctld times
# out): the job may well be alive, so nothing is submitted again and no CRAM counts as kept on its account.
# The whole list of the user's jobs is asked for, not the one id: squeue -j errors on an id it has purged.
job_gone() {
  local jid out
  jid="$(cat "$WORK_DIR/dispatched/$1" 2>/dev/null)" || return 0
  jid="${jid%%;*}"; [ -n "$jid" ] || return 0
  command -v squeue >/dev/null 2>&1 || return 0                 # no SLURM here: nothing can be queued
  out="$(squeue -h -u "${USER:-$(id -un)}" -o %i 2>/dev/null)" || return 2
  if grep -qx -- "$jid" <<< "$out"; then return 1; fi        # a here-string, not a pipe: grep -q quitting early must not read as "gone" under pipefail
  return 0
}
tries_of() { if [ -s "$WORK_DIR/dispatched/$1.tries" ]; then wc -l < "$WORK_DIR/dispatched/$1.tries" | tr -d ' '; else echo 0; fi; }
dispatch() {  # sample line [md5]: submit 01b_dose_sample.sh
  local id
  mkdir -p "$WORK_DIR/dispatched"
  id="$(sbatch --parsable "$EX_DIR/01b_dose_sample.sh" "$@")" && [ -n "$id" ] || { echo "[$1] sbatch failed" >&2; return 1; }
  echo "$id" > "$WORK_DIR/dispatched/$1"; echo "$id" >> "$WORK_DIR/dispatched/$1.tries"
}
# hold_lock FD FILE: an exclusive lock on FILE through file descriptor FD (flock, where it exists). Returns 1
# only when another process holds it; where locks cannot be taken (no flock, a file system without them) it
# says so and returns 0. `wait` blocks until the lock is free instead.
hold_lock() {
  local rc=0
  command -v flock >/dev/null 2>&1 || return 0
  eval "exec $1<>\"\$2\"" 2>/dev/null || { echo "note: cannot open $2; running without a lock" >&2; return 0; }
  if [ "${3:-}" = wait ]; then flock "$1" || rc=$?; else flock -n "$1" || rc=$?; fi
  [ "$rc" = 1 ] && [ "${3:-}" != wait ] && return 1
  [ "$rc" = 0 ] || echo "note: flock on $2 failed (status $rc); running without a lock" >&2
  return 0
}

# check_candidate_panels: the CANDIDATE_PANELS are panels (k = 31, class lines, k-mers) that share no k-mer and no
# class name with the bundle panel, EXTRA_PANELS and FETCH_PANELS, nor with each other. Standard library Python in
# the image, no data needed; a pass is remembered in $WORK_DIR/candidate_panels.checked until a panel file or the
# image changes.
check_candidate_panels() {
  local stamp="$WORK_DIR/candidate_panels.checked" key x shipped="" rc=0
  [ -n "${CANDIDATE_PANELS:-}" ] || return 0
  if [ -n "$CANDIDATE_PANELS_ERROR" ]; then echo "ERROR: CANDIDATE_PANELS: $CANDIDATE_PANELS_ERROR" >&2; return 1; fi
  for x in "$BUNDLE/panel.k31.tsv.gz" $SHIPPED_SCAN_PANELS ${FETCH_PANELS:-}; do
    case " $shipped " in *" $x "*) ;; *) shipped="$shipped $x" ;; esac
  done
  key="$(echo "image ${SIF:-$NGSDOSE_BIN} $( [ -z "$SIF" ] || wc -c < "$SIF" | tr -d ' ')"; echo "beside$shipped"
         for x in $CANDIDATE_PANEL_FILES; do echo "candidate $x $(md5_of "$x")"; done)"
  [ -s "$stamp" ] && [ "$(cat "$stamp")" = "$key" ] && return 0
  # shellcheck disable=SC2086
  ngsdose_py python3 - $CANDIDATE_PANEL_FILES -- $shipped <<'PY' || rc=$?
import gzip, os, sys
from collections import Counter
args = sys.argv[1:]
cands, beside = args[:args.index("--")], args[args.index("--") + 1:]

def opened(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)

def header(path):
    k, names = None, []
    with opened(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##k="):
                k = int(line[4:])
            elif line.startswith("##class\t"):
                names.append(dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:] if "=" in f).get("name", ""))
    return k, names

def kmers(path):
    with opened(path) as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                yield line.split("\t", 1)[0]

try:
    bad, owner, ks, classes = [], {}, set(), {}
    beside = [p for p in beside if os.path.exists(p) or print(f"note: {p} not found; the candidates are not checked against it", file=sys.stderr)]
    for p in beside:
        k, classes[p] = header(p)
        ks.add(k)
        for n in classes[p]:
            owner.setdefault(n, p)
    of, shared, total = {}, Counter(), 0
    for p in cands:
        k, classes[p] = header(p)
        if not classes[p]:
            bad.append(f"{p}: not an ngs-dose panel (no ##class lines)")
            continue
        if ks and k not in ks:
            bad.append(f"{p}: k={k}, the panels beside it have k={', '.join(map(str, sorted(ks)))}")
        for n in classes[p]:
            if n in owner:
                bad.append(f"{p}: class {n} is a class of {owner[n]} too")
            owner.setdefault(n, p)
        n0 = len(of)
        for km in kmers(p):
            q = of.setdefault(km, p)
            if q != p:
                shared[(q, p)] += 1
        if len(of) == n0:
            bad.append(f"{p}: no k-mers")
        total += len(of) - n0
    for p in beside:
        hit = Counter(of[km] for km in kmers(p) if km in of)
        shared.update({(q, p): n for q, n in hit.items()})
    for (q, p), n in sorted(shared.items()):
        bad.append(f"{q} shares {n:,} k-mers with {p} (classes {', '.join(classes[p])}): the engine would drop them from both")
    for b in bad:
        print(f"  {b}", file=sys.stderr)
    if not bad:
        print(f"candidate panels: {', '.join(n for p in cands for n in classes[p])} ({total:,} k-mers in {len(cands)} file(s)); no k-mer "
              f"and no class name shared with each other or with the {len(beside)} panel(s) loaded beside them", file=sys.stderr)
    sys.exit(1 if bad else 0)
except (OSError, EOFError, ValueError) as e:
    print(f"  {e}", file=sys.stderr)
    sys.exit(2)
PY
  case "$rc" in
    0) mkdir -p "$WORK_DIR" && printf '%s\n' "$key" > "$stamp"; return 0 ;;
    1) echo "ERROR: the candidate panels clash with the panels loaded beside them (above). A k-mer claimed by two loaded panels" >&2
       echo "  is dropped from both, so the classes the cohort already measures would be counted with other k-mers from here on." >&2
       echo "  Rebuild the candidates with those panels as background (ngs-dose panel), or rename the class." >&2 ;;
    *) echo "ERROR: could not check the candidate panels (above)" >&2 ;;
  esac
  return 1
}

# A fetch plan (FETCH_CLASSES, FETCH_PRESET, FETCH_BUDGET_MB, FETCH_CAPTURE): fetch_plan_key lists what the plan is
# built from - the fetchplan arguments and the SHA-256 of the menu, of every file it names and of every file among
# the arguments but the engine - so that a change of any of them is seen. use_fetch_plan builds the plan once into FETCH_PLAN_DIR
# (under a lock; the first run that needs it builds it, the others wait and reuse it), refuses a plan built from
# something else, and sets FETCH_SINKS_BED, FETCH_CONTROLS, FETCH_PANELS and FETCH_FLAGS to the plan's. The plan's
# count flags are checked against the engine that will fetch: --classes and --allow-missing-sinks need an engine newer than fae1124.
fetch_plan_key() {
  ngsdose_py python3 - "$@" <<'PY'
import hashlib, os, sys
args = sys.argv[1:]
menu = args[args.index("--menu") + 1]

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

# the engine's file is not hashed: an engine rebuilt from the same source must not refuse the plan (EXPECTED_ENGINE_BUILD
# watches the build, and the count flags are checked against the engine itself)
files, head = [a for i, a in enumerate(args) if os.path.isfile(a) and (i == 0 or args[i - 1] != "--engine")], None
with open(menu) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        row = line.rstrip("\n").split("\t")
        if head is None:
            head = row
            continue
        for col in ("panel", "sinks", "stats"):
            v = dict(zip(head, row)).get(col, "-").strip()
            if v not in ("", "-", "."):
                v = v if os.path.isabs(v) else os.path.join(os.path.dirname(menu), v)
                if os.path.isfile(v):
                    files.append(v)
print("fetchplan " + " ".join(args))
for f in dict.fromkeys(files):
    print(f"sha256 {sha(f)} {f}")
PY
}
use_fetch_plan() {
  local dir="${FETCH_PLAN_DIR%/}" args=() cost=() crais="" key x s url rc=0 n=0 tried=0 tmp help ctrl
  [ -z "$FETCH_PLAN" ] || return 0                                   # loaded already in this shell
  if [ -z "$FETCH_CLASSES$FETCH_PRESET$FETCH_BUDGET_MB" ]; then
    echo "ERROR: FETCH_CAPTURE and FETCH_PLAN_ARGS tune a fetch plan; choose its options with FETCH_CLASSES, FETCH_PRESET or" >&2
    echo "  FETCH_BUDGET_MB (FETCH_PRESET=core is 45S, 5S and DJ, as the fae1124 image fetches; core_tel adds the telomeric repeat)" >&2
    return 1
  fi
  if [ "${FETCH_PANELS_DEFAULT:-0}" != 1 ] && [ -n "${FETCH_PANELS:-}" ]; then
    echo "ERROR: FETCH_PANELS and a fetch plan (FETCH_CLASSES, FETCH_PRESET, FETCH_BUDGET_MB) both say what a fetch loads, and" >&2
    echo "  the plan takes its panels from the menu: unset FETCH_PANELS" >&2
    return 1
  fi
  if ! help="$(ngsdose_py ngsdose fetchplan --help 2>&1)"; then
    echo "ERROR: FETCH_CLASSES, FETCH_PRESET, FETCH_BUDGET_MB and FETCH_CAPTURE need \`ngsdose fetchplan\`, which the ngsdose of" >&2
    echo "  this image does not have (the fae1124 image, which counted the cohort's first 1,748 genomes, predates it). Unset" >&2
    echo "  them to fetch as before, or use a newer image (NGSDOSE_IMAGE, SIF) with EXPECTED_ENGINE_BUILD set to match." >&2
    return 1
  fi
  args=(--menu "$FETCH_MENU")
  # shellcheck disable=SC2206
  { [ -z "$FETCH_CLASSES" ] || args+=(--classes $FETCH_CLASSES); [ -z "$FETCH_PRESET" ] || args+=(--preset $FETCH_PRESET); }
  [ -z "$FETCH_BUDGET_MB" ] || args+=(--budget-mb "$FETCH_BUDGET_MB")
  [ -z "$FETCH_CAPTURE" ] || args+=(--capture "$FETCH_CAPTURE")
  # shellcheck disable=SC2206
  [ -z "$FETCH_PLAN_ARGS" ] || args+=($FETCH_PLAN_ARGS)
  # the plan's count flags depend on what the engine takes (--classes): fetchplan asks `ngs-dose` on the image's PATH,
  # which is the engine with SIF set; without it, the engine is NGSDOSE_BIN
  case " ${args[*]} " in *" --engine "*|*" --engine="*) ;; *) if [ -z "$SIF" ]; then case "$help" in *--engine*) args+=(--engine "$NGSDOSE_BIN") ;; esac; fi ;; esac
  key="$(fetch_plan_key "${args[@]}")" || { echo "ERROR: could not read the fetch menu $FETCH_MENU, or a file it names (above)" >&2; return 1; }
  mkdir -p "$(dirname "$dir")"
  hold_lock 7 "$dir.lock" wait
  if [ ! -d "$dir" ]; then
    if [ -n "$FETCH_CRAI" ]; then
      crais="$FETCH_CRAI"
    elif [ -r "$MANIFEST" ]; then                                    # the first samples' indexes; a failed download ends the search
      while IFS=$'\t' read -r s url _ <&3; do
        [ "$n" -lt "$FETCH_COST_SAMPLES" ] && [ "$tried" -lt $(( 3 * FETCH_COST_SAMPLES )) ] || break
        tried=$((tried + 1))
        case "$url" in
          *://*) x="$WORK_DIR/crai/$s.crai"; mkdir -p "$WORK_DIR/crai"
                 crai_ok "$x" || download "$url.crai" "$x" gzip || break ;;
          *)     x="$url.crai"; [ -s "$x" ] || continue ;;
        esac
        crais="$crais $x"; n=$((n + 1))
      done 3< "$MANIFEST"
    fi
    if [ -n "$crais" ] && have_file "$REF_FASTA.fai"; then
      # shellcheck disable=SC2206
      cost=(--crai $crais --contigs "$REF_FASTA.fai")
    elif [ -n "$FETCH_BUDGET_MB" ]; then
      echo "ERROR: FETCH_BUDGET_MB needs CRAM indexes to cost the options on (FETCH_CRAI, or the manifest's first samples' indexes)" >&2
      echo "  and $REF_FASTA.fai; neither could be had" >&2
      exec 7>&-; return 1
    else
      echo "note: the fetch plan is not costed (no CRAM index, or no $REF_FASTA.fai)" >&2
    fi
    tmp="$dir.tmp.$$"; rm -rf "$tmp"; mkdir -p "$tmp"
    ngsdose_py ngsdose fetchplan "${args[@]}" ${cost[@]+"${cost[@]}"} -o "$tmp/plan" > "$tmp/plan.table.txt" 2> "$tmp/plan.notes.txt" || rc=$?
    cat "$tmp/plan.notes.txt" >&2
    if [ "$rc" != 0 ] || [ ! -s "$tmp/plan.sinks.bed" ] || [ ! -s "$tmp/plan.panels.txt" ] || [ ! -e "$tmp/plan.count_flags.txt" ] \
       || [ ! -s "$tmp/plan.controls.txt" ]; then
      echo "ERROR: ngsdose fetchplan ${args[*]} failed (above); no plan was written" >&2
      rm -rf "$tmp"; exec 7>&-; return 1
    fi
    printf '%s\n' "$key" > "$tmp/settings.txt"
    echo "${crais# }" > "$tmp/costed_on.txt"
    mv "$tmp" "$dir"
    echo "fetch plan built in $dir:" >&2; sed 's/^/  /' "$dir/plan.table.txt" >&2
  elif [ "$(cat "$dir/settings.txt" 2>/dev/null)" != "$key" ]; then
    echo "ERROR: the fetch plan in $dir was built from other settings or files:" >&2
    diff "$dir/settings.txt" <(printf '%s\n' "$key") | sed 's/^/  /' >&2
    echo "  A cohort's fetches should share one plan (each counts file records the SHA-256 of the sinks BED it used). To fetch" >&2
    echo "  by the new one, remove $dir or set FETCH_PLAN_DIR to a new directory; the fetches made so far keep the old plan." >&2
    exec 7>&-; return 1
  fi
  exec 7>&-
  # the controls FASTA the plan was costed on (fetchplan --controls: the bundle's controls.lite200, say), which `ngsdose
  # estimate` takes only when it is the bundle's controls.fa.gz or a subset the bundle names (controls.NAME.bed beside it)
  ctrl="$(head -n 1 "$dir/plan.controls.txt" 2>/dev/null)"
  if [ -z "$ctrl" ] || ! have_file "$ctrl"; then
    echo "ERROR: the fetch plan in $dir names no controls FASTA, or one that is not there (${ctrl:-plan.controls.txt is empty}):" >&2
    echo "  build it (ngs-dose controls -b BED -T REFERENCE --flank 1000 -o FASTA) or remove $dir" >&2
    return 1
  fi
  x="$(basename "$ctrl")"; x="${x#controls.}"; x="${x%.fa.gz}"
  if ! ngsdose_py python3 -c 'import os, sys; sys.exit(not os.path.samefile(os.path.dirname(sys.argv[1]), sys.argv[2]))' "$ctrl" "$BUNDLE" 2>/dev/null \
     || { [ "$(basename "$ctrl")" != controls.fa.gz ] && { [ "$(basename "$ctrl")" != "controls.$x.fa.gz" ] || ! have_file "$BUNDLE/controls.$x.bed"; }; }; then
    echo "ERROR: the fetch plan in $dir fetches the control regions of $ctrl, which is neither $BUNDLE/controls.fa.gz nor" >&2
    echo "  the FASTA of a subset the bundle names ($BUNDLE/controls.NAME.fa.gz beside controls.NAME.bed): \`ngsdose estimate\`" >&2
    echo "  would refuse every counts file it made. Use --controls $BUNDLE/controls.lite200.bed (or none) and remove $dir" >&2
    return 1
  fi
  FETCH_PLAN="$dir/plan"; FETCH_SINKS_BED="$dir/plan.sinks.bed"; FETCH_PANELS_DEFAULT=0
  FETCH_PANELS="$(tr '\n' ' ' < "$dir/plan.panels.txt")"; FETCH_PANELS="${FETCH_PANELS% }"
  FETCH_FLAGS="$(tr '\n' ' ' < "$dir/plan.count_flags.txt")"; FETCH_FLAGS="${FETCH_FLAGS% }"
  FETCH_CONTROLS="$ctrl"
  if [ -n "$FETCH_FLAGS" ]; then                     # the flags older engines lack must be ones this engine takes
    help="$(ngsdose_engine count --help 2>&1)" || true
    for x in $FETCH_FLAGS; do
      case "${x%%=*}" in --classes|--allow-missing-sinks) ;; *) continue ;; esac
      case "$help" in *"${x%%=*}"*) ;; *)
        echo "ERROR: the fetch plan in $dir needs \`count ${x%%=*}\` (plan.count_flags.txt), which this engine does not take." >&2
        echo "  The fae1124 engine takes neither --classes nor --allow-missing-sinks, so a plan that loads only part of a panel" >&2
        echo "  file's classes needs a newer image (NGSDOSE_IMAGE, SIF, EXPECTED_ENGINE_BUILD); a plan of whole panel files" >&2
        echo "  (FETCH_PRESET=core or core_tel) needs neither. Remove $dir after changing either." >&2
        return 1 ;;
      esac
    done
  fi
  for x in $FETCH_PANELS; do have_file "$x" || { echo "ERROR: $x, a panel of the fetch plan in $dir, is not there" >&2; return 1; }; done
  echo "note: the fetch follows the plan in $dir: $(grep -vc '^#' "$FETCH_SINKS_BED") sink intervals; controls $(basename "$FETCH_CONTROLS");" \
    "panels $FETCH_PANELS; flags ${FETCH_FLAGS:-none}" >&2
}

# estimate_sinks_args (call it in the script's own shell): ESTIMATE_SINKS_ARGS, the `--fetch-sinks BED...` of `ngsdose
# estimate`, which knows the sinks BED of a fetch by its SHA-256 (the counts file records no more) and reports the named
# sub-options of a satellite family (DXZ1, DYZ3, DYZ1, DYZ2) only for fetches whose BED it knows. The BEDs are
# FETCH_SINKS_KNOWN (space-separated, e.g. the plans of earlier fetches), or else the plan in FETCH_PLAN_DIR, if there is
# one. Without either, and with an ngsdose that predates the option (the fae1124 image), it is empty: the bundle's own
# sinks BEDs are known without it.
FETCH_SINKS_KNOWN="${FETCH_SINKS_KNOWN:-}"
estimate_sinks_args() {
  local given="$FETCH_SINKS_KNOWN" beds="" x d help
  ESTIMATE_SINKS_ARGS=()
  [ -n "$given" ] || [ ! -s "${FETCH_PLAN_DIR%/}/plan.sinks.bed" ] || given="${FETCH_PLAN_DIR%/}/plan.sinks.bed"
  for x in $given; do                                        # as absolute paths, bound, the same inside the container
    [ -f "$x" ] || { echo "ERROR: FETCH_SINKS_KNOWN: $x is not a file on this host" >&2; return 1; }
    d="$(cd "$(dirname "$x")" && pwd)"; beds="$beds $d/$(basename "$x")"
    [ -z "$SIF" ] || bind_dir "$d"
  done
  beds="${beds# }"
  [ -n "$beds" ] || return 0
  help="$(ngsdose_py ngsdose estimate --help 2>&1)" || true
  case "$help" in *--fetch-sinks*) ;; *)
    echo "note: this image's ngsdose estimate has no --fetch-sinks; the sub-options of fetches made with $beds are not reported" >&2
    return 0 ;;
  esac
  # shellcheck disable=SC2034,SC2206
  ESTIMATE_SINKS_ARGS=(--fetch-sinks $beds)
}

# The image and the scripts have to agree. check_fetch_panels: every class of FETCH_PANELS has sinks in
# the image's bundle (an engine older than NGS-DOSE 645ae55 would fetch it anyway and count only the reads
# that happen to fall in the controls and in other classes' sinks: about 1% of the telomeric repeat). The
# default FETCH_PANELS drops such a panel, with a note; an explicit one stops. With a fetch plan asked for, the
# plan says what is loaded instead (use_fetch_plan). The candidate panels are checked first (check_candidate_panels),
# so that the stager stops before it stages anything. Call it in the script's own
# shell (not in $(...)): it sets FETCH_PANELS. check_build: a counts file was written by EXPECTED_ENGINE_BUILD.
check_fetch_panels() {
  local x rc keep=""
  check_candidate_panels || return 1
  if [ -n "$FETCH_CLASSES$FETCH_PRESET$FETCH_BUDGET_MB$FETCH_CAPTURE$FETCH_PLAN_ARGS" ]; then use_fetch_plan; return; fi
  [ ! -d "$FETCH_PLAN_DIR" ] || echo "WARNING: $FETCH_PLAN_DIR holds a fetch plan, but FETCH_CLASSES, FETCH_PRESET and FETCH_BUDGET_MB are unset" \
    "here (not exported?): this fetch reads the bundle's sinks, as without a plan" >&2
  [ -n "${FETCH_PANELS:-}" ] || return 0
  if [ "${FETCH_PANELS_DEFAULT:-0}" = 1 ]; then
    for x in $FETCH_PANELS; do
      if ! have_file "$x"; then echo "note: fetch panel $x is not in this image; not loaded" >&2; continue; fi
      rc=0; ngsdose_py python3 "$EX_DIR/check_counts.py" fetch-panels "$BUNDLE/sinks.bed" "$x" || rc=$?
      case "$rc" in
        0) keep="$keep $x" ;;
        1) echo "note: the image's bundle predates the sinks of $(basename "$x") (above), so the fetch leaves that panel" >&2
           echo "  out and loads the bundle's classes alone, as the fae1124 image did. Set FETCH_PANELS to choose." >&2 ;;
        *) echo "ERROR: could not check $x against $BUNDLE/sinks.bed (above)" >&2; return 1 ;;
      esac
    done
    FETCH_PANELS="${keep# }"; FETCH_PANELS_DEFAULT=0
    return 0
  fi
  for x in $FETCH_PANELS; do have_file "$x" || { echo "ERROR: fetch panel $x not found" >&2; return 1; }; done
  # shellcheck disable=SC2086
  ngsdose_py python3 "$EX_DIR/check_counts.py" fetch-panels "$BUNDLE/sinks.bed" $FETCH_PANELS && return 0
  echo "ERROR: FETCH_PANELS loads a class the bundle has no sinks for (above), so a fetch would silently miss" >&2
  echo "  nearly all of its reads. The image predates those sinks (the telomeric repeat's arrived in NGS-DOSE" >&2
  echo "  45c9f9a; the fae1124 image has none). Unset FETCH_PANELS (the default leaves such a panel out), set" >&2
  echo "  FETCH_PANELS=\"\" to fetch the bundle's classes only, or use an image whose bundle has them (NGSDOSE_IMAGE," >&2
  echo "  SIF), and set EXPECTED_ENGINE_BUILD to match." >&2
  return 1
}
engine_build_of() { ngsdose_py python3 "$EX_DIR/check_counts.py" engine-build "$1"; }
check_build() {  # counts file
  local b
  [ -n "${EXPECTED_ENGINE_BUILD:-}" ] || return 0
  b="$(engine_build_of "$1")" || b=""
  case "$b" in "$EXPECTED_ENGINE_BUILD"*) return 0 ;; esac
  echo "ERROR: $1 was written by engine build ${b:-(unreadable)}, not EXPECTED_ENGINE_BUILD=$EXPECTED_ENGINE_BUILD:" >&2
  echo "  the image is not the one this cohort is pinned to (SIF, NGSDOSE_IMAGE)" >&2
  return 1
}
note_build() {  # what the counts so far were made with, against EXPECTED_ENGINE_BUILD; a note, never a stop
  local f first="" b
  for f in "$WORK_DIR"/counts_scan/*.json.gz; do [ -s "$f" ] && { first="$f"; break; }; done
  [ -n "$first" ] || return 0
  b="$(engine_build_of "$first" 2>/dev/null)" || return 0
  if [ -z "${EXPECTED_ENGINE_BUILD:-}" ]; then
    echo "note: the counts so far carry engine_build ${b:0:7}; EXPECTED_ENGINE_BUILD=${b:0:7} keeps new ones to it"
  else
    case "$b" in "$EXPECTED_ENGINE_BUILD"*) ;; *) echo "WARNING: the counts so far carry engine_build ${b:0:7}, new ones will carry $EXPECTED_ENGINE_BUILD: the cohort mixes engine builds" >&2 ;; esac
  fi
}
