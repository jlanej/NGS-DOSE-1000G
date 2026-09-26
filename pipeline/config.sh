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
# staged CRAMs for the full-accuracy path (01b_dose_sample.sh): $CRAM_DIR/<sample>.cram(.crai)
CRAM_DIR="${CRAM_DIR:-$WORK_DIR/crams}"

# CRAM decoding needs the reference: a local FASTA (recommended) or htslib's REF_PATH/REF_CACHE
REF_FASTA="${REF_FASTA:-$WORK_DIR/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa}"

# what the container has to see (bind sources must exist: 00_setup.sh creates them)
APPTAINER_BINDS="${APPTAINER_BINDS:-$WORK_DIR,$REPO,$(dirname "$REF_FASTA"),$CRAM_DIR}"
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
    if [ "$aria" = 0 ] && [ "$rc" = 22 ]; then rm -f "$part"; fi       # an HTTP error, e.g. a range past the end of a whole file
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

# The image and the scripts have to agree. check_fetch_panels: every class of FETCH_PANELS has sinks in
# the image's bundle (an engine older than NGS-DOSE 645ae55 would fetch it anyway and count only the reads
# that happen to fall in the controls and in other classes' sinks: about 1% of the telomeric repeat). The
# default FETCH_PANELS drops such a panel, with a note; an explicit one stops. Call it in the script's own
# shell (not in $(...)): it sets FETCH_PANELS. check_build: a counts file was written by EXPECTED_ENGINE_BUILD.
check_fetch_panels() {
  local x rc keep=""
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
