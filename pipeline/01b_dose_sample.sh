#!/usr/bin/env bash
# =============================================================================
# 01b_dose_sample.sh - NGS-DOSE for one sample whose CRAM is already local
# =============================================================================
# The full-accuracy path: a whole-file scan of a staged CRAM (every class, the
# satellite and telomere panels, where every class read was aligned and what else
# sits in those bins: 1 kb for positional classes, 10 kb for compositional ones),
# then - about five seconds more - a targeted fetch of the same file, so that every
# sample of the cohort says exactly what the fast mode would have returned for it.
# The CRAM is removed only once the outputs are verified: both counts files whole,
# both from EXPECTED_ENGINE_BUILD when that is set (config.sh), and the fetch's
# control and region reads identical to the scan's, as they are by construction,
# with every class read the scan placed inside the sinks; a fetch that differs
# (e.g. through a cut index) is removed and the CRAM kept. A counts file is removed
# only while the CRAM is here to count it again: without the CRAM (a requeue after
# the tidy, a second job) a bad fetch is moved aside to <file>.rejected, so that the
# stager stages the CRAM again, and a file of another build is left where it is.
# Before any counting: the index must be whole (gzip -t), and FETCH_PANELS is
# checked against the bundle's sinks (config.sh). One job per sample at a time: a second one
# waits for the first (flock on $WORK_DIR/dispatched/<sample>.lock, where flock exists).
#
# Same calling convention as NGS-PCA's 01b_mosdepth_sample.sh, so whatever stages
# CRAMs for that pipeline (its aria2 download manager, its stage watcher, a Globus
# batch) can dispatch this job instead:
#   sbatch 01b_dose_sample.sh SAMPLE_ID MANIFEST_LINE_NUM [CRAM_MD5]
# Expects $CRAM_DIR/SAMPLE_ID.cram and .cram.crai. KEEP_CRAMS=1 leaves them in place.
# =============================================================================
#SBATCH --job-name=ngsdose_sample
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/dose_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
SAMPLE="${1:?usage: 01b_dose_sample.sh SAMPLE_ID MANIFEST_LINE_NUM [CRAM_MD5]}"
LINE_NUM="${2:-0}"; CRAM_MD5="${3:-}"
CRAM="$CRAM_DIR/$SAMPLE.cram"; CRAI="$CRAM.crai"
SCAN_OUT="$WORK_DIR/counts_scan/$SAMPLE.json.gz"; FETCH_OUT="$WORK_DIR/counts_fetch/$SAMPLE.json.gz"
mkdir -p "$WORK_DIR/counts_scan" "$WORK_DIR/counts_fetch"
echo "== NGS-DOSE: $SAMPLE (manifest line $LINE_NUM), started $(date)"

tidy() { [ "${KEEP_CRAMS:-0}" = 1 ] && return 0; [ -e "$CRAM" ] && echo "   removed $CRAM"; rm -f "$CRAM" "$CRAI"; }
valid() { [ -s "$1" ] && gzip -t "$1" 2>/dev/null; }          # a counts file that exists and is whole
same_reads() { ngsdose_py python3 "$EX_DIR/check_counts.py" same-reads "$SCAN_OUT" "$FETCH_OUT" "$BUNDLE/sinks.bed"; }
agree() {  # the fetch saw exactly the scan's control and region reads; a fetch that did not is removed
  local rc=0
  same_reads || rc=$?
  [ "$rc" = 0 ] && return 0
  if [ "$rc" = 1 ]; then rm -f "$FETCH_OUT"; echo "   the fetch does not match the scan (above): $FETCH_OUT removed" >&2; fi
  return "$rc"                                           # 1: they differ (fetch removed); 2: they could not be compared
}
own_build() {  # a counts file from another engine build is removed
  check_build "$1" && return 0
  rm -f "$1"; echo "   removed $1" >&2; return 1
}

check_fetch_panels || exit 1
mkdir -p "$WORK_DIR/dispatched"
hold_lock 9 "$WORK_DIR/dispatched/$SAMPLE.lock" wait     # a second job for this sample waits for the first
# counts are removed and made again only while the CRAM and a whole index are here to make them from
if [ -s "$CRAM" ] && [ -s "$CRAI" ] && gzip -t "$CRAI" 2>/dev/null; then have_cram=1; else have_cram=0; fi
if [ "$have_cram" = 1 ]; then
  for f in "$SCAN_OUT" "$FETCH_OUT"; do valid "$f" && { own_build "$f" || true; }; done
  if valid "$SCAN_OUT" && valid "$FETCH_OUT" && agree; then echo "   already done"; tidy; exit 0; fi   # idempotent under a requeue
elif valid "$SCAN_OUT" && valid "$FETCH_OUT"; then
  # nothing here to count again with: judge the counts, remove none of them
  ok=1; for f in "$SCAN_OUT" "$FETCH_OUT"; do check_build "$f" || ok=0; done
  [ "$ok" = 1 ] || { echo "ERROR: counts of another engine build (above), left in place: without $CRAM and a whole index they can only be redone by removing them and staging it again" >&2; exit 1; }
  rc=0; same_reads || rc=$?
  [ "$rc" = 0 ] && { echo "   already done"; tidy; exit 0; }
  [ "$rc" = 1 ] && { mv -f "$FETCH_OUT" "$FETCH_OUT.rejected"; echo "ERROR: the fetch does not match the scan (above) and there is no CRAM with a whole index to fetch again from: moved to $FETCH_OUT.rejected; the stager stages the CRAM again for a new fetch" >&2; exit 1; }
  echo "ERROR: could not compare the fetch of $SAMPLE with its scan (above); both left in place" >&2; exit 1
fi
[ -s "$CRAM" ] && [ -s "$CRAI" ] || { echo "ERROR: $CRAM or its index is missing" >&2; exit 1; }
# htslib reads a cut .crai without an error, and a fetch through it silently loses every contig past the cut
gzip -t "$CRAI" 2>/dev/null || { echo "ERROR: $CRAI is not a whole index (a download cut short?); the CRAM is kept" >&2; exit 1; }
if [ -n "$CRAM_MD5" ]; then
  got="$(md5_of "$CRAM")"
  [ "$got" = "$CRAM_MD5" ] || { echo "ERROR: MD5 mismatch for $CRAM ($got, expected $CRAM_MD5)" >&2; exit 1; }
  echo "   MD5 verified"
fi

common=(-i "$CRAM" --index "$CRAI" -T "$REF_FASTA" -c "$BUNDLE/controls.fa.gz" -p "$BUNDLE/panel.k31.tsv.gz" -s "$SAMPLE" -@ "${SLURM_CPUS_PER_TASK:-$THREADS}")
if ! valid "$SCAN_OUT"; then
  scan=("${common[@]}" -m scan)
  # a scan without a panel that was asked for is a scan to be done again: stop rather than skip it
  for x in $EXTRA_PANELS; do have_file "$x" || { echo "ERROR: extra panel $x not found" >&2; exit 1; }; scan+=(-p "$x"); done
  ngsdose_engine count "${scan[@]}" -o "$SCAN_OUT"
  own_build "$SCAN_OUT" || exit 1
fi
fetch=("${common[@]}" -m fetch --sinks "$BUNDLE/sinks.bed")
for x in ${FETCH_PANELS:-}; do fetch+=(-p "$x"); done
if ! valid "$FETCH_OUT"; then
  ngsdose_engine count "${fetch[@]}" -o "$FETCH_OUT"
  own_build "$FETCH_OUT" || exit 1
fi
valid "$SCAN_OUT" && valid "$FETCH_OUT" || { echo "ERROR: counts missing or damaged for $SAMPLE; the CRAM is kept" >&2; exit 1; }
rc=0; agree || rc=$?
if [ "$rc" = 1 ]; then echo "ERROR: the fetch of $SAMPLE does not match its scan; the CRAM is kept" >&2; exit 1; fi
if [ "$rc" != 0 ]; then echo "ERROR: could not compare the fetch of $SAMPLE with its scan (above); both counts files and the CRAM are kept" >&2; exit 1; fi
tidy
echo "== $SAMPLE complete, $(date)"
