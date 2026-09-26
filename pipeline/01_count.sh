#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_count
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/count_%A_%a.out
# Per-sample counting, one block of the manifest per array task. Submit from this directory
# (config.sh is read from the submit directory, and logs/ is relative to it):
#   N=$(wc -l < $WORK_DIR/manifest.tsv); sbatch --array=0-$(( (N - 1) / ${SAMPLES_PER_TASK:-10} ))%25 01_count.sh
#   N=$(wc -l < $WORK_DIR/manifest.hprc.tsv)
#   MODE=scan MANIFEST=$WORK_DIR/manifest.hprc.tsv sbatch --array=0-$(( (N - 1) / ${SAMPLES_PER_TASK:-10} ))%10 01_count.sh
# (%25 caps concurrent tasks: the bottleneck is the network, and a few hundred parallel range
#  readers against one S3 prefix earn throttling, not speed.)
# Without SLURM:  for i in $(seq 0 $(( (N - 1) / ${SAMPLES_PER_TASK:-10} ))); do SLURM_ARRAY_TASK_ID=$i bash 01_count.sh; done
# Idempotent: finished samples are skipped, failures are reported and retried by re-submitting. A task
# exits non-zero when any of its samples failed, when the manifest cannot be read, or when an array
# from 0 does not reach the end of the manifest; each ends by saying how many manifest samples lack
# counts.
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
mkdir -p "$COUNTS_DIR" "$WORK_DIR/crai" "$LOG_DIR"
task="${SLURM_ARRAY_TASK_ID:-0}"
[ -r "$MANIFEST" ] || { echo "ERROR: manifest $MANIFEST not readable" >&2; exit 1; }
N=$(( $(wc -l < "$MANIFEST") ))
# an array from 0 without gaps (the submit lines above) has to reach the end of the manifest; one that
# re-runs a few tasks (--array=17,40) is not checked
if [ "${SLURM_ARRAY_TASK_MIN:-}" = 0 ] && [ "${SLURM_ARRAY_TASK_COUNT:-}" = "$(( ${SLURM_ARRAY_TASK_MAX:-0} + 1 ))" ] \
   && [ $(( (SLURM_ARRAY_TASK_MAX + 1) * SAMPLES_PER_TASK )) -lt "$N" ]; then
  echo "ERROR: --array=0-$SLURM_ARRAY_TASK_MAX at SAMPLES_PER_TASK=$SAMPLES_PER_TASK covers $(( (SLURM_ARRAY_TASK_MAX + 1) * SAMPLES_PER_TASK )) of the manifest's $N lines; submit --array=0-$(( (N - 1) / SAMPLES_PER_TASK ))" >&2
  exit 1
fi
[ "$MODE" = scan ] || check_fetch_panels || exit 1
# a stalled HTTPS connection can hang forever: bound every sample where `timeout` exists (GNU coreutils)
TMO=(); command -v timeout >/dev/null 2>&1 && TMO=(timeout "$SAMPLE_TIMEOUT")
run() { if [ -n "$SIF" ]; then "${TMO[@]}" apptainer exec --bind "$APPTAINER_BINDS" "$SIF" ngs-dose "$@"; else "${TMO[@]}" "$NGSDOSE_BIN" "$@"; fi; }

fail=0
while IFS=$'\t' read -r sample cram _; do
  out="$COUNTS_DIR/$sample.json.gz"
  [ -s "$out" ] && continue
  args=(-i "$cram" -T "$REF_FASTA" -c "$BUNDLE/controls.fa.gz" -p "$BUNDLE/panel.k31.tsv.gz" -@ "$THREADS" -s "$sample" -o "$out.tmp.gz")
  if [ "$MODE" = scan ]; then
    for x in $EXTRA_PANELS; do have_file "$x" || { echo "ERROR: extra panel $x not found" >&2; exit 1; }; args+=(-p "$x"); done
    args+=(-m scan)
  else
    args+=(-m fetch --sinks "$BUNDLE/sinks.bed")
    for x in ${FETCH_PANELS:-}; do have_file "$x" || { echo "ERROR: fetch panel $x not found" >&2; exit 1; }; args+=(-p "$x"); done
    case "$cram" in
      *://*) crai="$WORK_DIR/crai/$sample.crai"
             crai_ok "$crai" || download "$cram.crai" "$crai" gzip || { echo "FAILED $sample (index download)" >&2; fail=1; continue; }
             args+=(--index "$crai") ;;
      *)     # a local file: name the index if it is a .crai beside the CRAM, otherwise htslib finds it (.bai, .csi)
             for crai in "$cram.crai" "${cram%.cram}.crai"; do [ -s "$crai" ] && { args+=(--index "$crai"); break; }; done ;;
    esac
  fi
  # transient failures are the norm at this scale: a refused connection fails at once, a dead one
  # makes the engine give up after --stall-timeout seconds (exit 75); either way, try again
  ok=0
  for attempt in 1 2 3; do
    if run count "${args[@]}" && [ -s "$out.tmp.gz" ]; then
      # a counts file from another engine than the cohort's: every sample of this task would be the same
      check_build "$out.tmp.gz" || { rm -f "$out.tmp.gz"; exit 1; }
      ok=1; break
    fi
    echo "attempt $attempt failed: $sample" >&2; rm -f "$out.tmp.gz" "$out.tmp.gz.partial"; sleep $(( attempt * 30 ))
  done
  if [ "$ok" = 1 ]; then mv "$out.tmp.gz" "$out"; else echo "FAILED $sample" >&2; fail=1; fi
done < <(sed -n "$(( task * SAMPLES_PER_TASK + 1 )),$(( (task + 1) * SAMPLES_PER_TASK ))p" "$MANIFEST")
left=0; while IFS=$'\t' read -r sample _; do [ -s "$COUNTS_DIR/$sample.json.gz" ] || left=$((left + 1)); done < "$MANIFEST"
echo "task $task done; $left of $N manifest samples without counts in $COUNTS_DIR yet"
exit $fail
