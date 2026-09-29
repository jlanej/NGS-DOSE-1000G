#!/usr/bin/env bash
# Dispatch 01b_dose_sample.sh for every manifest sample whose CRAM and index have landed in
# $CRAM_DIR and whose counts are not there yet. It never downloads and never deletes, so it can
# run beside whatever is staging the files (NGS-PCA's aria2 download manager, Globus, rsync);
# a file that aria2 is still writing (it keeps a .aria2 control file beside it), or an index that is
# not yet whole (gzip -t), is left alone. A sample whose job is queued or running is left alone too, and
# so is every sample with a job while squeue does not answer; one whose job has ended without counts is
# submitted again, up to MAX_TRIES jobs (config.sh).
# Run it again, or from cron, as files arrive; it exits non-zero when sbatch refuses a job (the rest
# wait for the next run).   RUN_LOCAL=1 runs the jobs here instead of sbatch.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
cd "$EX_DIR"; mkdir -p logs "$WORK_DIR/dispatched"
check_fetch_panels || exit 1
note_build
n=0; line=0; failed=0; held=0; unknown=0
while IFS=$'\t' read -r sample _ md5; do
  line=$((line + 1))
  cram="$CRAM_DIR/$sample.cram"
  [ -s "$cram" ] && [ ! -e "$cram.aria2" ] && crai_ok "$cram.crai" || continue
  [ -s "$WORK_DIR/counts_scan/$sample.json.gz" ] && [ -s "$WORK_DIR/counts_fetch/$sample.json.gz" ] && continue
  if [ "${RUN_LOCAL:-0}" = 1 ]; then
    bash "$EX_DIR/01b_dose_sample.sh" "$sample" "$line" "${md5:-}"
  else
    job_gone "$sample" || { [ $? = 2 ] && unknown=$((unknown + 1)); continue; }     # queued, running, or unknown
    if [ "$(tries_of "$sample")" -ge "$MAX_TRIES" ]; then
      echo "[$sample] $(tries_of "$sample") jobs ended without counts (last: logs/dose_$(cat "$WORK_DIR/dispatched/$sample" 2>/dev/null).out); not submitted again - remove $WORK_DIR/dispatched/$sample.tries to retry" >&2
      held=$((held + 1)); continue
    fi
    dispatch "$sample" "$line" "${md5:-}" || { failed=1; break; }     # a refusal (a submit limit) would refuse the rest too
  fi
  n=$((n + 1))
done < "$MANIFEST"
msg=""; [ "$held" = 0 ] || msg="; $held held back after $MAX_TRIES jobs without counts"
[ "$unknown" = 0 ] || msg="$msg; $unknown with a job left for the next run, as squeue did not answer"
echo "dispatched $n sample(s)$msg"
[ "$failed" = 0 ] || { echo "ERROR: sbatch refused a job; the remaining samples are submitted by the next run" >&2; exit 1; }
