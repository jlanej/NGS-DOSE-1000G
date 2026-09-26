#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_stage
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=96:00:00
#SBATCH --output=logs/stage_%j.out
# Stage the cohort a few CRAMs at a time and hand each one to 01b_dose_sample.sh (whole-file scan +
# targeted fetch, then the CRAM is removed) the moment it has landed and its MD5 has been checked.
# The WAN is the scarce, failure-prone resource, so it gets a few gentle multi-connection
# transfers; counting is ordinary compute and gets whatever the scheduler grants. At most
# MAX_LOCAL_CRAMS files are on disk at any time (~15 GB each).
#
#   sbatch 01_stage_and_dose.sh            # or, on a login / transfer node, in tmux:  bash 01_stage_and_dose.sh
#   LIMIT=5 bash 01_stage_and_dose.sh      # the first five samples only: a smoke test
#   RUN_LOCAL=1 ...                        # count here instead of submitting jobs
#
# Idempotent: samples with both counts files are skipped; a sample whose job is still queued or
# running is left alone, and so is any sample with a job while squeue does not answer. Each run first goes through the CRAMs already on disk - verified, not
# downloaded again, and submitted again if their last job ended without counts (at most MAX_TRIES
# jobs per sample, config.sh) - and then downloads the rest. A download goes to <file>.part and
# takes its name only once whole (the CRAM against the MD5 of the sequence index, the index by
# gzip -t), so an interrupted one is resumed. aria2c comes from the container image when SIF is set
# and from the host otherwise; the script falls back to curl (run on the host) only when that
# aria2c cannot be run.
#
# Back-pressure waits only before a download. A CRAM a job has kept after a failure frees no room,
# so the stager stops, naming those CRAMs, when MAX_KEPT_CRAMS of them are on disk, or when the disk
# is full and nothing on it has a job that could free it - each on five checks in a row, a minute
# apart (a job whose state squeue cannot tell counts as alive). Fix the cause (logs/dose_<job>.out)
# and run it again.
#
# One stager at a time. Stop a stager started from an older version of these scripts before starting
# this one: this one resumes the partial downloads the older one left under the CRAM's own name
# (X.cram beside X.cram.aria2) by renaming them to X.cram.part, and leaves alone any modified in the
# last 15 minutes, which may still be being written. Two stagers of this version refuse to run
# together (a lock on $WORK_DIR/stager.lock, where flock exists).
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
DOWNLOAD_SLOTS="${DOWNLOAD_SLOTS:-4}"                # concurrent transfers
DOWNLOAD_CONNECTIONS="${DOWNLOAD_CONNECTIONS:-8}"    # connections per transfer
MAX_LOCAL_CRAMS="${MAX_LOCAL_CRAMS:-40}"
MAX_KEPT_CRAMS="${MAX_KEPT_CRAMS:-10}"               # CRAMs kept by failed jobs before the stager stops
LIMIT="${LIMIT:-0}"
cd "$EX_DIR" || exit 1; mkdir -p logs "$CRAM_DIR" "$WORK_DIR/dispatched" "$WORK_DIR/counts_scan" "$WORK_DIR/counts_fetch"

done_sample() { [ -s "$WORK_DIR/counts_scan/$1.json.gz" ] && [ -s "$WORK_DIR/counts_fetch/$1.json.gz" ]; }
staged() { [ -s "$CRAM_DIR/$1.cram" ] || [ -e "$CRAM_DIR/$1.cram.part" ]; }      # whole or partly downloaded
n_local() { find "$CRAM_DIR" -maxdepth 1 \( -name '*.cram' -o -name '*.cram.part' \) | wc -l | tr -d ' '; }
kept_crams() {  # whole CRAMs whose latest job has ended without counts: "sample(job id)" per line
  local f s
  for f in "$CRAM_DIR"/*.cram; do
    [ -s "$f" ] && [ ! -e "$f.aria2" ] || continue
    s="$(basename "$f" .cram)"
    [ -s "$WORK_DIR/dispatched/$s" ] || continue
    done_sample "$s" || ! job_gone "$s" || echo "$s(job $(cat "$WORK_DIR/dispatched/$s"))"
  done
}
room_can_free() {  # a transfer of this run, or a job whose CRAM is on disk, is still going (or may be)
  local f s
  [ -n "$(jobs -rp)" ] && return 0
  for f in "$CRAM_DIR"/*.cram "$CRAM_DIR"/*.cram.part; do
    [ -e "$f" ] || continue
    s="$(basename "$f")"; s="${s%.part}"; s="${s%.cram}"
    job_gone "$s" || return 0
  done
  return 1
}
stop=""
STOP_CHECKS=5                                          # checks in a row, a minute apart, before a stop
wait_for_room() {  # back-pressure before a download: 01b_dose_sample.sh removes each CRAM once its counts are verified
  local waited=0 kept n_kept full=0 stuck=0
  while [ "$(n_local)" -ge "$MAX_LOCAL_CRAMS" ]; do
    kept="$(kept_crams)"; n_kept="$(printf '%s' "$kept" | grep -c .)"
    if [ "$n_kept" -ge "$MAX_KEPT_CRAMS" ]; then full=$((full + 1)); else full=0; fi
    if room_can_free; then stuck=0; else stuck=$((stuck + 1)); fi
    if [ "$full" -ge "$STOP_CHECKS" ]; then
      stop="$n_kept CRAMs kept on disk by jobs that ended without counts (MAX_KEPT_CRAMS=$MAX_KEPT_CRAMS): $(echo $kept)"
      return 1
    fi
    if [ "$stuck" -ge "$STOP_CHECKS" ]; then
      stop="$(n_local) CRAMs on disk (MAX_LOCAL_CRAMS=$MAX_LOCAL_CRAMS) and none has a job that would remove it${kept:+; kept by jobs that ended without counts: $(echo $kept)}"
      return 1
    fi
    [ $((waited % 10)) = 0 ] && echo "waiting for room: $(n_local) CRAMs on disk (MAX_LOCAL_CRAMS=$MAX_LOCAL_CRAMS)${kept:+; kept by jobs that ended without counts: $(echo $kept)}"
    sleep 60; waited=$((waited + 1))
  done
}

stage_one() {  # sample url md5 line
  local sample="$1" url="$2" md5="$3" line="$4" cram="$CRAM_DIR/$1.cram"
  # a partial download from earlier versions of this script (under its final name, beside aria2c's
  # control file): resumed under the .part name
  if [ -e "$cram.aria2" ]; then mv -f "$cram" "$cram.part" 2>/dev/null; mv -f "$cram.aria2" "$cram.part.aria2"; fi
  if [ -s "$cram" ] && [ "$(md5_of "$cram")" = "$md5" ]; then
    echo "[$sample] already on disk, MD5 verified"
  else
    rm -f "$cram"
    download "$url" "$cram" "$md5" || { echo "[$sample] FAILED to download" >&2; return 1; }
  fi
  # a cut index loses every contig past the cut without an error, so it has to be whole, not just there
  if ! crai_ok "$cram.crai"; then
    rm -f "$cram.crai" "$cram.crai.aria2"
    download "$url.crai" "$cram.crai" gzip || { echo "[$sample] FAILED to download the index" >&2; return 1; }
  fi
  if [ "${RUN_LOCAL:-0}" = 1 ]; then
    bash "$EX_DIR/01b_dose_sample.sh" "$sample" "$line"            # MD5 checked above
  else
    dispatch "$sample" "$line" && echo "[$sample] staged, job $(cat "$WORK_DIR/dispatched/$sample")"
  fi
}

hold_lock 8 "$WORK_DIR/stager.lock" || { echo "ERROR: another stager is running on $WORK_DIR (it holds $WORK_DIR/stager.lock); stop it first" >&2; exit 1; }
check_fetch_panels || exit 1                  # every job would fail (or, with an old engine, fetch the wrong thing)
note_build
have_aria2                                    # decided once, before the transfers run in the background
started=0
# pass 1: what is already on disk, whole or partial (no room needed); pass 2: downloads
for pass in 1 2; do
  line=0
  while IFS=$'\t' read -r sample url md5; do
    line=$((line + 1))
    done_sample "$sample" && continue
    case "$url" in *://*) ;; *) [ "$pass" = 1 ] && echo "[$sample] the manifest gives a local path, nothing to stage" >&2; continue ;; esac
    if staged "$sample"; then [ "$pass" = 1 ] || continue; else [ "$pass" = 2 ] || continue; fi
    # a partial download of an older stager under the CRAM's own name, still being written perhaps
    if [ "$pass" = 1 ] && { [ -e "$CRAM_DIR/$sample.cram.aria2" ] || [ -e "$CRAM_DIR/$sample.cram.crai.aria2" ]; } \
       && [ -n "$(find "$CRAM_DIR" -maxdepth 1 -name "$sample.cram*" -mmin -15 2>/dev/null)" ]; then
      echo "[$sample] a partial download modified in the last 15 minutes (another stager's?): left alone this run" >&2
      continue
    fi
    job_gone "$sample" || { [ $? = 2 ] && echo "[$sample] squeue did not answer: its job $(cat "$WORK_DIR/dispatched/$sample") is taken as alive" >&2; continue; }
    if [ "${RUN_LOCAL:-0}" != 1 ] && [ "$(tries_of "$sample")" -ge "$MAX_TRIES" ]; then
      echo "[$sample] $(tries_of "$sample") jobs ended without counts (last: logs/dose_$(cat "$WORK_DIR/dispatched/$sample" 2>/dev/null).out); not submitted again - remove $WORK_DIR/dispatched/$sample.tries to retry" >&2
      continue
    fi
    [ "$LIMIT" -gt 0 ] && [ "$started" -ge "$LIMIT" ] && break 2
    [ "$pass" = 1 ] || wait_for_room || break 2
    while [ "$(jobs -rp | wc -l)" -ge "$DOWNLOAD_SLOTS" ]; do sleep 5; done
    stage_one "$sample" "$url" "$md5" "$line" &
    started=$((started + 1))
  done < "$MANIFEST"
done
wait
left=0; while IFS=$'\t' read -r sample _; do done_sample "$sample" || left=$((left + 1)); done < "$MANIFEST"
echo "staging pass finished: $started sample(s) started; $left of $(( $(wc -l < "$MANIFEST") )) without counts yet (queued jobs included). Run again to continue."
kept="$(kept_crams)"
[ -z "$kept" ] || echo "CRAMs kept on disk by jobs that ended without counts (see logs/dose_<job>.out): $(echo $kept)"
[ -z "$stop" ] || { echo "STOPPED: $stop. Fix the cause and run again." >&2; exit 1; }
