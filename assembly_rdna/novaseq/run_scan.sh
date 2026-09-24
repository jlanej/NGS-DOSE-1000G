#!/usr/bin/env bash
# Whole-file scan of Google NovaSeq BAMs: checks that the NYGC-learned sinks capture this pipeline's class reads
set -uo pipefail
cd "$(dirname "$0")"
ENGINE=${ENGINE:-/Users/Kitty/git/NGS-DOSE/target/release/ngs-dose}; B=${B:-/Users/Kitty/git/NGS-DOSE/resources/GRCh38}; REF=${REF:-/Users/Kitty/git/NGS-DOSE/work/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa}
for run in "$@"; do
  url=$(awk -F'\t' -v r=$run '$1==r{print $4}' manifest.tsv); s=${run%%.*}
  [ -s scan/$run.json.gz ] && continue
  "$ENGINE" count -m scan -@ 6 -i "$url" -T "$REF" -p $B/panel.k31.tsv.gz -c $B/controls.fa.gz -s "$s" -o scan/$run.tmp.json.gz > scan/$run.log 2>&1 && mv scan/$run.tmp.json.gz scan/$run.json.gz
  echo "$run $(date +%T) $([ -s scan/$run.json.gz ] && echo ok || echo FAIL)"
done
