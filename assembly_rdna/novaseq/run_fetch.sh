#!/usr/bin/env bash
# NGS-DOSE fetch-mode counts on a second NovaSeq pipeline (Google Health / GIAB, Baid et al. 2020), plus the
# NYGC 1000 Genomes CRAMs of the three people sequenced by both (the CEPH trio), with the cohort engine and bundle.
set -uo pipefail
cd "$(dirname "$0")"
ENGINE=${ENGINE:-/Users/Kitty/git/NGS-DOSE/target/release/ngs-dose}; B=${B:-/Users/Kitty/git/NGS-DOSE/resources/GRCh38}; REF=${REF:-/Users/Kitty/git/NGS-DOSE/work/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa}
tail -n +2 manifest.tsv | while IFS=$'\t' read -r run s plat url; do
  [ -s counts/$run.json.gz ] && continue
  for try in 1 2 3; do
    "$ENGINE" count -m fetch -@ 8 -i "$url" -T "$REF" -p $B/panel.k31.tsv.gz -c $B/controls.fa.gz --sinks $B/sinks.bed -s "$s" -o counts/$run.tmp.json.gz > counts/$run.log 2>&1 && mv counts/$run.tmp.json.gz counts/$run.json.gz && break
  done
  echo "$run $(date +%T) $([ -s counts/$run.json.gz ] && echo ok || echo FAIL)"
done
