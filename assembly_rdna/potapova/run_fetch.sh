#!/usr/bin/env bash
# NGS-DOSE fetch-mode counts for the 1000 Genomes samples of Potapova et al. 2025 (Table S1) not yet in the cohort run
set -uo pipefail
cd "$(dirname "$0")"
ENGINE=${ENGINE:-/Users/Kitty/git/NGS-DOSE/target/release/ngs-dose}; B=${B:-/Users/Kitty/git/NGS-DOSE/resources/GRCh38}; REF=${REF:-/Users/Kitty/git/NGS-DOSE/work/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa}
while IFS=$'\t' read -r s url md5; do
  [ -s counts/$s.json.gz ] && continue
  for try in 1 2 3; do
    "$ENGINE" count -m fetch -@ 8 -i "$url" -T "$REF" -p $B/panel.k31.tsv.gz -c $B/controls.fa.gz --sinks $B/sinks.bed -s "$s" -o counts/$s.tmp.json.gz > counts/$s.log 2>&1 && mv counts/$s.tmp.json.gz counts/$s.json.gz && break
  done
done < manifest.tsv
