#!/usr/bin/env bash
# Heterozygous sites in one genome's fetched reads, for allele_balance.py: SAMPLE CRAM REFERENCE.fa REGIONS.bed OUT_DIR
# (the CRAM: the small file of the containers that hold the bundle's single-copy regions, as slice_fetch.py writes it;
# REGIONS.bed: the bundle's controls.bed). Needs samtools and bcftools.
set -euo pipefail
S="$1"; CRAM="$2"; REF="$3"; BED="$4"; OUT="$5"
mkdir -p "$OUT"
[ -s "$OUT/$S.vcf.gz" ] && exit 0
[ -s "$CRAM.crai" ] || samtools index "$CRAM"
bcftools mpileup -f "$REF" -R <(cut -f1-3 "$BED") -a AD,DP -q 20 -Q 20 -d 500 -Ou "$CRAM" 2> "$OUT/$S.log" |
  bcftools call -mv -Oz -o "$OUT/$S.vcf.gz.part" 2>> "$OUT/$S.log"
mv "$OUT/$S.vcf.gz.part" "$OUT/$S.vcf.gz"
