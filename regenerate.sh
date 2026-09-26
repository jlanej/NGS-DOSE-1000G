#!/usr/bin/env bash
# Rebuild docs/ from the counts files with this repository's `report` package. Needs the `ngsdose` library
# (pip install git+https://github.com/jlanej/NGS-DOSE, or a checkout at ../NGS-DOSE installed with `pip install -e`)
# and its GRCh38 bundle: NGSDOSE_RESOURCES, or the checkout's resources/GRCh38 when the package is installed editable.
# The figure and the PDF need matplotlib and are skipped without it.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
# Downloads go to a temporary name and take their final one only once whole and checked, so an error page or a cut
# transfer is never kept (a file under its final name is not fetched again).
PED=meta/20130606_g1k_3202_samples_ped_population.txt
if [ ! -s "$PED" ]; then
  if curl -sSfL --retry 5 -o "$PED.part" \
       https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt &&
     head -1 "$PED.part" | grep -q SampleID && [ "$(wc -l < "$PED.part")" -gt 3202 ]; then
    mv "$PED.part" "$PED"
  else
    rm -f "$PED.part"; echo "ERROR: could not download the pedigree ($PED)" >&2; exit 1
  fi
fi
# HPRC release-2 CenSat annotations of the cohort's members (open access, HPRC S3 bucket): each fetched once, kept
# gzipped (-n: no name or time in the header, so the bytes are the same every time). A download that fails, or that is
# not a BED (every line other than track/browser/# with at least 4 tab-separated fields), warns and moves on.
KEYS=meta/hprc_r2_censat.keys.txt
if [ -s "$KEYS" ]; then
  mkdir -p hprc_censat
  awk 'NR == FNR { if (FNR > 1) m[$2] = 1; next } { split($0, p, "/"); if (p[3] in m) print }' "$PED" "$KEYS" |
  while read -r key; do
    f="hprc_censat/$(basename "$key").gz"
    [ -s "$f" ] && continue
    if curl -sSfL --retry 5 -o "$f.tmp" "https://s3-us-west-2.amazonaws.com/human-pangenomics/$key" && [ -s "$f.tmp" ] &&
       awk -F'\t' '/^(track|browser|#)/ || NF == 0 { next } NF < 4 { bad = 1; exit } END { exit bad }' "$f.tmp" &&
       gzip -9n < "$f.tmp" > "$f.part"; then
      mv "$f.part" "$f"
    else
      echo "HPRC annotation not fetched (download failed or not a BED): $key" >&2
    fi
    rm -f "$f.tmp" "$f.part"
  done
fi
args=(-p "$PED" --hall meta/hall2021_MOESM1.txt --pilot pilot --pcs meta/ngspca/svd.pcs.txt -o docs --cache cache -j "${JOBS:-4}")
[ -n "$(ls counts_scan/*.json.gz 2>/dev/null)" ] && args+=(--scan counts_scan)
[ -n "$(ls counts_fetch/*.json.gz 2>/dev/null)" ] && args+=(--fetch counts_fetch)
[ -d hprc_censat ] && args+=(--censat hprc_censat)
[ -s meta/ngspca_sample_qc.tsv ] && args+=(--qc meta/ngspca_sample_qc.tsv)
[ -s assembly_rdna/tables/potapova_comparison.tsv ] && args+=(--ddpcr assembly_rdna/tables/potapova_comparison.tsv)
python3 -m report "${args[@]}"
if python3 -c 'import matplotlib' 2>/dev/null; then
  python3 -m report.evidence_figure --report docs/report.json --pilot pilot -o docs/evidence.png
  python3 -m report.trio_report --report docs/report.json --data docs/data -o docs/trio_report.pdf
else
  echo "figure and PDF skipped: matplotlib is not installed" >&2
fi
