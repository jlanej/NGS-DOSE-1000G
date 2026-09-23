#!/usr/bin/env bash
# Rebuild docs/ from the counts files. Needs `ngsdose` (pip install git+https://github.com/jlanej/NGS-DOSE)
# and, for the coverage-PC section, a checkout of NGS-DOSE (its example carries NGS-PCA's PCs for this cohort).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NGSDOSE="${NGSDOSE:-../NGS-DOSE}"
[ -s meta/20130606_g1k_3202_samples_ped_population.txt ] || curl -sSL -o meta/20130606_g1k_3202_samples_ped_population.txt \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
[ -s meta/hall2021_MOESM1.txt ] || cp "$NGSDOSE/example/1000G/pilot/hall2021_MOESM1.txt" meta/
# HPRC release-2 CenSat annotations of the cohort's members (open access, HPRC S3 bucket): each fetched once, kept
# gzipped (-n: no name or time in the header, so the bytes are the same every time); a failed fetch warns and moves on
KEYS="$NGSDOSE/example/1000G/hprc_r2_censat.keys.txt"
if [ -s "$KEYS" ]; then
  mkdir -p hprc_censat
  awk 'NR == FNR { if (FNR > 1) m[$2] = 1; next } { split($0, p, "/"); if (p[3] in m) print }' meta/20130606_g1k_3202_samples_ped_population.txt "$KEYS" |
  while read -r key; do
    f="hprc_censat/$(basename "$key").gz"
    [ -s "$f" ] && continue
    if curl -sSfL --retry 5 -o "$f.tmp" "https://s3-us-west-2.amazonaws.com/human-pangenomics/$key"; then gzip -9n < "$f.tmp" > "$f"; else echo "HPRC annotation not fetched: $key" >&2; fi
    rm -f "$f.tmp"
  done
fi
args=(-p meta/20130606_g1k_3202_samples_ped_population.txt --hall meta/hall2021_MOESM1.txt -o docs --cache cache -j "${JOBS:-4}")
[ -n "$(ls counts_scan/*.json.gz 2>/dev/null)" ] && args+=(--scan counts_scan)
[ -n "$(ls counts_fetch/*.json.gz 2>/dev/null)" ] && args+=(--fetch counts_fetch)
[ -s "$NGSDOSE/example/1000G/ngspca/svd.pcs.txt" ] && args+=(--pcs "$NGSDOSE/example/1000G/ngspca/svd.pcs.txt")
[ -s "$NGSDOSE/example/1000G/pilot/pilot_heldout.tsv" ] && args+=(--pilot "$NGSDOSE/example/1000G/pilot")
[ -d hprc_censat ] && args+=(--censat hprc_censat)
[ -s meta/ngspca_sample_qc.tsv ] && args+=(--qc meta/ngspca_sample_qc.tsv)
ngsdose report "${args[@]}"
# the evidence figure: eight panels from report.json and the pilot tables (needs matplotlib and the checkout)
if [ -s "$NGSDOSE/example/1000G/evidence_figure.py" ] && python3 -c 'import matplotlib' 2>/dev/null; then
  python3 "$NGSDOSE/example/1000G/evidence_figure.py" --report docs/report.json --pilot "$NGSDOSE/example/1000G/pilot" -o docs/evidence.png
  python3 "$NGSDOSE/example/1000G/trio_report.py" --report docs/report.json --data docs/data -o docs/trio_report.pdf
else
  echo "evidence figure skipped: needs $NGSDOSE and matplotlib" >&2
fi
