#!/usr/bin/env bash
# Rebuild docs/ from the counts files. Needs `ngsdose` (pip install git+https://github.com/jlanej/NGS-DOSE)
# and, for the coverage-PC section, a checkout of NGS-DOSE (its example carries NGS-PCA's PCs for this cohort).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NGSDOSE="${NGSDOSE:-../NGS-DOSE}"
[ -s meta/20130606_g1k_3202_samples_ped_population.txt ] || curl -sSL -o meta/20130606_g1k_3202_samples_ped_population.txt \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
[ -s meta/hall2021_MOESM1.txt ] || cp "$NGSDOSE/example/1000G/pilot/hall2021_MOESM1.txt" meta/
args=(-p meta/20130606_g1k_3202_samples_ped_population.txt --hall meta/hall2021_MOESM1.txt -o docs --cache cache -j "${JOBS:-4}")
[ -n "$(ls counts_scan/*.json.gz 2>/dev/null)" ] && args+=(--scan counts_scan)
[ -n "$(ls counts_fetch/*.json.gz 2>/dev/null)" ] && args+=(--fetch counts_fetch)
[ -s "$NGSDOSE/example/1000G/ngspca/svd.pcs.txt" ] && args+=(--pcs "$NGSDOSE/example/1000G/ngspca/svd.pcs.txt")
[ -d hprc_censat ] && args+=(--censat hprc_censat)
ngsdose report "${args[@]}"
# the evidence figure: eight panels from report.json and the pilot tables (needs matplotlib and the checkout)
if [ -s "$NGSDOSE/example/1000G/evidence_figure.py" ] && python3 -c 'import matplotlib' 2>/dev/null; then
  python3 "$NGSDOSE/example/1000G/evidence_figure.py" --report docs/report.json --pilot "$NGSDOSE/example/1000G/pilot" -o docs/evidence.png
else
  echo "evidence figure skipped: needs $NGSDOSE and matplotlib" >&2
fi
