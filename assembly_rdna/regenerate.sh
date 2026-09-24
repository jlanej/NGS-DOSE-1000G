#!/usr/bin/env bash
# Rebuild the HPRC r2 assembly-rDNA analysis and report. Run from anywhere.
#   bash assembly_rdna/regenerate.sh            # analysis + PDF from the current docs/data/cohort.tsv (seconds)
#   FULL=1 bash assembly_rdna/regenerate.sh     # also re-extract (~10 GB over HTTPS) and re-annotate all 466 assemblies (~1 h)
# Needs samtools, minimap2, pysam, pandas, scipy, matplotlib, reportlab. After new NGS-DOSE counts land, run the
# repository's regenerate.sh first (it rewrites docs/data/cohort.tsv), then this.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [ -n "${FULL:-}" ]; then
  python scripts/fetch_metadata.py
  python scripts/extract_regions.py
  ls data/seq/*.fa.gz | xargs -n1 basename | sed 's/\.fa\.gz$//' > work_units/names.txt
  rm -f work_units/*.paf.gz work_units/*.summary.json work_units/*.arrays.tsv work_units/*.copies.tsv.gz
  xargs -P "${JOBS:-5}" -n 4 python scripts/annotate_units.py --threads 2 < work_units/names.txt > work_units/all.jsonl
fi
python scripts/analyze.py
python scripts/method_accuracy.py
python scripts/build_report.py
