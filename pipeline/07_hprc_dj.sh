#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_hprc_dj
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/hprc_dj_%j.out
# The distal junction in the HPRC release-2 assemblies of cohort members: each haplotype FASTA (~900 MB)
# is fetched, screened for the bundle's DJ panel k-mers (ngs-dose panel --report), aligned to the
# core-masked DJ unit (minimap2), and removed. SAMPLES="HG00097 HG01891" restricts the run; otherwise
# every cohort member with a CenSat annotation (meta/hprc_r2_censat.keys.txt) is screened. Then
#   python3 pipeline/hprc_dj.py screen $WORK_DIR/hprc_dj -o meta/dj_hprc
# writes the two tables the report reads (--dj-assemblies meta/dj_hprc).
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
OUT="$WORK_DIR/hprc_dj"; mkdir -p "$OUT"; cd "$OUT"
command -v minimap2 >/dev/null || ngsdose_py minimap2 --version >/dev/null 2>&1 || { echo "minimap2 is needed (host or image)" >&2; exit 1; }
mm2() { if command -v minimap2 >/dev/null; then minimap2 "$@"; else ngsdose_py minimap2 "$@"; fi; }

[ -s DJ.core_masked.fa ] || ngsdose_py python3 "$EX_DIR/hprc_dj.py" masked-unit -o DJ.core_masked.fa -r "$BUNDLE"
printf 'DJ\tpositional\t%s\t0\n' "$BUNDLE/units/DJ.CHM13v2_chr21_2708299_3108298.fa" > dj.manifest.tsv

# sample and haplotype names, from the CenSat keys (…/<sample>_<hap>_hprc_r2_v1.cenSat.bed)
sed -E 's#.*/([A-Z]+[0-9]+)_([a-z0-9]+)_hprc_r2.*#\1 \2#' "$REPO/meta/hprc_r2_censat.keys.txt" | sort -u |
awk -v keep="${SAMPLES:-}" 'BEGIN { n = split(keep, k, " "); for (i = 1; i <= n; i++) want[k[i]] = 1 } n == 0 || ($1 in want)' |
while read -r sample hap; do
  tag="${sample}_${hap}"
  [ -s "$tag.report.tsv" ] && [ -s "$tag.masked.paf" ] && continue
  url="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/$sample/assemblies/release2/${sample}_${hap}_hprc_r2_v1.0.1.fa.gz"
  want_bytes="$(curl -sI "$url" | tr -d '\r' | awk 'tolower($1) == "content-length:" { print $2 }')"
  if ! curl -sfL --retry 5 -o "$tag.fa.gz" "$url" || [ -z "$want_bytes" ] || [ "$(stat -c %s "$tag.fa.gz" 2>/dev/null || stat -f %z "$tag.fa.gz")" != "$want_bytes" ]; then
    echo "[$tag] download failed or incomplete" >&2; rm -f "$tag.fa.gz"; continue
  fi
  ngsdose_engine panel -m dj.manifest.tsv -k 31 -b "$tag.fa.gz" --max-bg 100000 --report "$tag.report.tsv" -o "$tag.panel.tsv.gz" > /dev/null
  mm2 -x asm20 -c --secondary=no -t "${SLURM_CPUS_PER_TASK:-4}" DJ.core_masked.fa "$tag.fa.gz" > "$tag.masked.paf" 2> "$tag.mm2.log"
  rm -f "$tag.fa.gz" "$tag.panel.tsv.gz"
  echo "[$tag] screened"
done
echo "screens in $OUT; now: python3 $EX_DIR/hprc_dj.py screen $OUT -o $REPO/meta/dj_hprc"
