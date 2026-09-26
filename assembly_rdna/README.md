# Do HPRC release-2 assemblies contain the rDNA, and is it reliably measured?

**[HPRC_r2_rDNA_report.pdf](HPRC_r2_rDNA_report.pdf)**: background on assembly methods and rDNA, direct sequence tests on all
466 rows of the HPRC r2 assembly index (460 HPRC/HPP haplotypes plus CHM13, T2T-HG002 v1.1, HG06807 and GRCh38), and a
comparison with NGS-DOSE short-read copy numbers of the same people and their parents.

In short: every hifiasm haplotype holds 45S units (median 119 18S genes). They come as dozens of short pieces ending at contig
ends, about half of the true copy number, correlating with it at r ~ 0.55-0.6. That is haplotype-specific signal, not a
measurement. The 5S array is assembled whole and matches NGS-DOSE's 5S estimate at r = 0.99 (ratio 1.02; 57 people).

Which method is more accurate? For 45S, NGS-DOSE: the child's NGS-DOSE value tracks the parents' at r = 0.82 against 0.58 for
the assembly (36 trios; the same with parents measured by Hall et al. 2021), and it correlates with ddPCR (Potapova et al. 2025)
at r = 0.97 over 9 cell lines, ahead of CONKORD (0.89 over the same 9) and the 18S depth ratio (0.93). For 5S the two methods tie.
NGS-DOSE's weakness is absolute level, about 3% low against ddPCR and against the 10 copies of the distal junction.

These numbers, like the committed tables, figures and PDF, come from the cohort table of 2026-09-24 with 750 genomes (57 of them
people with an assembly, 36 trios). Rerun on the 1,259-genome `docs/data/cohort.tsv` (102 people, 72 trios), the conclusions hold:
the child tracks the parents at r = 0.78 with NGS-DOSE against 0.57 with the assembly (0.78 against 0.55 over 62 trios with
parents measured by Hall et al.), the assembly's 45S correlates with NGS-DOSE at r = 0.54, the 5S at r = 0.97 (ratio 1.02), and
the ddPCR comparison is unchanged (r = 0.97 over the 9 lines). `regenerate.sh` rebuilds everything from the current table; quote its numbers with the
cohort size they came from.

| path | what |
| --- | --- |
| `background.md` | the cited background section in full |
| `scripts/` | `fetch_metadata.py`, `extract_regions.py` (remote faidx of rDNA, short-arm and 5S regions), `tile_map.py` and `validate_full.py` (completeness check on whole assemblies), `annotate_units.py` (18S/28S/unit/5S/DJ copies, arrays, array ends, unit periods), `analyze.py` (tests, figures 1-3), `method_accuracy.py` (inheritance, distal junction, ddPCR, FISH and second-pipeline comparisons; figures 4-7 and the scorecard), `build_report.py`, `common.py` (shared helpers) |
| `tables/` | `haplotype_rdna.tsv` (one row per assembly), `arrays_all.tsv.gz`, `person_vs_ngsdose.tsv`, `trio_haplotypes.tsv`, `tests.json`, `validation.tsv`, `assemblies.tsv`, `extraction.tsv`, `acro_contigs.tsv`; from `method_accuracy.py`: `method_accuracy.json` and `scorecard.tsv` (the accuracy tests and the overview table), `potapova_comparison.tsv` (also read by the repository's report, `--ddpcr`), `potapova_arrays.tsv`, `potapova_fish_totals.tsv`, `novaseq_bridge.tsv`, `novaseq_estimates.tsv` |
| `figures/` | the seven report figures (fig1-3 from `analyze.py`, fig4-7 from `method_accuracy.py`; fig6 = second NovaSeq pipeline, fig7 = overview) |
| `data/` | index, per-assembly metadata, extracted sequence (7.9 GB; see `data/README.md`) |
| `potapova/` | Potapova et al. 2025 Tables S1 (ddPCR, CONKORD) and S2 (FISH per array) as TSV; NGS-DOSE fetch counts for the 7 of their 1000 Genomes lines that were not yet counted on 2026-09-24 (`run_fetch.sh`) and estimates for all 9 of their 1000 Genomes lines. All 9 are now in the cohort run: `method_accuracy.py` uses a line's cohort value whenever `docs/data/cohort.tsv` has it, and falls back to `potapova/est` otherwise |
| `novaseq/` | second NovaSeq pipeline: `manifest.tsv`, fetch counts (`run_fetch.sh`) and estimates calibrated with the NYGC cohort's window efficiencies for Google Health's GIAB PCR-free BAMs (HG002-4, CEPH trio NA12878/NA12891/NA12892) and for the NYGC 1000 Genomes CRAMs of the same CEPH trio; one whole-file scan of HG002 (`run_scan.sh`) with its sink-capture check `scan_check.json` |
| `regenerate.sh` | reruns the analysis and PDF from the current `docs/data/cohort.tsv`; `FULL=1` also re-extracts and re-annotates |

Committed: scripts, tables, figures, the PDF, `data/index.csv`, `data/regions/` (what was extracted from each assembly),
`data/validation/*.summary.json` and the per-assembly annotations in `work_units/` (`*.summary.json`, `*.arrays.tsv`,
`*.copies.tsv.gz`), which, together with the rest of `potapova/` (Tables S1 and S2; `est/` if present), `novaseq/est/*/cohort.tsv`
(and `novaseq/scan_check.json` if present) and the repository's `docs/data/cohort.tsv`, `docs/report.json` and `meta/` (pedigree,
Hall 2021), is all `regenerate.sh` needs. Not committed (`.gitignore`): the extracted sequence and assembly metadata (about
9 GB), the minimap2 alignments and the extraction tile screens (`data/regions/*.hits.tsv.gz`), which `FULL=1 bash regenerate.sh`
rebuilds. That run also needs KY962518.1, X12811.1, the CHM13 DJ core and the CHM13 CenSat BED at the local paths set in
`scripts/tile_map.py`, `scripts/annotate_units.py` and `scripts/common.py`. Also not committed: the validation tile hits
(`data/validation/*.hits*.tsv.gz`), which only `python scripts/validate_full.py` rebuilds, and `work_tiles/` (a one-off tile run
over each assembly's extracted sequence) and `logs/`, which no committed script writes.
`tables/validation_outside_hits.tsv` is committed gzipped.
