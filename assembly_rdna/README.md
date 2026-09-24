# Do HPRC release-2 assemblies contain the rDNA, and is it reliably measured?

**[HPRC_r2_rDNA_report.pdf](HPRC_r2_rDNA_report.pdf)**: background on assembly methods and rDNA, direct sequence tests on all
466 rows of the HPRC r2 assembly index (460 HPRC/HPP haplotypes plus CHM13, T2T-HG002 v1.1, HG06807 and GRCh38), and a
comparison with NGS-DOSE short-read copy numbers of the same people and their parents.

In short: every hifiasm haplotype holds 45S units (median 119 18S genes). They come as dozens of short pieces ending at contig
ends, about half of the true copy number, correlating with it at r ~ 0.6. That is haplotype-specific signal, not a measurement.
The 5S array is assembled whole and matches NGS-DOSE's 5S estimate at r = 0.99 (ratio 1.02).

Which method is more accurate? For 45S, NGS-DOSE: the child's NGS-DOSE value tracks the parents' at r = 0.82 against 0.58 for
the assembly (36 trios; the same with parents measured by Hall et al. 2021), and it correlates with ddPCR (Potapova et al. 2025)
at r = 0.97 over 9 cell lines, ahead of CONKORD (0.89) and the 18S depth ratio (0.93). For 5S the two methods tie. NGS-DOSE's weakness
is absolute level, about 3% low against ddPCR and against the 10 copies of the distal junction.

| path | what |
| --- | --- |
| `background.md` | the cited background section in full |
| `scripts/` | `fetch_metadata.py`, `extract_regions.py` (remote faidx of rDNA, short-arm and 5S regions), `tile_map.py` and `validate_full.py` (completeness check on whole assemblies), `annotate_units.py` (18S/28S/unit/5S/DJ copies, arrays, array ends, unit periods), `analyze.py` (tests and figures), `method_accuracy.py` (inheritance, distal junction, ddPCR and FISH comparisons), `build_report.py` |
| `tables/` | `haplotype_rdna.tsv` (one row per assembly), `arrays_all.tsv.gz`, `person_vs_ngsdose.tsv`, `trio_haplotypes.tsv`, `tests.json`, `validation.tsv`, `assemblies.tsv`, `extraction.tsv`, `acro_contigs.tsv` |
| `figures/` | the five report figures |
| `data/` | index, per-assembly metadata, extracted sequence (7.5 GB; see `data/README.md`) |
| `potapova/` | Potapova et al. 2025 Tables S1 (ddPCR, CONKORD) and S2 (FISH per array) as TSV; NGS-DOSE fetch counts and estimates for the 7 of their 1000 Genomes lines not yet in the cohort run (`run_fetch.sh`) |
| `regenerate.sh` | reruns the analysis and PDF from the current `docs/data/cohort.tsv`; `FULL=1` also re-extracts and re-annotates |

Committed: scripts, tables, figures, the PDF, `data/index.csv`, `data/regions/` (what was extracted from each assembly),
`data/validation/*.summary.json` and the per-assembly annotations in `work_units/` (`*.summary.json`, `*.arrays.tsv`,
`*.copies.tsv.gz`), which is all `regenerate.sh` needs. Not committed (`.gitignore`): the extracted sequence and assembly
metadata (about 9 GB), minimap2 alignments, tile hits and logs, all rebuilt by `FULL=1 bash regenerate.sh`.
`tables/validation_outside_hits.tsv` is committed gzipped.
