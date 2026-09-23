# NGS-DOSE on the 1000 Genomes 30× cohort

Results of running [NGS-DOSE](https://github.com/jlanej/NGS-DOSE) — rDNA copy number and other
multi-copy sequence from short-read genomes — on the expanded 1000 Genomes cohort (3,202 samples,
602 trios; NYGC 30× CRAMs, GRCh38). Published as the run proceeds, partial results included.

**The page: [jlanej.github.io/NGS-DOSE-1000G](https://jlanej.github.io/NGS-DOSE-1000G/)** —
what is measured and why, and the evidence that it works, recomputed from the files here.

![the evidence](docs/evidence.png)

The eight panels are the case that the method works, each drawn from `docs/report.json` and the
pilot's tables by `regenerate.sh`; the write-up — the number, what it rules out, and what is not
yet shown — is [docs/EVIDENCE.md in NGS-DOSE](https://github.com/jlanej/NGS-DOSE/blob/main/docs/EVIDENCE.md).

| path | what |
| --- | --- |
| `counts_scan/<sample>.json.gz` | whole-file scan of one CRAM: fragment-end counts per class and unit position, control-region counts, GC tables, where class reads were aligned and what else sits in those bins (~240 kB) |
| `counts_fetch/<sample>.json.gz` | the targeted fetch of the same CRAM: what a biobank-scale run would return (~70 kB) |
| `docs/` | the page (`index.html`), every number behind it (`report.json`) and every table (`data/*.tsv`); served by GitHub Pages |
| `docs/evidence.png` | the evidence figure, rendered from `report.json` |
| `meta/` | the pedigree, the published table and the NGS-PCA QC table the page compares against |
| `regenerate.sh` | rebuilds `docs/` from the counts |

The counts files are the primary data: no reads, no genotypes, nothing a 1000 Genomes
consent does not cover. Everything else is derived from them by `ngsdose report`, so anyone can
re-run the analysis — or a different one — without touching a CRAM.

## Regenerate

```bash
pip install git+https://github.com/jlanej/NGS-DOSE     # or the container image
bash regenerate.sh                                     # -> docs/
```

## Data sources and attribution

1000 Genomes 30× CRAMs: Byrska-Bishop et al., *Cell* 185:3426 (2022), AWS Open Data mirror;
1000 Genomes data are open access (https://www.internationalgenome.org/IGSR_disclaimer).
`meta/hall2021_MOESM1.txt`: Supplementary Data 1 of Hall, Turner & Queitsch, *Sci Rep* 11:449
(2021), doi:10.1038/s41598-020-80049-y, CC BY 4.0, unmodified.
`meta/ngspca_sample_qc.tsv`: per-sample QC of the same cohort from
[NGS-PCA-Manuscript](https://github.com/jlanej/NGS-PCA-Manuscript) (`1000G/qc_output/sample_qc.tsv`;
mosdepth coverage, mitochondrial copy number, X/Y coverage ratios, inferred sex, release batch), MIT, unmodified.
