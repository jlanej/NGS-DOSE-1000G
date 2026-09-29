# The direct tests behind docs/DJ.md

Run from the repository root with the `ngsdose` library installed and `docs/` regenerated (the scripts read the
published tables and the cached estimates). `selection.tsv` is the 600-genome set the tests were first run on
(193 whole-copy carriers, 257 between steps, 150 at the mode, by `DJ.step` on the fetch-mode estimates);
`HG00097_hap1.report.tsv.gz` is the `ngs-dose panel --report` screen of that haplotype (every k-mer of the unit and
its count in the assembly), the one per-k-mer screen kept here.

| script | question | answer on this cohort |
| --- | --- | --- |
| `fetch_capture.py` | does the targeted fetch miss junction reads? | the fetch holds 99.76% of the scan's DJ reads (range 99.64–99.83%, 600 genomes) |
| `presence_vs_gc.py` | is the GC-dependent deficit divergence among the copies? | exact-k-mer presence in HG00097's first haplotype is 0.95 per window and flat in GC (r = +0.08), while the raw estimate of 150 mode genomes falls from 10.08 to 9.26 between the lowest and highest GC bins (r = −0.26); not divergence |
| `acro_dosage.py` | are the values between steps mosaic aneuploidies of an acrocentric? | per-chromosome dosage SD 0.02–0.06 copies, r with the step −0.06; four of 600 genomes off by 0.2–0.9 copies, one (HG00142, chr13 +0.5) with a matching junction step |
| `trio_integers.py` | are integer states Mendelian, and transmitted? | whole unit: 468 of 468 fully called trios consistent; core: 477 of 477; single-carrier parents transmit 21 of 60 (fathers' losses 5 of 23, fathers' gains 5 of 14, mothers' losses 8 of 16, mothers' gains 3 of 7) |
| `between_profiles.py` | what do between-step profiles look like? | of the 257 selected: 221 partial (a run of blocks at another integer), 27 flat shifts, 9 noisy; of the 254 in the cohort more than 0.3 from an integer on the core: 210, 29, 15 |

The assembly screens themselves (`pipeline/07_hprc_dj.sh`, `pipeline/hprc_dj.py`) and the comparison the page
recomputes (`report/dj.py`) are described in `docs/DJ.md`.
