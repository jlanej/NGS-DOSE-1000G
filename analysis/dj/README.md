# The direct tests behind docs/DJ.md

Run from the repository root with the `ngsdose` library installed (0.2.0 or later) and `docs/` regenerated: the
scripts read the published tables, and those that need window-level profiles run the cohort layer on the cached
estimates (`cache/scan`, a few minutes; `mosaic_power.py` twenty). Settled calls are those on their whole numbers:
a fractional or an uncertain call is counted in none of the families' tallies. `selection.tsv` is the 600-genome set the first tests were run on (193
whole-copy carriers, 257 between steps, 150 at the mode, by the level before the calls existed).
`HG00097_hap1.report.tsv.gz` is the `ngs-dose panel --report` screen of that haplotype (every k-mer of the unit and
its count in the assembly), the one per-k-mer screen kept here.

| script | question | answer on this cohort |
| --- | --- | --- |
| `fetch_capture.py` | does the targeted fetch miss junction reads? | the fetch holds 99.76% of the scan's DJ reads (range 99.64–99.84%, 3,202 genomes) |
| `presence_vs_gc.py` | is the GC-dependent deficit divergence among the copies? | exact-k-mer presence in HG00097's first haplotype is 0.95 per window and flat in GC (r = +0.08), while the raw estimate of 150 mode genomes falls from 10.09 below 35% GC to 9.28 at 50–60% GC (r = −0.26); not divergence |
| `acro_dosage.py` | are the values between steps mosaic aneuploidies of an acrocentric? | per-chromosome dosage SD 0.02–0.05 copies, r with the step −0.08; ten of 600 genomes off by more than 0.15 copies, one (HG00142, chr13 +0.5) with a junction step to match |
| `polymorphic_offsets.py` | what factor do the polymorphic intervals' efficiencies need, by the assemblies? | +0.084 (5–15 kb), +0.034 (15–23 kb) and +0.062 (195–217 kb) in the logarithm, over 18 resolved genomes; the cohort's own comb gives +0.073, +0.028 and +0.045; a core stretch (300–340 kb) gives +0.004 |
| `trio_integers.py` | are whole numbers Mendelian, by the level and by the calls? | no trio is Mendelian-impossible by either; 452 of 602 trios have every member within 0.3 of a whole number on the whole unit, 458 on the core, and 556 have three settled calls |
| `between_profiles.py` | what are the genomes whose level lies between whole numbers? | of 309 more than 0.3 from a whole number on the core, the calls make 112 a genome with a copy that holds or lacks an end, 108 a ten-copy genome with a scale a few percent off, 31 fractional, 26 a ten-copy genome leaning more than 5%, 19 uncertain, 7 another whole number throughout, 6 another event of 40 kb or more |
| `partial_transmission.py` | did the calls miss the copies of children counted as lacking their parent's? | no: across the parent's breakpoint the children whose call lacks it step by +0.02 copies (median; quartiles −0.14 to +0.13), those whose call has it by +0.96; by the profile 24 of 55 pairs are passed on, by the calls 23 |
| `carriers_by_sex.py` | is the fathers' deficit a property of reading a man's genome? | no: 57 of 1,565 men and 49 of 1,559 women are at nine throughout, and the genomes at ten read 9.990 and 9.995; a father's loss passes to 2 of 8 sons and 2 of 14 daughters; nine throughout in 4.0% of fathers, 2.8% of mothers, 2.4% of children and 3.8% of the others |
| `intervals_by_batch.py` | the children read the distal intervals lower than their parents: batch or generation? | in 15–23 kb the later release batch reads 0.15 ± 0.03 copies lower than the earlier, and the 103 parents sequenced in it 0.14 ± 0.07 lower than the other parents: the batch; in 5–15 kb the batches differ by 0.11 ± 0.05 and the parents' two groups do not (0.00 ± 0.13), which settles nothing; the core stretches differ by 0.02 or less |

| `mosaic_power.py` | what do the calls make of a change in part of the cells? | put into forty genomes at ten copies and judged with the cohort: a junction lost in 0.4, 0.5 and 0.6 of the cells is fractional or uncertain in 60%, 88% and 58% of them, and from 0.6 it is called nine in nine tenths; a copy of the first 316 kb in 50%, 80% and 52%; 240–320 kb lost in 65%, 80% and 65%; the first 110 kb lost in 20%, 33% and 35%; at 0.3 of the cells 3% to 28% are flagged, at 0.1 none; the level reads the change at every share (0.31 below ten for a junction lost in 0.3 of the cells) |

The assembly screens themselves (`pipeline/07_hprc_dj.sh`, `pipeline/hprc_dj.py`) and everything the page
recomputes (`report/dj.py`, `report/report.py`) are described in `docs/DJ.md`.
