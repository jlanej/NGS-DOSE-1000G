# The case that it works

What would convince someone who doubts that rDNA copy number can be measured from short-read
whole-genome sequencing? Not a model, and not agreement with ourselves. This page collects the
evidence that exists so far, in the order a sceptic would ask for it, with the number, what it
rules out, and where it comes from. It is a dated snapshot: its numbers were computed by scripts
from the data committed on 2026-09-23, the [cohort page](https://jlanej.github.io/NGS-DOSE-1000G/)
(`python -m report`, 735 of 3,202 genomes and 149 of 602 trios, commit 4a276d4) and the
[pilot report](../pilot/pilot_report.md) (`evaluate_pilot.py`, 12 genomes each sequenced twice),
and they are not updated as the cohort grows. Where a later review found a number of that date
wrong, it is corrected here on the same data and the correction says so.

**For the current numbers, see the [cohort page](https://jlanej.github.io/NGS-DOSE-1000G/) and
its figure, [evidence.png](evidence.png), which `regenerate.sh` rewrites at every run.**

![the evidence, as of 2026-09-23](evidence_2026-09-23.png)

The figure is the one drawn on that date. Two of its labels are out of date: panel b counts
transmitted steps by the rule of that date (7 of 22; section 2 gives the corrected count), and
panel d labels the Platinum reads 2×100 (they are 2×101).

## 1. It reads known copy numbers correctly, in every genome

Every sample carries sequence whose copy number is not in question, measured under the same
fragment-GC model and the same 2 × observed / expected ratio as the rDNA: 80 held-out autosomal
regions (two copies), 60 chrX and 40 chrY regions (one or two, one or none, by sex), counted by
their alignment position, which tests the GC model and the control denominator; and the *distal
junction*, a 400-kb sequence present once on each of the ten acrocentric short arms (multi-copy,
paralogous, acrocentric, i.e. the kind of sequence the rDNA is), counted by the same k-mer
classification and placement on the unit as the rDNA, which also tests the k-mer assignment.

| known quantity | reads | n |
| --- | --- | --- |
| held-out autosomal sequence (2) | 1.999 ± 0.008 | 735 |
| chrX in men (1) / in women with an intact culture (2) | 0.991 ± 0.006 / 1.933 ± 0.025 | 350 / 370 |
| chrY in men with an intact Y (1) / in women (0) | 0.975 ± 0.014 / 0.007 at most | 343 / 385 |
| distal junction (10) | 9.66, robust SD 0.15 | 614 |

Sex read from the X and Y agrees with the pedigree in 735 of 735: the highest X in a man is
1.00, the lowest in a woman with an intact culture 1.85. Fifteen women and seven men read below
those lines because their cell line has lost an X or a Y in part of its cells, and an independent
coverage pipeline reads the same people the same way (finding 6). **What it rules out:** that the
fragment-GC model, the k-mer assignment or the control regions are wrong in a way that shows.
*Panel a.*

## 2. A ten-copy paralog steps in whole copies, and the steps are inherited

Relative to the cohort's level, the distal junction sits at whole numbers: 4 people at −2, 16 at
−1, 614 at 0, 21 at +1, with a robust SD of 0.15 copies around each (51 people sit between
steps). A step is a structural variant of an acrocentric short arm; a two-copy loss is the
signature of a Robertsonian translocation (about one person in a thousand carries one). HG00651
and her daughter HG00652 both read −2, and both also carry 10–30% less of every satellite family
of the acrocentric short arms — ACRO, SST1, β-satellite, HSat3, CER — while the pan-centromeric
α-satellite that every chromosome carries is unchanged: the arms are missing, not just the
junction. A second family, HG01204 and his child HG01206, also reads −2 in both, but only the
child carries less of the acrocentric satellites (ACRO 0.77 and SST1 0.71 of the cohort's
median, against the father's 1.11 and 1.12): not every two-copy loss of the junction takes whole
arms with it. Where a carrier parent and a child were both counted, the step was transmitted in
6 of 23 pairs (one pair unclassified) — fewer than the half a heterozygous variant passes on,
though not clearly below it at 22 classified pairs (two-sided binomial p = 0.05). A pair counts
as transmitted when the child is nearer the parent's step than zero and within half a copy of
it, and as not transmitted when nearer zero and within half a copy of zero; HG01607 (+1.49) and
her child HG01609 (+0.58) fit neither. (The count first published, 7 of 22, came from
overlapping windows that counted a child within half a copy of both as transmitted and left out
the pair that fits neither; a later review corrected the rule.) One child, HG01517, carries a
step that neither parent has (+0.77 copies): a new structural variant, or a change in part of
its cell line. **What it rules out:** that multi-copy acrocentric sequence cannot be resolved to
a single copy by this path. *Panel b.*

## 3. The measured rDNA variation is inherited

A child's dosage is the average of the parents' plus segregation; measurement error is not
inherited. In 149 trios the midparent slope for the calibrated 45S estimate is 1.05 ± 0.07,
i.e. a reliability of 1 (95% interval 0.93–1.21); the interval allows a measurement error of
at most 6.4% of a person's value. Held-out autosomal sequence, which has no true variance, reads
−0.07 (−0.32 to 0.19); the cell line's mitochondrial content, which is not in the nuclear genome,
0.06 (−0.19 to 0.29). **What it rules out:** that the between-person variation is a property of
the sample preparation or the culture rather than the genome. *Panel c.*

Three honest notes. First, the design helps: every one of the 149 children is among the 698
related genomes the 1000 Genomes 30× release added to its original 2,504, and 287 of their 298
parents are among the 2,504, so no batch is shared within a family to imitate inheritance. The
same design means a difference between the batches is a difference between the generations:
the children read 4.8% more 45S than their parents, which is generation or batch, and a batch
that reads the children on a slightly different scale moves the midparent slope with it. With
the children first put on their parents' scale, the 45S reads 1.03 (0.94–1.13); for HSat2 and
ACRO, whose children spread 1.15× and 0.86× as much as their parents, the two readings part
(1.23 against 1.07, and 0.86 against 0.99). The "Mendelian" reliability, 1.5 − Var(child −
midparent) / Var(parent), turns out to be R + 1 + ρ/2 − s² with s that ratio of spreads, so on
this design it measures the batch; the page no longer shows it. Second, the spousal correlation
after centring within population is 0.17 (−0.01 to 0.34), down from 0.36 at 42 trios: spouses
share no DNA, so it is either structure shared within couples or chance, and it barely moves a
reliability near 1. Third, because people differ in 45S copy number by 22% (CV) and any
competent estimator errs by a few percent, *every* estimator has a reliability near 1 within one
pipeline: the trios establish that the variation is real, not which estimator measures it best
(the paired difference between the calibrated estimate and the published 18S ratio is +0.03,
−0.01 to +0.07). Ranking estimators takes the next two findings. The page's section 3.6 ("Class
by class") gives every metric in these terms.

## 4. The same person, sequenced twice, years apart, on different instruments

Twelve pilot genomes have an independent older library (HiSeq 2500 2×126 or HiSeq 2000 2×101,
2012–15) beside their NovaSeq 2×150 library (2019): different chemistry, read length, insert
size, depth and alignment pipeline (bwa-mem, ALT-aware, against the same GRCh38 analysis set,
but an older version with bwakit/postalt post-processing), and a GC response that is the reverse
of NovaSeq's. The aligner itself is the same family, so this does not test another aligner.
Test–retest reliability (ICC) across the two technologies:

| estimator | ICC | within-person CV | offset between technologies |
| --- | --- | --- | --- |
| calibrated, anchors chosen out of sample | **0.978** | 2.8% | +2% |
| 18S depth ratio, as the literature computes it | 0.19 | 22.6% | −27% |
| the same, after removing its offset (a batch correction) | 0.87 | 7.3% | — |

The two DNA batches come from different cultures of each cell line, so these are upper bounds
on the measurement error. **What it rules out:** that the number is a property of the library.
This is the finding that separates the calibrated estimate from a depth ratio. *Panel d.*

## 5. What the model removes is the library, not the person

Even within one chemistry, libraries differ in GC bias: across the cohort the rate at which
65%-GC fragments were sequenced, relative to each library's mean, runs from 1.18 to 1.29
(middle 80%). The 18S depth ratio divided by the calibrated estimate of the same sample follows
that bias with **r = 0.84** (0.82–0.86); the same 18S region under the fragment-GC model,
r = 0.08 (0.01–0.15). Within the cohort the effect is a few percent (SD of the log ratio 0.033)
— small next to how people differ, which is why finding 3 cannot see it — but it is systematic,
and across technologies (finding 4) it is 27%. **What it rules out:** that the GC model is
decoration. *Panel e.*

## 6. Another pipeline, the same files

Hall, Turner & Queitsch (2021) published rDNA copy number for 2,419 of these genomes, from the
same CRAMs: the mean depth over a 145-bp segment of the 18S (reads re-aligned from FASTQ to the
rDNA unit) divided by their chromosome-1 depth, which matches our duplicate-excluded depth (an
inference from their published values, not their stated method). On the 551 samples shared so far,
their value against a comparable 18S depth ratio computed from our counts (the callable 50-bp bins
of the 18S, no GC model; halved to their per-haploid scale): r = 0.984, theirs 1.08× ours. Dividing ours by the
unflagged fraction of control reads, the denominator their values imply, brings this to 0.98 (SD
0.035); dropping duplicate-flagged reads on both sides gives 1.03. Most of the offset therefore
fits a denominator without duplicate-flagged reads (the flag is set for 9% of single-copy reads);
their Methods describe re-alignment from FASTQ without duplicate marking, so their 18S numerator
most likely keeps duplicates. Against NGS-DOSE's calibrated estimate (halved) the correlation is
lower (0.96) and their values run higher (1.22×), because their ratio has no GC model (the cohort
page's section 3.7 gives both). **What it rules out:** that the number is idiosyncratic to this
code — and it gives a likely explanation of the one difference. *Panel f* (drawn with the both-sides correction,
1.031).

A second independent pipeline covers the dosage path: NGS-PCA's per-sample QC of the same CRAMs
(mosdepth coverage, duplicate-flagged reads excluded) gives mitochondrial copies per cell that
agree with ours at r = 0.995 on 735 genomes (theirs 0.92× ours, again the duplicate flag: a
16.6-kb genome at thousands-fold depth saturates the positions a duplicate marker distinguishes,
and the ratio falls with depth), a chrX coverage ratio that agrees at r = 0.9998, sex in 735 of
735, and the same fifteen women and seven men with partial loss of an X or a Y. The page's
section 3.7 carries it (`--qc`).

## 7. Against long-read assemblies

Assemblies collapse the rDNA, so they are no truth for it; but 54 of the 735 genomes belong to
people with an HPRC release-2 assembly, whose CenSat annotation of both haplotypes gives the
size of every satellite array — the same kind of sequence, measured by the same k-mer machinery.
A family is compared only in the assemblies whose arrays of it have no material gap. (The
numbers here count the stand-alone gap records of the annotation, which the comparison of that
date missed; a later review corrected it on the same data, and the conclusions did not change.)
Per genome the two measurements agree closely for most families: robust SD of the log ratio 3.0%
for α-satellite HORs, 3.2% CER, 3.5% ACRO, 5.5% β-satellite, 5.8% HSat1B, 6.3% HSat3 and 6.9%
HSat1A. β-satellite and CER sit on a line below equality because the panel sees only part of
them (their k-mer recall); SST1 and SATR sit further below, because the HPRC annotation labels
several times more sequence as those families than the CHM13 annotation the panels were built
from. Across people only HSat1B (r = 0.99) and ACRO (0.96) keep r ≥ 0.95, because r also depends
on how much people differ: α-satellite HOR mass differs between people by 4%, so its r is 0.71.
When only six people could be compared, seven families passed that mark; six people flattered
the correlations. All 15 outlying comparisons (in HSat1A, HSat1B, HSat2, HSat3 and ACRO) are
genomes whose assembly holds less than the reads show, as where part of an array is missing from
an assembly without a marked gap. HSat2 compares poorly (r = 0.33 in the 27 assemblies that
close its arrays), though it is inherited as faithfully as the other classes (R = 1.07 with the
children on their parents' scale): what it measures is heritable, but the assemblies do not yet
confirm that it is the mass of HSat2. For SST1 and SATR the two measurements disagree per genome
by more than people differ. **What it rules out:** that the k-mer path only appears to work
because nothing independent has been held against it. *Panel g.* The page's section 3.8 carries
the figures, one per family.

## 8. The targeted fetch equals the whole-file scan

Fetch mode reads the control regions and the few intervals where the aligner puts class reads:
here the 80 sink intervals (3.3 Mb) of the 45S, the 5S and the distal junction, learned from
whole-file scans of two genomes of this pipeline (NYGC, bwa-mem 0.7.15). By the pipeline's
estimate that is about 0.5 GB of a 15-GB CRAM and a minute over the network; the fetches here
read a staged local copy, in about six seconds. For all 735 genomes counted both ways it returns
0.9997 of the scan's 45S estimate (range 0.9992–0.9999), 0.9999 of the 5S and 0.998 of the
distal junction; the lowest 45S sink capture in any scan is 99.91%. Run separately on the fetch
counts alone, the whole pipeline — calibration, cohort layer and trio test — gives the same 45S
reliability. **What it rules out:** that 45S, 5S and distal-junction dosage needs every whole
file. On this pipeline, sinks learned from two whole-file scans give the scan's answer from a
targeted fetch. A biobank on another pipeline (DRAGEN untested at that date; one genome has since been
checked, [NGS-DOSE DESIGN §12](https://github.com/jlanej/NGS-DOSE/blob/main/docs/DESIGN.md)) still needs whole-file scans
of a small subset of its own CRAMs to learn its sinks, because sinks are specific to the aligner
and the reference. The satellite families of finding 7 and the telomeric repeat were not part of
this comparison: no fetch here carried them, and they were measured by whole-file scans. (A test
after this date found that, in this pipeline's CRAMs, sinks learned from 30 scans hold at least
99.8% of nine satellite families in each of 200 other genomes (99.85% in two of three random draws), and 96.6% of HSat1B; no
satellite sinks ship yet. The page's [section 3.3](https://jlanej.github.io/NGS-DOSE-1000G/#paths) lists which classes a fetch measures now.)
*Panel h.*

## And PCs?

Coverage principal components (NGS-PCA's, or the internal ones from the control regions) are
regressed out of log(estimate), and the adjusted values are written back on the original scale
to separate columns (`.adj`, `.adj_ngspca`); the estimates reported above are unadjusted. The
number of PCs is chosen at the Marchenko–Pastur edge of the noise bulk and checked against the
known truths and the trios. On this cohort they matter little for the rDNA, and the reason is
worth stating: measurement error is at most a few percent of a person's value and the variation
between people is 22%, so the technical share of the 45S variance is about 2% — the 16 control
PCs above the edge remove 4% against 2.2% expected by chance. Where there *is* technical
variance they find it: the same 16 PCs remove 34% of the variance of the held-out autosomal
estimate (which has nothing but error to remove), 22% of the mitochondrial and 20% of the EBV
dosage. The cross-validated sweep picks 5 PCs for the autosomal control and none for the rDNA.
PCs are insurance for small effects on a noisy trait; for rDNA copy number the GC model and the
calibration do the work.

## What is not yet shown

- The absolute scale has one external check: ddPCR (Potapova et al. 2025) on twelve
  lymphoblastoid lines with a NovaSeq genome, nine of them 1000 Genomes lines in this pipeline,
  where NGS-DOSE reads about 0.97× the assay (r ≈ 0.97; this repository's
  [`assembly_rdna`](../assembly_rdna/) study and the page's section 3.7). Beyond those lines
  the scale rests on unit windows where three Illumina chemistries agree.
- 5S copy number: reliability 0.78 with an interval from 0.46 to 1.11 at 149 trios — undecided.
- In these trios generation and batch go together (finding 3): whether children differ from
  their parents in level or spread cannot be told from the sequencing batch.
- All 735 genomes are one chemistry and one pipeline, and the fetch's sinks were learned on
  that pipeline; DRAGEN alignments and other chemistries are untested (one DRAGEN genome
  per version has since been checked, [NGS-DOSE DESIGN §12](https://github.com/jlanej/NGS-DOSE/blob/main/docs/DESIGN.md)).
- Every sample is a lymphoblastoid cell line.

## Reproduce

The page, its numbers and the current figure are rebuilt from the counts at hand, so they give
today's numbers, not this document's:

```bash
bash regenerate.sh   # runs: python -m report -p meta/20130606_g1k_3202_samples_ped_population.txt --hall meta/hall2021_MOESM1.txt --pilot pilot --pcs meta/ngspca/svd.pcs.txt -o docs --cache cache -j 4 --scan counts_scan --fetch counts_fetch --censat hprc_censat --qc meta/ngspca_sample_qc.tsv --ddpcr assembly_rdna/tables/potapova_comparison.tsv; then, with matplotlib, report.evidence_figure (docs/evidence.png) and report.trio_report (docs/trio_report.pdf)
```

The numbers of 2026-09-23 are in `docs/report.json` and `docs/data/` at commit 4a276d4
(`git show 4a276d4:docs/report.json`), built from the 735 genomes committed then.
