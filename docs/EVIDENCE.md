# The case that it works

What would convince someone who doubts that rDNA copy number can be measured from short-read
whole-genome sequencing? Not a model, and not agreement with ourselves. This page collects the
evidence, in the order a sceptic would ask for it, with the number, what it rules out, and where
it comes from. Everything here is recomputed by scripts from committed data: the
[cohort page](https://jlanej.github.io/NGS-DOSE-1000G/) (`python -m report`; the run is complete:
all 3,202 genomes of the expanded 1000 Genomes cohort, scanned and fetched, and all 602 trios the
release sequenced, the pedigree's 603rd naming a parent it never sequenced; as of 2026-09-28) and the [pilot report](../pilot/pilot_report.md) (`evaluate_pilot.py`,
12 genomes each sequenced twice). An earlier version of this write-up, dated 2026-09-23, was made
from the first 735 genomes and 149 trios; where the complete cohort changed a conclusion, the
change is said.

![the evidence](evidence.png)

## 1. It reads known copy numbers correctly, in every genome

Every sample carries sequence whose copy number is not in question, measured by exactly the
code that measures the rDNA: 80 held-out autosomal regions (two copies), 60 chrX and 40 chrY
regions (one or two, one or none, by sex), and the *distal junction*, a 400-kb sequence present
once on each of the ten acrocentric short arms — multi-copy, paralogous, acrocentric, i.e.
the kind of sequence the rDNA is.

| known quantity | reads | n |
| --- | --- | --- |
| held-out autosomal sequence (2) | 1.997 ± 0.008 | 3,202 |
| chrX in men with one X (1) / in women with an intact culture (2) | 0.992 ± 0.006 / 1.939 ± 0.038 | 1,596 / 1,543 |
| chrY in men with an intact Y (1) / in women (0) | 0.983 ± 0.044 / 0.012 at most | 1,585 / 1,605 |
| distal junction (10) | 9.99 ± 0.33, its scale pinned to the cohort's mode (9.78 on the fragment-GC model's scale); robust SD 0.13 among the genomes at ten | 3,202 |

Sex read from the X and Y agrees with the pedigree in 3,200 of 3,202: the highest X in a man
with one X is 1.05, the lowest in a woman with an intact culture 1.85. Sixty-two women and
twelve men read below those lines because their cell line has lost an X or a Y in part of its
cells, and an independent coverage pipeline reads the same people the same way (finding 6). The
two mismatches are of two kinds, and the reads say which: HG02300 carries two X chromosomes and
no Y under a man's record (another person's DNA, a swapped sample), and NA19226 one X and no Y
under a man's record (a line that has lost its Y in culture). One man, HG01683, reads two X
chromosomes and a Y, as a 47,XXY karyotype would, by both routes. **What it rules out:** that the
fragment-GC model, the k-mer assignment or the control regions are wrong in a way that shows.
*Panel a.*

## 2. A ten-copy paralog is read in whole copies, along its length, and the readings are Mendelian

The distal junction is ten copies of a 400-kb sequence, and the method reads it twice over: as a
level (9.99 ± 0.33 copies once its scale is pinned to the cohort's mode; the genomes at ten scatter
with a robust SD of 0.13), and as whole numbers of copies called along the unit for every genome
([DJ.md](DJ.md); the page's section 3.2). Of 3,202 calls 3,173 are settled. 2,775 genomes hold ten
copies throughout, 109 nine, 7 eight, 27 eleven, 3 twelve and 1 thirteen; 220 carry a copy that
holds or lacks an end of the unit, the commonest a copy of the first 316 kb (63 genomes) and one
of the first 262 kb (49); 29 lie between two whole numbers throughout and are called uncertain.
A loss of a whole junction is a structural variant of an acrocentric short arm, and a two-copy
loss is the signature of a Robertsonian translocation (about one person in a thousand carries
one). Six of the seven two-copy carriers also read the ACRO1 composites of the acrocentric short
arms at 0.75–0.84 of the cohort's median, and most of them less of SST1, β-satellite, HSat3 and
CER as well, while the pan-centromeric α-satellite that every chromosome carries is unchanged
(0.91–1.10): the arms are missing, not just the junction. Two of them are parent and child:
HG01204 and HG01206 (PUR) both hold eight, as do HG00651 and her daughter HG00652 (CHS).

Two independent truths test the calls. **Assemblies:** in the HPRC release-2 assemblies of 28
cohort members, screened copy by copy, all four carriers of a lost junction hold nine copies
(HG01891, HG00621 and HG03521 in a resolved assembly, HG00658 in one whose paternal haplotype is
fragmented; the two Rhie et al. 2026 name as entire losses are the first two), all 7 partial copies
of the 19 resolved assemblies are called with the breakpoint within 4 kb, and the called state
equals the assembly's in 90.6% of the 5-kb blocks of the 18 genomes whose level the assembly shares
(in HG00673 the reads hold one copy more than the assembly throughout). The other 9 assemblies are
fragmented at the junction (copies cut by contig ends, fragments assembled twice, the ten copies
phased 2 + 8), and there the reads are the more coherent reading of the locus. **Trios:** of the
585 trios whose three calls are settled, 580 have no 5-kb block of the core that the parents'
states do not allow (99.5% of 34,852 blocks), and a parent explains 94% of the blocks in which a
child deviates from ten. One child (HG01517) holds a whole copy more than either parent and one
(NA18497) a partial copy that neither has. The junction's common deletions, which are in the germ
line beyond doubt, pass as a germline variant should: where one parent lacks one copy of
197–217 kb and the other none, 73 of 147 children lack one, and the child's value regressed on the
mean of its parents' has a slope of 0.95 (0.85–1.04) against a reliability of 0.95. The
measurement passes on what the germ line does.

What the calls do not explain is how often a whole copy is passed on. Where a carrier holds one
state throughout and the other parent ten, a whole-copy step passes to 15 of 46 children
(two-sided binomial p = 0.03 against one half). A partial copy, looked for in the child by its
breakpoint, passes to 33 of 81 (p = 0.12), and the two together to 48 of 127 (p = 0.008). A
father's loss of a junction passes to 4 of 23 children, a mother's to 10 of 17 (Fisher exact
p = 0.009, a comparison made after the fact). The readings leave little room: where a step was not
passed on, the parent's level lies within 0.3 copies of its whole number in 29 of 31 pairs and the
child's within 0.3 of ten in 29. Nine-copy calls are as common in men as in women (57 and 52).
Nor did the events arise in culture at a rate the children's lines would show: had the excess
arisen there, the 443 children of two parents at ten throughout would hold about 14 new events,
and they hold 2. What remains is a change in a donor's blood that a clonal line makes whole, the
likelier the older the donor, or a variant passed on less often than chance; neither is tested
here, and the deficit rests on 23 fathers. **What it rules out:** that multi-copy acrocentric
sequence cannot be resolved to a single copy by this path, that the values between whole numbers
are measurement noise, and that the junction's 2.2% deficit on the fragment-GC model's scale is a
loss of copies in the cell lines (nine genomes that read 9.37–9.79 on that scale hold ten complete
copies in their assemblies).
*Panel b.*

## 3. The measured rDNA variation is inherited

A child's dosage is the average of the parents' plus segregation; measurement error is not
inherited. In 602 trios the midparent slope for the calibrated 45S estimate is 0.955 ± 0.045: a
reliability of 0.95 (95% family-bootstrap interval 0.86–1.04), 0.94 (0.87–1.01) with the children
first put on their parents' scale. The interval allows a measurement error of at most 8.1% of a
person's value; the pilot's replicates across technologies put the measurement's own error at
2.8% (finding 4), so at the point estimate the few percent not inherited would lie in the cell
lines or in transmission rather than in the measurement. The negative controls read near zero:
held-out autosomal sequence, which has no true variance, −0.04 (−0.18 to 0.10); the cell line's
mitochondrial content, which is not in the nuclear genome, 0.07 (−0.08 to 0.22); its EBV load
0.15 (−0.14 to 0.54); sequencing depth −0.07 (−0.18 to 0.02). The 5S array, undecided at 149
trios (0.78, 0.46–1.11), is decided at 602: 1.02 (0.88–1.15), 0.95 (0.86–1.03) on the parents'
scale. **What it rules out:** that the between-person variation is a property of the sample
preparation or the culture rather than the genome. *Panel c.*

Four honest notes. First, the design helps: 597 of the 602 children are among the 698 related
genomes the 1000 Genomes 30× release added to its original 2,504, and 1,097 of their 1,204
parents are among the 2,504, so a batch shared within a family, which could imitate inheritance,
is all but absent (85 families have a parent in the child's batch). The same design means a
difference between the batches is a difference between the generations: the children read 3.5%
more 45S than their parents (+2.3% to +4.8%), which is generation or batch, and a batch that
reads the children on a slightly different scale moves the midparent slope with it. Where the
children's spread departs from their parents' the two readings part: HSat2, whose children
spread 1.08× as much, reads 1.03 from the slope and 0.96 on the parents' scale; the 5S, 1.08×,
1.02 and 0.95. The two bracket the reliability. Second, one negative control has an interval
wholly below zero, insert size (−0.31, −0.45 to −0.18), and the reason is the parents' own
correlation: father and mother agree at ρ = 0.20, both having been sequenced in the earlier
batch, whose libraries share a fragment-size distribution, and the correction R = b − ρ(1 − b)
subtracts what a batch shared by the parents would otherwise pass for transmission; the slope
itself, −0.09 ± 0.04, is no different from zero. Third, the spousal correlation of the 45S after
centring within population is 0.06 (−0.02 to 0.14), down from 0.17 at 149 trios and 0.36 at 42:
spouses share no DNA, and it barely moves the reliability. Fourth, because people differ in 45S
copy number by 22% (CV) and any competent estimator errs by a few percent, *every* estimator
has a reliability near 1 within one pipeline. At 602 trios the difference between them is
nevertheless resolved: the published 18S depth ratio, computed from the same reads, reads 0.90
(0.80–1.00), and the paired family bootstrap puts the calibrated estimate above it by +0.05
(+0.02 to +0.07; 999 of 1,000 draws), the same reads with the library model and the calibration
being the only difference. Ranking estimators still rests mainly on the next two findings. The
page's section 3.6 ("Class by class") gives every metric in these terms, and its split by the
sex of parent and child shows the one Y-linked class as such: the HSat1B array passes from
father to son at r = 0.80 and to daughters at 0.14, where every autosomal class reads about 0.5 in
all four pairings.

## 4. The same person, sequenced twice, years apart, on different instruments

Twelve pilot genomes have an independent older library (HiSeq 2500 2×126 or HiSeq 2000 2×100,
2012–15) beside their NovaSeq 2×150 library (2019): different chemistry, read length, insert
size, depth, aligner, and a GC response that is the reverse of NovaSeq's. Test–retest
reliability (ICC) across the two technologies:

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
65%-GC fragments were sequenced, relative to each library's mean, runs from 1.17 to 1.30
(middle 80%). The 18S depth ratio divided by the calibrated estimate of the same sample follows
that bias with **r = 0.85** (0.85–0.86; n = 3,202); the same 18S region under the fragment-GC model,
r = 0.05 (0.01–0.08). Within the cohort the effect is a few percent (SD of the log ratio 0.035)
— small next to how people differ, which is why finding 3 can only just see it — but it is
systematic, and across technologies (finding 4) it is 27%. **What it rules out:** that the GC
model is decoration. *Panel e.*

## 6. Another pipeline, the same files

Hall, Turner & Queitsch (2021) published rDNA copy number for 2,419 of these genomes, from the
same CRAMs, as the 18S depth relative to chromosome 1 with duplicate-flagged reads excluded.
On all 2,419: r = 0.984. Their values are 1.07× ours; re-applying their duplicate exclusion to
our counts brings this to 1.03 (SD 0.04) — most of the offset is the duplicate flag, which is
set for 5.2% of rDNA reads but 9.0% of single-copy reads (the collapsed rDNA hides duplicates
from the marker), by an amount that differs between samples (0.43–0.85 of the control rate).
Against the calibrated estimate their values run 1.22× (r = 0.96), the GC model and the
calibration on top of the duplicate flag. **What it rules out:** that the number is idiosyncratic
to this code — and it explains the one difference. *Panel f.*

A second independent pipeline covers the dosage path: NGS-PCA's per-sample QC of the same CRAMs
(mosdepth coverage, duplicate-flagged reads excluded) gives mitochondrial copies per cell that
agree with ours at r = 0.996 on 3,199 genomes (theirs 0.92× ours, again the duplicate flag: a
16.6-kb genome at thousands-fold depth saturates the positions a duplicate marker distinguishes,
and the ratio falls with depth), a chrX coverage ratio that agrees at r = 0.9997, sex in 3,198
of 3,199, and the same 62 women and 12 men with partial loss of an X or a Y. The page's
section 3.7 carries it (`--qc`).

## 7. Against long-read assemblies

Assemblies collapse the rDNA, so they are no truth for it: the HPRC release-2 assemblies of 200
cohort members hold about a third of the rDNA their measured copy number implies, in pieces no
longer than 1.7 Mb. But their CenSat annotation of both haplotypes gives the size of every
satellite array — the same kind of sequence, measured by the same k-mer machinery. Per genome
the two measurements agree closely for most families: robust SD of the log ratio 2.2% for CER,
2.8% ACRO, 3.1% α-satellite HORs, 5.0% β-satellite, 5.3% SATR, 6.9% HSat3, 7.7% HSat1A and 8.0%
HSat1B. β-satellite and CER sit on a line below equality because the panel sees only part of
them (their k-mer recall); SST1 and SATR sit further below, because the HPRC annotation labels
several times more sequence as those families than the CHM13 annotation the panels were built
from. Across people, r depends on how much people differ as well as on the agreement: HSat1B
(r = 1.00, people differ by 83%), ACRO (0.95), β-satellite (0.93) and CER (0.91) track the
assemblies; HSat1A (0.79) and HSat3 (0.76) less closely; the α-satellite HORs, on which the two
measurements agree to 3% while people differ by only 5%, reach 0.70. Of the 62 comparisons more
than three robust SDs from their family's median, 46 are genomes whose assembly holds less than
the reads show, as where part of an array is missing from an assembly without a marked gap. HSat2
compares poorly (r = 0.39 in the 79 assemblies that close its arrays; 121 leave part of an array
as a gap), though it is inherited as faithfully as the other classes (R = 0.96, 0.90 to 1.03,
with the children on their parents' scale): what it measures is heritable, but the assemblies do
not confirm that it is the mass of HSat2. For SST1 and SATR the two measurements disagree per
genome by more than people differ. **What it rules out:** that the k-mer path only appears to
work because nothing independent has been held against it. *Panel g.* The page's section 3.8
carries the figures, one per family.

## 8. The one-minute fetch equals the whole-file scan

Fetch mode reads the control regions and the few intervals where the aligner puts class reads —
about 0.5 GB of a 15-GB CRAM, a minute over the network, six seconds from disk. For all 3,202
genomes counted both ways it returns 0.9996 of the scan's 45S estimate (range 0.9988–0.9999),
0.9999 of the 5S and 0.998 of the distal junction; the lowest 45S sink capture in any scan is
99.85%. Run separately on the fetch counts alone, the whole pipeline — calibration, cohort layer
and trio test — gives the same 45S reliability, 0.95 (0.86–1.04). The telomeric repeat's sinks
entered the bundle part-way through the run, so 1,454 of the fetches carry it, at 0.9988 of the
scan. **What it rules out:** that a biobank would need the whole files. *Panel h.*

## And PCs?

Coverage principal components (NGS-PCA's, or the internal ones from the control regions) are
regressed out of the estimates, with the number chosen at the Marchenko–Pastur edge of the
noise bulk and checked against the known truths and the trios. On this cohort they matter
little for the rDNA, and the reason is worth stating: measurement error is at most a few
percent of a person's value and the variation between people is 22%, so the technical share of
the 45S variance is small — the 22 control PCs above the edge remove 2.4% of it against 0.7%
expected by chance, NGS-PCA's 40 genome-wide PCs 8.4% against 1.3%. Where there *is* technical
variance they find it: the same control PCs remove 26% of the variance of the held-out autosomal
estimate (which has nothing but error to remove), 31% of the mitochondrial and 15% of the EBV
dosage. The cross-validated sweep picks 3 PCs for the autosomal control, 6 for chrX, 2 for chrY,
13 for the distal junction and none for the rDNA. PCs are insurance for small effects on a noisy
trait; for rDNA copy number the GC model and the calibration do the work.

The one place they change a reading is between populations. Unadjusted, the median 45S copy
number runs from 434 in EUR to 524 in AFR (21% apart); with NGS-PCA's coverage PCs regressed out
the medians are AFR 486, EAS 488, AMR 476, SAS 469 and EUR 452 (8% apart). Those PCs are built
from genome-wide coverage and carry ancestry along with technique, so the adjusted spread is a
floor and the unadjusted one a ceiling on what is genomic in the difference.

## What is not yet shown

- The absolute scale has one external check: ddPCR (Potapova et al. 2025) on twelve
  lymphoblastoid lines with a NovaSeq genome, nine of them 1000 Genomes lines in this pipeline,
  where NGS-DOSE reads 0.96× the assay (r = 0.94, Spearman 0.97; the results repository's
  `assembly_rdna` study and the page's section 3.7). Beyond those lines the scale rests on unit
  windows where three Illumina chemistries agree.
- In these trios generation and batch go together (finding 3): whether children differ from
  their parents in level or spread cannot be told from the sequencing batch.
- All 3,202 genomes are one chemistry and one pipeline; DRAGEN alignments and other chemistries
  are untested.
- Every sample is a lymphoblastoid cell line.
- What the HSat2 panel measures is heritable but not yet confirmed to be HSat2 mass (finding 7),
  and the telomeric-repeat class is a relative measure of (TTAGGG)n content.
- Two engine builds counted the cohort (1,748 genomes with 0.1.0+fae1124, 1,454 with
  0.1.0+7772e32); NGS-DOSE's tests assert byte-identical counts across them, and the page's
  provenance records which build counted each genome.

## Reproduce

```bash
python -m report --scan counts_scan/ --fetch counts_fetch/ -p meta/20130606_g1k_3202_samples_ped_population.txt --hall meta/hall2021_MOESM1.txt --pilot pilot --pcs meta/ngspca/svd.pcs.txt --censat hprc_censat --qc meta/ngspca_sample_qc.tsv --ddpcr assembly_rdna/tables/potapova_comparison.tsv -o docs/
python -m report.evidence_figure --report docs/report.json --pilot pilot -o docs/evidence.png
python -m report.trio_report --report docs/report.json --data docs/data -o docs/trio_report.pdf
```
