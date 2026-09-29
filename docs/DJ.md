# The distal junction: whole numbers of copies along the unit, held against assemblies and trios

*NGS-DOSE-1000G, 2026-09-29. The distal junction (DJ) is the 400-kb sequence beside every rDNA array on the five
acrocentric short arms, ten copies in a diploid genome. It is the cohort's ten-copy truth class. Its estimate used to
be one number that read 9.7 in most people and sat between whole numbers in one genome in five. This document
records why, what the method now does about it (NGS-DOSE 0.2.0), and how the result stands against two independent
truths: the HPRC release-2 assemblies of 28 cohort members and the 602 trios. It also records what the assemblies
themselves get wrong at this locus. The page ([index.html, section 3.2](index.html#djsteps)) recomputes every number
here from the committed tables, except the direct tests of `analysis/dj/`, which are named where they are used.*

## Summary

- **The junction is read in whole numbers of copies along its length.** Every genome's profile over the unit is
  read as a chain of whole numbers with its breakpoints. Of 3,202 genomes 3,173 calls are settled: 2,775 hold ten
  copies throughout, 109 nine, 7 eight, 27 eleven, 3 twelve and 1 thirteen, and 220 carry a copy that holds or lacks
  an end of the unit. The other 29 lie between two whole numbers throughout and are called uncertain.
- **Three corrections put the level where the copies are.** The level is set on the core of the unit, without the
  four stretches where copies differ. Its scale is pinned to the cohort's mode, which the fragment-GC model puts at
  9.78 copies and the assemblies at ten. The common deletions' windows are put on the cohort's own comb of whole
  numbers. The level is then 9.99 ± 0.33 copies over the cohort, and the genomes at ten scatter with a robust SD of
  0.13.
- **The values between whole numbers are partial copies, and they recur.** A copy of the first 316 kb of the unit is
  carried by 63 genomes (2.0%), a copy of the first 262 kb by 49 (1.5%), and a copy that lacks the first 122 kb by 19.
  A whole-unit median reads each as a fraction of a copy; the chain reads eleven copies over 316 kb and ten beyond.
- **The assemblies agree where they are resolved.** In 19 of 28 genomes the assembly resolves the junction. There
  the reads' core level lies −0.09 ± 0.20 copies from the assembly's, the called state equals the assembly's in 90.6%
  of 5-kb blocks, and all 7 partial copies are called with the breakpoint within 4 kb. All four carriers of a lost
  junction hold nine copies in their assemblies.
- **The trios agree.** In the 585 trios with three settled calls the parents' states allow the child's in 99.54% of
  34,852 core blocks, and 580 trios have no block out of place. The common deletions pass to children at the rate a
  germline variant should: one copy in one parent is found in 73 of 147 children at 197–217 kb.
- **Whole copies pass to fewer than half of the children, and the reason is open.** A whole-copy step passes to 15
  of 46 children and a partial copy to 33 of 81, together 48 of 127 (p = 0.008 against one half). The readings
  leave little room for error, and the events did not arise in culture at a rate the children's lines would show.
  The deficit rests on a father's loss of a junction, which passes to 4 of 23 children where a mother's passes to
  10 of 17.
- **What the assemblies get wrong**, each visible in the reads: copies cut by a contig end (29 of 300 copies, in 9 of
  28 genomes), fragments of one copy assembled twice, the ten copies phased 2 + 8, and a haplotype's whole-unit
  k-mer median under-reading seven copies.

## 1. Why the junction, and what "ten" means

Each acrocentric short arm (13, 14, 15, 21, 22) carries an rDNA array flanked distally by the distal junction, a
sequence assembled in T2T-CHM13 on all five arms and present nowhere else. The panel is the 169,808 31-mers that
occur exactly once in each of CHM13's five junctions and nowhere in GRCh38 outside them. A read is assigned to the
class by any four of its k-mers, and its 5′ end is counted in a 250-bp window of the unit (1,182 of the 1,600 windows
hold panel k-mers). The estimator predicts each window's count from the library's fragment-GC response measured on
single-copy control regions, and the cohort layer (a median polish over genomes × windows) removes each window's
shared efficiency.

Rhie et al. 2026 (bioRxiv, DJCounter) typed the same 3,202 genomes from the multiplicity of once-per-copy DJ 31-mers
against each library's two-copy k-mer peak, binned to integers with a Gaussian mixture, and examined the atypical
genomes in the HPRC release-2 assemblies. They report 9 copies in 2.8–3.4% of people and 11 or more in 8.4–9.3%, seven
genomes at 8 (one, HG01204, a G-banded Robertsonian carrier), and in the assemblies partial duplications for most of
the gains and, for the losses, three of six assembled cleanly. Section 4.7 holds the two studies together.

## 2. How the junction is read

The rules below are the class's own, in the bundle's `calibration.json` (`ngsdose` 0.2.0; NGS-DOSE `docs/DESIGN.md`,
section 7). A class without rules is read as before.

**The profile.** Within a genome the calibrated window estimates read the same copy number in every window that a
junction copy holds whole. A copy that lacks part of the unit lowers the windows it lacks, and an extra partial copy
raises the windows it holds. Across the cohort the SD of the calibrated estimate per 5-kb sub-block has a floor of
0.48 copies, which is counting noise, and rises far above it over 0–20 kb, 160–165 kb and 200–215 kb, where copies
differ between people.

**The core.** The level is the median over the windows outside 0–30 kb, 128–137 kb, 155–170 kb and 190–232 kb: 925
windows, 12 of the 20 blocks of 20 kb in full. A copy that lacks the distal 22 kb does not move it. The level over
the whole unit is kept beside it as `DJ.cn_unit`.

**The scale.** On the scale the fragment-GC model sets in the anchor windows (40–60% GC), the core of the cohort's
main mode reads 9.78 copies (2,920 genomes). The junction is ten copies in nearly everyone and the assemblies
confirm it (section 4.4), so every estimate is multiplied by 1.0221. A cohort of fewer than fifty genomes is not
pinned, and a saved efficiency table carries the pin to genomes counted later.

**The polymorphic intervals.** Where a deletion is common the cohort's median genome lacks part of a copy, and a
median over genomes gives the interval's windows an efficiency that is too low. Every genome then reads the interval
too high by one factor, and its comb of whole numbers sits half a copy off. The cohort's own comb gives the factor,
with the reference state (every copy holding the interval) the highest that many genomes share:

| interval | factor | every copy holds it | one lacks it | two lack it | three or more | one or more extra |
| --- | --- | --- | --- | --- | --- | --- |
| 5–15 kb | ×1.075 | 35% | 32% | 18% | 7% | 8% |
| 15–23 kb | ×1.028 | 65% | 27% | 4% | 0.4% | 4% |
| 197–217 kb | ×1.046 | 59% | 30% | 9% | 1.5% | 0.5% |

The shares are of the 2,755 genomes whose level lies within a quarter of a copy of a whole number, the ones the comb
is fitted on. The resolved assemblies give the same factors to 0.02 in the logarithm
(`analysis/dj/polymorphic_offsets.py`: +0.084, +0.034 and +0.062 against the cohort's +0.073, +0.028 and +0.045).

**The chain of whole numbers** (`ngsdose.segments`). A window in state *k* is expected at *k* copies times the
genome's own scale. An event, with its two changes of state, costs as much as 27 windows a copy off, and the most
probable chain is found for each scale of a grid. A prior on the scale (SD 1.5%, the level's SD among the genomes at
ten) keeps it from explaining a copy away. One copy over about 10 kb, or two over 4 kb, can be found. The genomes'
window noise is 0.65 copies per 250 bp.

**The lean.** Some genomes' profiles rise or fall smoothly along the unit: robustly an SD of 1.5% across it, more
than 5% in 5.6% of genomes and more than 8% in 1.1%. The lean goes weakly with the release batch (more than 5% in
5.1% of the earlier batch's genomes and 7.4% of the later's) and with the library's GC response (|r| at most 0.08),
and not with depth (r = −0.02). None of them accounts for it. A chain of whole numbers would break a lean into a
step, so it is a parameter of the chain, like the scale (prior SD 3%). A smooth rise costs less as a lean than as a
change of state, and a step, which a lean fits badly on both sides, stays a step.

**The description.** A genome is described against ten copies where it holds ten over 40 kb or more of the core,
otherwise against the state that holds most of it. A gain that reaches an end of the unit is a *partial copy*, a loss
that does is a *partial loss* (a copy that lacks that end), and anything else is a local gain or loss. A genome whose
level lies between two whole numbers throughout has two readings, and the call says how far behind the second is:
below three log units it is *uncertain*.

**What is written.** The sample table carries `DJ.cn` (the level on the core, pinned), `DJ.cn_unit`, `DJ.copies`,
`DJ.partial`, `DJ.variants`, `DJ.scale_f`, `DJ.tilt`, `DJ.call` and `DJ.call_gap`. Every segment of every genome is in
`data/dj_calls.tsv`, and every genome's profile in 20-kb blocks in `data/dj_blocks.tsv`.

## 3. What the calls were held against

**Assemblies.** 28 cohort members with an HPRC release-2 assembly, chosen to cover the whole-copy steps, the values
between steps and the mode: 56 haplotypes (hap1/hap2 for Hi-C-phased lines, mat/pat for trio children). Each
haplotype FASTA was screened for the panel's k-mers with `ngs-dose panel --report`, which counts every k-mer of the
unit in the assembly. A complete junction copy contributes one, so five complete copies read 5 at every k-mer, a copy
lacking a segment reads 4 across it, and a nucleotide difference in one copy lowers single k-mers only. The median
over the k-mers of a 5-kb sub-block is the haplotype's copies there, and the sum over the two haplotypes is the
genome's. Independently, each assembly was aligned (minimap2, asm20) to the unit with everything outside a panel
k-mer masked to N, so that each junction copy appears as a run of alignments on one contig. The 300 copies found are
classed as complete (231), complete with internal gaps (8), distal-start (19: the copies that begin 22 kb in),
truncated at a contig end (29) or partial with an internal breakpoint (13). An assembly counts as *resolved* at the
junction when no copy is cut by a contig end inside the unit and at most one copy is partial: 19 of 28.

**Trios.** The 602 trios of the cohort; 585 have three settled calls. The test is made position by position: a
parent whose state at a position is ten plus *d* carries *d* on its two haplotypes between them and passes on one, so
a child's deviation must be a part of its father's plus a part of its mother's.

**Reproducing it.** `pipeline/07_hprc_dj.sh` fetches, screens and aligns each haplotype (about four minutes each on
a laptop, most of it the download) and `python3 pipeline/hprc_dj.py screen DIR -o meta/dj_hprc` writes the two
committed tables. `regenerate.sh` reads them with `--dj-assemblies meta/dj_hprc`.

## 4. Findings

### 4.1 The cohort in whole numbers

| state | genomes | |
| --- | --- | --- |
| eight copies throughout | 7 | six of them read the ACRO1 composites of the acrocentric short arms at 0.75–0.84 of the cohort's median: the arms are missing, as in a Robertsonian translocation |
| nine throughout | 109 | 3.4% of settled calls |
| ten throughout | 2,775 | 87.5% |
| eleven, twelve, thirteen throughout | 27, 3, 1 | |
| a copy that holds or lacks an end of the unit | 220 | 146 with a partial copy, 75 with a partial loss |
| uncertain | 29 | the level between two whole numbers throughout |

The breakpoints recur. In the core they fall most often at 315 kb (68 genomes), 260 kb (51), 120 kb (17) and 390 kb
(17), and the copies that hold or lack an end are, by their extent:

| partial copy or loss | genomes | of the cohort | passed on (section 4.3) |
| --- | --- | --- | --- |
| copy of 0–316 kb | 63 | 2.0% | 8 of 17 |
| copy of 0–262 kb | 49 | 1.5% | 8 of 21 |
| loss of 0–122 kb | 19 | 0.6% | 2 of 6 |
| copy of 217–300 kb | 9 | 0.3% | |
| copy of 300–400 kb | 9 | 0.3% | |
| loss of 0–80 kb | 9 | 0.3% | 0 of 1 |
| loss of 0–60 kb | 8 | 0.2% | 0 of 2 |
| loss of 290–400 kb | 6 | 0.2% | 0 of 2 |

Of the 309 genomes whose core level lies more than 0.3 copies from a whole number, the calls make 118 a genome with
a copy that holds or lacks an end, 119 a ten-copy genome whose scale is a few percent off, 35 a ten-copy genome with
a lean of more than 5%, 19 uncertain, 11 another whole number throughout and 7 another event of 40 kb or more
(`analysis/dj/between_profiles.py`).

### 4.2 Against the assemblies

| measure | value |
| --- | --- |
| core level, reads less assembly, 19 resolved genomes | −0.09 ± 0.20 copies (mean ± SD); within half a copy in 19 of 19 |
| the same over all 28 | −0.09 ± 0.24; within half a copy in 26 of 28 |
| 20-kb blocks that round to the same copy number | 80% of 560 |
| blocks the assembly puts off ten that the reads confirm | 237 of 286 |
| blocks the assembly puts at ten that the reads put half a copy or more away | 31 of 274 |
| called state equal to the assembly's, 5-kb blocks, 18 resolved genomes that share the assembly's level | 90.6% of 1,404 |
| the same in the 9 fragmented assemblies | 61% of 702 |
| partial copies of the resolved assemblies that are called | 7 of 7, the breakpoint within 4 kb in each |
| partial copies of all 28 assemblies that are called | 10 of 13 |

*The whole-copy losses.* HG01891, HG00621 and HG03521 hold nine copies in a resolved assembly, and HG00658 holds
nine in an assembly whose paternal haplotype is fragmented. The reads call nine throughout in all four. The two that
Rhie et al. name as entire losses are the same two here (maternal in HG01891, paternal in HG00621).

*The recurrent gain.* In HG00146, HG00232, HG01786, HG03942, NA20752 and NA20805 the assembly holds, on one contig, a
partial copy spanning the first 316 kb of the unit followed by a complete copy. The calls read eleven copies up to
313–316 kb and ten beyond in all six. HG03654's partial copy spans 3–287 kb and is called over 22–291 kb. HG00320's
spans 187–400 kb and is called over 185–400 kb.

*Partial losses.* HG00673 and HG02523 each hold a copy that begins at 122 kb, and the calls read one copy fewer over
0–121 kb and 0–123 kb. HG01981's maternal copy ends at 277 kb; the call reads one copy fewer over 282–376 kb and ten
again beyond, so it is listed as a local loss and not as a partial loss. Of the three partial copies that are not
called, two are in HG00642, whose two copies are in pieces at contig ends, and the third is HG01981's.

*Where the two disagree on the level.* HG00673 reads one copy more than its resolved assembly throughout (ten with a
partial loss, against nine with one). Its level is 9.49 on the core, and the call takes its scale to be 4.5% low.
HG02280 reads 9.58 where its assembly holds ten, and is called ten. Neither has a chromosome flagged in the control
regions. A copy the assembly does not hold, a junction changed in part of the cell line (the assembly's DNA and the
reads' DNA are different cultures), and a genome whose level the model misplaces are all possible.

### 4.3 In the trios: what is consistent, and what is passed on

*Consistency.* In the 585 trios with three settled calls, leaving out the blocks within 10 kb of a breakpoint of any
of the three, the parents' states allow the child's in 99.54% of 34,852 core blocks and in 98.5% of 5,640 blocks of
the polymorphic intervals. 580 trios have no core block out of place. Of the blocks in which a child deviates from
ten, a parent explains 93.8% in the core. HG01517 holds a whole copy more than either parent, and NA18497 a partial
copy (56–400 kb) that neither has.

*The common deletions pass as a germline variant should.* These deletions are old and common, in the germ line
beyond doubt, so they test the measurement: a method that lost or invented copies would pass them to fewer than half
of a carrier's children. Three readings were made (the table is the page's; `report/dj.py`, `inheritance`). Without
a call, the child's value in the interval is regressed on the mean of its parents', and under Mendel the slope equals
the value's reliability. By the values, a genome within 0.3 of a whole number is given it. By the calls, a genome has
the called state in the middle of the interval less the state beside it.

| interval | windows | slope of child on midparent (95%) | reliability | one copy passed on, by the values | by the calls | called in the child of two parents without it | child less the mean of its parents |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 5–15 kb | 39 | 0.94 (0.86–1.02) | 0.98 | 22 of 45 | 57 of 123 | 20 of 121 | −0.12 ± 0.03 |
| 15–23 kb | 28 | 0.81 (0.68–0.95) | 0.89 | 35 of 72 | 67 of 149 | 27 of 320 | −0.11 ± 0.02 |
| 197–217 kb | 49 | 0.95 (0.85–1.04) | 0.95 | 62 of 110 | 73 of 147 | 23 of 249 | −0.03 ± 0.02 |

Over the three intervals one copy passes to 119 of 227 children by the values and to 197 of 419 by the calls, and
each slope's interval holds its reliability. The measurement passes on what the germ line does. Two things qualify the
calls in these short intervals. A deletion is called in 8–17% of the children of two parents without it, which is the
calls' error where an event is near the limit of what the chain can find. And in the two distal intervals the
children read 0.11–0.12 copies lower than the mean of their parents, which Mendel does not allow. Most children were
sequenced in the later release batch, and in 15–23 kb the parents of that batch read lower by the same amount
(`analysis/dj/intervals_by_batch.py`), so there the difference is the libraries'.

*Whole copies and partial copies pass to fewer than half.* Where a carrier holds one state throughout and the other
parent ten, a whole-copy step passes to 15 of 46 children (two-sided binomial p = 0.03 against one half). A partial
copy is looked for in the child by its breakpoint and not by its name, because a copy that holds the first 316 kb
and a copy that lacks the last 84 kb are the same step down at 316 kb, described against different tens. Where a
parent's call has one breakpoint in the core and the other parent is at ten with none near it, the child's call has
the breakpoint in 33 of 81 pairs (p = 0.12). Together that is 48 of 127 (p = 0.008). Four things were tested.

| question | test | result |
| --- | --- | --- |
| are the carriers' or the children's readings in doubt? | the level of each, where a whole-copy step was not passed on | the parent's level lies within 0.3 copies of its whole number in 29 of 31 pairs (−0.07 ± 0.15) and the child's within 0.3 of ten in 29 (0.00 ± 0.14); the two parents further off are fathers called at eleven |
| did the calls miss the child's copy? | the child's profile across the parent's breakpoint (`analysis/dj/partial_transmission.py`) | the children whose call lacks the breakpoint step by +0.02 copies there (median), those whose call has it by +0.96; by the profile 28 of 65 pairs are passed on, by the calls 27 |
| is it a property of reading a man's genome? | nine-copy calls by sex (`analysis/dj/carriers_by_sex.py`) | 57 of 1,587 men and 52 of 1,586 women are at nine throughout, and the level of the genomes at ten is 9.990 and 9.993; a father's loss passes to 2 of 8 sons and 2 of 14 daughters |
| did the events arise in the cell lines? | new events in the children of two parents at ten throughout | 147 of 1,160 parents carry an event; had the excess over half-transmission arisen in culture, the 443 such children would hold about 14 new ones; they hold 2 |

By the parent, a father's loss of a junction passes to 4 of 23 children and a mother's to 10 of 17 (Fisher exact
p = 0.009, a comparison made after the fact). A father's gain passes to 0 of 4 and a mother's to 1 of 2. A partial
copy passes from a father in 21 of 49 pairs and from a mother in 12 of 32. Children carry fewer of these events than
their parents: 62 of 585 against 147 of 1,160.

What remains is a difference between the generations or in the germ line. A change in a donor's blood that a clonal
cell line makes whole would be commoner in older donors, and would be passed on by no one. A variant that is passed
on less often than chance would show the same counts. Neither is tested here. The partial copies alone are within
chance of one half, and the deficit rests on the whole-copy steps, a father's losses most of all: 23 pairs.

### 4.4 Why the fragment-GC model's scale reads 9.78

Nine genomes whose resolved assembly holds ten complete copies (HG00097, HG00438, HG00735, HG01255, HG02155, HG02258,
HG02280, HG02615, HG03239) read 9.37–9.79 on the fragment-GC model's scale (median 9.65), and 9.58–10.01 once it is
pinned. The deficit is in the measurement's scale: the assemblies hold the copies the reads do not fully count.

| hypothesis | test | result |
| --- | --- | --- |
| the targeted fetch misses junction reads | DJ reads in fetch against whole-file scan, 3,202 genomes (`fetch_capture.py`) | the fetch holds 99.76% (range 99.64–99.84%); not it |
| duplicate flags remove piled-up reads | engine source | duplicate-flagged reads are counted (tallied apart); not it |
| divergence among the ten copies loses k-mers | exact-k-mer presence per window in HG00097's assembly against window GC (`presence_vs_gc.py`) | presence 0.95 per haplotype and flat in GC (r = +0.08), while the raw estimate of the mode genomes falls from 10.09 below 35% GC to 9.28 at 50–60% GC (r = −0.26); and a read is classified by any four of its k-mers; not it |
| mosaic loss of an acrocentric chromosome in the line | per-chromosome dosage from the 800 control regions, 600 genomes (`acro_dosage.py`) | dosage SD 0.02–0.05 copies; r with the junction's step −0.08; ten genomes off by more than 0.15 copies; not the bulk |
| replication timing of the culture's DNA | the same twelve people on two chemistries (the pilot) | the older libraries read the junction 2.1% higher than the NYGC ones; a property of the library, not only of the DNA |
| the GC model's residual in repeat context | the cohort's window efficiencies against window GC | they fall with GC: −0.31 in the logarithm per unit GC for the junction (r = −0.26, 1,182 windows) and −0.53 for the 45S unit (r = −0.32, 167 windows), whose scale the shipped anchors correct; the surviving explanation |

The 45S unit's scale is pinned by anchor windows chosen on replicate pairs across chemistries. The junction has no
such windows, so its scale is pinned to the mode of the cohort it is counted in.

### 4.5 Which is more accurate?

For a genome whose junction assembled cleanly the two agree to the noise of the reads: 0.13 copies on the level, and
one 5-kb block in ten on the called state. Where the junction did not assemble cleanly, in 9 of 28 genomes, the reads
are the more coherent reading of the locus, and the profile says which blocks the assembly broke. The reads are
weakest in three places. Events shorter than about 10 kb at one copy are not called. In the short polymorphic
intervals a deletion is called in the child of two parents without it in 8–17% of trios, which is the calls' error
there. And 29 genomes cannot be given a whole number at all.

### 4.6 What the assemblies get wrong, and what the reads show

The junction sits beside the rDNA array, where long-read assemblies break, and its copies are more than 99%
identical. The catalogue of copies (`data/dj_hprc_copies.tsv`) and the profiles show:

1. **Copies cut by a contig end** (29 of 300 copies, in HG00146, HG00320, HG00642, HG00658, HG01943, HG01981, HG02040,
   HG02392, HG02523). The k-mers of the missing part are absent, and the assembly under-counts those blocks. In
   HG02523 two paternal copies begin at 86 and 109 kb at the starts of 0.36-Mb and 0.37-Mb contigs. The k-mer sums
   read 6–8 copies over 0–120 kb where the reads read 9, the one real partial loss. In HG00642, a ten-copy genome by
   the reads, the maternal 22–155 kb and paternal 3–134 kb pieces sit at contig ends and the 111–400 kb pieces on
   three contigs of 0.3–0.4 Mb. The k-mer sums read 8–9 over 0–120 kb and 11–12 over 120–160 kb, where the
   fragments overlap.
2. **Fragments assembled twice.** HG00658's paternal haplotype has a copy ending at 343 kb at a contig's end and two
   small contigs (0.17 and 0.23 Mb) each carrying the unit's last 80 kb. The k-mer sums read 11 over 320–400 kb, the
   reads 9.0, as everywhere else in this nine-copy genome. HG00146's first haplotype has a 0.05-Mb contig holding
   66–117 kb, and HG00320's second a 0.18-Mb contig with 368–400 kb.
3. **Phasing.** The ten copies split 2 + 8 (HG00146), 3 + 7 (HG02040), 2 + 7 (HG03521), 6 + 4 (HG00320, HG02155,
   HG02392) and 4 + 6 (HG01786, NA20752). The acrocentric short arms are phased by Hi-C or trio k-mers no better than
   their near-identical sequence allows. The sum is what the reads measure, and it is right.
4. **A haplotype's whole-unit k-mer median under-reads a high copy number.** HG01786's second haplotype holds six
   complete copies and the 316-kb partial. Its k-mer median over the unit is 6 while its block medians over the core
   are 7, because with seven copies the nucleotide differences among them leave fewer than half the k-mers at the
   full count. An exact-k-mer count, from an assembly or from reads, needs the median per block. Nine haplotypes
   show it.
5. **A locus the assembler could not lay out.** HG01943 has its junction in ten pieces of 0.1–1.1 Mb across both
   haplotypes, and the k-mer sums swing between 8 and 13 along the unit. The reads call twelve copies over 0–261 kb
   and ten beyond, with a lean of 12% that is the largest among the compared genomes: two extra copies of the distal
   261 kb. Its mother HG01942 carries the same two copies.
6. **Inverted segments within copies** appear as mixed strands in the alignments of many complete copies (the unit
   holds palindromes) and are not counted as breaks.

### 4.7 The comparison with Rhie et al. 2026

Their eight-copy genomes are the seven at eight throughout here. Their entire losses in HG00621 and HG01891 are
nine-copy calls, HG01981's partial loss is a loss over 282–376 kb, and the duplications in NA20752, NA20805, HG03654
and HG00320 are partial copies with the duplicated stretch placed. Their 9 copies in 2.8–3.4% of people are 3.4% at
nine throughout here. Their 11 or more in 8.4–9.3% compare with 5.9% that carry a whole or partial extra copy here
and 7.3% whose level is 10.3 or more; the two typings bin the partial copies differently, and the difference is not
resolved here. They report parents typed higher than their children (p = 3.7 × 10⁻⁶). On the pinned level here the
children's median is 10.03 and the parents' 9.99, and the children carry fewer large events than their parents.

## 5. What changed

**In NGS-DOSE (0.2.0, branch `claude/dj-calls`).** Class rules in the bundle (`calibration.json`): the windows left
out of the level, the pin of the scale to the cohort's mode, the polymorphic intervals and their comb, and the
chain's parameters. A segmentation module (`ngsdose/segments.py`). New columns in the cohort table and
`ngsdose cohort --segments`. `--no-class-rules` reads a class as before.

**In this repository.** The page's section 3.2 is written on the calls: the states, the partial copies by extent,
the position-by-position test in the trios, the transmission of the common deletions and of the large events, and the
assemblies against the calls. The sample table carries the calls, and `data/dj_calls.tsv` every segment. The figure
shows the assembly, the profile and the call of every compared genome. `analysis/dj/` holds the direct tests.

## 6. What remains

1. **The transmission of the large events.** The donors' ages, or a second tissue of the same people, would separate
   a change in the blood from a variant that is passed on less often. Long reads of a father at nine and his children
   would say which arm lacks the junction.
2. **A truth panel of 200.** A screen of all assembled cohort members costs about four minutes a haplotype. It would
   test the calls on genomes chosen without regard to what they read, which the 28 were not.
3. **The short intervals.** The polymorphisms at 130–135, 160–165 and 225–230 kb show in the segment map, are left
   out of the level and are not genotyped. A comb per interval, as the three common deletions have, would genotype
   them.
4. **A genome counted alone.** Without a cohort or a saved efficiency table its level rests on the GC model, a few
   percent low, and it has no calls.
5. **The two release batches.** In 15–23 kb the later release batch reads 0.15 copies lower than the earlier one,
   in parents as in children. Calibration by batch is a decision for the whole method and is not made here.

## 7. Data and reproducibility

| file | what |
| --- | --- |
| `meta/dj_hprc/haplotypes.tsv`, `meta/dj_hprc/copies.tsv` | the 56 screens: k-mer medians per 5-kb sub-block, and the alignment catalogue of 300 copies (`pipeline/hprc_dj.py screen`) |
| `data/dj_hprc.tsv`, `data/dj_hprc_blocks.tsv`, `data/dj_hprc_copies.tsv` | the page's comparison: per genome, per genome and block, and the copies |
| `data/dj_blocks.tsv`, `data/dj_calls.tsv`, `data/dj_segments.tsv` | every genome's profile in 20-kb blocks with its call; every segment of every call; the cohort's segment map |
| `dj_assemblies.png` | the figure (`python -m report.dj_figure`) |
| `pipeline/07_hprc_dj.sh`, `pipeline/hprc_dj.py` | fetch, screen and align the haplotypes; write the tables |
| `report/dj.py`, `report/dj_page.py`, `report/dj_figure.py` | the analysis, the page section, the figure |
| `analysis/dj/` | the direct tests outside the page ([README](../analysis/dj/README.md)) |

## 8. Limitations

The 28 genomes were chosen for what their estimates showed, so they over-represent steps and values between steps.
The assemblies are release-2 drafts whose acrocentric short arms are among the least resolved sequence they contain,
and the assembly's DNA is a different culture of each line from the reads' DNA. The exact-k-mer screen counts a copy
at every k-mer it holds verbatim, so a copy diverged by more than a few percent from CHM13 would be under-counted;
none was seen. The pin assumes that the cohort's mode is ten copies, which holds for a human cohort and would not for
a set of genomes chosen for their junctions. The chain's parameters were set on this cohort, with the trios and the
assemblies in view, so the agreement with both is not an out-of-sample test. The comparisons of transmission by the
parent's sex were made after the counts were seen.
