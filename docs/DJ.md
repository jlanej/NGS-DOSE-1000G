# The distal junction: NGS-DOSE against long-read assemblies

*NGS-DOSE-1000G, 2026-09-28. The distal junction (DJ) is the 400-kb sequence beside every rDNA array on the five
acrocentric short arms, ten copies in a diploid genome. It is the cohort's ten-copy truth class, and its
estimate reads 9.72 rather than 10 in most people while 20% of people sit off a whole number. This document
compares the estimate with the HPRC release-2 assemblies of 28 cohort members, explains where the deviations from
ten come from, records what the assemblies themselves get wrong at this locus, and says what follows for the method.
The page ([index.html, section 3.2](index.html#djsteps)) recomputes the assembly comparison from the committed
tables; the numbers here are the page's, with the direct tests from this analysis added.*

## Summary

- **Whole-copy steps are real.** Every carrier of a −1 step with a resolved assembly (HG01891, HG00621, HG00658,
  HG03521) has nine junction copies in it; the two Rhie et al. name as entire losses are the same two (maternal in
  HG01891, paternal in HG00621). Gains of +0.85 to +0.87 (HG01786, HG03942) are an extra *partial* copy.
- **The values between steps are partial variants, and one of them recurs.** Six gain genomes carry the same
  structure: a partial copy holding the first 316 kb of the unit in tandem with a complete copy, which reads 11 over
  the first 320 kb and 10 beyond; a whole-unit median puts it at +0.5 to +0.9. Partial losses (HG01981: one copy ends
  at 277 kb; HG02523 and HG00673: one copy begins at 122 kb) read −0.3 to −0.7 the same way. Across the cohort,
  breakpoints of half a copy or more fall at 20 kb in 1,410 genomes, at 200 kb in 904 and at 220 kb in 1,018: the
  distal 22 kb and the 200–215 kb segment are deletion polymorphisms of single copies (19 of 300 assembled copies begin
  22 kb in; one copy of HG00097 lacks 193–222 kb), and the 0–20 kb block alone reads 9 copies in 24% of genomes and
  11 in 21%.
- **The 2.8% deficit is in the scale, not in the copies.** Nine genomes whose resolved assembly holds ten complete
  copies read 9.32–9.77 on the cohort's scale (median 9.52), and the cohort's core sits at 9.72. The level rests on
  the fragment-GC model in windows of 40–60% GC; the cohort's window efficiencies fall with GC (−3.4% per 10% GC for
  the junction, −5.3% for the 45S unit), a residual of the model in repeat context. What was ruled out: the fetch
  (99.76% of the junction's reads), duplicate flags (counted), nucleotide divergence among the copies (the
  assemblies' exact-k-mer presence is 0.95 per haplotype but uncorrelated with GC, r = +0.08, while the raw estimate
  falls with GC, and a read is classified by any four of its k-mers), mosaic loss of an acrocentric (rare: 4 of 600 genomes examined), and the
  cell line's replication timing (the deficit differs between chemistries: the pilot's older libraries read the
  junction 2.1% higher than the NYGC libraries of the same people).
- **Which is more accurate?** For the 19 genomes whose assembly resolves the junction, NGS-DOSE's core level lies
  −0.08 ± 0.19 copies from the assembly's (mean ± SD; within half a copy in 19 of 19) and 78% of 20-kb blocks round
  to the same copy number, with 230 of 286 blocks the assembly puts off ten confirmed by the reads. For the other 9
  genomes the assembly is fragmented at the junction, and the reads are the more coherent reading of the locus. The
  reads' own noise is 0.14 copies per genome (the mode cluster's SD) and about 0.3 per 20-kb block; two resolved genomes
  (HG00673, HG02280) disagree with their assembly by 0.35–0.45 copies on the core for reasons neither side settles.
- **What the assemblies get wrong**, each visible in the reads: junction copies cut by a contig end (29 of 300
  copies, in 9 of 28 genomes), fragments of one copy assembled twice on small contigs (HG00658), the two haplotypes
  of a person splitting the ten copies 2 + 8, 3 + 7 or 2 + 7, and a haplotype's whole-unit k-mer median under-reading
  seven copies (HG01786). One zero-step genome (HG00642) has a fragmented junction the k-mer counts render as 8–9
  copies over its first 120 kb and 11–12 at 120–160 kb, where the reads are flat at ten.
- **What follows.** Pin the junction's scale to its core's mode (or to anchor windows chosen on replicate pairs, as
  the 45S unit's are); set the level on the core (17 of 20 blocks) and report the hyper-variable segments apart; call
  partial variants by segment with the recurrent breakpoints (22, 200, 220, 316 kb). None of this changes what the
  junction already shows about the method, that ten copies of a 400-kb acrocentric sequence are resolved to a
  single copy in every genome; it makes the number honest where the biology is not a whole number.

## 1. Why the junction, and what "ten" means

Each acrocentric short arm (13, 14, 15, 21, 22) carries an rDNA array flanked distally by the distal junction, a
sequence assembled in T2T-CHM13 on all five arms and present nowhere else. The panel is the 169,808 31-mers that
occur exactly once in each of CHM13's five junctions and nowhere in GRCh38 outside them; a read is assigned to the
class by any four of its k-mers, and its 5′ end is counted in a 250-bp window of the unit. The estimator predicts
each window's count from the library's fragment-GC response measured on single-copy control regions, and the cohort
layer (a median polish over genomes × windows) removes each window's shared efficiency and sets the scale on anchor
windows of 40–60% GC. The result is one number per genome, `DJ.cn`, expected to be 10.

Rhie et al. 2026 (bioRxiv, DJCounter) typed the same 3,202 genomes from the multiplicity of once-per-copy DJ 31-mers
against each library's two-copy k-mer peak, binned to integers with a Gaussian mixture, and examined the atypical
genomes in the HPRC release-2 assemblies. They report 9 copies in 2.8–3.4% of people and 11 or more in 8.4–9.3%, seven
genomes at 8 (one, HG01204, a G-banded Robertsonian carrier), and in the assemblies partial duplications for most of
the gains and, for the losses, three of six assembled cleanly. Their frequencies and ours agree once their integer
bins are read against our continuous values: at or below −0.7 of the level, 3.9% of the cohort here; at or above +0.3,
8.5%; and they added a mixture component between 10 and 11 for the partial duplications this document describes.

## 2. What was compared

**Genomes.** 28 cohort members with an HPRC release-2 assembly, chosen to cover the whole-copy steps (7), the values
between steps (13), the mode (7) and HG00097 as a first test: 56 haplotypes (hap1/hap2 for Hi-C-phased lines, mat/pat
for trio children).

**The assembly's reading.** Each haplotype FASTA (~900 MB) was screened for the panel's k-mers with
`ngs-dose panel --report`, which counts every k-mer of the unit in the assembly: a complete junction copy contributes
one, so five complete copies read 5 at every k-mer, a copy lacking a segment reads 4 across it, and a nucleotide
difference in one copy lowers single k-mers only. The median over the k-mers of a 5-kb sub-block is the haplotype's
copies there; the sum over the two haplotypes is the genome's. Independently, each assembly was aligned (minimap2,
asm20) to the unit with everything outside a panel k-mer masked to N, so that the repeat elements the unit shares
with the rest of the genome seed nothing and each junction copy appears as a run of alignments on one contig; every
copy's extent, contig and distance from the contig's ends are in `data/dj_hprc_copies.tsv`, classed as complete (231),
complete with internal gaps (8), distal-start (19: the copies that begin 22 kb in), truncated at a contig end (29)
or partial with an internal breakpoint (13). An assembly counts as *resolved* at the junction when no copy is cut by a
contig end inside the unit and at most one copy is partial: 19 of 28.

**The reads' reading.** The same cohort-calibrated window estimates, `C_w / exp(a_w)`, in the same 20-kb blocks;
their median over the unit is the page's `DJ.cn`. Because the cohort's scale reads 9.72 for ten copies (section 3.3),
the profiles are shown on the *ten-copy scale*, multiplied by 10 / 9.721 = 1.0287, the factor that puts the cohort's
core level at ten. The *core* is the unit without the segments that vary between people (0–20, 160–165 and 200–215 kb,
found from the cohort's own profiles: the sub-blocks whose cohort SD exceeds 1.4 times the median sub-block's).

**Reproducing it.** `pipeline/07_hprc_dj.sh` fetches, screens and aligns each haplotype (about four minutes each on a
laptop, dominated by the download) and `python3 pipeline/hprc_dj.py screen DIR -o meta/dj_hprc` writes the two
committed tables (`haplotypes.tsv`: the k-mer medians per sub-block; `copies.tsv`: the alignment catalogue). The page
reads them with `--dj-assemblies meta/dj_hprc` and writes `data/dj_hprc.tsv` (per genome), `data/dj_hprc_blocks.tsv`
(per genome and block), `data/dj_hprc_copies.tsv`, `data/dj_blocks.tsv` (every cohort genome's profile on the ten-copy
scale, with its core level) and `data/dj_segments.tsv` (the segment map), and draws `dj_assemblies.png`. The direct
tests of section 3.3 that lie outside the page are in `analysis/dj/`.

## 3. Findings

### 3.1 Whole-copy steps are real, uniform along the unit, and inherited

| genome | NGS-DOSE (step) | assembly, haplotypes | assembly core | note |
| --- | --- | --- | --- | --- |
| HG01891 | 8.62 (−1.10) | mat 4 + pat 5 | 9 | four complete copies on four maternal contigs; Rhie: DJ lost entirely (maternal) |
| HG00621 | 8.76 (−0.96) | mat 5 + pat 4 | 9 | Rhie: DJ lost entirely (paternal); a second copy lacks 240–300 kb, and the reads read 7.9–8.3 there |
| HG00658 | 8.69 (−1.03) | mat 5 + pat 4 | 9 | paternal haplotype fragmented; the reads flat at 9 |
| HG03521 | 8.89 (−0.83) | hap1 2 + hap2 7 | 9 | nine copies phased 2 + 7 |

The steps are flat along the unit: a lost junction lowers every block by one. The trios say the same. With integer
states called on the ten-copy scale (members within 0.3 of an integer; `analysis/dj/trio_integers.py`), every fully
called trio is Mendelian-consistent: 468 of 468 on the whole-unit level, 477 of 477 on the core, which calls nine
more trios because the polymorphic segments no longer pull members off an integer. Single-carrier parents transmitted
their step to 21 of 60 children (35%; fathers' losses 5 of 23, fathers' gains 5 of 14, mothers' losses 8 of 16,
mothers' gains 3 of 7), the deficit against one half the page also reports (27 of 72 by its own pairing). A parent
whose step arose in the cell line rather than the germ line would explain part of it.

### 3.2 The values between steps are partial variants

On the whole unit 8.8% of genomes lie more than 0.3 copies from an integer. Their profiles are structured, not
flat: of the 254 genomes more than 0.3 from an integer on the core, 210 have a run of blocks at one integer and the
rest at another, 29 are shifted as a whole and 15 are noisy (`analysis/dj/between_profiles.py`; 221, 27 and 9 of the
257 first selected on the whole unit). The assemblies show what the runs are.

*The recurrent gain.* In HG00146, HG00232, HG01786, HG03942, NA20752 and NA20805 the assembly holds, on one contig, a
partial copy spanning exactly the first 316 kb of the unit followed by a complete copy. The reads read 11 over blocks
0–320 kb and 10 over 320–400 kb in all six; their whole-unit values are +0.52 to +0.87. HG03654's partial copy spans
3–287 kb (Rhie: a smaller portion of the flank), and reads 11 over 20–280 kb. HG00320's spans 187–400 kb (Rhie: one arm
of the palindrome, flank intact to the rDNA) and reads 10 over 0–180 kb, 11 beyond. These are the "11 copies" of an
integer typing; they are duplications of most of a junction, and their frequency (6 of 8 of Rhie's assembled gains
share the structure) makes the 316-kb breakpoint a feature of the locus.

*Partial losses.* HG01981's maternal copy ends at 277 kb on a chromosome-scale contig, and its paternal copy at
328 kb at a contig's end; the reads read 10 over 0–280 kb and 8.4–9.1 beyond (−0.49 on the whole unit; Rhie: 9.1–9.4,
a partial loss of the flank). HG02523 has one maternal copy beginning at 122 kb; its reads read 9 over 0–120 kb.
HG00673 has a paternal copy beginning at 122 kb *and* one junction fewer; its reads read 8.6–9.1 over 0–120 kb and
9.3–9.6 beyond.

*The polymorphic segments.* Across the 3,202 profiles the cohort SD per 5-kb sub-block is 0.47 copies where only
counting noise acts and 0.8–1.3 over 0–20 kb, 0.76 at 160–165 kb and 0.9–1.0 over 200–215 kb. The first is the
distal 22 kb: 19 of the 300 assembled copies begin at 22–23 kb, in 12 of the 28 genomes, so the block reads 9 copies
in 24% of the cohort, 11 in 21%, 8 in 6% and 12 in 4%. The second is a deletion of 193–222 kb seen in one copy of
HG00097 and in one or more copies of HG01891, HG03521, NA20752 and NA20805 among others: the block reads 9 in 20% and 11 in 17%. These segments
contribute little to a whole-unit median (5% of the unit each) but everything to the block-level noise, and they are
left out of the core.

### 3.3 The level is 2.8% low, and the assemblies place the deficit in the scale

Nine genomes whose resolved assembly holds ten complete copies (HG00097, HG00438, HG00735, HG01255, HG02155, HG02258,
HG02280, HG02615, HG03239) read 9.32–9.77 on the cohort's scale (median 9.52; three of them were chosen for reading
between steps, which their partial variants outside the core explain). The cohort's core level has its median at 9.721
and its mode there. So the junction's deficit against ten is in the measurement's scale: the assemblies hold the
copies the reads do not fully count.

The scale is set by the fragment-GC model in the anchor windows (40–60% GC). Direct tests of where the residual comes
from (`analysis/dj/`):

| hypothesis | test | result |
| --- | --- | --- |
| the targeted fetch misses junction reads | DJ reads in fetch vs whole-file scan, 600 genomes | fetch holds 99.76% (range 99.64–99.83%); not it |
| duplicate flags remove piled-up reads | engine source | duplicate-flagged reads are counted (tallied apart); not it |
| divergence among the ten copies loses k-mers | per-window exact-k-mer presence in HG00097's assembly vs window GC | presence 0.95 per haplotype and flat in GC (r = +0.08; 0.93–0.96 from 35% to 60% GC), while the raw estimate of the mode genomes falls from 10.08 to 9.26 over the same bins (r = −0.26); and a read is classified by any four of its k-mers; not it |
| mosaic loss of an acrocentric chromosome in the line | per-chromosome dosage from the 800 control regions, 600 genomes | acrocentric dosage SD 0.02–0.06 copies; r with the DJ step −0.06; four genomes off by 0.2–0.9 copies (one, HG00142, with a DJ step of the same size); not the bulk |
| replication timing of the culture's DNA | the same twelve people on two chemistries | the older libraries read the junction 2.1% higher than the NYGC ones; a property of the library, not only of the DNA |
| the GC model's residual in repeat context | the cohort's window efficiencies against window GC | −0.34 log units per unit GC (r = −0.28, 1,182 windows) for the junction; −0.53 (r = −0.32) for the 45S unit, whose scale the shipped anchors correct; the surviving explanation |

The 45S unit's scale is pinned by anchor windows chosen on replicate pairs across chemistries
(`resources/GRCh38/anchors.json`); the junction has no such windows and its scale rests on the GC rule alone. Its
mode is ten copies in 87% of people, and the assemblies confirm ten in the genomes that read 9.5, so the scale can be
pinned to the core's mode, as the profiles on this page are.

### 3.4 Which is more accurate?

*Sample level.* Over the 19 resolved genomes the reads' core level lies −0.08 ± 0.19 copies from the assembly's
(within half a copy in every one; block means +0.04 ± 0.22; r = 0.95 between the assembly's mean and `DJ.cn`). Over
all 28 the SD doubles (0.38) because the fragmented assemblies scatter. The reads' precision is set by the mode
cluster's SD, 0.14 copies per genome.

*Block level.* 78% of 560 20-kb blocks round to the same copy number; the reads confirm 230 of the 286 blocks the
assembly puts off ten, and read 38 of the 274 blocks the assembly puts at ten as off by half a copy or more, most of
them in the fragmented genomes. The reads' noise per 20-kb block is about 0.3 copies (robust SD of the difference); at
5 kb it is 0.3–0.5, and a real one-copy step is visible in a single 5-kb sub-block.

*Where the reads may be wrong.* HG00673 reads 9.46 on the core where its resolved assembly holds 9, and HG02280 reads
9.64 against 10; neither has a chromosome flagged in the control regions, so a mosaic loss or gain of a whole
acrocentric is not the reason. A copy the assembly does not hold, a junction changed in part of the cell line (the
assembly's DNA and the reads' DNA are different cultures), or a genome whose level the cohort model misplaces are all
possible; the 29 flat-shifted genomes among the between-step ones are the same question at cohort scale.

*Where the assemblies are wrong* is the next section. On balance: for a genome whose junction assembled cleanly the
two agree to the noise of the reads; where it did not, the reads are the more coherent reading of the locus, and the
profile says which blocks the assembly broke. Neither gives an integer everywhere, because the biology is not an
integer everywhere.

### 3.5 Oddities in the assemblies, illuminated by NGS-DOSE

The junction sits beside the rDNA array, where long-read assemblies break, and its copies are more than 99% identical.
The catalogue of copies (`data/dj_hprc_copies.tsv`) and the profiles show:

1. **Copies cut by a contig end** (29 of 300 copies, in HG00146, HG00320, HG00642, HG00658, HG01943, HG01981, HG02040,
   HG02392, HG02523). The k-mers of the missing part are absent, and the assembly under-counts those blocks. In
   HG02523 two paternal copies begin at 86 and 109 kb at the starts of 0.36-Mb and 0.37-Mb contigs; the k-mer sums
   read 6–8 copies over 0–120 kb where the reads read 9 (the one real partial loss, a maternal copy beginning at
   122 kb inside a 2.6-Mb contig). In HG00642, a zero-step genome, the maternal 22–155 kb and paternal 3–134 kb pieces
   sit at contig ends and the 111–400 kb pieces on three 0.3–0.4-Mb contigs; the k-mer sums read 8–9 over 0–120 kb and
   11–12 over 120–160 kb where the fragments overlap, and the reads are flat at ten.
2. **Fragments assembled twice.** HG00658's paternal haplotype has a copy ending at 343 kb at a contig's end and two
   small contigs (0.17 and 0.23 Mb) each carrying the unit's last 80 kb: the k-mer sums read 11 over 320–400 kb, the
   reads 9.0, as everywhere else in this nine-copy genome. HG00146's first haplotype has a 0.05-Mb contig holding
   66–117 kb (the sums read 12 there; the reads 10.4–10.9 like the neighbouring blocks); HG00320's second, a 0.18-Mb
   contig with 368–400 kb.
3. **Phasing.** The ten copies split 2 + 8 (HG00146), 3 + 7 (HG02040), 2 + 7 (HG03521), 6 + 4 (HG00320, HG02155,
   HG02392), 4 + 6 (HG01786, NA20752): the acrocentric short arms are phased by Hi-C or trio k-mers no better than
   their near-identical sequence allows. The sum is what the reads measure, and it is right.
4. **A haplotype's whole-unit k-mer median under-reads a high copy number.** HG01786's second haplotype holds six
   complete copies and the 316-kb partial; its k-mer median over the unit is 6 while its block medians over the core
   are 7, because with seven copies nucleotide differences among them leave fewer than half the k-mers at the full
   count. An exact-k-mer count, from an assembly or from reads, needs the median per block.
5. **An unresolvable locus.** HG01943 (+1.18 on the whole unit) has its junction in ten pieces of 0.1–1.1 Mb across
   both haplotypes (3–244, 3–214, 327–400, 247–400 twice, 135–238, 133–235, 318–400, 3–85, 3–264 kb); the k-mer sums
   swing between 8 and 13 along the unit, and the reads read 11.4–12.5 over 0–260 kb (apart from the polymorphic 200–220 kb block) and 10.2–10.7 beyond: at least
   one and a half extra copies' worth of the distal 260 kb, a structure the assembler could not lay out but the
   profile describes.
6. **Inverted segments within copies** appear as mixed strands in the alignments of many complete copies (the unit
   holds palindromes) and are not counted as breaks.

### 3.6 The comparison with Rhie et al. 2026

Their eight-copy genomes (Robertsonian carriers) are the cohort's −2 steps here (7 genomes); their entire losses in
HG00621 and HG01891 are −0.96 and −1.10; HG01981's partial loss is −0.49 with the loss placed at 280–400 kb; the
duplications in NA20752, NA20805, HG03654 and HG00320 are +0.53 to +0.69 with the duplicated segment placed. Their 9.3%
"11 or more" corresponds to the 8.5% of genomes at or above +0.3 here, most of them the 316-kb duplication; their
integer bins (with an added intermediate component) and our continuous values describe the same structures. Two
things differ. They report parents typed higher than their children (p = 3.7 × 10⁻⁶); on the core level here parents
and children have the same median (9.98 and 10.01 on the ten-copy scale). And their scale is each library's two-copy
k-mer peak, which puts the mode at ten by construction; ours is the fragment-GC model's, which puts it at 9.72 and lets
the assemblies show that the deficit is the scale's.

## 4. What follows for the method

In the order of their effect on the number a user reads:

1. **Pin the junction's scale.** The level is 2.8% low on this cohort and 0.8% low on the pilot's older libraries; the
   junction is ten copies in nearly everyone and the assemblies confirm it. Either pin the class to its core's mode
   (the profiles on the page do this; it needs a cohort) or ship anchor windows for the junction chosen on replicate
   pairs across chemistries, as `anchors.json` does for the 45S unit (`pilot/evaluate_pilot.py --write-anchors` has
   the machinery), which needs no cohort.
2. **Report the core level, and the profile.** `DJ.cn_core`, the calibrated median over the 17 core blocks, is now
   in the sample table; the whole-copy steps are called on it in more trios (477 against 468 fully called at integers). The polymorphic segments (0–22 kb, 160–165 kb, 200–215 kb) belong in a table of their own states.
3. **Call partial variants by segment.** A whole-unit median reports the 316-kb duplication as +0.6 to +0.9 and a
   120-kb loss as −0.3. A segmentation of the profile with the recurrent breakpoints (22, 200, 220, 316 kb) gives
   each segment an integer and names the variant, as an assembly does; the cohort has enough carriers of each to
   define the segments once.
4. **Use the assemblies as truth, with the profile as the check.** Nineteen of 28 assemblies resolve the junction;
   for the rest the block profile says which blocks the assembly broke. A screen of all 200 assembled cohort members
   (about four minutes per haplotype, dominated by the download) would give a truth panel of 200 for this class, and
   the same screen applies to any class with a positional unit.
5. **Leave read counting alone.** The fetch reads 99.8% of the junction's reads; duplicates are counted; the
   classification tolerates the copies' divergence. Nothing in the counting path needs to change for the junction.

## 5. Data and reproducibility

| file | what |
| --- | --- |
| `meta/dj_hprc/haplotypes.tsv`, `meta/dj_hprc/copies.tsv` | the 56 screens: k-mer medians per 5-kb sub-block, and the alignment catalogue of 300 copies (`pipeline/hprc_dj.py screen`) |
| `data/dj_hprc.tsv`, `data/dj_hprc_blocks.tsv`, `data/dj_hprc_copies.tsv` | the page's comparison: per genome, per genome and block, and the copies |
| `data/dj_blocks.tsv`, `data/dj_segments.tsv` | every cohort genome's profile on the ten-copy scale with its core level; the cohort's segment map |
| `dj_assemblies.png` | the figure (`python -m report.dj_figure`) |
| `pipeline/07_hprc_dj.sh`, `pipeline/hprc_dj.py` | fetch, screen and align the haplotypes; write the tables |
| `report/dj.py`, `report/dj_page.py`, `report/dj_figure.py` | the analysis, the page section, the figure |
| `analysis/dj/` | the direct tests outside the page: fetch capture, exact-k-mer presence against GC, acrocentric dosage, trio consistency of integer calls, the between-step profile classes, the Rhie et al. text |

## 6. Limitations

The 28 genomes were chosen for what their estimates showed, so they over-represent steps and between-step values;
the comparison tests the method where it is interesting, not on a random draw. The assemblies are release-2 drafts
whose acrocentric short arms are among the least resolved sequence they contain (9 of 28 fragmented at the junction),
and the assembly's DNA is a different culture of each line from the reads' DNA, so a change in part of a line is
possible in either. The exact-k-mer screen counts a copy at every k-mer it holds verbatim; a copy diverged by more
than a few percent from CHM13 would be under-counted, and none was seen. The recurrent 316-kb breakpoint is placed to
the resolution of the masked alignment (a few kb); its sequence context is not examined here.
