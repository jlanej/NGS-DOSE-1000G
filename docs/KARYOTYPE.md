# Chromosomes in copies: a karyotype for every genome

*NGS-DOSE-1000G, 2026-10-01. Every counts file holds single-copy regions, counted by where the reads align: the
controls that every class is a ratio against, the held-out known-truth regions and, in counts made with NGS-DOSE
0.3.0's bundle, the karyotype windows, laid along every chromosome arm. NGS-DOSE 0.3.0 reads every chromosome of
every genome from them in copies, with the arms and stretches that hold another level, and keeps every level as
measured beside the whole number it is called. This document records what that finds in the 3,202 genomes of the
1000 Genomes cohort, all lymphoblastoid lines, and how far the readings can be trusted. The page
([index.html, section 3.10](index.html#karyotype)) recomputes every number here from the committed counts, except
the direct tests of `analysis/karyotype/`, which are named where they are used. The method is NGS-DOSE's
DESIGN.md, section 8.*

## Summary

- **Every genome gets a karyotype, written like one and measured by depth.** `46,XX`; `47,XY,+21`; `45,X`; a
  change in part of the cells with its share in brackets, `46,XX,-X[0.20]`; an arm, `+5q`, or a stretch with its
  span, `+3(165-194Mb)`. Of the 3,202 genomes, 2,951 are 46,XX or 46,XY throughout, 30 carry another whole number
  in every cell (a sex-chromosome complement other than XX or XY, or a chromosome or stretch at one, three or four
  copies), and 241 a chromosome or a stretch between two whole numbers: a change in part of the cells. None is
  uncertain.
- **The sex chromosomes are counted on scales of their own**, a second X reading 0.960 of the first as the inactive
  X of a growing culture does: 1,591 XY, 1,586 XX, 19 X, 3 XYY, 2 XXX and 1 XXY. The two published karyotypes in
  the cohort read as published: HG01683 47,XXY, HG03456 47,XYY. The X agrees with NGS-PCA's coverage ratio from
  the same files at r = 0.9995.
- **The windows make every chromosome as well known as every other.** In the 186 genomes counted with them, every
  autosome's level has a standard error of 0.0046 to 0.0092 copies, where the controls alone leave chromosome 19
  at 0.028 and chromosome 22 at 0.035; a gain of a whole chromosome in 5% of the cells is then found in 84 to 100%
  of genomes (66% for chromosome 22), and in 8% of the cells in all of them.
- **The errors are what they say.** z, a level's distance from its whole number in standard errors, has a robust
  SD of 1.00 over the cohort's autosomes; with the windows learned on half the genomes and read in the other half,
  1.07. Places along a chromosome put in random order give no stretch in 3,000 readings, and noise added as at
  lower depth, down to 4x, widens the errors and makes no call.
- **The calls stand against what does not depend on them.** Every whole chromosome called off shows in both arms
  alike (72 of 72, r = 0.994); the fetch of the same files gives the same karyotype for 3,202 of 3,202;
  a model learned on the first release reads the 698 genomes added later as the whole cohort's model does (691 karyotypes of 698 written alike, the rest at a threshold); and the alleles of the same reads give the same shares of the cells (r = 0.983 over 51 changes in 68 genomes), with two partial gains they do not bear out.
- **What the cohort holds is what cultures do.** Chromosomes 12 and 9 gained in part of the cells (48 and 22
  genomes), an X lost in part of the cells of 69 women and a Y in 26 men, clones that carry several chromosomes at
  one share, and stretches that recur. Two changes are shared by a child and a mother, and so were in the germ line
  or arose early.

## 1. How a chromosome is read

A region's depth, against what the genome's own fragment-GC curve expects, says how many copies of that sequence
the genome holds. A model learned on the cohort carries what stands between a region and a number of copies: the
region's efficiency (its median), the libraries' shared modes (11 components above the noise edge, learned and
scored on what the regions of a chromosome do relative to each other, so that a gained chromosome, which moves
all its regions alike, cannot be taken for a mode and learned away however common it is), the region's spread,
which is its weight, and a line across the chromosomes of each genome: in some libraries the chromosomes rich in
GC (19, 22, 17, 16) read low or high together, and each chromosome is set against what the other autosomes of
the same genome say about it. Along each chromosome a chain finds the levels; one level throughout is a whole
chromosome, a change at the centromere an arm, anything else a stretch of at least four places (a window is one
place, so one copy-number variant cannot make a stretch). A level 5 standard errors and 0.03 copies from its
whole number is fractional. The cohort says how far a level's error is from its regions' noise: 1.00 times, with
a floor of 0.0013 copies that more regions do not lower.

The windows did not exist when the cohort was counted. 186 genomes were counted again with them, by fetching from
the public CRAMs exactly the containers that hold the regions (`analysis/karyotype/slice_fetch.py`, about 450 MB
each) and counting those locally with `ngs-dose count -m fetch` (`counts_karyotype/`). The other 3,016 are read
from the 980 regions their counts hold. One model, learned on all 3,202, reads both: a region is learned on the
genomes that hold it. The 186 are 52 genomes whose earlier regions read most clearly off two copies or off XX and
XY, 120 in which nothing was called, drawn at random (30 of each sex in each of the two release batches), and 14
more chosen for a stretch or a clone that the first reading showed.

## 2. The sex chromosomes

| complement | genomes | of which |
| --- | ---: | --- |
| XY | 1,591 | |
| XX | 1,586 | |
| X | 19 | 4 with one X in every cell (three women, and NA19226, a man with no Y); 13 women whose line has lost an X in most cells, read as one X with a second in part of them (`45,X,+X[0.39]`); 2 men whose line has lost the Y in most cells (`45,X,+Y[0.31]`, `45,X,+Y[0.40]`) |
| XYY | 3 | HG03456 (published), HG01967 (a Y lost again in 30% of the cells), HG02869 |
| XXX | 2 | NA20821, and HG02805 with an X lost in 35% of the cells |
| XXY | 1 | HG01683 (published) |

In part of the cells, over all the genomes: an X lost in 69 and gained in 21, a Y lost in 26 and gained in 17
(the share of the cells from the 10th to the 90th percentile: 6 to 38% for an X lost, median 13%; 5 to 27% for a Y
lost, median 11%). These are cell lines: an X or a Y lost in culture is what most of them are.

HG01683 was found to be XXY by Richmond and colleagues in Illumina's Polaris sequencing of the same line
([PLoS Comput Biol 2021](https://doi.org/10.1371/journal.pcbi.1008815)), and HG03456 is described as XYY in the
long-read assemblies of Logsdon and colleagues ([Nature 2025](https://doi.org/10.1038/s41586-025-09140-6)); both read
so here. The pedigree's sex is not what the chromosomes say in two genomes: HG02300 holds two X and no Y under a
man's record, and NA19226 one X and no Y under a man's record (a line that lost its Y, or another person's DNA;
depth cannot tell which). NGS-PCA's coverage ratio of the whole X chromosome, from the same files by another
method, agrees with the X read here at r = 0.9995 over 3,199 genomes, and at r = 0.934 among the genomes with two
X, where the spread is mostly loss in part of the cells. A man's X is known to 0.0042 copies with the windows
(0.0059 without), a woman's to 0.0097 (0.0115), a Y to 0.0061 (0.0083).

A sample whose X and Y were off their whole numbers by the same amount in opposite directions would be marked as
two kinds of cells, or two people's DNA; none is.

## 3. The autosomes

122 genomes hold an autosome, an arm or a stretch off two copies. Whole chromosomes: 98 gained and 4 lost, nearly
all in part of the cells: chromosome 12 in 49 genomes (one, NA12739, in every cell: 47,XY,+12), chromosome 9 in
22, chromosome 11 in 9, chromosome 15 in 5, chromosome 6 in 3. Gains of 12 and 9 are what lymphoblastoid lines
acquire in culture. Where several chromosomes are gained in the same share of the cells, one clone carries them
all, and the shares agree without anything from outside: NA10843 reads chromosomes 5, 9, 10, 12 and 15 each in 0.20
to 0.22 of the cells, NA12248 four chromosomes in 0.15 to 0.22, NA21143 three in 0.15 to 0.16.

49 arms and stretches are called; the largest departures, in standard errors:

| sample | where | copies | cells | z |
| --- | --- | ---: | --- | ---: |
| NA12236 | chromosome 14 from 35 Mb to the end | 3.01 | every cell | 117 |
| NA20533 | chromosome 13 from 51 Mb to the end | 2.78 | 78% | 90 |
| HG01103 | chromosome 2, 140–179 Mb | 1.32 | 68% | −82 |
| NA12342 | chromosome 11, 102–115 Mb | 1.00 | every cell | −77 |
| NA12274 | chromosome 11 from 63 Mb to the end | 2.62 | 62% | 67 |
| NA06984 | chromosome 14 from 80 Mb to the end | 2.97 | every cell | 66 |
| HG03363 | chromosome 11, 4–41 Mb (from the end of the short arm) | 2.79 | 79% | 66 |
| NA19454 | chromosome 20, short arm | 2.78 | 78% | 65 |
| NA19355 | chromosome 2, 140–165 Mb | 1.44 | 56% | −52 |
| NA19750 | chromosome 2, 140–154 Mb | 1.42 | 58% | −50 |
| HG03061 | chromosome 3 from 165 Mb to the end | 2.91 | 91% | 49 |

Stretches recur. The long arm of chromosome 2 is lost in part of the cells of 7 genomes, three of them from the
same point (between the regions either side of 139.9 Mb) to 154, 165 and 179 Mb, in 56 to 68% of the cells; NA07055
loses 149–165 Mb in 78%. The end of the long arm of chromosome 4 is gained in 6, from 66 to 154 Mb on; chromosome
11 from 45 Mb in 5, chromosome 14 from 35 Mb in 4, chromosome 5 from 104 Mb in 3.

The distal junction, one on every acrocentric short arm, should follow a whole acrocentric chromosome. For 11
genomes with a whole acrocentric chromosome off, its level lies within 0.35 of ten plus the chromosome's excess in
nine. NA11891 reads chromosome 15 gained in 94% of the cells and the junction at 8.53 where 10.94 would follow: the
extra 15 brings no junction, and one or two are missing elsewhere; NA12058 reads 0.37 above its line.

## 4. Checks from inside the measurement

- **The two arms.** Of the 72 autosomes called off as a whole whose arms can each be read, both arms lie on the
  same side of two copies in all 72, their departures agree at r = 0.994, and the arms' difference spreads by 1.16
  of its own standard errors. Of the 32 within a tenth of a copy of two, the smallest calls made, all 32 agree in
  direction, and in 26 each arm alone is two standard errors off.
- **The quiet side.** A culture gains chromosomes and seldom loses an autosome, so the levels below two copies
  show the measurement's own scatter: of 70,430 autosomal levels at two copies, 73 lie more than 3 standard errors
  below, where a normal scatter gives 95; 284 lie as far above.
- **What cannot be real** (`checks.py`). With the places of every chromosome in random order (a window's pieces kept together), a real stretch is scattered and cannot be found: 1,500 readings of 125 clean genomes with the windows and 1,500 clean genomes without give no stretch and no whole chromosome.
- **Lower depth** (`checks.py`). Counting noise added to the 125 clean genomes with the windows as at 30x down to 4x, read against the model learned at 37x: the genomes' noise factor rises from 1.04 to 2.74, z keeps a robust SD of 0.97 to 1.02, and one chromosome is called in 750 readings; a typical autosome's standard error goes from 0.0064 copies to 0.0103 at 10x and 0.0159 at 4x.
- **The fetch.** Read from the targeted fetch of the same files against the scan's model, 3,202 of 3,202
  karyotypes are written exactly as from the scan's counts of the same regions (`checks.py`: the same, 3,202 of 3,202).
- **A model that is fixed** (`fixed_model.py`). A model learned on the 2,504 genomes of the first release and applied unchanged to the 698 added later (another sequencing batch) writes 691 of their 698 karyotypes as the whole cohort's own model does; the other seven differ in a share's last digit (HG01053 -X[0.50] against [0.49]) or by a stretch at the threshold (NA12274's 7p, NA10856's 8 Mb of chromosome 12), and the levels differ by 0.0004 to 0.0021 copies, 0.05 to 0.12 of their standard error. A genome's reading is its own, not the cohort's.
- **Read back** (`spike_in.py`). A chromosome gained or lost in a share of the cells, put into the 125 clean genomes with the windows and read back against the cohort's model, is found in 49 to 60% of them at 3% of the cells for chromosomes 1, 13 and 18, in 84 to 100% at 5% (chromosome 22, with 75 regions: 66%) and in all of them at 8%; an arm of 41 to 84 regions from 8 to 12% of the cells (18p, 28 regions: 55% at 8%, 99% at 12%); an X lost in a woman from 8% (58% at 5%), an X gained or a Y lost or gained in a man from 5% (98 to 100%). The shares read back are the shares put in. With the `screen` set a whole autosome is found from 8 to 12% of the cells (chromosome 22: 26% at 8%, 96% at 12%) and the sex chromosomes from 5 to 8%; arms of few regions are not read on their own.
- **The estimator's own flag.** Of the 72 chromosomes the estimator sets aside from its denominator (a pooled
  depth 4% off), 63 are whole chromosomes here and 9 a stretch and not the whole (HG03363 +11(4-41Mb)[0.79], …).

## 5. Against the alleles of the same reads

Depth says how many copies a chromosome is held in; the alleles of the same reads say it again, independently. Two
copies carry a heterozygous site's alleles in equal shares; a copy gained in a share f of the cells makes them
(1 + f) : 1, a copy lost 1 : (1 − f). In the reads fetched for 68 of the genomes (54 with a change or a
sex-chromosome complement other than XX and XY, 14 in which nothing was called) heterozygous sites were called on the
regions of the bundle (`analysis/karyotype/allele_calls.sh`, bcftools), and the spread of their allele fractions,
beyond what each site's depth gives by chance and beyond the genome's own chromosomes at two copies (d = 0.018 in
the median genome), was turned into a share of the cells (`allele_balance.py`; `meta/karyotype_allele_balance.tsv`).

| | changes | agreement |
| --- | ---: | --- |
| two copies expected, 5% of the cells or more, at least 30 sites | 51 | r = 0.983; robust SD of the difference 0.036; 47 within a tenth of each other |
| of which in a tenth of the cells or more | 43 | r = 0.981; robust SD 0.031 |

The changes in every cell agree: NA12236's 14q (1.01 by depth, 0.96 by the alleles), NA12739's chromosome 12 (0.99,
0.94), NA06984's 14q (0.97, 0.96). Where depth reads one copy, the sites vanish or sit at the minority's share:
NA12342's 11q loss in every cell holds no heterozygous site where the genome's density gives 91; NA07055's 2q,
lost in 78% of the cells, holds 2 where 46 are expected; NA20533's 17p, lost in 80%, holds its sites at a minor
allele of 0.20 in the median, where a loss in 80% of the cells gives 0.17. HG01258's 9q34 gain (0.13 by depth)
reads 0.19 by the alleles.

Two partial gains are not borne out by the alleles: NA12341's long arm of chromosome 4 from 114 Mb (0.19 of the
cells by depth, 0.05 by the alleles, 454 sites) and HG00703's chromosome 9 from 104 Mb (0.15 and 0.01, 210 sites).
Either the extra copies hold both homologs, or depth reads something of those libraries there as a gain; the
reads cannot say which, and gains of the end of 4q recur in six genomes of the cohort. The calls are depth's, and
these two are not confirmed. HG00623's chromosome 12 goes the other way: 0.72 of the cells by depth, 0.94 by the
alleles, more imbalance than a third copy in 72% of the cells gives; part of the cells may hold one homolog only,
twice (a change that leaves the copies as they were, which depth cannot see), or two extra copies of it. Where one X is
expected and a second is in part of the cells (four women whose line has lost an X in most cells), heterozygous
sites exist only where the second X is, and the caller misses the rarest; they read 0.07 to 0.11 above depth.
Below a twentieth of the cells the alleles cannot confirm a call: their own scatter is larger.

## 6. The windows: what they add, and what a fetch pays

The 186 genomes were read twice against the same model: with every region, and with the regions their earlier
counts hold. A chromosome's standard error is 1.7 times smaller with the windows in the median chromosome (1.4 to
3.1 from the 10th to the 90th percentile). 77 chromosomes, arms and stretches are called in both readings, at
levels that agree at r = 0.9998 and differ by 0.008 copies (robust SD), 0.54 of what their standard errors
predict. 16 are called only with the windows: a stretch of chromosome 17 lost in 80% of NA20533's cells, the short
arm of 17 lost in 31% of NA20509's, chromosome 22 gained in 16% of NA12248's (read from 3 regions without the
windows, from 75 with them), chromosome 21 gained in 7% of HG00594's, and stretches of a few Mb. One is called
only without them, HG03615's chromosome 1 at 1.95 copies, which the windows read as a stretch of its short arm.

The GC line across a genome's chromosomes is measured best with the windows: its slope, in copies per ten points
of GC, has a robust SD of 0.0077 over the 186 genomes where the fit's own error is 0.0052, and 14 lie more than
three of their standard errors from zero where chance gives 0.5. In HG00246 and HG01112, the two furthest,
chromosome 19 read 1.966 and 1.958 copies without the fit and 1.995 and 1.980 with it.

What a fetch pays, as `ngs-dose count -m fetch` reads it (median of 13 NYGC CRAM indexes, NGS-DOSE
`docs/fetch_examples.md`), and how well each control set reads the chromosomes of the 186:

| control set | regions | MB of a CRAM | SE of an autosome, copies: median (range) |
| --- | ---: | ---: | --- |
| `all`: every region | 3,247 | 404.4 | 0.0060 (0.0046–0.0092) |
| `karyotype`: 200 controls and every window | 2,647 | 277.2 | 0.0068 (0.0059–0.0092) |
| `base`: the regions before the windows | 982 | 247.9 | 0.0103 (0.0065–0.0349) |
| `screen`: 200 controls, a quarter of the autosomal windows, every sex-chromosome window | 977 | 123.4 | 0.0111 (0.0091–0.0172) |
| `lite200` | 382 | 112.9 | 0.0173 (0.0116–0.0371) |

With every region a whole chromosome gained in 5% of the cells is found in 66 to 100% of the genomes it is put
into, and in 8% of the cells in all of them; with `screen`, from 8% (26 to 100%) and in 12% of the cells in all.
Of the 73 chromosomes, arms and stretches that every region finds in a tenth of the cells or more, `karyotype` and
`base` find 67 and `screen` and `lite200` 56: what the lighter sets miss are stretches of a few Mb and arms of few
regions. The sex-chromosome complement is the same in every set but for one genome in `screen`, HG02966, whose
isodicentric Y the set cannot resolve; its note says so.

## 7. What depth does not see, and what is not yet shown

- A change that leaves the copies as they were: a balanced translocation, an inversion, both copies of a
  chromosome from one parent. A whole genome in three copies, which moves every chromosome alike.
- The acrocentric short arms, which hold no single-copy sequence (their junctions are the page's section 3.2).
- A change in fewer cells than the standard error allows is reported as its level and not called.
- Whether a change was in the donor or arose in the culture, the reads of one sample cannot say; a relative can.
  These are cell lines, and nearly all of what is found is the culture's. These are not clinical karyotypes.
- 3,016 of the genomes are read without the windows. Fetching them for all 3,202 (`pipeline/README.md`, "Chromosomes
  in copies") would read every chromosome as the 186 are read here, and let the bundle's model of the windows be
  learned on the whole cohort.
- No genome of the cohort carries a constitutional trisomy 21 or 18; what such a genome reads is the spike-ins'
  answer, not an observation.

## Reproduce

```bash
bash regenerate.sh                                        # the page: section 3.10, data/karyotype_events.tsv, data/karyotype_model.json.gz
python3 analysis/karyotype/fixed_model.py                 # a model learned on the first release, applied to the 698 added later
python3 analysis/karyotype/checks.py                      # places in random order, lower depth, the fetch
python3 analysis/karyotype/spike_in.py                    # changes put into clean genomes and read back
python3 analysis/karyotype/allele_balance.py --vcfs VCF_DIR   # the alleles of the same reads (allele_calls.sh makes the VCFs)
```
