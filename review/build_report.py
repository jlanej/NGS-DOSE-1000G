"""Builds review/NGS-DOSE_review_2026-09-23.pdf from the tables written by review_analysis.py."""
import json, os
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak,
                                KeepTogether, ListFlowable, ListItem)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "NGS-DOSE_review_2026-09-23.pdf")

rep = json.load(open(f"{ROOT}/docs/report.json"))
summ = json.load(open(f"{HERE}/review_summary.json"))
bt = pd.read_csv(f"{HERE}/batch_vs_generation.tsv", sep="\t").set_index("quantity")
rel = pd.read_csv(f"{HERE}/transmission_reanalysis.tsv", sep="\t").set_index("metric")
sx = pd.read_csv(f"{HERE}/transmission_by_sex_reanalysis.tsv", sep="\t").set_index("pair")
kt, bio, hall, repl = rep["known_truth"], rep["biology"], rep["hall"], rep["replicates"]["table"]
meta = rep["meta"]

ss = getSampleStyleSheet()
body = ParagraphStyle("body", parent=ss["BodyText"], fontName="Helvetica", fontSize=9.2, leading=12.2, spaceAfter=5)
small = ParagraphStyle("small", parent=body, fontSize=7.8, leading=10, textColor=colors.HexColor("#333333"))
cap = ParagraphStyle("cap", parent=small, spaceBefore=2, spaceAfter=10)
h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontName="Helvetica-Bold", fontSize=13.5, spaceBefore=10, spaceAfter=5)
h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontName="Helvetica-Bold", fontSize=10.5, spaceBefore=8, spaceAfter=3)
title = ParagraphStyle("title", parent=ss["Title"], fontName="Helvetica-Bold", fontSize=16, leading=20, spaceAfter=4)
cell = ParagraphStyle("cell", parent=body, fontSize=7.8, leading=9.6, spaceAfter=0)
cellb = ParagraphStyle("cellb", parent=cell, fontName="Helvetica-Bold")

P = lambda t, s=body: Paragraph(t, s)
def bullets(items, style=body):
    return ListFlowable([ListItem(P(i, style), leftIndent=10) for i in items], bulletType="bullet", start="•",
                        leftIndent=12, bulletFontSize=8)
def table(rows, widths, header=True, zebra=True):
    data = [[P(str(c), cellb if (header and i == 0) else cell) for c in r] for i, r in enumerate(rows)]
    t = Table(data, colWidths=[w * inch for w in widths], repeatRows=1 if header else 0)
    st = [("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.black),
          ("LINEBELOW", (0, -1), (-1, -1), 0.6, colors.black), ("LINEABOVE", (0, 0), (-1, 0), 0.6, colors.black),
          ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    if zebra:
        st += [("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f3f3f3")) for i in range(2, len(rows), 2)]
    t.setStyle(TableStyle(st))
    return t
f1 = lambda x: f"{x:.1f}"; f2 = lambda x: f"{x:.2f}"; f3 = lambda x: f"{x:.3f}"
def ci(r, a="R", lo="R_lo", hi="R_hi"): return f"{r[a]:.2f} ({r[lo]:.2f} to {r[hi]:.2f})"

story = []
story += [P("NGS-DOSE and NGS-DOSE-1000G: state of the findings and an independent reading", title),
          P(f"Review of the repositories <i>NGS-DOSE</i> (commit fa67723) and <i>NGS-DOSE-1000G</i> (commit 4a276d4), "
            f"with a re-analysis of the committed cohort tables. Data as of {meta['as_of']}: {meta['n']} of "
            f"{meta['total']:,} genomes and {summ['n_trios']} of 602 trios.", small), Spacer(1, 6)]

# ---------------- summary ----------------
story += [P("Summary", h1)]
story += [P(
    "NGS-DOSE measures how many copies of the 45S and 5S ribosomal DNA units, and how much of several satellite "
    "families, a person carries, using existing GRCh38 short-read CRAMs. The first quarter of the 1000 Genomes 30× "
    "cohort has been counted. My reading of the evidence is below. The measurement is technically sound. It reads "
    "sequence of known copy number correctly in every genome. A one-minute targeted fetch reproduces the "
    "whole-file scan. Between-person variation in 45S copy number is real and inherited. Independent "
    "recomputation from the committed tables reproduces every headline trio statistic.")]
story += [P(
    "What has <b>not</b> been shown is narrower than the repository's own framing sometimes suggests. The trios "
    "cannot rank estimators: within this single-pipeline cohort, the published 18S depth ratio correlates with "
    f"the calibrated estimate at r = {summ['r_cn_flat']:.3f}, and it transmits just as well. The case for the "
    "calibration therefore rests on twelve cross-chemistry replicate pairs. It holds there, but with a wide "
    "interval. The absolute scale has no orthogonal anchor. The cohort so far contains no African-ancestry genomes "
    "and only four South Asian ones.")]
story += [P(
    "The re-analysis adds one new finding, set out in section 4.2. The later 1000 Genomes release batch (698 samples) "
    "carries a technical signature <i>inside</i> the rDNA unit. GC-extreme segments read 2–5% higher relative to "
    "18S, and insert sizes are 2.5% longer. The ten parents who were sequenced in that batch show the signature "
    "as strongly as the children do. The window-efficiency model assumes one set of efficiencies per chemistry. "
    "That assumption does not hold between batches of one chemistry at one centre. This matters more for biobank "
    "use (multiple centres and years) than for the trio conclusion, which survives it.")]

story += [P("Verdict by claim", h2)]
V = [["Claim", "Status", "Basis"],
     ["Known copy numbers are read correctly", "Established",
      f"Autosomal 1.999 ± 0.008 (n = 735); male X 0.991; female Y ≤ 0.007; sex matches pedigree in 735/735"],
     ["Fetch mode equals scan", "Established", "45S fetch/scan median 0.9997 (range 0.9992–0.9999) in 735 genomes"],
     ["45S variation is real and inherited", "Established",
      f"R = {ci(rel.loc['rDNA45S.cn'])}; negative controls (truth, mtDNA, EBV) are about 0"],
     ["Same numbers as an independent pipeline", "Established (implementation)",
      f"r = {hall['flat']['r']:.3f} with Hall et al. 2021 (n = {hall['n']}); offset explained by the duplicate flag"],
     ["Calibration beats the 18S ratio", "Supported across chemistries; not shown within one",
      f"Pilot within-person CV 2.6% vs 7.0% after batch-centring (ratio 2.6, 95% CI 1.3–4.3, n = 12); trio "
      f"paired ΔR = +{summ['paired_R_cn_minus_flat'][0]:.3f} ({summ['paired_R_cn_minus_flat'][1]:.3f} to "
      f"{summ['paired_R_cn_minus_flat'][2]:.3f})"],
     ["Absolute rDNA copy number", "Not established",
      "Anchors are chosen from where chemistries agree; the 18S ratio reads about 14% higher in the same library; no ddPCR on cohort samples"],
     ["5S copy number", "Undecided", f"R = {ci(rel.loc['rDNA5S.cn'])}; the 68%-GC unit has no anchor windows"],
     ["Distal-junction whole-copy steps", "Promising",
      "Integer clustering (−2: 4, −1: 16, +1: 21); 7 of 22 transmitted (fewer than the 50% expected); 51 genomes between steps"],
     ["Satellite masses (HSat1A/3, HORs, ACRO, CER, bSat)", "Tracks assemblies per genome",
      "Per-genome log-ratio robust SD 3–7.5% against HPRC; across-person r is limited by how little people differ"],
     ["HSat2, SST1, SATR, telomeric repeat", "Relative or unvalidated", "HSat2 r = 0.30 against assemblies; SST1/SATR annotation mismatch; TEL R = 0.33"]]
story += [table(V, [1.75, 1.45, 3.85])]

# ---------------- repos & methods ----------------
story += [PageBreak(), P("1  The two repositories", h1)]
story += [P(
    "<b>NGS-DOSE</b> is the method: a Rust counting engine (<font face='Courier'>ngs-dose count</font>, about 2,800 "
    "lines, htslib) and a Python modelling layer (<font face='Courier'>ngsdose</font>, about 5,000 lines: "
    "estimate, cohort, adjust, pcsweep, trios, sinks, report, selftest). It also holds the GRCh38 resource bundle, "
    "the experimental satellite and telomere panels, the 1000 Genomes SLURM pipeline and 12-sample pilot, and a CI "
    "test suite (simulated genome with known truth, NA12878 subsample, mock trio cohort, bundle integrity). "
    "<font face='Courier'>docs/DESIGN.md</font> is a 725-line specification that records measurements on real data "
    "behind each design decision. It includes an assumption audit (section 15) in which nine assumptions were "
    "found false and fixed before the cohort run. <font face='Courier'>docs/EVIDENCE.md</font> is the argument, "
    "addressed to a sceptic.")]
story += [P(
    "<b>NGS-DOSE-1000G</b> is the results repository. It holds per-sample counts files for every genome processed "
    "(735 scans and 735 fetches, about 240 kB and 70 kB each), the pedigree, the Hall et al. 2021 table, NGS-PCA's "
    "per-sample QC, and HPRC release-2 CenSat annotations for 200 cohort members. It also holds the page, "
    "tables and PDF generated from these by <font face='Courier'>regenerate.sh</font>. The counts files contain no reads "
    "or genotypes, so every downstream number can be recomputed without touching a CRAM. That design choice is "
    "sound, and it made this review possible.")]
story += [P("2  Methods in brief", h1)]
story += [bullets([
    "<b>Counting.</b> Each read is tested against class-diagnostic 31-mers (45S unit KY962518.1; 5S unit X12811.1; "
    "a 400-kb distal-junction core). A k-mer is kept only if it occurs nowhere else in GRCh38 or T2T-CHM13. Reads "
    "are placed on the unit by their k-mer hits, not by the aligner. Fragment 5′ ends are counted. The duplicate "
    "flag is ignored in numerator and denominator alike, because MarkDuplicates under-flags collapsed repeats "
    "(rDNA 5.2% flagged against 9.1% for single-copy sequence in this cohort).",
    "<b>Two modes.</b> <i>Scan</i> reads the whole CRAM and does not depend on read placement. <i>Fetch</i> reads only "
    "the control regions plus 80 learned sink intervals (3.3 Mb), where 99.9% or more of class reads land: about 0.5 GB "
    "of a 15.8-GB file, taking roughly a minute over HTTPS.",
    "<b>Library model.</b> Per sample, a Poisson spline of fragment-end density on fragment-scale GC (window equal to "
    "the insert size, following Benjamini &amp; Speed 2012) is fitted over 800 single-copy control regions (10.1 Mb). "
    "The copy number of a unit window is 2 × observed / expected.",
    "<b>Unit calibration.</b> Segments of the rDNA unit drop out beyond what any genome-wide GC curve predicts, and "
    "which segments drop out depends on the chemistry. The model log C<sub>iw</sub> = c<sub>i</sub> + a<sub>w</sub> + "
    "e<sub>iw</sub> is fitted by median polish, with the scale pinned on 13 anchor windows (3.25 kb) on which "
    "NovaSeq and two HiSeq chemistries agree in the same individuals.",
    "<b>Known truths in every sample.</b> These are 80 held-out autosomal regions (2 copies), 60 chrX and 40 chrY "
    "regions, and the distal junction (10 copies, one per acrocentric arm), each measured by the same code path. "
    "mtDNA and EBV copies per cell are recorded as culture covariates.",
    "<b>Cohort layer.</b> Coverage-PC adjustment uses a Marchenko–Pastur edge fitted from its square-root shape, "
    "with a cross-validated sweep against the truths and trios. Trio transmission reliability is estimated from "
    "the midparent slope, corrected for spousal correlation and centred within population, with family-bootstrap "
    "intervals, a permuted-family null and a paired bootstrap between estimators."])]

story += [PageBreak(), P("3  Results as reported", h1)]
story += [Image(f"{ROOT}/docs/evidence.png", width=7.0 * inch, height=7.0 * inch * 1216 / 2560),
          P("<b>Figure 1.</b> The repository's evidence figure (docs/evidence.png, reproduced unchanged), 735 genomes. "
            "(a) known copy numbers; (b) distal-junction steps; (c) child against midparent, 149 trios; (d) same person "
            "on two technologies, 12 pilot pairs; (e) residual GC bias of the 18S depth ratio; (f) against Hall et al. "
            "2021; (g) satellite mass against HPRC assemblies; (h) fetch against scan.", cap)]
T = [["Quantity", "Value", "n"],
     ["Held-out autosomal (truth 2)", f"{kt['auto']['mean']:.3f} ± {kt['auto']['sd']:.3f}", kt['auto']['n']],
     ["chrX in men (1) / women with an intact culture (2)", f"{kt['chrX']['M']['mean']:.3f} / {kt['sex']['women_intact']['mean']:.3f}",
      f"{kt['chrX']['M']['n']} / {kt['sex']['women_intact']['n']}"],
     ["chrY in men with an intact Y (1) / women (0)", f"{kt['sex']['men_intact_Y']['mean']:.3f} / {kt['chrY']['F']['mean']:.4f}",
      f"{kt['sex']['men_intact_Y']['n']} / {kt['chrY']['F']['n']}"],
     ["Distal junction (truth 10): median, robust SD", f"{kt['DJ']['median']:.2f}, {kt['DJ']['mad_sd']:.2f}", kt['DJ']['n']],
     ["Mosaic X loss / Y loss (LCL)", f"{kt['sex']['n_mosaic_X']} / {kt['sex']['n_mosaic_Y']}", 735],
     ["45S calibrated copy number: mean (CV); 10th–90th percentile", f"{rep['rdna']['rDNA45S.cn']['mean']:.0f} ({100*bio['cn45_cv']:.0f}%); "
      f"{rep['rdna']['rDNA45S.cn']['q10']:.0f}–{rep['rdna']['rDNA45S.cn']['q90']:.0f}", 735],
     ["18S depth ratio (published estimator), mean", f"{rep['rdna']['rDNA45S.18S.flat']['mean']:.0f} (slope against calibrated {bio['gc_bias']['flat_vs_cn']['slope']:.2f})", 735],
     ["5S calibrated copy number, mean", f"{rep['rdna']['rDNA5S.cn']['mean']:.0f}", 735],
     ["Trio reliability, 45S calibrated", ci(rel.loc['rDNA45S.cn']), 149],
     ["Pilot test–retest ICC: calibrated / 18S ratio / 18S ratio offset-removed",
      f"{repl['calibrated']['icc']:.3f} / {repl['flat']['icc']:.2f} / {repl['flat_centred']['icc']:.2f}", 12],
     ["18S ratio / calibrated vs library GC bias (r); same under GC model",
      f"{bio['gc_bias']['flat_vs_gc']['r']:.2f}; {bio['gc_bias']['modelled_vs_gc']['r']:.2f}", 735],
     ["Hall et al. 2021: r; their/our ratio; after their duplicate rule",
      f"{hall['flat']['r']:.3f}; {hall['flat_ratio']:.3f}; {hall['dup_corrected_ratio']:.3f}", hall['n']],
     ["mtDNA copies against NGS-PCA QC (r)", f"{rep['ngspca_qc']['mtdna']['r']:.3f}", 735],
     ["rDNA in HPRC r2 assemblies, as fraction of the implied array", "0.38 (0.22–0.62); longest piece ≤ 1.1 Mb", 54],
     ["45S vs 5S; 45S vs mtDNA (r)", f"{bio['cn45_vs_5S']['r']:.2f} ({bio['cn45_vs_5S']['r_lo']:.2f} to {bio['cn45_vs_5S']['r_hi']:.2f}); "
      f"{bio['cn45_vs_chrM']['r']:.2f} ({bio['cn45_vs_chrM']['r_lo']:.2f} to {bio['cn45_vs_chrM']['r_hi']:.2f})", 735]]
story += [P("<b>Table 1.</b> Headline numbers, read from docs/report.json and checked against the re-analysis.", small),
          table(T, [3.6, 2.7, 0.75])]
story += [Spacer(1, 4), P(
    f"Processing so far follows sample-ID order. The superpopulations done are AMR {meta['by_superpop']['AMR']['done']}/"
    f"{meta['by_superpop']['AMR']['total']}, EAS {meta['by_superpop']['EAS']['done']}/{meta['by_superpop']['EAS']['total']}, "
    f"EUR {meta['by_superpop']['EUR']['done']}/{meta['by_superpop']['EUR']['total']}, SAS {meta['by_superpop']['SAS']['done']}/"
    f"{meta['by_superpop']['SAS']['total']} and AFR {meta['by_superpop']['AFR']['done']}/{meta['by_superpop']['AFR']['total']}. "
    "All 735 genomes were produced by one engine build, one bundle and one sink set, and all have end-of-file markers.", small)]

# ---------------- independent reanalysis ----------------
story += [PageBreak(), P("4  Independent re-analysis", h1)]
story += [P(
    "Everything in this section was recomputed from docs/data/cohort.tsv, the pedigree and the NGS-PCA QC table "
    "(which carries the release batch). The code is review/review_analysis.py, and its tables are written alongside "
    "this PDF. Nothing was re-counted from CRAMs, and the engine and test suite were not re-run.")]
story += [P("4.1  The trio statistics reproduce", h2)]
TR = [["Metric", "R (95% CI)", "R, children rescaled", "child SD / parent SD"]]
for m in rel.index:
    r = rel.loc[m]
    TR.append([r.label, ci(r), ci(r, "R_rescaled", "R_rescaled_lo", "R_rescaled_hi"), f2(r.sd_child_over_parent)])
story += [table(TR, [2.5, 1.55, 1.55, 1.2])]
story += [Spacer(1, 4), Image(f"{HERE}/fig_transmission.png", width=7.0 * inch, height=7.0 * inch * 3.3 / 7.2),
          P("<b>Figure 2.</b> (a) Transmission reliability (midparent slope, spousal-corrected, values centred within "
            "population; 4,000 family bootstraps) for 149 complete trios. Squares rescale the children to their parents' "
            "spread, which removes a between-batch scale difference. (b) Single-parent slopes for calibrated 45S by the "
            "sex of parent and child. Under additive inheritance each is expected near (R + ρ)/2 ≈ 0.6.", cap)]
story += [P(
    f"My estimates match the repository's to the second decimal: R = {rel.loc['rDNA45S.cn','R']:.3f} for calibrated 45S "
    f"and a paired difference from the 18S ratio of +{summ['paired_R_cn_minus_flat'][0]:.3f} "
    f"({summ['paired_R_cn_minus_flat'][1]:.3f} to {summ['paired_R_cn_minus_flat'][2]:.3f}). All three negative controls "
    "sit at zero. The lower bound of about 0.91 means measurement error is at most a few percent of a person's value. "
    "It does not mean error is absent. As the repository says, with a between-person CV of 22%, every competent "
    "estimator reaches R ≈ 1. The trio test establishes that the variation is genetic. It does not say which "
    "estimator is best.")]
story += [P(
    "The ‘heritability is one by construction’ premise deserves a qualification. rDNA arrays are reported to "
    "rearrange frequently in meiosis (Stults et al., <i>Genome Res</i> 2008), and LCL cultures can change arrays "
    "somatically. Both would lower R without being measurement error. R ≈ 1 therefore bounds error and "
    "array instability together, and it suggests that net dosage changes in transmission are small relative to "
    "between-person variation. That is a biological finding in its own right, worth stating.")]

story += [P("4.2  New: the later release batch leaves a signature inside the rDNA unit", h2)]
story += [Image(f"{HERE}/fig_batch_signature.png", width=7.0 * inch, height=7.0 * inch * 2.9 / 7.2),
          P("<b>Figure 3.</b> Within-sample ratios of rDNA segments, and the library insert size, for parents "
            "sequenced in the original 2,504-sample batch (n = 289), the ten parents sequenced in the 698-sample batch, "
            "and children (all in the 698 batch, n = 151). Bars mark medians. Effects in the titles come from one "
            "regression per quantity on batch, generation and population (n = 735), so batch and generation are "
            "estimated jointly. The ten 698-batch parents separate the two.", cap)]
B = [["Quantity (log scale)", "Batch effect (698 vs 2504)", "p", "Generation effect", "p"]]
for q in ["28S / 18S (GC model)", "28S / 18S (no GC model)", "5'ETS / IGS (GC model)", "18S depth ratio / calibrated 45S",
          "insert-size median", "library rate at 65% GC", "calibrated 45S copy number", "18S depth ratio (published)",
          "held-out autosomal (truth 2)", "distal junction (truth 10)"]:
    r = bt.loc[q]
    B.append([q, f"{r.batch_effect_pct:+.1f}% ± {r.batch_se_pct:.1f}", f"{r.batch_p:.2g}",
              f"{r.generation_effect_pct:+.1f}% ± {r.generation_se_pct:.1f}", f"{r.generation_p:.2g}"])
story += [table(B, [2.45, 1.55, 0.6, 1.5, 0.6])]
story += [Spacer(1, 4), P(
    "The repository reports that children read 4.8% more 45S than their parents. It attributes this to "
    "‘generation or batch’ and handles it by rescaling. The within-sample ratios decide part of that question. "
    "Ratios between segments of the same unit cancel copy number, so a shift in them can only be technical. The "
    "28S/18S ratio (+5.4% without the GC model, +2.0% with it), the 5′ETS/IGS ratio (+3.3%) and the 18S-ratio-to-"
    "calibrated ratio (+2.5%) all move with the batch. None of these ratios moves with generation "
    "(generation p ≥ 0.17 for each). The ten parents sequenced in the later batch carry the full shift. Insert size shows the same "
    "pattern (+2.5%, p = 8 × 10<super>−5</super>). The known truths do not move, so the denominator is unaffected.")]
story += [P(
    "Three implications follow. First, window efficiencies a<sub>w</sub> differ between batches of one chemistry "
    "at one centre, and the ‘per library type’ calibration absorbs only their average. Second, part of the "
    "child–parent offset in total 45S is probably technical as well. With only ten later-batch parents the split "
    "cannot be made for the total: batch +3.0% ± 6.7, generation +1.7% ± 6.8. The published 18S ratio is more "
    "exposed to the shift than the calibrated estimate. Third, for UK Biobank, All of Us or any multi-centre "
    "cohort, efficiencies should be fitted per sequencing batch or centre, or batch should be carried as a "
    "covariate in the unit model. The trio result is unaffected: a batch offset moves the intercept, and the "
    "scale difference is small (child/parent SD 1.02).")]

story += [P("4.3  Where calibration earns its keep", h2)]
story += [P(
    f"Within the cohort the calibrated estimate and the 18S depth ratio correlate at r = {summ['r_cn_flat']:.3f}. "
    "The 18S ratio carries residual GC bias (its ratio to the calibrated estimate follows the library's 65%-GC rate "
    "at r = 0.84), but the SD of that log ratio is only 3.3%, small beside 22% between people. The pilot replicates "
    "are where the methods separate. Removing the 27% offset between technologies from the 18S ratio (a batch "
    "correction any analyst would apply) leaves a within-person CV of "
    f"{100*summ['replicate_within_cv_flat_centred']:.1f}%, against {100*summ['replicate_within_cv_cal']:.1f}% for the "
    f"calibrated estimate. The ratio is {summ['replicate_cv_ratio_ci'][1]:.1f}× (bootstrap 95% CI "
    f"{summ['replicate_cv_ratio_ci'][0]:.1f}–{summ['replicate_cv_ratio_ci'][2]:.1f}, n = 12 pairs). The advantage is "
    "real but rests on twelve people. It is the most important number to firm up, for example with more HGSVC/"
    "Platinum replicates or 1000 Genomes samples resequenced on other platforms.")]

story += [P("4.4  Leads to watch, not yet findings", h2)]
story += [bullets([
    f"<b>Mother-to-son transmission.</b> Slopes are {sx.loc['father to son','slope']:.2f}, "
    f"{sx.loc['father to daughter','slope']:.2f}, {sx.loc['mother to son','slope']:.2f} and "
    f"{sx.loc['mother to daughter','slope']:.2f} (father–son, father–daughter, mother–son, mother–daughter). The "
    "repository's permutation heterogeneity test gives p = 0.009 for 45S, but it was one of 23 metrics tested. rDNA "
    "is autosomal and no mechanism predicts this particular pattern, so it should be treated as a pre-specified "
    "test for 602 trios rather than as a result.",
    f"<b>Population differences.</b> In the 2,504 batch, population explains {100*summ['pop_R2']:.1f}% of log 45S "
    f"variance (p = {summ['pop_F_p']:.1g}, n = {summ['pop_n']}). European means are the lowest (about 430) and East Asian "
    "the highest (about 480). African genomes, which carry the most sequence diversity and are the hardest test of "
    "the k-mer panels and sinks, have not been counted.",
    "<b>The distal junction reads 9.66, not 10.</b> The repository proposes that late-replicating sequence is "
    "under-represented in cycling cultures (the S-phase hypothesis). The prediction that DJ and the female X "
    f"move together is not borne out (r = {bio['DJ_vs_chrX_female']['r']:.2f}, {bio['DJ_vs_chrX_female']['r_lo']:.2f} to "
    f"{bio['DJ_vs_chrX_female']['r_hi']:.2f}). DJ does correlate weakly with the first control PC (r ≈ 0.18 in the 2,504 batch). "
    "Sequence divergence of the DJ core from the CHM13-derived k-mers in some individuals is an alternative that "
    "has not been excluded. The older-library pilot reading of 9.93 argues for a culture or library effect.",
    "<b>Junction steps are transmitted less than half the time</b> (7 of 22, two-sided p = 0.13), and 51 genomes sit "
    "between integer states. Some steps may therefore be somatic changes in the cell lines. With 602 trios this becomes "
    "a clean Mendelian test.",
    "<b>Satellite correlations need context.</b> HSat1B's across-person r = 0.99 is largely a sex effect (the family "
    "is mostly on Yq; CV across people 0.90). α-satellite HOR mass differs so little between people (CV about 5%) "
    "that its transmission reliability is 0.43, and it is best treated as a QC quantity. HSat2 is inherited "
    "(R = 1.07 rescaled) but agrees poorly with assemblies (r = 0.30), so what it measures is heritable but not "
    "yet shown to be HSat2 mass."])]
story += [P("4.5  Biology that the cohort already speaks to", h2)]
story += [P(
    f"45S and 5S copy numbers are uncorrelated (r = {bio['cn45_vs_5S']['r']:.2f}) in the same people and libraries, "
    "and 45S is unrelated to mitochondrial copies per cell (r = 0.02). The first agrees with Hall et al. 2021 and with the "
    "UK Biobank analysis of Raj et al. (medRxiv, January 2026), who describe the two arrays as varying "
    "independently. The second does not support the inverse rDNA–mtDNA relation of Gibbons et al. 2014, at least "
    "in lymphoblastoid lines, where mtDNA content is a culture property with no transmission (R ≈ 0.06). These are "
    "among the more publishable observations so far, because the known truths and trios control them.")]

# ---------------- interpretation ----------------
story += [P("5  My take", h1)]
story += [P(
    "This is careful work. It checks its own claims more rigorously than most published rDNA copy-number methods. "
    "The strongest design decisions are to count from positions rather than summaries (counts files can be "
    "re-modelled without CRAMs), to measure known copy numbers in every sample by the same code path, and to "
    "make the fetch/scan equivalence empirical rather than assumed. The duplicate-flag finding (rDNA is "
    "under-flagged, so depth tools that drop duplicates over-read it by a library-specific factor) is simple, "
    "general and useful to anyone using mosdepth-type depth for repeats. The assumption audit in DESIGN.md "
    "section 15 is a good model for how to document a pipeline.")]
story += [P(
    "The main reservation concerns how strongly the conclusions are worded. ‘Reliability 1’ is a lower bound of "
    "0.91, and trios cannot distinguish the method from its simplest competitor within one batch. The decisive "
    "comparison is the cross-chemistry replicate, which rests on n = 12. The absolute scale depends on choosing "
    "windows where chemistries agree. That is a reasonable heuristic, but three Illumina chemistries could share a "
    "bias, and the 18S ratio and the calibrated estimate differ by about 14% in the same library. Until ddPCR (or a "
    "characterised genome on the same instrument) is compared, absolute values such as ‘460 copies’ should be "
    "reported as calibrated units rather than copies. The batch signature in section 4.2 shows that the "
    "calibration's key assumption, one set of efficiencies per chemistry, is already violated within this cohort. "
    "The fix is straightforward.")]
story += [P(
    "For its stated purpose, a biobank-scale measurement of relative rDNA dosage, the method is ready to be tried. "
    "The value it adds over the 18S ratio will scale with how heterogeneous the sequencing is. The value "
    "that is unique to it (known truths, the fetch mode, the 5S and junction classes, the satellites) does not depend "
    "on that comparison. One competitor, Raj et al. 2026, is described as using a closely related fragment-GC approach (not yet verified against their methods) on 490,383 "
    "UK Biobank genomes, with code announced at github.com/calico/rDNA_CN. A head-to-head comparison on these "
    "1000 Genomes samples, which they also analysed, would be the most informative single external check.")]
story += [P(
    "A note on independence. The repositories credit their method design and implementation to an AI system working "
    "with the user, and this review comes from the same kind of system. I recomputed the statistics from the committed "
    "tables and looked for weaknesses, but that is not independent human verification. The code, the statistical "
    "arguments and the interpretation still need to be checked by a person before publication, as the NGS-DOSE README "
    "itself states.", small)]

story += [P("6  Recommended next steps, in priority order", h1)]
story += [ListFlowable([ListItem(P(t), leftIndent=12) for t in [
    "<b>Finish the cohort, starting with the African and South Asian superpopulations</b> (893 and 597 genomes "
    "outstanding). Check sink capture, k-mer hit fractions and DJ level by ancestry before trusting cross-population "
    "comparisons.",
    "<b>Make the unit calibration batch-aware.</b> Fit a<sub>w</sub> per release batch, or add batch to the "
    "median polish, and report how the child–parent offset, R and the 18S comparison change. Use the full cohort's "
    "non-child members of the 698 batch as the bridge.",
    "<b>Obtain an orthogonal absolute anchor.</b> Run ddPCR for 18S/28S/IGS on a few dozen Coriell LCLs already in "
    "the cohort, ideally including the pilot trios, or sequence CHM13/HG002 on NovaSeq alongside published ddPCR.",
    "<b>Pre-register the 602-trio tests.</b> These are the paired estimator comparison, the 5S reliability, DJ step "
    "Mendelian transmission with de novo rate, the four sex pairings, and the S-phase prediction. That keeps the "
    "final analysis confirmatory.",
    "<b>Enlarge the replicate set.</b> Every 1000 Genomes sample with an independent older library (HGSVC, "
    "Platinum, the phase-3 low-coverage data at reduced precision) adds to the one comparison that ranks estimators.",
    "<b>Test DRAGEN alignments.</b> UK Biobank and All of Us are DRAGEN-aligned. Scanning a handful of DRAGEN CRAMs "
    "and re-learning the sinks is the gate to biobank use.",
    "<b>Compare against a blood-derived genome.</b> Even a few people with both blood and LCL sequencing would put a "
    "number on the culture effects (X/Y loss, S-phase, EBV) that LCL-only data cannot separate from biology.",
    "<b>Compare head-to-head with calico/rDNA_CN</b> once the code is public, and settle a licence for NGS-DOSE."]],
    bulletType="1", leftIndent=14)]

story += [P("Appendix: files in this directory", h2)]
story += [table([["File", "Contents"],
                 ["review_analysis.py", "Re-analysis: batch/generation regressions, trio reliabilities, sex-specific slopes, figures"],
                 ["build_report.py", "Builds this PDF from the tables below and docs/report.json"],
                 ["batch_vs_generation.tsv", "Section 4.2 regression table (population-adjusted)"],
                 ["transmission_reanalysis.tsv", "Section 4.1 reliabilities with bootstrap intervals"],
                 ["transmission_by_sex_reanalysis.tsv", "Section 4.4 single-parent slopes by sex"],
                 ["review_summary.json", "Scalars quoted in the text (paired ΔR, replicate CVs, population R²)"],
                 ["fig_transmission.png, fig_batch_signature.png", "Figures 2 and 3"]], [2.6, 4.4])]
story += [Spacer(1, 4), P("Regenerate from the repository root: <font face='Courier'>python review/review_analysis.py &amp;&amp; "
                          "python review/build_report.py</font> (needs pandas, statsmodels, matplotlib, reportlab).", small)]

def on_page(c, d):
    c.setFont("Helvetica", 7); c.setFillColor(colors.grey)
    c.drawString(0.75 * inch, 0.5 * inch, "NGS-DOSE review · data as of 2026-09-23 · 735 genomes, 149 trios")
    c.drawRightString(7.75 * inch, 0.5 * inch, str(d.page))
doc = SimpleDocTemplate(OUT, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                        topMargin=0.7 * inch, bottomMargin=0.75 * inch, title="NGS-DOSE review", author="")
doc.build(story, onFirstPage=on_page, onLaterPages=on_page)
print(OUT)
