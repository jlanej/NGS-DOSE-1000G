#!/usr/bin/env python
"""Builds assembly_rdna/HPRC_r2_rDNA_report.pdf from tables/tests.json, the tables and figures/."""
import json, os
import numpy as np
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak, ListFlowable, ListItem

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB, FIG = os.path.join(BASE, "tables"), os.path.join(BASE, "figures")
OUT = os.path.join(BASE, "HPRC_r2_rDNA_report.pdf")
T = json.load(open(f"{TAB}/tests.json"))
dip = pd.read_csv(f"{TAB}/person_vs_ngsdose.tsv", sep="\t", index_col=0)
ref = T["references"]
ph = T["per_haplotype"]
hp = T["hap_parent"]
asof = open(f"{BASE}/../docs/report.json").read().split('"as_of": "')[1][:10]
n_cohort = sum(1 for _ in open(f"{BASE}/../docs/data/cohort.tsv")) - 1

ss = getSampleStyleSheet()
body = ParagraphStyle("b", parent=ss["BodyText"], fontName="Helvetica", fontSize=9.2, leading=12.2, spaceAfter=5)
small = ParagraphStyle("s", parent=body, fontSize=7.8, leading=10, textColor=colors.HexColor("#333333"))
cap = ParagraphStyle("c", parent=small, spaceBefore=2, spaceAfter=10)
h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontName="Helvetica-Bold", fontSize=13.5, spaceBefore=10, spaceAfter=5, keepWithNext=1)
h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontName="Helvetica-Bold", fontSize=10.5, spaceBefore=8, spaceAfter=3)
title = ParagraphStyle("t", parent=ss["Title"], fontName="Helvetica-Bold", fontSize=16, leading=20, spaceAfter=4)
cell = ParagraphStyle("cell", parent=body, fontSize=7.8, leading=9.6, spaceAfter=0)
cellb = ParagraphStyle("cellb", parent=cell, fontName="Helvetica-Bold")
P = lambda t, s=body: Paragraph(t, s)
def bullets(items, style=body):
    return ListFlowable([ListItem(P(i, style), leftIndent=10) for i in items], bulletType="bullet", start="•", leftIndent=12, bulletFontSize=8)
def table(rows, widths):
    data = [[P(str(c), cellb if i == 0 else cell) for c in r] for i, r in enumerate(rows)]
    t = Table(data, colWidths=[w * inch for w in widths], repeatRows=1)
    st = [("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEABOVE", (0, 0), (-1, 0), 0.6, colors.black),
          ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.black), ("LINEBELOW", (0, -1), (-1, -1), 0.6, colors.black),
          ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    st += [("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f3f3f3")) for i in range(2, len(rows), 2)]
    t.setStyle(TableStyle(st))
    return t
def rr(d, k=2): return f"{d['r']:.{k}f} ({d['lo']:.{k}f} to {d['hi']:.{k}f})"
def fig(name, w, h):
    return Image(f"{FIG}/{name}", width=w * inch, height=h * inch)

n18, n5, ndj = ph["n_18S"], ph["n_5S_units"], ph["n_DJ_copies"]
sides = T["array_sides"]; tot_sides = sum(sides.values())
place = T["units_by_place"]; tot_units = sum(place.values())
r45, r5, rc, rdj = T["r_45S"], T["r_5S"], T["r_censat"], T["r_DJ"]
q45, q5 = T["ratio45"], T["ratio5"]
ceil = T["hap_parent_ceiling"]
c13 = ref["chm13v2.0_maskedY_rCRS"]; hgp, hgm = ref["hg002v1.1.pat"], ref["hg002v1.1.mat_MT"]

S = []
M = json.load(open(f"{TAB}/method_accuracy.json"))
inh = {(r["cls"], r["parents"]): r for r in M["inheritance"]}
i45, i45h, i5 = inh[("45S", "NGS-DOSE")], inh[("45S", "Hall 2021")], inh[("5S", "NGS-DOSE")]
DJ, DD, AR = M["distal_junction"], M["ddpcr"], M["arrays"]
dn, dc, dfl = DD["ngsdose"], DD["conkord_same9"], DD["ratio18S_flat"]
dg, da, dc12 = DD["ngsdose_google"], DD["ngsdose_any"], DD["conkord_same12"]
RC = DD["r_ci"]
PCMP = pd.read_csv(f"{TAB}/potapova_comparison.tsv", sep="\t").set_index("sample")
NV = M["novaseq"]
SC = json.load(open(f"{BASE}/novaseq/scan_check.json")) if os.path.exists(f"{BASE}/novaseq/scan_check.json") else None
S += [P("HPRC release-2 assemblies, the rDNA, and how NGS-DOSE compares", title),
      P(f"Direct sequence tests on all {T['n_haplotypes']} HPRC/HPP release-2 haplotype assemblies (CHM13, T2T-HG002 v1.1, HG06807 and "
        f"GRCh38 as references), set against NGS-DOSE short-read copy numbers of the same people and their parents, the independent "
        f"Hall et al. 2021 pipeline, and the ddPCR and FISH measurements of Potapova et al. 2025. NGS-DOSE counts as of {asof} "
        f"({n_cohort:,} genomes; {T['n_people_counted']} assembled people and {T['n_trios']} trios), plus 16 genomes counted for this report: 7 NYGC 1000 Genomes CRAMs of ddPCR-measured people, and 6 genomes from a second NovaSeq "
        f"pipeline (Google Health's GIAB set) with the NYGC CRAMs of the 3 people sequenced by both.", small),
      Spacer(1, 6)]
S += [P("Summary", h1)]
S += [bullets([
    f"<b>The assemblies contain rDNA, but not the arrays.</b> All {T['n_hifiasm']} hifiasm haplotypes hold 45S units (median "
    f"{n18['50%']:.0f} 18S genes), in a median of {ph['n_arrays']['50%']:.0f} pieces of 1–10 units. {100*sides.get('contig_end',0)/tot_sides:.0f}% "
    f"of piece ends are contig ends and {100*place.get('unplaced contig',0)/tot_units:.0f}% of units are on unplaced contigs. They hold "
    f"{100*q45['50%']:.0f}% of the copies short reads count ({100*DD['assembly_18S']['median_ratio']:.0f}% of the ddPCR copies), and their "
    "per-chromosome distribution does not follow FISH. The 5S array, by contrast, is assembled whole in "
    f"{100*T['frac_hap_5S_single_closed']:.0f}% of haplotypes.",
    f"<b>For 45S, NGS-DOSE is the more accurate method.</b> Measured against the parents, which no error in the child's measurement can "
    f"imitate, the child's NGS-DOSE value tracks inheritance at r = {i45['r_ngsdose']:.2f} and the assembly's at {i45['r_assembly']:.2f} "
    f"(paired difference +{i45['diff']:.2f}, 95% CI {i45['diff_lo']:.2f} to {i45['diff_hi']:.2f}; {i45['n']} trios). The result is the same "
    f"with the parents measured by Hall et al. (+{i45h['diff']:.2f}). Against ddPCR, an orthogonal assay, NGS-DOSE correlates at "
    f"r = {dn['r']:.2f} (n = {dn['n']}), ahead of the published 18S depth ratio ({dfl['r']:.2f}) and the CONKORD k-mer pipeline ({dc['r']:.2f}). "
    f"Adding HG002, HG003 and HG004 from a second NovaSeq pipeline gives r = {da['r']:.2f} over all {da['n']} ddPCR-measured people "
    f"(CONKORD {dc12['r']:.2f}).",
    f"<b>The estimate transfers between sequencing pipelines; the 18S depth ratio does not.</b> Three people sequenced by both NYGC and "
    f"Google Health on NovaSeq 6000, with different inserts and GC bias, get 45S values that differ by {NV['rDNA45S.cn']['mean_abs_pct']:.1f}% "
    f"on average with NGS-DOSE and by {NV['rDNA45S.18S.flat']['mean_abs_pct']:.0f}% with the 18S depth ratio computed from the same reads.",
    f"<b>For 5S the two methods agree, so both are validated.</b> The assembled 5S unit count equals NGS-DOSE's estimate (ratio "
    f"{q5['50%']:.2f} ± {q5['std']:.3f}, r = {r5['r']:.3f}), and the inheritance test finds no difference between them "
    f"({i5['diff']:+.2f}, {i5['diff_lo']:.2f} to {i5['diff_hi']:.2f}). This also shows the inheritance test does not favour short reads "
    "where the assembly is complete.",
    f"<b>NGS-DOSE's weakness is absolute scale, not precision.</b> It reads about {abs(100*(dn['median_ratio']-1)):.0f}% below ddPCR "
    f"(median ratio {dn['median_ratio']:.2f}) and {100*(1-DJ['ngsdose_mean']/10):.1f}% below the known 10 copies of the distal junction "
    f"({DJ['ngsdose_mean']:.2f}). Both shortfalls point the same way, which suggests one shared cause. For the distal junction the "
    "assemblies count whole copies exactly where their short arms are complete; NGS-DOSE detects single-copy steps, and all "
    f"{DJ['ngs_steps']} of its steps recur in the assemblies in the same direction ({DJ['ngs_steps_confirmed']} of the same size)."])]
S += [PageBreak(), P("The comparison at a glance", h2)]
S += [fig("fig7_at_a_glance.png", 6.1, 6.1 * 10.6 / 7.4),
      P("<b>Overview figure.</b> Every comparison in this report on one page, with the measures defined underneath. (a) Each method's "
        "45S estimate against ddPCR for the 12 ddPCR-measured lines: NGS-DOSE on the NYGC 1000 Genomes reads (filled squares) and on Google's "
        "GIAB NovaSeq reads (open), CONKORD, the 18S depth ratio, and the HPRC assembly (5 people). (b) The same as estimate ÷ ddPCR, one point "
        "per person. (c) Child's 45S against its parents' mean in 36 trios, the child measured by NGS-DOSE (squares) and by its assembly "
        "(circles), both as deviations from the superpopulation mean; lines are least-squares fits. (d) The inheritance correlations with "
        "bootstrap 95% CIs, for 45S and 5S. (e) Pearson r with ddPCR and 95% CI (Fisher z), over all 12 lines and over the 5 with an "
        "assembly. (f) Assembled copies against NGS-DOSE in the same 57 people, for 45S (filled) and 5S (open). (g) Distal-junction copies, "
        "known to be 10 per diploid genome. (h) The same three people counted from two sequencing pipelines. Figures 2 and 4–6 show each "
        "test in full.", cap)]
SCT = M["scorecard"]
S += [table([["Test", "n", "NGS-DOSE", "CONKORD", "18S ratio", "Assembly", "Best"]] +
            [[r['test'], r["n"], r["ngsdose"], r["conkord"], r["ratio18S"], r["assembly"], r["best"]] for r in SCT],
            [2.05, 0.7, 0.82, 0.82, 0.7, 0.75, 1.2]),
      P("<b>Overview table.</b> The numbers behind the overview figure (tables/scorecard.tsv). Level is the median of estimate ÷ ddPCR; "
        "1 is agreement. CONKORD exists only for Potapova et al.'s 12 lines. The 18S depth ratio ranks people as well as NGS-DOSE within one "
        "pipeline (inheritance), but it reads 9% high and moves by 4–25% when the same person is sequenced again elsewhere. The assembly's "
        "r with ddPCR rests on 5 people, one of them the curated T2T HG002, and its interval spans almost the whole range.", cap)]

# ---------------- background ----------------
S += [PageBreak(), P("1  Background: assembly methods and the rDNA", h1)]
S += [P("<b>The problem.</b> Each 45S unit is about 44.8 kb (the reference unit KY962518.1 is 44,838 bp). Units lie head to tail in five arrays "
        "on the short arms of chromosomes 13, 14, 15, 21 and 22, between the distal junction (DJ, toward the telomere) and the proximal "
        "junction. A diploid genome typically carries a few hundred units. The T2T consortium cites 315 ± 104 (Nurk et al. 2022, citing "
        "Parks et al. 2018), and ten HPRC lymphoblastoid lines ranged from over 400 to over 600 (Potapova et al. 2025). "
        "Individual arrays range from none or a few units to more than 150, and from 50 kb to more than 6 Mb (Stults et al. 2008; "
        "Potapova et al. 2025). Arrays are reported to rearrange frequently in meiosis (Stults et al. 2008), though recent pedigree work "
        "finds array sizes largely stable across a generation (Potapova et al. 2025; Cechova et al. 2025, preprint). The short arms also "
        "exchange sequence between different acrocentrics (Guarracino et al. 2023). The 5S genes form a separate array of 2.2-kb units at "
        "1q42 (Stults et al. 2008 measured an average of 98 units).")]
S += [P("<b>Why assemblers fail on it.</b> PacBio HiFi reads (13.5 kb on average and 99.8% accurate when first described; Wenger et al. 2019) are shorter than one "
        "unit. Ultra-long nanopore reads (N50 over 100 kb; Jain et al. 2018) span a few units but not a megabase-scale array of "
        "near-identical units. hifiasm (Cheng et al. 2021, 2022, 2024) builds phased string graphs from HiFi, optionally adding ultra-long "
        "reads, with phasing from parental k-mers (trio) or Hi-C. Verkko and Verkko2 (Rautiainen et al. 2023; Antipov et al. 2025) build "
        "de Bruijn graphs from HiFi, resolve them with ultra-long reads, and scaffold with Hi-C. In every such graph an rDNA array is a "
        "tangle. On linearisation it becomes a gap, a collapsed unit at a contig end, a small unplaced contig or a scattered set of "
        "copies. When the CHM13 reads were assembled with hifiasm, one unit variant with 12 true copies was placed in seven locations and "
        "another with 15 copies in one (Rautiainen 2024). Verkko2 can screen rDNA out of the main assembly, using KY962518 by default.")]
S += [P("<b>How complete genomes handled it.</b> T2T-CHM13 represents three of its five arrays as model sequence: blocks of each "
        "array's dominant unit variant, sized by estimated copy number. It holds 219 complete copies in 9.9 Mb, close to half the ddPCR "
        "diploid estimate of 409 for this effectively haploid genome (Nurk et al. 2022). T2T-HG002 v1.1 is complete except that nine of "
        "ten arrays are sized N-gaps; the one resolved array, paternal chr13, holds six full units plus a partial one (Hansen et al. "
        "2026; Potapova et al. 2025). Verkko2 scaffolded 55% of the acrocentrics it tested telomere to telomere, typically with the rDNA "
        "left as a gap (Antipov et al. 2025). In 156 short arms assembled from a pedigree the rDNA was a predictable collapse region "
        "(Lin et al. 2026). Only short arrays, or one array assembled with targeted ultra-long reads (Cechova et al. 2025, preprint), have "
        "been assembled completely. Ribotin reconstructs unit variants (morphs) and their relative abundance from long reads, but not "
        "array order or length (Rautiainen 2024).")]
S += [P(f"<b>HPRC release 2.</b> Release 2 (announced May 2025; 232 people) consists mostly of hifiasm assemblies with HiFi and "
        f"ultra-long input, trio- or Hi-C-phased; the index lists 452 hifiasm and 12 Verkko haplotype rows. The consortium reports that "
        "only 6.8% of acrocentric short arms were assembled contiguously without structural flags. Its T2T classification ignores the "
        "short arm up to and including the rDNA (Lucas et al. 2026, preprint), so a 'T2T' acrocentric says nothing about its rDNA. The "
        "CenSat annotation labels rDNA-like sequence found by a profile HMM, not arrays. <b>Prior expectation:</b> rDNA present as "
        "fragments and unplaced pieces, rarely as whole arrays, and not proportional to true copy number. The tests below check each "
        "part of that expectation.", body)]

# ---------------- methods ----------------
S += [P("2  Data and methods", h1)]
S += [bullets([
    "<b>Assemblies.</b> All 466 rows of the HPRC r2 index: 460 HPRC/HPP haplotypes (452 hifiasm 0.19.7–0.19.9, trio- or Hi-C-phased; "
    "8 Verkko 2.2.1, Hi-C), four extramural Verkko haplotypes (T2T-HG002 v1.1 maternal and paternal, HG06807) and two reference rows "
    "(CHM13 v2.0, GRCh38). Per haplotype the files used were the FASTA index, chromAlias, gap BED, T2T table, CenSat BED and the "
    "minigraph-cactus chain to GRCh38.",
    "<b>Targeted extraction</b> with samtools faidx over HTTPS (about 9.6 GB of bgzf blocks; 40 Gb of sequence). It took every CenSat rDNA "
    "interval ±200 kb; the short-arm segment of every contig assigned to an acrocentric (at least 12 Mb, and 1 Mb past the last rDNA "
    "or gap); every unplaced contig carrying CenSat rDNA or passing a 45S, 5S or DJ tile screen (505 contigs that CenSat leaves "
    "unlabelled were added this way); and the 5S locus lifted from GRCh38 chr1:228.6 Mb through the chain.",
    "<b>Completeness check.</b> Eight whole assemblies (hifiasm trio and Hi-C, Verkko, HG002) were tiled into 5-kb windows and mapped to "
    "the 45S and 5S units and the DJ core. In all eight, 100% of 45S- and 5S-matching sequence fell inside the extracted regions. "
    "CenSat alone missed 0–20% of 45S sequence per haplotype, all of it on unlabelled unplaced contigs.",
    "<b>Unit annotation</b> (minimap2 -x asm20, all secondary chains, repeat filter off). A <i>gene copy</i> is an 18S alignment over at "
    "least 95% of its length at at least 95% identity; a <i>5S unit</i> is a full X12811.1 alignment. An <i>array</i> is a run of gene "
    "copies on one contig with successive starts at most 120 kb apart. Each array end is a <i>contig end</i> (within 100 kb), a "
    "<i>gap</i> (an N-run within 100 kb) or <i>flank</i> (non-rDNA sequence); a DJ alignment of 50 kb or more within 1 Mb marks the "
    "distal side. A <i>unit period</i> is the sequence from one 18S start to the next; byte-identical periods mark templated or model "
    "sequence. DJ copies are DJ-core alignments of at least 200 kb.",
    f"<b>Annotator validation on known references.</b> CHM13: {c13['n_18S']} 18S and {c13['n_28S']} 28S copies (published: 219 complete "
    f"copies), 5 arrays, all 5 closed and DJ-flanked, and {100*c13['n_periods_identical']/c13['n_periods']:.0f}% identical periods, which "
    f"flags the model arrays correctly. HG002 paternal: {hgp['n_18S']} copies in {hgp['n_arrays']} pieces, one closed array of "
    f"{hgp['units_in_closed']} copies on chr13 (published: six full units plus a partial), the rest ending at N-gaps. HG002 maternal: "
    f"{hgm['n_18S']} copies, no closed array, 10 gap ends (published: gaps).",
    "<b>Orthogonal data.</b> Potapova et al. 2025 (Cell Genomics 5:101031), Table S1 (ddPCR and CONKORD totals for 12 LCLs) and Table S2 "
    "(FISH per array), transcribed from the supplementary PDF. Hall et al. 2021 Supplementary Data 1 (18S, 28S, 5S per haploid genome, "
    "from the same NYGC CRAMs by an independent pipeline).",
    f"<b>Short reads.</b> NGS-DOSE calibrated copy numbers (rDNA45S.cn, rDNA5S.cn, DJ.cn) from docs/data/cohort.tsv: {T['n_people_counted']} "
    f"hifiasm-assembled people so far, and {T['n_trios']} trio-phased children with both parents counted. Pearson r is reported with 5,000 "
    "bootstrap resamples. The perfect-measure range for a transmitted haplotype against its parent's diploid total comes from simulating "
    "five arrays per haplotype with log-normal sizes."])]

# ---------------- Q1 ----------------
S += [P("3  Do the assemblies contain rDNA?", h1)]
S += [fig("fig1_what_assemblies_hold.png", 7.0, 7.0 * 2.7 / 7.2),
      P(f"<b>Figure 1.</b> (a) Complete 18S genes per haplotype, by assembler, phasing and hifiasm version (bars: medians). The dashed line "
        f"is the NGS-DOSE median diploid copy number halved, for the {T['n_people_counted']} counted people; the dotted line is CHM13's "
        f"model arrays. (b) Units per array for {T['arrays_n']:,} hifiasm arrays (log scale). (c) What lies beyond each array end. (d) Where "
        "the units sit: a chromosome-length acrocentric scaffold, an unlocalised acrocentric piece, or a contig not assigned to a chromosome.", cap)]
G = T["by_group"]
S += [table([["Assembly group", "18S copies / hap", "arrays / hap", "closed arrays / hap", "5S units / hap", "DJ copies / hap"]] +
            [[g, f"{G['n_18S'][g]:.0f}", f"{G['n_arrays'][g]:.0f}", f"{G['n_closed'][g]:.1f}", f"{G['n_5S_units'][g]:.0f}", f"{G['n_DJ_copies'][g]:.0f}"]
             for g in ["hifiasm trio (0.19.7)", "hifiasm trio (0.19.9)", "hifiasm Hi-C (0.19.8-9)", "Verkko"]] +
            [["CHM13 v2.0 (model arrays)", c13["n_18S"], c13["n_arrays"], c13["n_closed"], c13["n_5S_units"], c13["n_DJ_copies"]],
             ["HG002 v1.1 paternal / maternal", f"{hgp['n_18S']} / {hgm['n_18S']}", f"{hgp['n_arrays']} / {hgm['n_arrays']}",
              f"{hgp['n_closed']} / {hgm['n_closed']}", f"{hgp['n_5S_units']} / {hgm['n_5S_units']}", "5 / 5"]],
            [2.2, 0.95, 0.85, 1.05, 0.95, 0.95]),
      P("<b>Table 1.</b> Medians per haplotype. hifiasm groups: trio 0.19.7 n = 122, trio 0.19.9 n = 130, Hi-C n = 200; Verkko 2.2.1 Hi-C n = 8.", cap)]
S += [P(
    f"All {T['n_hifiasm']} hifiasm haplotypes contain 45S units: between {n18['min']:.0f} and {n18['max']:.0f} complete 18S genes, median "
    f"{n18['50%']:.0f}, and 28S counts agree (r = 0.997; 28S/18S median 0.96). What they do not contain is arrays. The units come in a median "
    f"of {ph['n_arrays']['50%']:.0f} pieces per haplotype. Half of all pieces hold three units or fewer, and the largest piece in a "
    f"haplotype holds a median of {ph['largest_array']['50%']:.0f}. Of {tot_sides:,} piece ends, {100*sides.get('contig_end',0)/tot_sides:.0f}% are "
    f"contig ends, {100*sides.get('flank',0)/tot_sides:.0f}% run into non-rDNA sequence and {100*sides.get('gap',0)/tot_sides:.1f}% into N-gaps. "
    "This is the signature of an assembly graph broken inside a tangle, not of arrays sized and bridged by the scaffolder. Placement is "
    f"weak: {100*place.get('unplaced contig',0)/tot_units:.0f}% of units lie on unplaced contigs, "
    f"{100*place.get('acrocentric, unlocalised piece',0)/tot_units:.0f}% on unlocalised acrocentric pieces, and "
    f"{100*place.get('acrocentric, chromosome-length',0)/tot_units:.0f}% on chromosome-length acrocentric scaffolds.")]
S += [P(
    f"About one array per haplotype is flanked by non-rDNA sequence on both sides ({T['arrays_closed']} arrays in {T['n_hifiasm']} haplotypes; "
    f"median {T['closed_units']['50%']:.0f} units, maximum {T['closed_units']['max']:.0f}). {T['arrays_closed_DJ']} of them lie next to a distal "
    "junction. They may be genuinely short arrays, the kind FISH finds on a minority of chromosomes; HG002's resolved chr13 array "
    "(7 copies here) is one. Without read depth their sizes cannot be confirmed. The short arms themselves are largely present up to "
    f"the rDNA: a median of {ndj['50%']:.0f} distal junctions per haplotype, 82% of haplotypes within 4–6 (5 expected), with occasional "
    "extra copies that may be real or duplicated. The eight Verkko haplotypes hold a median of 15.5 gene copies, mostly beside N-gaps, "
    "consistent with Verkko's rDNA handling.")]

# ---------------- Q2 ----------------
S += [P("4  Is the assembled rDNA reliably measured?", h1)]
S += [fig("fig2_against_short_reads.png", 7.0, 7.0 * 3.0 / 7.2),
      P(f"<b>Figure 2.</b> (a) Assembled 18S genes (both haplotypes) against NGS-DOSE 45S copy number in the same person, n = {T['n_people_counted']}; "
        "solid line y = x, dashed line the median assembled fraction. (b) The same for the 5S array. (c) Transmitted haplotype against "
        f"each parent's short-read copy number in {T['n_trios']} trio-phased children. Filled markers: the parent who transmitted the "
        "haplotype; open markers: the other parent. The grey band is the 95% range for a perfect measurement at this n.", cap)]
S += [table([["Comparison", "n", "Pearson r (95% CI)", "Spearman ρ", "Assembly / short reads"],
             ["45S: assembled 18S (diploid) vs NGS-DOSE 45S", r45["n"], rr(r45), f"{r45['rho']:.2f}", f"{q45['50%']:.2f} ({q45['min']:.2f}–{q45['max']:.2f})"],
             ["45S: CenSat rDNA bp vs NGS-DOSE 45S", rc["n"], rr(rc), f"{rc['rho']:.2f}", "—"],
             ["5S: assembled units (diploid) vs NGS-DOSE 5S", r5["n"], rr(r5, 3), f"{r5['rho']:.3f}", f"{q5['50%']:.3f} ± {q5['std']:.3f}"],
             ["DJ: assembled copy-equivalents vs NGS-DOSE DJ", rdj["n"], rr(rdj), f"{rdj['rho']:.2f}", "—"],
             ["45S maternal haplotype vs mother / vs father", hp["45S_mat~mother"]["n"], f"{rr(hp['45S_mat~mother'])} / {hp['45S_mat~father']['r']:.2f}", "", ""],
             ["45S paternal haplotype vs father / vs mother", hp["45S_pat~father"]["n"], f"{rr(hp['45S_pat~father'])} / {hp['45S_pat~mother']['r']:.2f}", "", ""],
             ["5S paternal haplotype vs father / vs mother", hp["5S_pat~father"]["n"], f"{rr(hp['5S_pat~father'])} / {hp['5S_pat~mother']['r']:.2f}", "", ""],
             ["5S maternal haplotype vs mother / vs father", hp["5S_mat~mother"]["n"], f"{rr(hp['5S_mat~mother'])} / {hp['5S_mat~father']['r']:.2f}", "", ""],
             ["Perfect measurement, haplotype vs its parent (simulated)", T["n_trios"], f"{ceil['median']:.2f} ({ceil['lo']:.2f} to {ceil['hi']:.2f})", "", ""]],
            [2.75, 0.35, 1.9, 0.75, 1.25]),
      P("<b>Table 2.</b> Assembly content against NGS-DOSE short-read copy number.", cap)]
S += [P(
    f"<b>45S.</b> The assemblies hold {100*q45['50%']:.0f}% of the units the short reads count (interquartile range {100*q45['25%']:.0f}–"
    f"{100*q45['75%']:.0f}%). The fraction varies between people with a CV of about {100*q45['std']/q45['mean']:.0f}% and falls slightly as "
    f"copy number rises (Spearman ρ = {T['ratio45_vs_cn']:.2f}). It is higher for Hi-C-phased (hifiasm 0.19.8–9) than trio-phased assemblies "
    f"({T['ratio45_by_phasing'].get('hic', float('nan')):.2f} against {T['ratio45_by_phasing'].get('trio', float('nan')):.2f}; phasing and version "
    f"are confounded here). Across people the assembled count correlates with copy number at r = {r45['r']:.2f}. That is informative, but "
    "much weaker than short reads measured against themselves on two technologies (ICC 0.98). An assembly-based estimate of 45S copy "
    "number would therefore carry roughly 20% error per person and a two-fold bias. The CenSat-based figure used on the cohort page "
    f"tracks copy number less well (r = {rc['r']:.2f}), because CenSat leaves some rDNA contigs unlabelled.")]
S += [P(
    "<b>The haplotype test.</b> In trio-phased children, each haplotype's 18S count correlates with the copy number of the parent who "
    f"transmitted it (maternal {hp['45S_mat~mother']['r']:.2f}, paternal {hp['45S_pat~father']['r']:.2f}) and not with the other parent "
    f"({hp['45S_mat~father']['r']:.2f}, {hp['45S_pat~mother']['r']:.2f}). Both own-parent values lie inside the range a perfect measurement "
    "would give, but at 36 trios that range is wide enough to include a moderately noisy one too, so this test shows specificity, not "
    "precision. The partial content is not random. hifiasm keeps a share of each haplotype's rDNA that scales with how much that "
    f"haplotype carries, and phasing puts it on the right haplotype. The {T['n_trios']}-trio intervals are wide, and this should be rerun "
    "as counts accumulate.")]
S += [P(
    f"<b>5S: the positive control.</b> The 5S array is short enough (here 29–231 units per haplotype, about 0.06–0.5 Mb) to assemble. It comes out as a single "
    f"array flanked by unique sequence in {100*T['frac_hap_5S_single_closed']:.0f}% of haplotypes. Its assembled unit count matches NGS-DOSE's "
    f"5S estimate almost exactly: median ratio {q5['50%']:.3f}, SD {q5['std']:.3f}, r = {r5['r']:.3f}. This does two things. It shows that "
    "where an array is assembled, the unit count is right. And it validates NGS-DOSE's 5S estimate, which rests entirely on the "
    "fragment-GC model, to about 3.5% per person with about 2% absolute bias. The trio reliability of that estimate (0.78, interval "
    "0.43–1.10) had left it undecided. Haplotype-versus-parent correlations for 5S are noisier "
    f"({hp['5S_pat~father']['r']:.2f} and {hp['5S_mat~mother']['r']:.2f}), as expected with one array per haplotype and 36 trios.")]
S += [fig("fig3_unit_checks.png", 7.0, 7.0 * 2.4 / 7.2),
      P("<b>Figure 3.</b> (a) Share of unit periods byte-identical to another period in the same haplotype. CHM13's model arrays score 89%; "
        f"HPRC hifiasm arrays score {100*T['periods_identical']:.1f}%. (b) Median unit period per hifiasm array; dotted lines mark the five "
        "CHM13 arrays. (c) Distal-junction copies per hifiasm haplotype.", cap)]
S += [P(
    "<b>Unit-level checks.</b> The units hifiasm keeps look like real units rather than templated copies. Only "
    f"{100*T['periods_identical']:.1f}% of unit periods are byte-identical to another, against 89% for CHM13's model arrays and 5% for HG002. "
    "Unit periods cluster at 43–46 kb, as in CHM13's arrays. What the assemblies get wrong is therefore quantity and order, not unit "
    "sequence. That fits the literature: unit variants can be reconstructed, arrays cannot.")]

# ---------------- interpretation ----------------
S += [PageBreak(), P("5  Which method is more accurate?", h1)]
S += [P("Neither method can be the truth for the other, so four tests are used that do not require one: inheritance, a sequence of "
        "known copy number, an orthogonal laboratory assay (ddPCR), and per-chromosome FISH.")]
S += [P("5.1  Inheritance", h2)]
S += [P(
    "A child's rDNA dosage follows its parents'. Measurement error in the child does not, because the parents were sequenced as "
    "different people, from different libraries and in a different release batch. So if the parents are measured once, and the "
    "child is measured by two methods, the method whose child value tracks the parents' mean more closely carries less error. With the "
    "same parents for both, the squared ratio of the two correlations estimates the ratio of the two methods' reliabilities. Values were "
    "centred within superpopulation, and the two methods were bootstrapped as a pair over trios.")]
S += [fig("fig4_method_accuracy.png", 7.0, 7.0 * 2.8 / 7.4),
      P(f"<b>Figure 4.</b> (a, b) The same {i45['n']} children, measured by the HPRC assembly (18S genes, both haplotypes) and by NGS-DOSE, "
        "against their parents' mean NGS-DOSE 45S copy number (centred within superpopulation). (c) Child-versus-parents correlations with "
        "95% bootstrap intervals, for 45S and 5S, with the parents measured by NGS-DOSE or by the independent Hall et al. 2021 pipeline.", cap)]
S += [P(
    f"For 45S the NGS-DOSE child value tracks the parents at r = {i45['r_ngsdose']:.2f} ({i45['r_ngsdose_lo']:.2f}–{i45['r_ngsdose_hi']:.2f}), "
    f"the assembly at {i45['r_assembly']:.2f} ({i45['r_assembly_lo']:.2f}–{i45['r_assembly_hi']:.2f}). The paired difference is "
    f"+{i45['diff']:.2f} ({i45['diff_lo']:.2f} to {i45['diff_hi']:.2f}), and NGS-DOSE is ahead in {100*i45['p_ngsdose_better']:.0f}% of "
    f"resamples. With the parents measured by Hall et al. instead, the difference is +{i45h['diff']:.2f} ({i45h['diff_lo']:.2f} to "
    f"{i45h['diff_hi']:.2f}; {i45h['n']} trios). By the squared ratio, the assembly count carries about {100*i45['reliability_ratio']:.0f}% of "
    "the between-person signal NGS-DOSE does. With Hall-measured parents the assembly-only test can use "
    f"{M['assembly_hall_parents']['n']} trios, including African families; it gives r = {M['assembly_hall_parents']['r']:.2f} "
    f"({M['assembly_hall_parents']['lo']:.2f}–{M['assembly_hall_parents']['hi']:.2f}).")]
S += [P(
    f"<b>5S is the control.</b> There the assemblies are complete, and the two methods tie ({i5['r_assembly']:.2f} against "
    f"{i5['r_ngsdose']:.2f}; difference {i5['diff']:+.2f}, {i5['diff_lo']:.2f} to {i5['diff_hi']:.2f}). So the test has no built-in "
    "preference for short reads, and the 45S gap reflects the assemblies' incompleteness. One bias could remain: a heritable effect "
    "shared by short-read methods, such as rDNA sequence variants that change k-mer or read recovery, would favour them. The tie on 5S "
    "and the independent Hall parents argue against it, but do not exclude it.")]
S += [P("5.2  A known truth: the distal junction", h2)]
S += [fig("fig5_truths_and_assays.png", 7.0, 7.0 * 5.4 / 7.4),
      P(f"<b>Figure 5.</b> (a) Distal-junction copies (truth 10 per diploid genome) by NGS-DOSE and by assembly in the same {DJ['n']} people; "
        "assembly values jittered vertically; purple marks people NGS-DOSE places a whole copy from the cohort level. (b) Short-read "
        f"45S estimates against ddPCR (Potapova et al. 2025, Table S1; bars ±1 SD) for all {da['n']} ddPCR-measured lines. Filled: the "
        f"{dn['n']} in the NYGC 1000 Genomes 30× set; open: HG002, HG003 and HG004 on Google Health's GIAB NovaSeq pipeline (section 5.5). "
        "CONKORD values are Potapova et al.'s, from their own reads. (c) Assembled 18S genes against ddPCR for all "
        f"{DD['assembly_18S']['n']} people with both, with NGS-DOSE for each (HG002 from the Google pipeline, open); grey lines join each "
        "person's two values. (d) Totals summed over the ten FISH-measured arrays (Table S2) against NGS-DOSE (filled NYGC, open Google) and "
        "the assemblies; Potapova et al. scaled FISH "
        "to CONKORD totals, so these totals are CONKORD's. (e) FISH units per chromosome (both homologues) against assembled units placed "
        "on that chromosome. (f) The same as shares of each person's total.", cap)]
S += [P(
    f"The assemblies count distal junctions in whole copies: exactly 10 in {100*DJ['assembly_exact10']:.0f}% of these people, and 9–11 in "
    f"{100*DJ['assembly_all_9_11']:.0f}% of all {DJ['assembly_all_n']} assembled people. NGS-DOSE reads {DJ['ngsdose_mean']:.2f} ± "
    f"{DJ['ngsdose_sd']:.2f}, {100*(1-DJ['ngsdose_mean']/10):.1f}% below the truth, which is the method's known shortfall on this "
    f"sequence. On whole-copy steps the two agree. NGS-DOSE places {DJ['ngs_steps']} people one copy or more from the cohort level. The "
    f"assembly shows a step in the same direction in all {DJ['ngs_steps']}, and of the same size in {DJ['ngs_steps_confirmed']}. Of the "
    f"{DJ['ngs_normal']} people NGS-DOSE reads at the cohort level, the assembly gives 10 for {DJ['ngs_normal_asm_10']}, 9 for "
    f"{DJ['ngs_normal_asm_low']} and 11 for {DJ['ngs_normal_asm_high']}. Most of those discordant cases are probably short arms that failed "
    "to assemble or were duplicated, but without the reads that cannot be confirmed. <b>Verdict:</b> the assembly gets the absolute level "
    "right; NGS-DOSE is about 3% low but resolves single-copy changes reliably.")]
S += [P("5.3  An orthogonal assay: ddPCR", h2)]
S += [P(
    "Potapova et al. (2025) measured total rDNA copy number by droplet digital PCR in 12 lymphoblastoid lines. They also report their own "
    "short-read k-mer estimate (CONKORD), which they used to scale their FISH measurements. Nine of the lines are in the 1000 Genomes 30× "
    "set. Two were already in the cohort run, and the other seven were counted here in fetch mode with the same engine and bundle, "
    "then calibrated with the cohort's saved window efficiencies. This procedure reproduces the cohort's values for the two already "
    "counted to within 0.3%. HG002, HG003 and HG004 are not in 1000 Genomes; their values here come from a second NovaSeq pipeline and "
    "are reported separately in Table 3 and section 5.5.")]
S += [table([["Method (n = 9 unless noted)", "r with ddPCR", "Level (median ratio)", "Mean |error|", "SD of log ratio"],
             ["NGS-DOSE, calibrated", f"{dn['r']:.2f}", f"{dn['median_ratio']:.2f}", f"{dn['mean_abs_pct']:.1f}%", f"{DD['residual_sd_log']['ngsdose']:.3f}"],
             ["CONKORD (Potapova et al.)", f"{dc['r']:.2f}", f"{dc['median_ratio']:.2f}", f"{dc['mean_abs_pct']:.1f}%", f"{DD['residual_sd_log']['conkord']:.3f}"],
             ["18S depth ratio (published estimator, from NGS-DOSE counts)", f"{dfl['r']:.2f}", f"{dfl['median_ratio']:.2f}", f"{dfl['mean_abs_pct']:.1f}%", f"{DD['residual_sd_log']['ratio18S_flat']:.3f}"],
             [f"NGS-DOSE, Google GIAB NovaSeq (n = {dg['n']}: HG002–HG004)", f"{dg['r']:.2f}", f"{dg['median_ratio']:.2f}", f"{dg['mean_abs_pct']:.1f}%", f"{dg['sd_log']:.3f}"],
             [f"NGS-DOSE, both pipelines (n = {da['n']})", f"{da['r']:.2f}", f"{da['median_ratio']:.2f}", f"{da['mean_abs_pct']:.1f}%", f"{da['sd_log']:.3f}"],
             [f"CONKORD, all {dc12['n']}", f"{dc12['r']:.2f}", f"{DD['conkord']['median_ratio']:.2f}", f"{DD['conkord']['mean_abs_pct']:.1f}%", f"{DD['conkord']['sd_log']:.3f}"],
             [f"HPRC assembly, 18S genes (n = {DD['assembly_18S']['n']}, incl. HG002 v1.1)", f"{RC['assembly_5']['r']:.2f} ({RC['assembly_5']['lo']:.2f} to {RC['assembly_5']['hi']:.2f})", f"{DD['assembly_18S']['median_ratio']:.2f}",
              f"{DD['assembly_18S']['mean_abs_pct']:.0f}%", "—"],
             [f"ddPCR replicate CV (Potapova et al.)", "", "", f"median {100*DD['ddpcr_cv_median']:.0f}%", ""]],
            [2.6, 0.95, 1.25, 1.05, 1.15]),
      P("<b>Table 3.</b> Short-read and assembly estimates against ddPCR. Without HG02053, where ddPCR (713) exceeds both k-mer "
        f"pipelines (about 615), r is {DD['without_HG02053']['ngsdose']['r']:.2f} (NGS-DOSE), {DD['without_HG02053']['conkord']['r']:.2f} "
        f"(CONKORD) and {DD['without_HG02053']['ratio18S_flat']['r']:.2f} (18S ratio). Assembly r with 95% CI (Fisher z).", cap)]
S += [P(
    f"NGS-DOSE ranks these people most like ddPCR does (r = {dn['r']:.2f}), and its scatter around its own offset "
    f"({100*DD['residual_sd_log']['ngsdose']:.0f}%) is about the size of ddPCR's replicate variation. Its level is "
    f"{abs(100*(dn['median_ratio']-1)):.0f}% low (mean of log ratios {DD['bias_pct']['ngsdose']:+.1f}%), the same direction and size as its "
    f"distal-junction shortfall. CONKORD sits closest in level ({dc['median_ratio']:.2f}) but correlates less well, and the 18S depth ratio "
    f"is {100*(dfl['median_ratio']-1):.0f}% high. The assemblies hold about {100*DD['assembly_18S']['median_ratio']:.0f}% of the ddPCR copies. "
    f"Only 5 people have both ddPCR and an assembly, so the assembly's correlation with ddPCR is weakly determined: r = "
    f"{RC['assembly_5']['r']:.2f} (95% CI {RC['assembly_5']['lo']:.2f} to {RC['assembly_5']['hi']:.2f}). Most of the negative sign comes from HG002, "
    f"whose curated T2T v1.1 assembly holds only {PCMP.loc['HG002','assembly_18S']:.0f} 18S genes while ddPCR gives it the highest total of the five ({PCMP.loc['HG002','ddpcr']:.0f}); "
    f"without it r = {RC['assembly_4_noHG002']['r']:.2f} (n = 4, {RC['assembly_4_noHG002']['lo']:.2f} to {RC['assembly_4_noHG002']['hi']:.2f}). "
    f"NGS-DOSE on the same 5 people gives r = {RC['ngsdose_same5']['r']:.2f} ({RC['ngsdose_same5']['lo']:.2f} to {RC['ngsdose_same5']['hi']:.2f}) "
    f"and on the same 4, {RC['ngsdose_same4']['r']:.2f}. The inheritance test (section 5.1) is the better-powered comparison of the two. "
    "Nine samples over a narrow range (476–713 copies) give wide intervals on each r, and the DNA for ddPCR came from different cultures "
    "than the NYGC sequencing, which adds scatter that neither method can remove.")]
S += [P("5.4  Per array: FISH", h2)]
FT = M["fish_totals"]
S += [P(
    f"<b>Summed over the ten arrays.</b> Each person's FISH total equals the CONKORD total to within {FT['fish_minus_conkord_max']:.0f} "
    f"copies, because Potapova et al. scaled the fluorescence shares by CONKORD. Summed FISH is therefore not an independent total. "
    f"Against it, NGS-DOSE gives r = {FT['ngsdose']['r']:.2f} and ratio {FT['ngsdose']['median_ratio']:.2f} (n = {FT['ngsdose']['n']}), "
    f"and the assemblies ratio {FT['assembly']['median_ratio']:.2f} (n = {FT['assembly']['n']}; Figure 5d). The independent total "
    "is ddPCR (section 5.3).")]
S += [P(
    "<b>Per chromosome.</b> The FISH shares themselves are orthogonal to sequencing. "
    f"For the five people with both FISH and an assembly, between "
    f"{100*min(v for k,v in AR['placed_fraction'].items() if k!='HG002'):.0f}% and {100*max(v for k,v in AR['placed_fraction'].items() if k!='HG002'):.0f}% "
    "of each HPRC assembly's units are on contigs assigned to an acrocentric (all of HG002's). The share each chromosome gets in the "
    f"assembly is essentially unrelated to its FISH share (r = {AR['r_share']:.2f} over {AR['n_chrom']} chromosome totals), and so is the "
    f"number of units (r = {AR['r_units']:.2f}; Figure 5e, f). So the "
    "assemblies do not say how an individual's rDNA is divided among the five chromosomes. NGS-DOSE cannot say either, because it "
    "measures the total only.")]

S += [P("5.5  A second sequencing pipeline", h2)]
_scan = ""
if SC:
    _scan = (f" A whole-file scan of the HG002 Google BAM checks the fetch: the NYGC-learned sinks capture "
             f"{100*SC['frac_45S']:.2f}% of 45S, {100*SC['frac_5S']:.2f}% of 5S and {100*SC['frac_DJ']:.2f}% of distal-junction reads that the scan finds.")
else:
    _scan = " A whole-file scan to confirm that the NYGC-learned sinks capture this pipeline's reads was still running when this report was built."
gi = NV["google_insert"]
S += [P(
    "The NYGC 1000 Genomes CRAMs are one pipeline: one sequencing centre, one library protocol, one aligner build. Google Health's "
    "GIAB set (Baid et al. 2020) is a second: NovaSeq 6000 PCR-free 2×151 libraries made separately from the same cell lines, aligned "
    "with bwa-mem 0.7.17 to the same hs38DH reference. It covers HG002, HG003 and HG004, which have ddPCR but no 1000 Genomes reads, and "
    "the CEPH trio NA12878, NA12891 and NA12892, which NYGC also sequenced. All six were counted in fetch mode with the cohort engine, "
    "bundle and sinks, and calibrated with the NYGC cohort's window efficiencies, unchanged." + _scan)]
S += [fig("fig6_second_pipeline.png", 7.0, 7.0 * 2.9 / 7.4),
      P("<b>Figure 6.</b> (a) The CEPH trio on the two pipelines: percentage change, Google minus NYGC, for each estimate and for two library "
        "properties. Blue: NGS-DOSE classes; grey: the 18S depth ratio and the known-copy controls; orange: library properties. (b) NGS-DOSE "
        f"against ddPCR by pipeline (filled: NYGC, {dn['n']}; open: Google, {dg['n']}). (c) Median insert size and relative coverage of 65%-GC "
        f"sequence for the {len(json.load(open(f'{TAB}/method_accuracy.json'))['novaseq']['google_DJ'])} Google genomes against the NYGC cohort.", cap)]
p45, p18, p5 = NV["rDNA45S.cn"]["pct"], NV["rDNA45S.18S.flat"]["pct"], NV["rDNA5S.cn"]["pct"]
S += [P(
    f"<b>Same person, two pipelines.</b> The Google libraries have shorter inserts ({min(gi.values()):.0f}–{max(gi.values()):.0f} bp against "
    f"a cohort median of {NV['cohort_insert_median']:.0f}), and two of them, NA12891 and NA12892, have more GC bias than any NYGC genome "
    f"(coverage at 65% GC {NV['google_gc65']['NA12891']:.2f} and {NV['google_gc65']['NA12892']:.2f}; cohort maximum {NV['cohort_gc65_max']:.2f}). "
    f"NGS-DOSE's 45S values for the CEPH trio change by {p45[0]:+.1f}%, {p45[1]:+.1f}% and {p45[2]:+.1f}%. The 18S depth ratio from the same "
    f"reads changes by {p18[0]:+.0f}%, {p18[1]:+.0f}% and {p18[2]:+.0f}%, most in the two GC-biased libraries. The autosomal control and the "
    f"distal junction change by 1.5% or less. 5S reads {abs(np.mean(p5)):.0f}% lower on the Google pipeline in all three "
    f"({p5[0]:+.1f}%, {p5[1]:+.1f}%, {p5[2]:+.1f}%), so the 5S scale is not fully pipeline-independent.")]
S += [P(
    f"<b>Against ddPCR.</b> HG002, HG003 and HG004 read {dg['median_ratio']:.2f}× ddPCR (NYGC: {dn['median_ratio']:.2f}×), and the three "
    f"fall on the same line as the nine NYGC genomes. Across all {da['n']}, r = {da['r']:.2f}, against {dc12['r']:.2f} for CONKORD. The "
    "Google offset is within the spread of the NYGC one; three samples cannot show whether it is a real pipeline difference. "
    f"The distal junction reads {min(NV['google_DJ']):.2f}–{max(NV['google_DJ']):.2f} on the Google pipeline, the same range as in the NYGC cohort.")]

S += [P("6  NGS-DOSE: advantages and disadvantages, on this evidence", h1)]
S += [table([["", "Evidence"],
             ["<b>Advantages</b>", ""],
             ["Tracks true between-person variation in 45S better than any alternative tested",
              f"Inheritance r {i45['r_ngsdose']:.2f} vs assembly {i45['r_assembly']:.2f}; ddPCR r {dn['r']:.2f} vs CONKORD {dc['r']:.2f} and 18S ratio "
              f"{dfl['r']:.2f}; cross-technology ICC 0.98 vs 0.19 (pilot)"],
             ["Transfers to a second sequencing pipeline without re-learning",
              f"CEPH trio, NYGC vs Google NovaSeq: 45S changes {NV['rDNA45S.cn']['mean_abs_pct']:.1f}% on average (18S ratio {NV['rDNA45S.18S.flat']['mean_abs_pct']:.0f}%); "
              f"ddPCR r {da['r']:.2f} over both pipelines (n = {da['n']})"],
             ["Correct 5S copy number", f"Equals whole assembled 5S arrays: ratio {q5['50%']:.2f}, SD {100*q5['std']:.1f}%, r {r5['r']:.3f} (n = {r5['n']})"],
             ["Single-copy resolution in a paralogous sequence", f"All {DJ['ngs_steps']} distal-junction steps recur in the assemblies ({DJ['ngs_steps_confirmed']} with the same size)"],
             ["Whole-genome coverage of a cohort at low cost", "Every short-read genome; about 1 minute and 0.5 GB per genome in fetch mode. The assemblies cover 232 people and hold about half the 45S"],
             ["<b>Disadvantages</b>", ""],
             ["Absolute level a few percent low", f"{100*(dn['median_ratio']-1):+.0f}% vs ddPCR; distal junction {DJ['ngsdose_mean']:.2f} for 10; if the offset is constant, a ddPCR calibration would correct both"],
             ["Totals only", "No per-chromosome or per-haplotype array sizes; FISH provides these (and assemblies, as shown here, do not reliably)"],
             ["Scale depends on chemistry and batch", f"Window efficiencies differ between the two 1000 Genomes release batches (review of 2026-09-23); on the Google pipeline 5S reads "
              f"{abs(np.mean(NV['rDNA5S.cn']['pct'])):.0f}% lower and 45S {dg['median_ratio']:.2f}× ddPCR (NYGC {dn['median_ratio']:.2f}×); non-NovaSeq chemistries untested here"],
             ["Fetch mode depends on the aligner", "Sinks are learned from NYGC bwa-mem CRAMs; DRAGEN alignments untested"],
             ["Shares the cell-line limitation", "Culture changes to arrays (and S-phase effects) are measured as if they were genotype"]],
            [2.5, 4.5])]

S += [P("7  Interpretation", h1)]
S += [bullets([
    "<b>HPRC r2 assemblies are not a truth set for 45S copy number, per array or in total.</b> They hold about half of each person's "
    "units, in pieces that end at contig ends, distributed among chromosomes in a way FISH does not support. What they hold is "
    "haplotype-specific and scales with the true content, but carries about half the between-person signal NGS-DOSE does.",
    "<b>They are a truth set for 5S</b>, and there NGS-DOSE agrees to about 3.5%: an orthogonal confirmation of an NGS-DOSE rDNA class.",
    "<b>The best current picture of the absolute scale</b> is that NGS-DOSE reads about 3–5% low: the ddPCR offset and the distal-junction "
    "offset agree. A ddPCR-anchored correction on a few dozen cohort cell lines would settle it.",
    "<b>For the cohort page</b>, the statement that assemblies hold about a third of the rDNA should read about half, measured by gene "
    "copies. The CenSat label misses unplaced rDNA contigs."])]
S += [P("8  Limits", h1)]
S += [bullets([
    f"The short-read comparisons rest on {T['n_people_counted']} people and {T['n_trios']} trios, none of African ancestry, because only "
    f"{n_cohort:,} of 3,202 genomes have been counted so far. The scripts recompute every test from docs/data/cohort.tsv.",
    f"The ddPCR comparison has {da['n']} samples ({dn['n']} NYGC, {dg['n']} Google), and 5 people with both ddPCR and an assembly. ddPCR and "
    "sequencing used DNA from different cultures.",
    "The second-pipeline comparison rests on three people sequenced by both pipelines; both are NovaSeq 6000 PCR-free with bwa-mem on "
    "hs38DH, so it does not test other instruments, PCR-amplified libraries or other aligners.",
    "No read depth or read tiling was examined, so closed arrays and extra or missing DJ copies cannot be distinguished from collapses, "
    "false duplications or failed short arms.",
    "Standalone 18S copies (5% of hifiasm units) may include dispersed 18S-bearing fragments rather than array units.",
    "Hi-C-phased haplotypes cannot be assigned to a parent; phasing method and hifiasm version are confounded.",
    "Potapova et al.'s tables were transcribed from the supplementary PDF (checked: per-array values sum to the published totals "
    "within 2 copies); the seven additional NGS-DOSE counts used a sinks file that differs from the cohort's only by added telomere intervals."])]
S += [P("9  Files and regeneration", h1)]
S += [P("Everything lives in <font face='Courier'>assembly_rdna/</font>: background.md (the full cited background); data/ (index, regions; "
        "see data/README.md); potapova/ (Tables S1–S2 as TSV, the fetch script, counts files and estimates for the seven added genomes); novaseq/ (manifest, fetch and "
        "scan scripts, counts and estimates for the Google GIAB genomes and the NYGC CEPH trio); "
        "scripts/ (fetch_metadata.py, extract_regions.py, tile_map.py, validate_full.py, annotate_units.py, analyze.py, method_accuracy.py, "
        "build_report.py); tables/ (haplotype_rdna.tsv, arrays_all.tsv.gz, person_vs_ngsdose.tsv, trio_haplotypes.tsv, tests.json, "
        "method_accuracy.json, potapova_comparison.tsv, potapova_arrays.tsv, potapova_fish_totals.tsv, novaseq_bridge.tsv, novaseq_estimates.tsv, "
        "validation.tsv); figures/. "
        "<font face='Courier'>bash assembly_rdna/regenerate.sh</font> reruns the analysis and this PDF from new counts.", body)]
S += [P("References", h2)]
refs = ["Antipov D et al. (2025) Verkko2. Genome Res 35:1583. doi:10.1101/gr.280383.124",
        "Baid G et al. (2020) An extensive sequence dataset of gold-standard samples for benchmarking and development. bioRxiv. doi:10.1101/2020.12.11.422022",
        "Cechova M et al. (2025) Complete genomes of a multi-generational pedigree. bioRxiv. doi:10.64898/2025.12.14.693655",
        "Cheng H et al. (2021) hifiasm. Nat Methods 18:170; (2022) Nat Biotechnol 40:1332; (2024) Nat Methods 21:967. doi:10.1038/s41592-024-02269-8",
        "Guarracino A et al. (2023) Recombination between heterologous human acrocentric chromosomes. Nature 617:335. doi:10.1038/s41586-023-05976-y",
        "Hall AN, Turner TN, Queitsch C (2021) Thousands of high-quality sequencing samples fail to show meaningful correlation between 5S and 45S ribosomal DNA arrays in humans. Sci Rep 11:449. doi:10.1038/s41598-020-80049-y",
        "Hansen NF et al. (2026) A complete diploid human genome benchmark. Cell 189:4857. doi:10.1016/j.cell.2026.06.016",
        "Jain M et al. (2018) Nanopore sequencing with ultra-long reads. Nat Biotechnol 36:338. doi:10.1038/nbt.4060",
        "Kim J-H et al. (2018) Chromosome 21 rDNA by TAR cloning. Nucleic Acids Res 46:6712. doi:10.1093/nar/gky442",
        "Liao W-W et al. (2023) A draft human pangenome reference. Nature 617:312. doi:10.1038/s41586-023-05896-x",
        "Lin J et al. (2026) Human acrocentric chromosome short-arm de novo mutation and recombination. Cell 189:4876. doi:10.1016/j.cell.2026.05.035",
        "Lucas JK et al. (2026) HPRC2. bioRxiv. doi:10.64898/2026.07.21.739710",
        "Nurk S et al. (2022) The complete sequence of a human genome. Science 376:44. doi:10.1126/science.abj6987",
        "Parks MM et al. (2018) Variant ribosomal RNA alleles. Sci Adv 4:eaao0665 (cited via Nurk et al.)",
        "Potapova TA et al. (2025) Chromosome-specific epigenetic control and transmission of rDNA arrays. Cell Genomics 5:101031. doi:10.1016/j.xgen.2025.101031",
        "Rautiainen M et al. (2023) Verkko. Nat Biotechnol 41:1474. doi:10.1038/s41587-023-01662-6",
        "Rautiainen M (2024) Ribotin. Bioinformatics 40:btae124. doi:10.1093/bioinformatics/btae124",
        "Stults DM et al. (2008) Genomic architecture and inheritance of human rRNA gene clusters. Genome Res 18:13. doi:10.1101/gr.6858507",
        "Wenger AM et al. (2019) Accurate circular consensus long-read sequencing. Nat Biotechnol 37:1155. doi:10.1038/s41587-019-0217-9"]
S += [P(r, small) for r in refs]

def on_page(c, d):
    c.setFont("Helvetica", 7); c.setFillColor(colors.grey)
    c.drawString(0.75 * inch, 0.5 * inch, f"HPRC r2 assemblies and rDNA · NGS-DOSE counts as of {asof}")
    c.drawRightString(7.75 * inch, 0.5 * inch, str(d.page))
doc = SimpleDocTemplate(OUT, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch, topMargin=0.7 * inch,
                        bottomMargin=0.75 * inch, title="HPRC r2 assemblies, rDNA and NGS-DOSE")
doc.build(S, onFirstPage=on_page, onLaterPages=on_page)
print(OUT)
