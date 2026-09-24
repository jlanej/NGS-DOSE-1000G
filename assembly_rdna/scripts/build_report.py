#!/usr/bin/env python
"""Builds assembly_rdna/HPRC_r2_rDNA_report.pdf from tables/tests.json, the tables and figures/."""
import json, os
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
S += [P("Do HPRC release-2 assemblies contain the rDNA, and is it reliably measured?", title),
      P(f"Direct sequence tests on all {T['n_haplotypes']} HPRC/HPP release-2 haplotype assemblies, with CHM13, T2T-HG002 v1.1, "
        f"HG06807 and GRCh38 as references, against NGS-DOSE short-read copy numbers from the 1000 Genomes 30× CRAMs "
        f"(results as of {asof}, {n_cohort:,} genomes counted: {T['n_people_counted']} assembled people and {T['n_trios']} trios counted).", small), Spacer(1, 6)]
S += [P("Summary", h1)]
S += [bullets([
    f"<b>Do the assemblies contain rDNA? Yes, but not the arrays.</b> Every one of the {T['n_hifiasm']} hifiasm haplotypes holds "
    f"45S units: a median of {n18['50%']:.0f} complete 18S genes per haplotype (10th–90th percentile {n18['10%']:.0f}–{n18['90%']:.0f}). "
    f"They are spread over a median of {ph['n_arrays']['50%']:.0f} pieces of 1–10 units. "
    f"{100*sides.get('contig_end',0)/tot_sides:.0f}% of piece ends are contig ends, and "
    f"{100*place.get('unplaced contig',0)/tot_units:.0f}% of units sit on contigs not assigned to any chromosome. Only "
    f"{100*place.get('acrocentric, chromosome-length',0)/tot_units:.0f}% sit on a chromosome-length acrocentric. About one array per "
    f"haplotype is flanked on both sides (median {T['closed_units']['50%']:.0f} units, almost all next to the distal junction); these "
    "may be genuinely short arrays. The eight Verkko haplotypes hold far fewer units (median 15.5), with rDNA mostly left as gaps.",
    f"<b>Is it reliably measured? For 45S, no; for 5S, yes.</b> Summed over both haplotypes, the assemblies hold "
    f"{100*q45['50%']:.0f}% of the 45S copies NGS-DOSE measures in the same person (range {100*q45['min']:.0f}–{100*q45['max']:.0f}%). "
    f"Across people they correlate at r = {rr(r45)}. The 5S array, by contrast, is assembled whole and gap-free in "
    f"{100*T['frac_hap_5S_single_closed']:.0f}% of haplotypes. Its unit count matches NGS-DOSE at a ratio of {q5['50%']:.3f} (SD "
    f"{q5['std']:.3f}), with r = {rr(r5,3)}.",
    f"<b>The partial 45S content is not noise.</b> In {T['n_trios']} trio-phased children, the maternal haplotype's 18S count tracks the "
    f"mother's short-read copy number (r = {hp['45S_mat~mother']['r']:.2f}) and not the father's ({hp['45S_mat~father']['r']:.2f}). The paternal "
    f"haplotype tracks the father ({hp['45S_pat~father']['r']:.2f}) and not the mother ({hp['45S_pat~mother']['r']:.2f}). A perfect "
    f"measurement would give about {ceil['median']:.2f}. What hifiasm keeps scales with what the haplotype carries, but at roughly half the level "
    "and with about 20% person-to-person scatter in the assembled fraction.",
    "<b>What this means for NGS-DOSE.</b> Assemblies cannot validate 45S copy number, but they give a strong external check of NGS-DOSE's 5S "
    "estimate, which the trio test had left undecided. They also show that the CenSat-based figure on the cohort page (about 38% of "
    "rDNA assembled) understates what the assemblies hold, because CenSat leaves some rDNA contigs unlabelled."])]

S += [P("Verdict", h2)]
S += [table([["Question", "Answer", "Evidence"],
             ["Do HPRC r2 assemblies contain 45S rDNA?", "Yes, fragments, in every haplotype",
              f"{n18['50%']:.0f} 18S copies per haplotype (hifiasm); {100*sides.get('contig_end',0)/tot_sides:.0f}% of array ends are contig ends"],
             ["Are the 45S arrays complete?", "No, except about one short array per haplotype",
              f"{T['arrays_closed']} of {T['arrays_n']:,} arrays flanked on both sides; median {T['closed_units']['50%']:.0f} units"],
             ["Is assembled 45S a copy-number measurement?", "No", f"{100*q45['50%']:.0f}% of short-read copies; r = {r45['r']:.2f}"],
             ["Is assembled 45S haplotype-specific signal?", "Yes", f"own parent r ≈ {hp['45S_mat~mother']['r']:.2f} / {hp['45S_pat~father']['r']:.2f}; other parent ≈ 0"],
             ["Is the 5S array assembled and correct?", "Yes", f"{100*T['frac_hap_5S_single_closed']:.0f}% single gap-free arrays; ratio {q5['50%']:.2f}, r = {r5['r']:.3f}"],
             ["Are the short arms present up to the rDNA?", "Mostly", f"distal junctions per haplotype: median {ndj['50%']:.0f} (5 expected), 82% within 4–6"]],
            [2.1, 1.75, 3.15])]

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
S += [P("5  Interpretation", h1)]
S += [bullets([
    "<b>HPRC r2 assemblies are not a truth set for 45S copy number, per array or in total.</b> They hold about half of each person's "
    "units, broken into dozens of pieces that end at contig ends. The share is roughly proportional to the haplotype's true content, "
    "but the scatter is too large for per-person measurement. This agrees with the literature and with CHM13 and HG002, where the rDNA "
    "had to be modelled or left as gaps.",
    "<b>They are a truth set for 5S</b>, and on it NGS-DOSE's short-read estimate agrees to about 3.5%. That is the first orthogonal "
    "confirmation of an NGS-DOSE rDNA class. It suggests, but does not prove, that the fragment-GC model is well calibrated for a uniformly GC-rich (68%) unit.",
    "<b>For the cohort page</b>, the statement that assemblies hold about a third of the rDNA should be revised to about half, measured by "
    "gene copies. The CenSat label misses unplaced rDNA contigs and is not the right quantity to report.",
    "<b>Useful by-products.</b> The approximately one closed short array per haplotype, typically DJ-flanked, is a candidate set for "
    "per-array truth, but only after read-depth validation. The haplotype–parent signal shows that the assembled pieces are correctly "
    "phased, which matters for anyone using them to study rDNA unit variants."])]
S += [P("6  Limits", h1)]
S += [bullets([
    f"Short-read comparisons rest on {T['n_people_counted']} people and {T['n_trios']} trios, none of African ancestry, because only {n_cohort:,} of "
    "3,202 genomes have been counted so far. analyze.py recomputes every test from docs/data/cohort.tsv, so rerunning it after "
    "regenerate.sh updates the report.",
    "No read depth or read tiling was examined, so closed arrays and extra DJ copies cannot be distinguished from collapses or false "
    "duplications. Coverage-based flaggers are known to misbehave in rDNA.",
    "The gene-copy count includes 18S copies that stand alone (arrays of one: 5% of hifiasm units). Some may be dispersed "
    "18S-bearing fragments, such as the genuine 5′ETS/18S piece GRCh38 places outside the arrays, rather than array units.",
    "Hi-C-phased haplotypes cannot be assigned to a parent, so the haplotype test uses trio-phased assemblies only. Phasing method and "
    "hifiasm version are confounded.",
    "All DNA is from lymphoblastoid lines, and the assemblies and the NYGC short reads come from different DNA preparations. Culture "
    "changes to the arrays would appear as disagreement, although for 5S there is almost none."])]
S += [P("7  Files and regeneration", h1)]
S += [P("Everything lives in <font face='Courier'>assembly_rdna/</font>: background.md (the full cited background); data/ (index, metadata, "
        "extracted sequence, regions; see data/README.md); scripts/ (fetch_metadata.py, extract_regions.py, tile_map.py, validate_full.py, "
        "annotate_units.py, analyze.py, build_report.py); tables/ (haplotype_rdna.tsv, arrays_all.tsv.gz, person_vs_ngsdose.tsv, "
        "trio_haplotypes.tsv, tests.json, validation.tsv); figures/. <font face='Courier'>bash assembly_rdna/regenerate.sh</font> "
        "reruns the analysis and this PDF from new counts. Extraction and annotation need rerunning only if the assemblies change.", body)]
S += [P("References", h2)]
refs = ["Antipov D et al. (2025) Verkko2. Genome Res 35:1583. doi:10.1101/gr.280383.124",
        "Cechova M et al. (2025) Complete genomes of a multi-generational pedigree. bioRxiv. doi:10.64898/2025.12.14.693655",
        "Cheng H et al. (2021) hifiasm. Nat Methods 18:170; (2022) Nat Biotechnol 40:1332; (2024) Nat Methods 21:967. doi:10.1038/s41592-024-02269-8",
        "Guarracino A et al. (2023) Recombination between heterologous human acrocentric chromosomes. Nature 617:335. doi:10.1038/s41586-023-05976-y",
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
                        bottomMargin=0.75 * inch, title="HPRC r2 assemblies and rDNA")
doc.build(S, onFirstPage=on_page, onLaterPages=on_page)
print(OUT)
