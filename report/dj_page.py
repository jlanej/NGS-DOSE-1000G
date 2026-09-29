"""The distal-junction part of the page beyond the whole-copy steps: the profile along the unit, the
cohort's segment map and core level, the comparison with the HPRC assemblies, what the assemblies
get wrong, and what follows for the method. Rendered inside section 3.2."""
from __future__ import annotations

import html

import numpy as np

from .report import fmt

esc = lambda x: html.escape(str(x))


def _kb(a, b):
    return f"{a // 1000}–{b // 1000} kb"


def _pct(x, nd=0):
    return fmt(x, nd, pct=True) if x is not None else "–"


def gc_slope(eff: dict | None) -> tuple[float, float, int] | None:
    """The cohort's DJ window efficiencies against window GC: slope (log efficiency per unit GC) and r."""
    e = (eff or {}).get("DJ")
    if not e:
        return None
    a = np.array([np.nan if v is None else v for v in e["a"]], float)
    gc = np.array([np.nan if v is None else v for v in e["gc"]], float)
    m = np.isfinite(a) & np.isfinite(gc)
    if m.sum() < 30:
        return None
    return float(np.polyfit(gc[m], a[m], 1)[0]), float(np.corrcoef(gc[m], a[m])[0, 1]), int(m.sum())


def render(P, data: dict, rows: list[dict]):
    dj = data.get("dj")
    if not dj:
        return
    seg, st = dj["segments"], dj["stats"]
    asm = dj.get("assemblies")
    n = dj["n"]
    hyper = seg.get("hyper_intervals") or []
    top_bp = sorted(st["breakpoints"].items(), key=lambda kv: -kv[1])[:4]
    states = st["states"]

    def state_txt(kb):
        s = states.get(str(kb)) or states.get(kb) or {}
        tot = sum(s.values()) or 1
        return ", ".join(f"{k} copies in {100 * v / tot:.0f}%" for k, v in sorted(s.items(), key=lambda kv: int(kv[0])) if v / tot >= 0.03)

    # ------------------------------------------------------------------ the profile
    P.h("<h3>The junction in profile: what one number hides</h3>")
    P.h(f'''<p>The junction is not one number. Within a genome the cohort-calibrated window estimates read the same copy number in every window that
a junction copy holds whole, so the estimate along the 400-kb unit, in 20-kb blocks, shows what a single median hides: a copy that lacks part of
the unit lowers the blocks it lacks, an extra partial copy raises the blocks it holds. Across the {n:,} genomes the cohort SD of the calibrated
estimate per 5-kb sub-block has a floor of {fmt(seg["floor"], 2)} copies (counting noise) and rises far above it in
<strong>{", ".join(_kb(a, b) for a, b in hyper) if hyper else "no segment"}</strong>: those segments carry structural polymorphisms of the junction copies
(figure {"below" if dj.get("figure") else "in <code>data/dj_segments.tsv</code>"}). Steps of half a copy or more between adjacent blocks fall at the block boundaries at
{", ".join(f"{kb} kb in {v:,} genomes" for kb, v in top_bp)}. In the first block (0–20 kb) genomes read {state_txt(0)}; in the 200–220 kb block, {state_txt(200)}.
The rest of the unit is the <strong>core</strong>{f" ({len(seg.get('core_kb', []))} of 20 blocks)" if seg.get("core_kb") else ""}, and its median is a level the polymorphic segments cannot move.</p>''')
    w, c = st["whole"], st["core"]
    P.h(f'''<p>On the cohort's scale the main mode sits at {fmt(st["mode_level"], 2)} copies, not 10. Pinned so that the core's median reads 10 (a factor of
{fmt(st["pin"], 4)}), the whole-unit level lies within 0.2 copies of an integer in <strong>{_pct(w["within_0_2"], 1)}</strong> of genomes and the core level in
<strong>{_pct(c["within_0_2"], 1)}</strong>; {_pct(w["between"], 1)} and {_pct(c["between"], 1)} sit more than 0.3 from any integer. The mode cluster's SD is
{fmt(w.get("mode_sd"), 3)} copies on the whole unit and {fmt(c.get("mode_sd"), 3)} on the core. Every genome's profile on the pinned scale is in
<code>data/dj_blocks.tsv</code> (<code>DJ.cn_core</code> in the sample table is the core level on the cohort's scale).</p>''')

    if not asm:
        P.h('<p class="small">The comparison with HPRC assemblies appears when <code>--dj-assemblies</code> names the tables of <code>pipeline/hprc_dj.py</code>.</p>')
        return

    # ------------------------------------------------------------------ against the assemblies
    S, A = asm["stats"], asm["samples"]
    res = S["resolved"]
    k10 = S["known10"]
    reads10 = [x for x in k10["reads_cn"] if x is not None]
    P.h(f'<h3>Against the HPRC release-2 assemblies of {S["n"]} cohort members</h3>')
    P.h(f'''<p>For {S["n"]} genomes of the cohort with an HPRC release-2 assembly ({S["n_haplotypes"]} haplotypes, chosen to cover the whole-copy
steps, the values between steps and the mode; <code>pipeline/07_hprc_dj.sh</code>), every k-mer of the panel was counted in each haplotype FASTA
(<code>ngs-dose panel --report</code>: a complete junction copy contributes one; single-nucleotide differences lower single k-mers, so the
median per block is robust to them), and the assembly was aligned to the unit with everything outside a panel k-mer masked, so that each junction copy
appears as a run of alignments on one contig ({S["n_copies"]} copies: {", ".join(f"{v} {k}" for k, v in S["copy_classes"].items())}). The assembly's copies per
20-kb block are the two haplotypes' block medians summed; an assembly counts as <em>resolved</em> at the junction when no copy is cut by a contig end inside
the unit and at most one copy is partial. The reads' profile is the same calibrated estimate in blocks, on the ten-copy scale.</p>''')
    tiles = [("Genomes compared", f'{S["n"]}', f'{S["n_resolved"]} with the junction resolved in both haplotypes'),
             ("Core level, resolved assemblies", f'{fmt(res.get("core_diff_mean"), 2)} ± {fmt(res.get("core_diff_sd"), 2)}', f'reads − assembly, mean ± SD over {res.get("n", 0)} genomes; within half a copy in {res.get("core_within_0_5", 0)} of {res.get("n", 0)}'),
             ("Blocks", f'{_pct(S["blocks"]["concordance"])}', f'of {S["blocks"]["n"]:,} 20-kb blocks round to the same copy number; reads confirm {S["blocks"]["deviating_confirmed"]} of {S["blocks"]["deviating"]} blocks the assembly puts off 10')]
    if k10["n"]:
        tiles.append(("Ten-copy genomes by assembly", f'{fmt(float(np.median(reads10)), 2) if reads10 else "–"}', f'median of what NGS-DOSE reads on the cohort\'s scale in the {k10["n"]} genomes whose resolved assembly holds ten complete copies (range {fmt(min(reads10), 2)}–{fmt(max(reads10), 2)}); the cohort\'s core median is {fmt(st["core_median"], 2)}'))
    if S["partial_316"]:
        tiles.append(("The recurrent gain", f'{S["partial_316"]} of {S["n_partial"]}', "partial copies in the assemblies span exactly the first 316 kb of the unit, in tandem with a complete copy"))
    P.tiles(tiles)
    pts = [dict(x=t["assembly_mean"], y=t["reads_mean"], label=t["sample"], si=0 if t["resolved"] else 1,
                extra=[f"assembly {t['haplotypes']}", f"NGS-DOSE {fmt(t['reads_cn'], 2)} ({t['reads_step']:+.2f})" if t["reads_cn"] is not None else "", t["group"]]) for t in A]
    P.chart("djasm", dict(type="scatter", points=pts, legend=["assembly resolved at the junction", "assembly fragmented at the junction"], xlabel="HPRC assembly: copies per 20-kb block, mean over the unit",
                          ylabel="NGS-DOSE, ten-copy scale", identity=True),
            "The distal junction against long-read assemblies, genome by genome",
            f"Each dot is one of {S['n']} genomes of the 1000 Genomes 30× cohort with an HPRC release-2 assembly: the assembly's copies per 20-kb block (both haplotypes' "
            f"k-mer block medians summed) averaged over the unit, against NGS-DOSE's calibrated profile averaged the same way and pinned to the ten-copy scale; the line is identity. "
            f"Filled: the junction assembled without contig breaks in both haplotypes; open: fragmented. Over the resolved genomes the reads read {fmt(res.get('mean_diff_mean'), 2)} ± {fmt(res.get('mean_diff_sd'), 2)} "
            f"copies from the assembly (mean ± SD); the fragmented ones scatter because the assembly does, not the reads.")
    if dj.get("figure"):
        P.h(f'<figure><img src="{esc(dj["figure"])}" alt="the distal junction: segment map, per-genome agreement, and profiles against the assemblies" style="width:100%">'
            f'<figcaption>The junction in profile. A: cohort SD of the calibrated estimate per 5-kb sub-block over {n:,} genomes, the noise floor dashed, the hyper-variable segments shaded. '
            f'B: as the chart above. C: every compared genome\'s profile in 20-kb blocks, the assembly\'s row above the reads\' row, blue below ten copies and red above; hatched blocks are the '
            f'hyper-variable segments left out of the core. Redrawn from <code>report.json</code> and <code>data/dj_hprc_blocks.tsv</code> by <code>python -m report.dj_figure</code>.</figcaption></figure>')
    P.h("<p>Genome by genome (the assembly's haplotypes as the k-mer median over the unit; the core is the median over the core blocks, the mean over all twenty; "
        "NGS-DOSE on the cohort's scale, then its core and mean on the ten-copy scale; the last column is what Rhie et al. 2026 report for the same assemblies):</p>")
    hdr = ["sample", "class", "assembly (k-mer medians)", "assembly core", "assembly mean", "NGS-DOSE (step)", "reads core", "reads mean", "resolved", "partial / truncated copies", "Rhie et al. 2026"]
    P.table([[t["sample"], t["group"], t["haplotypes"], fmt(t["assembly_core"], 0), fmt(t["assembly_mean"], 2),
              f'{fmt(t["reads_cn"], 2)} ({t["reads_step"]:+.2f})' if t["reads_cn"] is not None else "–", fmt(t["reads_core"], 2), fmt(t["reads_mean"], 2),
              "yes" if t["resolved"] else "no", "; ".join(x for x in (t["partial_extents"], t["truncated_extents"]) if x) or "–", t["rhie"] or "–"]
             for t in sorted(A, key=lambda t: (t["assembly_mean"], t["sample"]))], hdr, numeric={3, 4, 5, 6, 7}, flagged=lambda r: r[8] == "no", wrap=True)
    P.h(f'<p class="small">Every block of every compared genome: <code>data/dj_hprc_blocks.tsv</code>; every junction copy the alignments found, with its contig, extent and class: <code>data/dj_hprc_copies.tsv</code>; the per-genome table: <code>data/dj_hprc.tsv</code>.</p>')

    # ------------------------------------------------------------------ what the assemblies get wrong
    odd = asm.get("oddities") or []
    P.h("<h3>What the assemblies get wrong at the junction, and what the reads show</h3>")
    trunc = [o for o in odd if o["truncated"]]
    over = [o for o in odd if o["over_blocks"] and not o["resolved"]]
    under = [o for o in odd if o["under_blocks"] and o["truncated"]]
    split = [o for o in odd if o["phasing_split"] >= 3]
    kva = [o for o in odd if o["unit_vs_blocks"]]
    unexplained = [t for t in A if t["resolved"] and abs(t["diff_core"]) >= 0.35]
    P.h(f'''<p>The junction sits beside the rDNA array, where long-read assemblies break, and its copies are {">"}99% identical, so an assembly's
count of it is only as good as its contigs there. The screens show four ways an assembly misreads the junction, each visible in the reads:</p><ul>
<li><strong>Copies cut by a contig end</strong> ({S["n_truncated"]} of {S["n_copies"]} copies, in {len(trunc)} genomes: {", ".join(esc(o["sample"]) for o in trunc) or "none"}). The k-mers of the missing part are
absent from the assembly, so it under-counts those blocks{f"; the reads read them at their full copy number in {', '.join(esc(o['sample']) for o in under)}" if under else ""}.</li>
<li><strong>Fragments assembled twice</strong>. Small contigs holding the same end of the junction (two of them of 0.2 Mb in HG00658's paternal haplotype, both with the unit's last 80 kb) add a copy the reads do not see
{f": the assembly exceeds the reads by three quarters of a copy or more in {', '.join(esc(o['sample']) + ' (' + ', '.join(str(b) for b in o['over_blocks']) + ' kb)' for o in over)}" if over else ""}.</li>
<li><strong>Phasing.</strong> The two haplotypes of a person hold five copies each, one per acrocentric homologue, but the assemblies split them {", ".join(f"{esc(o['sample'])} {o['haplotypes'].replace(' + ', ' / ')}" for o in split) if split else "unevenly in none"}:
the acrocentric short arms are phased by Hi-C or trio k-mers no better than their sequence allows. The sum over both haplotypes is what the reads measure, and it is right.</li>
<li><strong>A haplotype's k-mer median over the whole unit under-reads a high copy number</strong>{f" ({', '.join(esc(o['sample']) + ': ' + esc(o['unit_vs_blocks']) for o in kva)})" if kva else " (none here)"}: with six or more copies, nucleotide differences among them leave fewer than half the unit's k-mers at
the full count, so an exact-k-mer count from an assembly (as from reads, in a method that counts k-mer multiplicity) needs the median per block, not per unit. NGS-DOSE classifies a read by any four of its k-mers and is untouched by single differences.</li></ul>''')
    if unexplained:
        P.h(f'''<p>Where a resolved assembly and the reads disagree by a third of a copy or more on the core, {", ".join(f"{esc(t['sample'])} (assembly {fmt(t['assembly_core'], 0)}, reads {fmt(t['reads_core'], 2)})" for t in unexplained)},
neither side is confirmed: a copy the assembly does not hold, a junction lost or gained in part of the cell line (the assembly's DNA and the reads' DNA are different cultures of the same line), or a genome whose
level the cohort model misplaces are all possible. {"None of them has a chromosome flagged in the control regions (a chromosome 4% off its expected dosage is), so a mosaic loss or gain of a whole acrocentric is not the reason." if not any(t["flagged_chromosomes"] for t in unexplained) else "Flagged in the control regions: " + "; ".join(esc(t["sample"]) + " (" + esc(t["flagged_chromosomes"]) + ")" for t in unexplained if t["flagged_chromosomes"]) + "."}</p>''')
    n_ds = S["n_distal_start"]
    P.h(f'''<p>Two structures of the junction copies themselves recur across the assemblies and explain the cohort's segment map. <strong>Copies that begin about 22 kb into the unit</strong>
({n_ds} of {S["n_copies"]} copies, in {sum(1 for t in A if t["n_distal_start"])} of {S["n"]} genomes): the distal 22 kb is a deletion polymorphism, which is why the first block reads
{state_txt(0)} across the cohort. And <strong>a partial copy holding the first 316 kb of the unit, in tandem with a complete copy</strong>
({S["partial_316"]} of the {S["n_partial"]} partial copies; in every gain genome with a resolved assembly): the recurrent gain is not an eleventh junction but a duplication of most of one,
which reads 11 over the first 320 kb and 10 beyond, and a whole-unit median of it lands between steps. The 200–220 kb block is a deletion polymorphism of single copies (one copy of HG00097's first haplotype lacks 193–222 kb),
read at {state_txt(200)} across the cohort.</p>''')

    # ------------------------------------------------------------------ what follows
    gs = gc_slope(data.get("efficiencies"))
    cap = ((data.get("modes") or {}).get("capture") or {}).get("DJ") or {}
    P.h("<h3>What follows for the method</h3>")
    items = []
    if k10["n"] and reads10:
        items.append(f"<li><strong>The scale.</strong> The {k10['n']} genomes whose resolved assembly holds ten complete copies read {fmt(min(reads10), 2)}–{fmt(max(reads10), 2)} on the cohort's scale (median {fmt(float(np.median(reads10)), 2)}; some were chosen for reading between steps), and the cohort's core sits at {fmt(st['core_median'], 2)}: "
                     f"the level, which rests on the fragment-GC model in windows of 40–60% GC, is {_pct(1 - st['core_median'] / 10, 1)} low for the junction"
                     + (f". The cohort's window efficiencies fall with window GC (slope {fmt(gs[0], 2)} in log efficiency per unit GC, r = {fmt(gs[1], 2)} over {gs[2]:,} windows), as the 45S unit's do: a residual of the GC model in repeat context that the anchors do not remove" if gs else "")
                     + ". The junction has ten copies in nearly everyone, and the assemblies confirm it in the genomes that read low, so its scale can be pinned to the core's mode (or to anchor windows chosen on replicate pairs, as the 45S unit's are) rather than left on the GC model.</li>")
    items.append(f"<li><strong>The core, not the whole unit, as the level.</strong> The hyper-variable segments ({', '.join(_kb(a, b) for a, b in hyper) or 'none found'}) are polymorphisms of single copies; left out of the level, "
                 f"the whole-copy steps are called on the {len(seg.get('core_kb', []))} core blocks and the segments reported apart, with their own states per genome (<code>data/dj_blocks.tsv</code>).</li>")
    items.append("<li><strong>Partial variants by segment.</strong> A whole-unit median puts a 316-kb duplication at +0.6 to +0.9 and a 120-kb partial loss at −0.3 to −0.5 (the values between steps, "
                 f"{_pct(w['between'], 1)} of the cohort on the whole unit); a segmentation of the profile with the recurrent breakpoints (22, 200, 220, 316 kb) calls each segment's integer and names the variant, as the assemblies do.</li>")
    items.append("<li><strong>Read counting is not the limit.</strong> " + (f"The targeted fetch reads {_pct(cap.get('median'), 2)} of the junction's reads (median over {cap.get('n', 0):,} genomes; section 3.3); " if cap.get("n") else "")
                 + "duplicate-flagged reads are counted, so the ten-fold pile-up of junction reads on GRCh38 costs nothing; a read is classified by any four of its k-mers, so sequence differences among the copies, "
                 "which lower an assembly's exact-k-mer counts by 5% per haplotype, do not lower the reads' count.</li>")
    items.append("<li><strong>Assemblies as truth, with care.</strong> Of the compared genomes the junction was fully resolved in " + f"{S['n_resolved']} of {S['n']}; for the rest the reads are the more coherent reading of the locus, "
                 "and the profile identifies which blocks the assembly broke. A screen of all 200 assembled cohort members (two hours of downloads per twenty) would give the panel of truths a second class needs.</li>")
    P.h("<ul>" + "".join(items) + "</ul>")
