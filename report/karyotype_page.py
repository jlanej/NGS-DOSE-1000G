"""Section 3.10 of the page: chromosomes in copies. What is read and how, every chromosome's precision, the sex
chromosomes, the autosomes, the checks from the inside and against what else is known, and what depth cannot see."""
from __future__ import annotations

import html

from .report import fmt

esc = lambda x: html.escape(str(x))
WORDS = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}


def _n(k: int) -> str:
    return WORDS.get(k, f"{k:,}")


def _pct(x, nd=0) -> str:
    return fmt(x, nd, pct=True) if x is not None else "–"


def _mb(a, b) -> str:
    return f"{a / 1e6:.0f}–{b / 1e6:.0f} Mb" if b - a >= 5e6 else f"{a / 1e6:.1f}–{b / 1e6:.1f} Mb"


def _signed(x, nd=2) -> str:
    s = f"{x:+.{nd}f}"
    return s.replace("-", "−")


def _and(xs: list[str]) -> str:
    return ", ".join(xs[:-1]) + (" and " if len(xs) > 1 else "") + xs[-1] if xs else ""


def render(P, data: dict, rows: list[dict], sex_levels=None):
    k = data.get("karyotype")
    if not k:
        return
    m, ru, st, sx, au, ck = k["model"], k["rules"], k["status"], k["sex"], k["autosomes"], k["checks"]
    P.section("karyotype", "3.10 Chromosomes in copies: a karyotype for every genome, from the single-copy regions", "Chromosomes")
    single = m.get("single_values") or []
    lay = k.get("layouts")
    lw, lb = (lay or {}).get("with_windows"), (lay or {}).get("without")
    if lw and lb:
        read_from = (f"""{lw["n"]:,} genomes are read here from counts that hold the windows ({lw["regions"]:,} regions each; fetched for this page from the public files, the containers that hold them and no more),
the other {lb["n"]:,} from the {lb["regions"]:,} regions that the cohort's counts hold (the windows did not exist when it was scanned).""")
    elif lw:
        read_from = f"""All {lw["n"]:,} genomes are read from counts that hold the windows ({lw["regions"]:,} regions each)."""
    else:
        read_from = (f"""{m["regions"]:,} regions are read here ({m["autosomal"]:,} autosomal, {m["x"]:,} on chrX, {m["y"]:,} on chrY"""
                     + ("; " + _and(sorted(set(single))) + " as one pooled value each, in estimates written before the regions were kept singly" if single else "") + ").")
    err = m.get("error") or {}
    tl = k.get("tilt")
    P.h(f"""<p>Every counts file holds single-copy regions: the controls that every class is a ratio against, the held-out known-truth regions on the autosomes, chrX and chrY (section 3.1), and, in counts made with
NGS-DOSE 0.3.0's bundle, the <em>karyotype windows</em>: pieces of clean single-copy sequence, a few to a window, spread along every arm. A region's depth, relative to what the genome's own GC curve expects, says how many
copies of that sequence the genome holds; a chromosome is its regions read together, by alignment position, in a whole-file scan and in a targeted fetch alike. {read_from}</p>
<p>{f"A model learned on the cohort ({m['n']:,} genomes; a region that only part of the cohort holds is learned on that part)" if not m.get("saved") else f"The bundle's model, learned on {m['n']:,} genomes of the 1000 Genomes Project (a cohort of fifty or more learns its own),"} carries what stands between a region and a number of copies.
The region's <em>efficiency</em>: its median over the cohort. The libraries' <em>shared modes</em>: {m["components"]} components that stand above the noise edge (late-replicating sequence reads low in DNA from cycling cells, GC residue),
learned and scored on what the regions of a chromosome do <em>relative to each other</em>, so that a chromosome gained or lost, which moves all its regions alike, cannot be taken for a mode and learned away, however common it is in the cohort.
The region's <em>spread</em>, which is its weight. {"And one thing the contrasts cannot see: in some libraries the chromosomes rich in GC (19, 22, 17, 16) read low or high together, by the chromosome's GC and with no pattern inside a chromosome. Each chromosome is therefore set against the line that the <em>other</em> autosomes of the same genome give (a common shift and a slope on the chromosomes' GC, fitted so that a chromosome far off the line has no say); a chromosome takes no part in its own correction. " if tl else ""}Along each chromosome a chain then finds the levels: one level throughout is a whole chromosome,
a change at the centromere an arm, anything else a stretch (at least {ru["min_loci"]} places; a window is one place, and one region that a copy-number variant moves cannot buy a stretch). Every level is kept as measured, with its distance from the
nearest whole number: {fmt(ru["level_z"], 0)} standard errors and {fmt(ru["min_off"], 2)} copies from it, a level is called <em>fractional</em>, a change in part of the cells.
{f"A level's error is its regions' noise in that genome, and the cohort says how far that is from the whole of it: {fmt(err['factor'], 2)} times, with a floor of {fmt(err['floor'], 4)} copies that more regions do not lower. " if err else ""}A second X is read on a scale of its own:
it reads {fmt(m["u"], 3) if m.get("u") else "–"} of the first here, as an inactive X that replicates late does in a growing culture.
Columns: <code>karyotype</code> (written like one: <code>47,XY,+21</code>; a change in part of the cells with its share in brackets, <code>46,XX,-X[0.20]</code>; a stretch with its span),
<code>karyotype.status</code>, <code>sex_chromosomes</code>, <code>chr1.copies</code> … <code>chrY.copies</code> with their <code>.z</code>, <code>karyotype.events</code> and <code>karyotype.gc_tilt</code>. The method: NGS-DOSE's DESIGN.md, section 8; the evidence: <a href="KARYOTYPE.md">KARYOTYPE.md</a>.</p>""")
    n_part = st.get("fractional", 0)
    P.tiles([("Genomes read", f'{k["n"]:,}', (f'{lw["n"]:,} with the karyotype windows ({lw["regions"]:,} regions), {lb["n"]:,} without ({lb["regions"]:,})' if lw and lb
                                              else f'{m["regions"]:,} regions each; region noise {fmt(k["noise"][1], 2) if k.get("noise") else "–"} of the cohort\'s in the median genome')),
             ("46,XX or 46,XY", f'{k["plain"]:,}', "two copies of every autosome throughout, and the sex chromosomes of a woman or a man"),
             ("In every cell", f'{k["whole_number"]:,}', "another whole number: a sex-chromosome complement other than XX or XY, or a chromosome or stretch at one, three or four copies"),
             ("In part of the cells", f'{n_part:,}', "a chromosome or a stretch between two whole numbers" + (f'; {st.get("uncertain", 0):,} read too coarsely to settle' if st.get("uncertain") else ""))])

    # the sex chromosomes
    comp = sx["by_complement"]
    P.chart("karyo_xy", dict(type="scatter", x="chrX.copies", y="chrY.copies", group=dict(col="sex", levels=sex_levels or [["M", "male"], ["F", "female"]]), xlabel="chrX, copies",
                             ylabel="chrY, copies"),
            "The sex chromosomes of every genome, in copies",
            f"Each point is one of {k['n']:,} genomes: its chrX and chrY in copies, on the scales that put the cohort's one-X and two-X genomes at one and two and its one-Y genomes at one; coloured by the sex the pedigree gives. "
            "Genomes on whole numbers are a complement (XX, XY, X, XXY, XYY, XXX); a genome between two of them has lost or gained the chromosome in part of its cells. "
            "A sample that held two people's DNA would lie on the line from XX to XY.")
    txt = ("<p><strong>The sex chromosomes.</strong> By their whole numbers: " + _and([f"{n:,} {esc(c)}" for c, n in sorted(comp.items(), key=lambda kv: -kv[1])]) + ". ")
    wn = [t for t in sx["not_plain"] if t["complement"] not in ("XX", "XY")]
    if wn:
        txt += ("Other than XX and XY: " + "; ".join(f"{esc(t['sample'])} {esc(t['karyotype'])} (X {fmt(t['x'], 2)}, Y {fmt(t['y'], 2)}"
                                                    + (f", pedigree {esc(t['reported'])}" if t["reported"] else "") + ")" for t in wn[:14]) + (f"; and {len(wn) - 14} more" if len(wn) > 14 else "") + ". ")
    if sx["mismatch"]:
        txt += (f"The pedigree's sex is not what the chromosomes say in {_n(len(sx['mismatch']))}: " + ", ".join(esc(s) for s in sx["mismatch"][:10]) + ". ")
    ym = sx.get("y_lost_in_most") or []
    if ym:
        txt += (f"{_n(len(ym)).capitalize()} {'man reads' if len(ym) == 1 else 'men read'} as one X with a Y left in part of the cells ("
                + ", ".join(f"{esc(d['sample'])}, Y {fmt(d['y'], 2)}" for d in ym[:6]) + "): lines that lost the Y in most of their cells, which is written "
                "<code>45,X,+Y[share]</code>. ")
    pt = sx["part"]
    share = lambda d, first=False: ((f" ({_pct(d['q'][0])} to {_pct(d['q'][2])} of the cells{' from the 10th to the 90th percentile' if first else ''}, "
                                      f"median {_pct(d['q'][1])})") if d["n"] >= 3 else "")
    txt += (f"In part of the cells: an X lost in {pt['x_lost']['n']:,} genomes{share(pt['x_lost'], True)} and gained in {pt['x_gained']['n']:,}{share(pt['x_gained'])}, "
            f"a Y lost in {pt['y_lost']['n']:,}{share(pt['y_lost'])} and gained in {pt['y_gained']['n']:,}{share(pt['y_gained'])}. "
            "These are cell lines: an X or a Y lost in culture is what most of them are.")
    P.h(txt + "</p>")
    pub = [t for t in sx["not_plain"] if t["sample"] in ("HG01683", "HG03456")]
    if pub or sx.get("ngspca"):
        txt = "<p><strong>Against what else is known of them.</strong> "
        got = {t["sample"]: t for t in pub}
        if "HG01683" in got:
            txt += (f"HG01683 reads {esc(got['HG01683']['karyotype'])} here and was found to be XXY by Richmond and colleagues in Illumina's Polaris sequencing of the same line "
                    '(<a href="https://doi.org/10.1371/journal.pcbi.1008815">PLoS Comput Biol 2021</a>). ')
        if "HG03456" in got:
            txt += (f"HG03456 reads {esc(got['HG03456']['karyotype'])} and is described as XYY in the long-read assemblies of Logsdon and colleagues "
                    '(<a href="https://doi.org/10.1038/s41586-025-09140-6">Nature 2025</a>). ')
        ng = sx.get("ngspca")
        if ng:
            txt += (f"NGS-PCA's coverage ratio of the whole X chromosome, from the same files by another method, agrees with the X read here at r = {fmt(ng['r'], 4)} over {ng['n']:,} genomes"
                    + (f" and at r = {fmt(ng['r_two'], 3)} among the genomes with two X" if ng.get("r_two") is not None else "") + ". ")
        P.h(txt + "</p>")
    if sx["not_plain"]:
        shown = sx["not_plain"][:40]
        P.h(f"<p>The {len(shown)} genomes furthest from XX and XY of {len(sx['not_plain']):,} that are not plainly one of them (every genome's values are in <code>data/cohort.tsv</code>):</p>")
        P.table([[t["sample"], t["pop"], t["reported"], t["karyotype"], fmt(t["x"], 3), fmt(t["y"], 3), t["note"] or "–"] for t in shown],
                ["sample", "population", "pedigree", "karyotype", "chrX, copies", "chrY, copies", "note"], numeric={4, 5}, wrap=True, filter_box=True)

    # every chromosome's precision and what was found
    seen = f"({fmt(ru['level_z'], 0)} standard errors, and never below {_pct(ru['min_off'])})"
    if lw and lb:
        P.h(f"""<p><strong>Every chromosome.</strong> How well a chromosome's level is known is set by the regions counted on it: the standard error of a typical genome's level, in copies, and the share of the cells from which
a gain or a loss of one copy is called {seen}. The controls were chosen for another purpose and lie unevenly: five on chromosome 19, three on chromosome 22. The windows are laid along every arm, and with them every autosome
is known to about half a percent of a copy.</p>""")
        bw = {c["label"]: c for c in lw["chromosomes"]}
        bb = {c["label"]: c for c in lb["chromosomes"]}
        P.table([[c["label"], bw[c["label"]]["regions"] if c["label"] in bw else "–", fmt(bw[c["label"]]["se"], 4) if c["label"] in bw else "–", _pct(bw[c["label"]]["seen_from"], 1) if c["label"] in bw else "–",
                  bb[c["label"]]["regions"] if c["label"] in bb else "–", fmt(bb[c["label"]]["se"], 4) if c["label"] in bb else "–", _pct(bb[c["label"]]["seen_from"], 1) if c["label"] in bb else "–",
                  c["whole"] or "–", c["part"] or "–", c["stretches"] or "–"] for c in k["chromosomes"]],
                ["chromosome", "regions, with the windows", "SE of the level, copies", "one copy seen from this share of the cells", "regions, without", "SE of the level, copies", "one copy seen from",
                 "more or less in every cell", "in part of the cells", "arms and stretches"], numeric={1, 2, 3, 4, 5, 6, 7, 8, 9}, wrap=True)
    else:
        P.h(f"""<p><strong>Every chromosome.</strong> How well a chromosome's level is known is set by the regions counted on it: the standard error of a typical genome's level, in copies, and the share of the cells from which
a gain or a loss of one copy is called {seen}.</p>""")
        P.table([[c["label"], c["regions"], c["places"], fmt(c["se"], 4), _pct(c["seen_from"], 1), c["whole"] or "–", c["part"] or "–", c["stretches"] or "–"] for c in k["chromosomes"]],
                ["chromosome", "regions", "places", "SE of the level, copies", "one copy seen from this share of the cells", "more or less in every cell", "in part of the cells", "arms and stretches"],
                numeric={1, 2, 3, 4, 5, 6, 7}, wrap=True)

    # the genomes counted with the windows, read with and without them
    wc = k.get("windows")
    if wc:
        txt = (f"<p><strong>What the windows add.</strong> The {wc['n']:,} genomes counted with the windows were also read from the regions their earlier counts hold, against the same model. "
               f"A chromosome's standard error is {fmt(wc['se_ratio'][1], 1)} times smaller with the windows in the median chromosome ({fmt(wc['se_ratio'][0], 1)} to {fmt(wc['se_ratio'][2], 1)} from the 10th to the 90th percentile). ")
        if wc["both"]:
            txt += (f"{wc['both']:,} chromosomes, arms and stretches are called off their whole number in both readings" + (f", at levels that agree at r = {fmt(wc['r'], 4)}" if wc.get("r") is not None else "")
                    + (f"; the two differ by {fmt(wc['diff_sd'], 3)} copies (robust SD), {fmt(wc['diff_z'], 2)} of what their standard errors predict" if wc.get("diff_sd") is not None else "") + ". ")
        ow, ob = wc["only_windows"], wc["only_base"]
        txt += (f"{_n(len(ow)).capitalize()} are called only with the windows" + (" (below)" if ow else "") + f" and {_n(len(ob))} only without"
                + (": " + "; ".join(f"{esc(o['sample'])} {esc(o['label'])}" + (f" ({fmt(o['with_windows'], 2)} copies with the windows)" if o.get("with_windows") is not None else "") for o in ob[:6]) if ob else "") + ". ")
        P.h(txt + "</p>")
        if ow:
            P.table([[o["sample"], o["label"], fmt(o["copies"], 2), fmt(o["z"], 0), o["regions"], fmt(o["base_copies"], 2) if o.get("base_copies") is not None else "–", o["base_regions"]] for o in ow[:25]],
                    ["sample", "called only with the windows", "copies", "z", "regions", "the chromosome without them, copies", "its regions without"], numeric={2, 3, 4, 5, 6}, wrap=True)
    if tl and tl.get("groups"):
        txt = ("<p><strong>Chromosomes that follow their GC.</strong> The slope of the line across a genome's chromosomes, in copies per ten points of GC "
               "(chromosome 4 is 38% GC, chromosome 19 48%): "
               + "; ".join(f"over the {g['n']:,} genomes read {esc(g['label'])}, a robust SD of {fmt(g['sd'], 4)} where the fit's own error is {fmt(g['se'], 4)} "
                           f"in the median genome, and {g['beyond3']:,} beyond three of their standard errors (chance gives {fmt(g['expected3'], 0)})" for g in tl["groups"]) + ". "
               "Libraries differ in this, and the fit measures it best where the windows give the GC-rich chromosomes their regions. At the ends of the line, the GC-rich "
               "chromosomes would have read as lost or gained in a few percent of the cells without it: ")
        txt += "; ".join(f"{esc(e['sample'])} (slope {_signed(e['tilt'], 3)}, z {_signed(e['z'], 1)}): " + ", ".join(
            f"chromosome {c['chrom'][3:]} {fmt(c['without'], 3)} without the fit and {fmt(c['copies'], 3)} with it" for c in e["chromosomes"][:2]) for e in tl["ends"][:3]) + ". "
        txt += ("A chromosome takes no part in its own correction, and the fit's error is carried into each chromosome's. The slope and its z are in "
                "<code>karyotype.gc_tilt</code> and <code>karyotype.gc_tilt_z</code>.")
        P.h(txt + "</p>")
    al = k.get("alleles")
    if al and al.get("rows"):
        oc = al.get("one_copy") or []
        one_copy = ("" if not oc else
                    " Where one copy is expected and a second is in part of the cells (a line read as one X), heterozygous sites exist only where the second is, "
                    "and the alleles stand at its share to one: " + "; ".join(f"{esc(o['sample'])} {esc(o['event'])}, {fmt(o['by_alleles'], 2)} by the alleles" for o in oc[:6])
                    + (f"; and {len(oc) - 6} more" if len(oc) > 6 else "") + " (the caller misses the sites whose second allele is rarest, so these read high where the share is small).")
        P.h(f"""<p><strong>Against the alleles.</strong> Depth says how many copies; the alleles of the same reads say it again, independently. Two copies carry a heterozygous site's alleles in equal shares; a third copy in a share
f of the cells makes them (1 + f) : 1, a lost one 1 : (1 − f). In the reads fetched for the windows, heterozygous sites were called and the spread of their allele fractions, beyond what each site's depth gives by chance, turned into a share of the cells
(<code>analysis/karyotype/allele_balance.py</code>). {al["n"]:,} chromosomes, arms and stretches called here in a twentieth of the cells or more have enough sites{f", in {al['genomes']} genomes" if al.get("genomes") else ""}:
the two shares agree at r = {fmt(al["r"], 3)} and differ by {fmt(al["diff_sd"], 3)} (robust SD), {al.get("near", 0)} of them within a tenth of each other{f"; on the chromosomes that depth reads at two copies in the same genomes the alleles are balanced (d = {fmt(al['quiet_d'], 3)} in the median)" if al.get("quiet_d") is not None else ""}.
{("The alleles do not bear out " + _and([f"{esc(o['sample'])}'s {esc(o['event'])} ({fmt(o['by_alleles'], 2)} of the cells by the alleles)" for o in al.get("not_borne") or []]) + ": the extra copies may hold both homologs, or depth reads something of those libraries there as a gain; the calls are depth's.") if al.get("not_borne") else ""}
A change that depth reads in every cell and the alleles in none would be a copy with both alleles alike; none is seen. Below a twentieth of the cells the alleles cannot confirm a call: their own scatter is larger.{one_copy}</p>""")
        P.table([[r_["sample"], r_["event"], r_["sites"], fmt(r_["by_depth"], 2), fmt(r_["by_alleles"], 2)] for r_ in al["rows"][:30]],
                ["sample", "called by depth", "heterozygous sites", "share of the cells, by depth", "by the alleles"], numeric={2, 3, 4}, wrap=True, filter_box=True)

    # the autosomes
    byc = [c for c in au["by_chromosome"] if c["whole"] or c["part"]]
    if byc or au["stretches"]:
        top = sorted(byc, key=lambda c: -(c["whole"] + c["part"]))
        txt = (f"<p><strong>The autosomes.</strong> {au['genomes']:,} genomes hold a chromosome, an arm or a stretch off two copies. Whole chromosomes: {au['gained']:,} gained and {au['lost']:,} lost, "
               + ("most often " + _and([f"chromosome {c['chrom'][3:]} ({c['whole'] + c['part']})" for c in top[:5]]) if top else "none") + ". "
               "Gains of chromosomes 12 and 9 are what lymphoblastoid lines acquire in culture; nearly all are in part of the cells. ")
        if au["clones"]:
            c0 = au["clones"][0]
            txt += (f"Where several chromosomes are gained in the same share of the cells, one clone carries them all: {esc(c0['sample'])} reads {esc(', '.join(c0['chromosomes']))}, "
                    f"{len(c0['chromosomes'])} chromosomes within {fmt(c0['spread'], 2)} of each other" + (f"; {len(au['clones']) - 1} more genomes are like it" if len(au["clones"]) > 1 else "") + ". "
                    "That the shares agree is a check of the measurement that needs nothing from outside it. ")
        if au["recurrent"]:
            txt += ("Stretches that recur: " + "; ".join(f"chromosome {r['chrom'][3:]} {_mb(r['start'], r['end'])}, {'gained' if r['gain'] else 'lost'} in {r['n']} genomes" for r in au["recurrent"][:4]) + ". ")
        P.h(txt + "</p>")
        if au["stretches"]:
            shown = au["stretches"][:30]
            P.h(f"<p>Arms and stretches, the {len(shown)} furthest from two copies of {len(au['stretches']):,} (all events: <code>data/karyotype_events.tsv</code>):</p>")
            P.table([[s_["sample"], s_["chrom"][3:], {"p": "short arm", "q": "long arm", "pter": "to the end of the short arm", "qter": "to the end of the long arm", "inner": "inside"}.get(s_["span"], s_["span"]),
                      _mb(s_["start"], s_["end"]), s_["regions"], fmt(s_["copies"], 2), s_["cells"], fmt(s_["z"], 0)] for s_ in shown],
                    ["sample", "chromosome", "where", "span", "regions", "copies", "in", "z"], numeric={4, 5, 7}, wrap=True, filter_box=True)

    # the checks
    txt = "<p><strong>Checks.</strong> "
    a = ck.get("arms")
    if a:
        txt += (f"<em>The two arms.</em> A whole chromosome gained or lost shows in both arms alike. Of {a['n']} autosomes called off as a whole whose arms can each be read, both arms lie on the same side of two copies in {a['same_side']}"
                + (f", and their departures agree at r = {fmt(a['r'], 3)}" if a.get("r") is not None else "") + f"; the arms' difference, in its own standard errors, has a spread of {fmt(a['z_sd'], 2)}"
                + (f". Of the {a['small']} within a tenth of a copy of two, the smallest calls made, {a['small_same']} agree in direction and in {a['small_both']} each arm alone is two standard errors off" if a.get("small") else "") + ". ")
    t = ck.get("tails") or {}
    if t.get("rows"):
        r3 = t["rows"][0]
        txt += (f"<em>The quiet side.</em> A culture gains chromosomes and seldom loses an autosome, so the levels below two copies show the measurement's own scatter: of {t['n']:,} autosomes at two copies, "
                f"{r3['below']} lie more than 3 standard errors below, where a normal scatter gives {fmt(r3['expected'], 0)}; {r3['above']} lie as far above. ")
    f_ = ck.get("flags") or {}
    if f_.get("n"):
        txt += (f"<em>The estimator's flag.</em> The estimator sets a chromosome aside from its denominator when its pooled depth departs by 4%. Of its {f_['n']} flags, {f_['whole']} are whole chromosomes here"
                + (f", {f_['stretch']} are a stretch of the chromosome and not the whole of it ({esc(', '.join(f_['stretches'][:4]))}{' …' if len(f_['stretches']) > 4 else ''})" if f_.get("stretch") else "")
                + (f", {f_['neither']} are nothing here" if f_.get("neither") else "") + ". ")
    fe = k.get("fetch")
    if fe and fe.get("n"):
        txt += (f"<em>The fetch.</em> Read from the targeted fetch of the same file against the same model, {fe['identical']:,} of {fe['n']:,} genomes are written exactly as from the scan"
                + ("" if fe["identical"] == fe["n"] else "; the others: " + "; ".join(f"{esc(d['sample'])} {esc(d['scan'])} | {esc(d['fetch'])}" for d in fe["differ"][:5])) + ". ")
    P.h(txt + "</p>")
    dj = ck.get("dj") or []
    whole_dj = [d for d in dj if d["whole_chromosome"]]
    if dj:
        txt = ("<p><strong>The distal junction.</strong> Every acrocentric short arm carries one junction, so a whole chromosome 13, 14, 15, 21 or 22 more is a junction more, and a stretch of its long arm is not. "
               + "; ".join(f"{esc(d['sample'])}: {esc(d['events'])}, junction {fmt(d['dj'], 2)}" for d in dj[:8]) + ". ")
        if whole_dj:
            n_ok = sum(1 for d in whole_dj if abs((d["dj"] - 10) - d["extra"]) < 0.35)
            txt += f"For the {_n(len(whole_dj))} with a whole acrocentric chromosome off, the junction's level lies within 0.35 of ten plus the chromosome's excess in {_n(n_ok)}. "
        P.h(txt + "</p>")
    rel = ck.get("relatives") or {}
    if rel.get("pairs"):
        sh = rel["shared"]
        P.h("<p><strong>In the families.</strong> A change made in a culture is not in a relative's culture; one a child shares with a parent was in the germ line. "
            + (f"Of {rel['pairs']:,} pairs of a child with an event and a parent in the cohort, {_n(len(sh))} share one: "
               + "; ".join(f"{esc(s_['child'])} {esc(s_['child_event'])} and {esc(s_['role'])} {esc(s_['parent'])} {esc(s_['parent_event'])}" for s_ in sh[:6]) + "."
               if sh else f"None of the {rel['pairs']:,} pairs of a child with an event and a parent in the cohort shares one.") + "</p>")
    P.h("""<p><strong>What depth does not see.</strong> A change that leaves the number of copies as it was: a balanced translocation, an inversion, both copies of a chromosome from one parent. A whole genome in three copies,
which moves every chromosome alike. The short arms of the acrocentric chromosomes, which hold no single-copy sequence (their junctions are counted in section 3.2). A chromosome with fewer than three places is not read,
and an arm with fewer than five is not read on its own. A change in fewer cells than the table above gives is reported as its level and not called. Whether a change was in the donor or arose in the culture, the reads of one sample cannot say:
in cell lines most are the culture's, and the relatives are the test. These are not clinical karyotypes.</p>""")
    P.end()
