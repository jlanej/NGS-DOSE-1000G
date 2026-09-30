"""The distal-junction part of the page beyond the whole-copy steps: the profile along the unit (its scale, its core,
its polymorphic intervals), the integer copy states called along it, the comparison with the HPRC assemblies, what
the assemblies get wrong, and what the method now does and does not. Rendered inside section 3.2."""
from __future__ import annotations

import html

import numpy as np

from .report import fmt

esc = lambda x: html.escape(str(x))


def _kb(a, b):
    return f"{a // 1000}–{b // 1000} kb"


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _unexplained(t: dict) -> str:
    """One compared genome whose reads and assembly differ: what each holds, and what the call says of the reads."""
    bits = [f"assembly {fmt(t['assembly_core'], 0)}", f"reads {fmt(t['reads_core'], 2)}"]
    if t.get("status") in ("fractional", "uncertain"):
        bits.append(str(t["status"]))
    if _f(t.get("off")) is not None:
        bits.append(f"level {_signed(_f(t['off']))} from {t.get('copies')}")
    if t.get("fractional") not in (None, "", "none"):
        bits.append(str(t["fractional"]))
    return f"{esc(t['sample'])} ({esc(', '.join(bits))})"


def _signed(x, nd=2):
    if x is None:
        return "–"
    s = f"{x:+.{nd}f}"
    return s[1:] if float(s) == 0 else s.replace("-", "−")


def _pct(x, nd=0):
    return fmt(x, nd, pct=True) if x is not None else "–"


def _states(d: dict, floor=0.005) -> str:
    """'0: 59%, −1: 30%, …' for a share-per-state dict, the reference state first."""
    items = sorted(((int(k), float(v)) for k, v in d.items() if float(v) >= floor), key=lambda kv: (kv[0] != 0, -kv[1]))
    name = lambda k: "every copy holds it" if k == 0 else (f"{abs(k)} cop{'y lacks' if abs(k) == 1 else 'ies lack'} it" if k < 0 else f"{k} more")
    return ", ".join(f"{name(k)} in {100 * v:.0f}%" for k, v in items)


def gc_slope(eff: dict | None) -> tuple[float, float, int] | None:
    """The cohort's DJ window efficiencies against window GC: slope (log efficiency per unit GC) and r."""
    e = (eff or {}).get("DJ")
    if not e:
        return None
    a = np.array([np.nan if v is None else v for v in e["a"]], float)
    gc = np.array([np.nan if v is None else v for v in e["gc"]], float)
    m = np.isfinite(a) & np.isfinite(gc)
    if e.get("level"):
        m &= np.array(e["level"], bool)
    if m.sum() < 30:
        return None
    return float(np.polyfit(gc[m], a[m], 1)[0]), float(np.corrcoef(gc[m], a[m])[0, 1]), int(m.sum())


def lean_relations(rows: list[dict]) -> dict:
    """What the lean of a genome's profile goes with: the release batch (the share of genomes leaning more than 5%, per
    batch of fifty genomes or more), the depth and the library's GC response (correlations)."""
    def val(r, k):
        try:
            x = float(r.get(k))
        except (TypeError, ValueError):
            return np.nan
        return x
    tilt = np.array([val(r, "DJ.tilt") for r in rows])
    ok = np.isfinite(tilt)
    out = dict(n=int(ok.sum()), batches=[], r={})
    if ok.sum() < 50:
        return out
    batch = np.array([str(r.get("ngspca.batch") or "") for r in rows])
    for b in sorted({x for x in batch[ok] if x and x != "unknown"}):
        m = ok & (batch == b)
        if m.sum() >= 50:
            out["batches"].append(dict(batch=b, n=int(m.sum()), beyond_5=float(np.mean(np.abs(tilt[m]) > 0.05)), mean=float(tilt[m].mean())))
    for k, name in (("depth", "depth"), ("gc_rel_35", "GC response at 35%"), ("gc_rel_65", "GC response at 65%")):
        v = np.array([val(r, k) for r in rows])
        m = ok & np.isfinite(v)
        if m.sum() >= 50 and v[m].std() > 0:
            out["r"][name] = float(np.corrcoef(tilt[m], v[m])[0, 1])
    return out


def off_whole_numbers(P, dj: dict, data: dict):
    """The genomes that sit off whole numbers: what is measured, how many there are against chance, whom they are found
    in, what speaks for a change in part of the cells and what for the library, and the genomes themselves."""
    o = dj.get("off")
    if not o:
        return
    fr, lv, st = o.get("rules") or {}, o["level"], o["steps"]
    s = o["status"]
    lz, ez, mh = float(fr.get("level_z", 3.0)), float(fr.get("event_z", 4.0)), float(fr.get("min_height", 0.25))
    spread = fr.get("spread")
    P.h('<h3 id="djfractions">Off the whole numbers: a change in part of the cells?</h3>')
    P.h(f"""<p>Whole numbers are what a germ line holds. A junction lost or gained in part of the cells (a culture that is changing, a rearrangement in some of a donor's
blood) leaves the profile a fraction of a copy off them, as a lost X leaves a woman's X dosage between one and two. A chain of whole numbers would hide it: the nearest whole number is
called and the difference goes into the genome's scale. So the difference is measured and kept, in two ways.
<strong>The level.</strong> <code>DJ.off</code> is a genome's level less the whole numbers called, in copies, and <code>DJ.off_z</code> the same in robust SDs of the cohort's scales
({_pct(spread, 2) if spread else "–"} here{f", {fmt(o['spread_copies'], 2)} copies at ten" if o.get("spread_copies") else ""}; from {esc(str(fr.get("spread_from", "the class's rules")))}).
<strong>A stretch of the unit.</strong> In what the whole numbers leave, a level per stretch is fitted together with the genome's lean and its slope on the windows' GC
(a step is partly like a lean, and the profile of a library whose GC response the model left in follows GC); a step is proposed where it improves the fit as much as a change of state must, and it is kept where it
stands {fmt(ez, 0)} robust SDs from what the same fit finds in the cohort's other genomes and is {fmt(mh, 2)} copies high. The cohort is the measure because its profiles wander more than counting alone allows.
A call that is not uncertain is then <em>fractional</em> where the level lies {fmt(lz, 0)} SDs or more from its whole number or a stretch stands off; <code>DJ.fractional</code> gives the stretches with their height
against the copies the genome is described against, as <code>+0.53:106-400kb</code>.{"" if st.get("looked_for") else " In a cohort this small no step is looked for, and the level is judged against the spread the rules record."}</p>""")
    P.tiles([("Settled", f'{s.get("settled", 0):,}', f'of {o["n"]:,} calls: on their whole numbers'),
             ("Fractional", f'{s.get("fractional", 0):,}', f'{lv["n"]:,} by the level ({lv["below"]} below their whole number, {lv["above"]} above), {st["n"]:,} by a stretch of the unit'),
             ("Uncertain", f'{s.get("uncertain", 0):,}', f'between two whole numbers throughout; the level of {o["uncertain_off"]["beyond"]} of them lies {fmt(lz, 0)} SDs or more from the nearer'),
             ("By chance", fmt(lv["by_chance"], 0), f'levels {fmt(lz, 0)} SDs off that a normal scatter of the scales would give among {lv["of"]:,} genomes; {lv["n"]:,} are seen')])
    if lv.get("beyond"):
        P.table([[f'more than {fmt(b["copies"], 2)}', b["n"], b["below"], b["above"], fmt(b["by_chance"], 0)] for b in lv["beyond"]],
                ["level's distance from its whole number, copies", "genomes", "below", "above", "by chance"], numeric={1, 2, 3, 4})
        P.h(f"""<p class="small">Among the {lv["of"]:,} calls that are not uncertain; by chance: what a normal scatter with the robust SD of the offsets ({fmt(o["spread_copies"], 3)} copies) would give.
The tail is heavier than the scatter's, and heavier below the whole numbers than above them: a junction lost in part of the cells lies below.</p>""")
    who = lambda g: f'{esc(str(g["group"]))} {g["fractional"] + g["uncertain"]} of {g["n"]:,} ({_pct((g["fractional"] + g["uncertain"]) / g["n"], 1) if g["n"] else "–"})'
    sh, bp, osn = o.get("shared"), o.get("by_population"), o.get("other_signs") or {}
    txt = "<p><strong>Whom they are found in.</strong> Fractional or uncertain: " + "; ".join(who(g) for g in o["by_role"] if g["n"]) + ". "
    txt += "By sex: " + "; ".join(who(g) for g in o["by_sex"] if g["n"]) + ". "
    if len(o.get("by_batch") or []) > 1:
        txt += "By release batch: " + "; ".join(who(g) for g in o["by_batch"] if g["n"]) + ". "
    if osn.get("off") and osn.get("p") is not None:
        txt += (("They do not come with the other signs of a culture that has changed: " if osn["p"] > 0.05 else "Of the other signs of a culture that has changed: ")
                + f"an X or a Y lost in part of the cells, or a chromosome off its dosage, is seen in {osn['off_with']} of {osn['off']} of them and in "
                f"{osn['settled_with']:,} of {osn['settled']:,} settled genomes (Fisher exact p = {fmt(osn['p'], 2)}). ")
    P.h(txt + "</p>")
    if sh or bp:
        txt = "<p><strong>What the libraries make of it.</strong> A part of a genome's offset is not the genome's. "
        if sh:
            txt += (f"Spouses share no genes, and their offsets go together: r = {fmt(sh['spouses_r'], 2)} ({fmt(sh['spouses_lo'], 2)}–{fmt(sh['spouses_hi'], 2)}; {sh['n_trios']:,} trios), as a child's goes with its "
                    f"father's ({fmt(sh['child_father_r'], 2)}) and its mother's ({fmt(sh['child_mother_r'], 2)}). What a family shares here is how its samples were handled. ")
        if bp:
            ext = lambda g: f"{esc(g['pop'])} ({esc(g['superpop'])}) {_signed(g['mean'])}"
            txt += (f"Populations differ in their mean offset by more than their ancestry accounts for: {_pct(bp['share'], 0)} of the offset's variance lies between the {bp['n']} populations and "
                    f"{_pct(bp['share_superpop'], 0)} between the continental groups, from " + ", ".join(ext(g) for g in bp["lowest"]) + " to " + ", ".join(ext(g) for g in bp["highest"])
                    + " copies. The populations were collected and sequenced in groups. ")
        txt += ("A level judged against the whole cohort's spread therefore carries its group's offset with it, a tenth of a copy at most here, and a genome whose relatives are off the same way "
                "is more likely a matter of the libraries than of its cells: the table marks them. Calibration by batch is not done here.")
        P.h(txt + "</p>")
    if lv.get("with_relatives") or st.get("children_looked_at"):
        txt = "<p><strong>In the families.</strong> "
        if lv.get("with_relatives"):
            txt += (f"Of the {lv['n']} genomes whose level is off, {lv['with_relatives']} have a relative in the cohort, and in {lv['relatives_same_way']} of them a relative lies two SDs or more off the same way. ")
        if st.get("children_looked_at"):
            k = st["children_with_the_breakpoint"]
            txt += (("A step of fractional height is not passed on as a whole copy is: " if 2 * k < st["children_looked_at"] else "A step of fractional height in a parent, and its children: ")
                    + f"of the {st['children_looked_at']} children of a parent that carries one, {'none' if k == 0 else k} "
                    f"{'have' if k > 1 else 'has'} a breakpoint within 15 kb of the parent's, at a whole number or at a fraction. "
                    "A variant of the germ line would be in half of them; a change in part of a parent's cells, or a matter of the parent's library, in none. ")
        P.h(txt + "</p>")
    if st.get("breakpoints"):
        top = sorted(st["breakpoints"], key=lambda b: -b["n"])[:6]
        P.h(f"""<p><strong>Where the stretches end.</strong> {st["events"]} stretches in {st["n"]} genomes; their ends inside the unit fall most often near {", ".join(f"{b['kb']} kb ({b['n']})" for b in top)}.</p>""")
    rows_t = o.get("table") or []
    if rows_t:
        shown = rows_t[:40]
        ev = lambda t_: "–" if t_["fractional"] in (None, "", "none") else t_["fractional"].replace(";", "; ") + f' (z {str(t_["fractional_z"]).replace(";", "; ")})'
        P.h(f"""<p>The genomes, the {len(shown)} furthest off of {len(rows_t)} (all of them in <code>data/dj_fractional.tsv</code>): the level, the copies called, the level's distance from them, the stretches of
fractional height, what the chain of whole numbers called, the relatives in the cohort with their own offsets, and where a relative is off the same way.</p>""")
        P.table([[t_["sample"], t_["pop"] or "", t_["role"], t_["call"], fmt(t_["level"], 2), t_["copies"], f'{_signed(t_["off"])} ({_signed(t_["off_z"], 1)})', ev(t_),
                  (t_["whole_numbers"] or "").replace("none", "–").replace(";", "; "), (t_["relatives"] or "–"), t_["relatives_off_the_same_way"] or "–", t_["assembly"] or "–"] for t_ in shown],
                ["sample", "population", "in the pedigree", "call", "level", "copies", "level − copies (z)", "stretches of fractional height", "whole numbers called", "relatives and their offsets",
                 "off the same way", "assembly"], numeric={4, 5}, wrap=True, filter_box=True)
    P.h(f"""<p><strong>What a fraction is, and what is not seen.</strong> The measurement cannot say what a fraction is: a change in part of the cells, a library unlike the cohort's, or a copy of the germ line
that the panel reads in part. It says where the whole numbers do not hold. A level is found from {fmt(lz * (o.get("spread_copies") or 0), 2)} copies off its whole number, a stretch from about
0.4 copies over 30 kb or more, which for one copy is a change in two fifths to three fifths of the cells; nearer a whole number the whole number is called and nothing is said, so a change in a quarter of the cells
passes as none and one in three quarters as a whole copy. A step near the middle of the unit is partly a lean and is found less often than one near an end or a stretch inside the unit.</p>""")


def render(P, data: dict, rows: list[dict]):
    dj = data.get("dj")
    if not dj:
        return
    seg, st = dj["segments"], dj["stats"]
    asm = dj.get("assemblies")
    n = dj["n"]
    sc = st.get("scale") or {}
    pinned = sc.get("rule") == "mode" or sc.get("table_rule") == "mode"
    calls = st.get("calls")
    steps = (data.get("known_truth") or {}).get("DJ_steps") or {}
    hyper = seg.get("hyper_intervals") or []
    excl = seg.get("excluded") or []

    # ------------------------------------------------------------------ the profile: scale, core, polymorphic intervals
    P.h("<h3>The junction in profile: its scale, its core, and the intervals that vary</h3>")
    P.h(f'''<p>The junction is not one number. Within a genome the cohort-calibrated window estimates read the same copy number in every window that
a junction copy holds whole: a copy that lacks part of the unit lowers the windows it lacks, an extra partial copy raises the windows it holds. Across the
{n:,} genomes the cohort SD of the calibrated estimate per 5-kb sub-block has a floor of {fmt(seg["floor"], 2)} copies (counting noise) and rises far above it in
<strong>{", ".join(_kb(a, b) for a, b in hyper) if hyper else "no segment"}</strong> (figure {"below" if dj.get("figure") else "in <code>data/dj_segments.tsv</code>"}): there, junction copies differ.
Three things follow for the estimate, and the method's rules for this class (the bundle's <code>calibration.json</code>) carry them{"" if dj.get("rules") else "; this run had no such rules"}.</p>''')
    items = []
    if excl:
        items.append(f'''<li><strong>The level is set on the core.</strong> A genome's level is the median over the windows outside {", ".join(_kb(a, b) for a, b in excl)}
({len(seg.get("core_kb", []))} of the 20 blocks of 20 kb in full), so a copy that lacks the distal 22 kb does not move it. The core level lies within 0.2 copies of a whole number in
{_pct(st["core"]["within_0_2"], 1)} of genomes (the level over the whole unit, <code>DJ.cn_unit</code>, in {_pct(st["whole"]["within_0_2"], 1)}); the genomes near a whole number scatter about it with an SD of
{fmt(st["core"].get("mode_sd"), 3)} copies.</li>''')
    if pinned:
        items.append(f'''<li><strong>The scale is pinned to the cohort's mode.</strong> On the scale the fragment-GC model sets in the anchor windows the core of the main mode reads
{fmt(sc.get("mode_on_anchors"), 2)} copies ({sc.get("n_main", 0):,} genomes), not {fmt(sc.get("expected"), 0)}; every estimate is multiplied by {fmt(sc.get("factor"), 4)}.
The assemblies show that the deficit is the scale's and not the genomes' (below), and a saved efficiency table carries the pin to genomes counted later.</li>''')
    elif sc.get("rule") == "anchors" and dj.get("rules"):
        items.append("<li><strong>The scale is the anchors'.</strong> The cohort is too small to pin the scale to its mode; the level rests on the fragment-GC model in the anchor windows, "
                     "which reads this class a few percent low in the full cohort.</li>")
    offs = st.get("offsets") or []
    if offs:
        items.append("<li><strong>Polymorphic intervals are put on whole numbers.</strong> Where a deletion is common the cohort's median genome lacks part of a copy, and a median over genomes gives "
                     "the interval's windows an efficiency too low: every genome then reads the interval too high by the same factor, and a comb of whole numbers half a copy off. "
                     "The cohort's own comb gives the factor, with the reference state (every copy holding the interval) the highest that many genomes share: "
                     + "; ".join(f"{esc(o['name'])}, ×{fmt(float(np.exp(o['offset'])), 3)}" + (f" ({_states(o['states'])})" if o.get("states") else "") for o in offs)
                     + ". The assemblies agree with these factors to 0.02 (<code>analysis/dj/polymorphic_offsets.py</code>).</li>")
    if items:
        P.h("<ul>" + "".join(items) + "</ul>")

    # ------------------------------------------------------------------ the calls
    if calls:
        P.h("<h3>Whole numbers of copies along the unit</h3>")
        wh = calls.get("whole") or {}
        bp = sorted(calls["breakpoints"], key=lambda b: -b["n"])
        core_bp = [b for b in bp if not any(a - 5000 <= b["kb"] * 1000 <= e + 5000 for a, e in excl) and 5 < b["kb"] < 395][:6]
        holds = lambda k, v: (f"<strong>{v:,}</strong> hold{'s' if v == 1 else ''} {k} copies throughout" if k != 10 else f"{v:,} hold ten throughout")
        part = steps.get("partial") or {}
        md = steps.get("mendel") or {}
        chain = dj.get("chain") or {}
        tau, sig = chain.get("tau"), calls.get("sigma_median")
        event_cost = (f"an event, with its two changes of state, costs as much as {4 * tau * sig * sig:.0f} windows a copy off" if tau and sig
                      else "a change of state has a cost, so that a few windows off do not make an event")
        lr = lean_relations(rows)
        lean_txt = ""
        if lr["batches"] or lr["r"]:
            bits = []
            if len(lr["batches"]) > 1:
                bits.append("the release batch (more than 5% in " + " and ".join(f"{_pct(b['beyond_5'], 1)} of the {b['n']:,} genomes of batch {esc(b['batch'])}" for b in lr["batches"]) + ")")
            gc = [v for k, v in lr["r"].items() if k.startswith("GC")]
            if gc:
                bits.append(f"the library's GC response (|r| at most {fmt(max(abs(x) for x in gc), 2)})")
            if "depth" in lr["r"]:
                bits.append(f"depth (r = {_signed(lr['r']['depth'])})")
            lean_txt = "It goes weakly with " + ", ".join(bits[:-1]) + (" and " if len(bits) > 1 else "") + bits[-1] + ", and none of them accounts for it: a property of the library or of the culture. "
        P.h(f'''<p>Every genome's profile is read as a chain of whole numbers (<code>ngsdose.segments</code>): a window in state <em>k</em> is expected at <em>k</em> copies times the
genome's own scale, {event_cost}, and the most probable chain is found for each scale of a grid, a prior
(SD {_pct(chain.get("scale_sd", 0.015), 1)}, the level's SD among the genomes at ten) keeping the scale from explaining a copy away. One copy over about 10 kb, or two over 4 kb, can be found; the genomes' window noise is {fmt(calls["sigma_median"], 2)} copies per 250 bp.
Some genomes' profiles lean, rising or falling smoothly along the unit (SD {_pct(calls.get("tilt_mad_sd"), 1)} across it, robustly; more than 5% in {_pct(calls.get("tilt_beyond_5"), 1)} of genomes, more than 8% in
{_pct(calls.get("tilt_beyond_8"), 1)}). {lean_txt}A chain of whole numbers would break a lean into a step, so the lean is a parameter of the
chain, like the scale (prior SD {_pct(chain.get("tilt_sd", 0.03), 0)}): a smooth rise costs less as a lean than as a change of state, and a step, which a lean fits badly on both sides, stays a step.
A genome is described against ten copies where it holds ten over 40 kb or more of the core, otherwise against the state that holds most of it; a gain that reaches an end of the unit is a
<strong>partial copy</strong>, a loss that does a <strong>partial loss</strong> (a copy that lacks that end), anything else a local gain or loss. A genome whose level lies between two whole numbers
throughout has two readings, and the call says how far behind the second is: below three log units it is <em>uncertain</em>. A call that is not
uncertain is <em>settled</em> where the genome sits on its whole numbers and <em>fractional</em> where its level or a stretch of its unit sits off them (below).</p>
<p>Of {calls["n"]:,} genomes, {(calls.get("status") or {}).get("settled", calls["n"] - calls["scale_uncertain"]):,} are settled, {(calls.get("status") or {}).get("fractional", 0):,} fractional and {calls["scale_uncertain"]:,} uncertain. Of the settled, {", ".join(holds(k, v) for k, v in sorted(wh.items(), key=lambda kv: kv[0]) if v)}
(the polymorphic intervals aside): the whole-junction losses and gains. {calls["with_partial_copy"]:,} genomes carry a partial copy and {calls.get("with_partial_loss", 0):,} a partial loss;
{calls["with_large_event"]:,} hold an event of 40 kb or more of any kind, and {calls["off_integer"]} a segment of that length more than a third of a copy from its state. The genomes' scales have an SD of
{_pct(calls.get("scale_sd"), 1)}. Breakpoints recur: in the core, most at {", ".join(f"{b['kb']} kb ({b['n']})" for b in core_bp)}.
Every segment of every genome is in <code>data/dj_calls.tsv</code>, with its mean as the reads give it (<code>raw</code>) beside its mean on the genome's scale; <code>DJ.copies</code>, <code>DJ.partial</code>,
<code>DJ.variants</code>, <code>DJ.scale_f</code>, <code>DJ.call</code>, <code>DJ.off</code> and <code>DJ.fractional</code> are in the sample table.</p>''')
        if md.get("n_trios"):
            c_, p_ = md["core"], md["polymorphic"]
            P.h(f'''<p><strong>The calls are Mendelian.</strong> In the {md["n_trios"]:,} trios whose three calls are settled, position by position in {md["size_kb"]}-kb blocks (those within {md["edge_kb"]} kb of a
breakpoint of any of the three left out): a child's deviation from ten is a part of its father's plus a part of its mother's in <strong>{_pct(c_["share"], 2)}</strong> of {c_["blocks"]:,} blocks of the core
and {_pct(p_["share"], 1)} of {p_["blocks"]:,} blocks of the polymorphic intervals; {md["clean"]:,} trios have no block of the core out of place. Of the blocks in which a child deviates, a parent explains
{_pct(c_["explained"], 1)} in the core ({c_["child_dev"]:,} blocks) and {_pct(p_["explained"], 1)} in the polymorphic intervals. Where one parent deviates by one copy and the other not at all, the child has it in
{_pct(c_["passed_share"], 1)} of the blocks of the core ({c_["informative"]:,}) and {_pct(p_["passed_share"], 1)} of those of the polymorphic intervals ({p_["informative"]:,}), where one half is expected. With the blocks beside a
breakpoint left out, these are the blocks of events longer than 20 kb, inside the polymorphic intervals as outside them; the short deletions of those intervals are tested apart, below.</p>''')
        if part.get("kinds"):
            kinds = [k for k in part["kinds"] if k["n"] >= 5][:14]
            P.h(f'''<p>The copies that hold or lack an end of the unit, by their extent (ends snapped to the breakpoints that recur), with what the trios show. A parent's copy is looked for in the child by its
breakpoint, not by its name: a copy that holds the first 316 kb and a copy that lacks the last 84 kb are the same step down at 316 kb, described against different tens. Where a parent's call has one breakpoint in the core
and the other parent is at ten with none near it, the child's call has the breakpoint (within {part.get("tol_kb", 12)} kb, stepping the same way) in <strong>{part["transmitted"]} of {part["transmitted"] + part["not_transmitted"]}</strong>
pairs{f" (two-sided binomial p = {fmt(part['p_half'], 2)} against one half)" if part.get("p_half") is not None else ""}{f"; {part['several']} pairs whose parent has several breakpoints are not read" if part.get("several") else ""}{f"; {len(part['de_novo'])} child{'ren carry' if len(part['de_novo']) > 1 else ' carries'} a copy of this kind while both parents are at ten throughout ({', '.join(esc(x) for x in part['de_novo'][:8])}{' and others' if len(part['de_novo']) > 8 else ''})" if part.get("de_novo") else "; no child carries one while both parents are at ten throughout"}.</p>''')
            P.table([[k["kind"], k["n"], fmt(k["n"] / max(calls["n"], 1), 1, pct=True), k["transmitted"], k["not_transmitted"]] for k in kinds],
                    ["partial copy or loss (interval of the unit)", "genomes", "of the cohort", "transmitted", "not transmitted"], numeric={1, 2, 3, 4})
        inh, ln, npd = dj.get("inheritance") or [], steps.get("lines") or {}, steps.get("not_passed") or {}
        if inh:
            P.h(f"""<p><strong>What is passed on: the common deletions.</strong> The deletions of the polymorphic intervals are old and common, in the germ line beyond doubt, so they test the
measurement: a method that lost or invented copies would pass them to fewer than half of a carrier's children. Three readings, in the {inh[0]["n_trios"]:,} trios. Without a call, the child's value in the
interval (the median of its windows, less the genome's level) is regressed on the mean of its parents': under Mendel the slope equals the value's reliability, the share of its variance that is not counting noise
(which the core gives). By the values, a genome within 0.3 of a whole number of copies is given it; by the calls, a genome has the called state in the middle of the interval less the state beside it. In both, where
one parent lacks one copy and the other none, the child lacks one or none, and one half should lack one.</p>""")
            cnt = lambda c: "–" if not c else f'{c["passed"]} of {c["pairs"]} ({_pct(c["passed"] / c["pairs"])})'
            new_ = lambda c: "–" if not c or not c["both_at_zero"] else f'{c["new"]} of {c["both_at_zero"]}'
            ps = [c["p"] for r in inh for c in (r.get("by_value"), r.get("by_calls")) if c]
            short = lambda name: str(name).replace("distal ", "").replace("-", "–").replace(" ", "\u00a0")
            P.table([[short(r["name"]), r["windows"], f'{fmt(r["slope"], 2)} ({fmt(r["slope_lo"], 2)}–{fmt(r["slope_hi"], 2)})', fmt(r["reliability"], 2), cnt(r.get("by_value")), cnt(r.get("by_calls")), new_(r.get("by_calls")),
                      f'{_signed(r["child_minus_midparent"])} ± {fmt(r["child_minus_midparent_se"], 2)}'] for r in inh],
                    ["interval", "windows", "slope (95%)", "reliability", "passed on (values)", "passed on (calls)", "new in child", "child − parents"],
                    numeric={1, 2, 3, 4, 5, 6, 7}, wrap=True)
            shifted = [r for r in inh if abs(r["child_minus_midparent"]) > 2.5 * r["child_minus_midparent_se"]]
            P.h(f"""<p class="small">Slope: of the child's value on the mean of its parents', with its 95% interval by resampling families. Passed on: the children that lack one copy, of those with one parent that lacks
one and one that lacks none, by the values and by the calls{f"; no count differs from one half (two-sided binomial p from {fmt(min(ps), 2)} to {fmt(max(ps), 2)})" if ps and min(ps) > 0.05 else (f"; two-sided binomial p against one half from {fmt(min(ps), 3)} to {fmt(max(ps), 2)}" if ps else "")}.
New in child: a deletion called in the child of two parents without it, which is the calls' error in these short intervals, where an event is near the limit of what a chain of whole numbers can find.
Child − parents: the child's value less the mean of its parents', in copies (mean ± SE over the trios), zero under Mendel{"; in " + " and ".join(esc(short(r["name"])) for r in shifted) + " it is not: the libraries of the two generations differ in these windows (most children were sequenced in a later release batch than their parents), and a count by whole numbers carries the difference as error" if shifted else ""}.</p>""")
        if ln.get("pairs"):
            wt, wn = steps.get("transmitted", 0), steps.get("transmitted", 0) + steps.get("not_transmitted", 0)
            txt = (f"<p><strong>What is passed on: whole copies and partial copies.</strong> These pass to fewer children than one half. A whole-copy step passes to {wt} of {wn}, a partial copy (one breakpoint) to "
                   f"{part.get('transmitted', 0)} of {part.get('transmitted', 0) + part.get('not_transmitted', 0)}, together <strong>{ln['transmitted']} of {ln['pairs']}</strong>"
                   + (f" (two-sided binomial p = {fmt(ln['p_half'], 3)} against one half)" if ln.get("p_half") is not None else "") + ". "
                   + (f"The readings leave little room: where a whole-copy step was not passed on, the parent's level lies within 0.3 copies of its whole number in {npd['parent']['within_0_3']} of {npd['parent']['n']} pairs "
                      f"({_signed(npd['parent']['mean'])} ± {fmt(npd['parent']['sd'], 2)}, mean ± SD) and the child's within 0.3 of ten in {npd['child']['within_0_3']} ({_signed(npd['child']['mean'])} ± {fmt(npd['child']['sd'], 2)}). " if npd else ""))
            if ln.get("expected_new") is not None:
                txt += (f"Nor did the events arise in the cell lines, as far as the children can show: {ln['parents_carrying']:,} of the {ln['parents']:,} parents carry one, and had the excess over half-transmission "
                        f"({_pct(ln['excess_share'])} of their events) arisen in culture, the lines of the {ln['both_plain']:,} children of two parents at ten throughout would hold about {fmt(ln['expected_new'], 0)} new ones; "
                        f"they hold {len(ln['new'])}" + (f" ({', '.join(esc(x) for x in ln['new'][:6])})" if ln["new"] else "") + ". "
                        "What remains is a difference between the generations or in the germ line: a change in a donor's blood that a clonal line makes whole, the likelier the older the donor, "
                        "or a variant passed on less often than chance. Neither is tested here, and the deficit rests on the whole-copy steps, a father's losses most of all.")
            P.h(txt + "</p>")

    off_whole_numbers(P, dj, data)
    if not asm:
        P.h('<p class="small">The comparison with HPRC assemblies appears when <code>--dj-assemblies</code> names the tables of <code>pipeline/hprc_dj.py</code>.</p>')
        return

    # ------------------------------------------------------------------ against the assemblies
    S, A = asm["stats"], asm["samples"]
    res = S["resolved"]
    k10 = S["known10"]
    un10 = [x for x in k10.get("reads_unpinned") or [] if x is not None]
    sc_calls = S.get("calls") or {}
    P.h(f'<h3>Against the HPRC release-2 assemblies of {S["n"]} cohort members</h3>')
    P.h(f'''<p>For {S["n"]} genomes of the cohort with an HPRC release-2 assembly ({S["n_haplotypes"]} haplotypes, chosen to cover the whole-copy
steps, the values between steps and the mode; <code>pipeline/07_hprc_dj.sh</code>), every k-mer of the panel was counted in each haplotype FASTA
(<code>ngs-dose panel --report</code>: a complete junction copy contributes one; single-nucleotide differences lower single k-mers, so the
median per block is robust to them), and the assembly was aligned to the unit with everything outside a panel k-mer masked, so that each junction copy
appears as a run of alignments on one contig ({S["n_copies"]} copies: {", ".join(f"{v} {k}" for k, v in S["copy_classes"].items())}). The assembly's copies per
block are the two haplotypes' block medians summed; an assembly counts as <em>resolved</em> at the junction when no copy is cut by a contig end inside
the unit and at most one copy is partial.</p>''')
    tiles = [("Genomes compared", f'{S["n"]}', f'{S["n_resolved"]} with the junction resolved in both haplotypes'),
             ("Core level, resolved assemblies", f'{fmt(res.get("core_diff_mean"), 2)} ± {fmt(res.get("core_diff_sd"), 2)}', f'reads − assembly, copies, mean ± SD over {res.get("n", 0)} genomes; within half a copy in {res.get("core_within_0_5", 0)} of {res.get("n", 0)}')]
    rc = sc_calls.get("resolved_certain") or {}
    if rc.get("concordance") is not None:
        tiles.append(("Called states", _pct(rc.get("concordance_same_level"), 1), f'of {rc.get("sub_blocks_same_level", 0):,} 5-kb sub-blocks, the call equals the assembly, in the {rc.get("bulk_equal", 0)} resolved genomes with a settled call whose level the assembly shares'
                      + (f'; in {len(rc["bulk_differs"])} ({", ".join(rc["bulk_differs"])}) the reads hold a copy more or less than the assembly throughout' if rc.get("bulk_differs") else "")))
    if sc_calls.get("partial_in_assemblies"):
        off = sc_calls.get("partial_off_kb")
        tiles.append(("Partial copies", f'{sc_calls["partial_found"]} of {sc_calls["partial_in_assemblies"]}',
                      "partial copies of the resolved assemblies called"
                      + (f", the breakpoint inside the unit within {off} kb of the assembly's in each" if off is not None else "")
                      + " (an end in the distal 30 kb, where copies themselves begin at 3 or at 22 kb, is read as the unit's start)"))
    if k10["n"] and un10 and pinned:
        tiles.append(("Ten-copy genomes by assembly", f'{fmt(float(np.median(un10)), 2)}', f'median on the fragment-GC model\'s scale in the {k10["n"]} genomes whose resolved assembly holds ten complete copies (range {fmt(min(un10), 2)}–{fmt(max(un10), 2)}): what the pin corrects'))
    if S["partial_316"]:
        tiles.append(("The recurrent gain", f'{S["partial_316"]} of {S["n_partial"]}', "partial copies in the assemblies span exactly the first 316 kb of the unit, in tandem with a complete copy"))
    P.tiles(tiles)
    pts = [dict(x=t["assembly_mean"], y=t["reads_mean"], label=t["sample"], si=0 if t["resolved"] else 1,
                extra=[f"assembly {t['haplotypes']}", f"NGS-DOSE {fmt(t['reads_cn'], 2)}; {t.get('copies')} complete copies, partial {t.get('partial') or 'none'}" if t["reads_cn"] is not None else "", t["group"]]) for t in A]
    P.chart("djasm", dict(type="scatter", points=pts, legend=["assembly resolved at the junction", "assembly fragmented at the junction"], xlabel="HPRC assembly: copies per 20-kb block, mean over the unit",
                          ylabel="NGS-DOSE: calibrated profile, mean over the unit", identity=True),
            "The distal junction against long-read assemblies, genome by genome",
            f"Each dot is one of {S['n']} genomes of the 1000 Genomes 30× cohort with an HPRC release-2 assembly: the assembly's copies per 20-kb block (both haplotypes' "
            f"k-mer block medians summed) averaged over the unit, against NGS-DOSE's calibrated profile averaged the same way; the line is identity. "
            f"Filled: the junction assembled without contig breaks in both haplotypes; open: fragmented. Over the resolved genomes the reads read {fmt(res.get('mean_diff_mean'), 2)} ± {fmt(res.get('mean_diff_sd'), 2)} "
            f"copies from the assembly (mean ± SD); the fragmented ones scatter because the assembly does, not the reads.")
    if dj.get("figure"):
        P.h(f'<figure><img src="{esc(dj["figure"])}" alt="the distal junction: segment map, per-genome agreement, and profiles and calls against the assemblies" style="width:100%">'
            f'<figcaption>The junction in profile. A: cohort SD of the calibrated estimate per 5-kb sub-block over {n:,} genomes, the noise floor dashed, the hyper-variable segments shaded. '
            f'B: as the chart above. C: every compared genome in 20-kb blocks: the assembly, NGS-DOSE\'s calibrated profile, and the integer states called from it; blue below ten copies and red above; hatched blocks are '
            f'left out of the level. Redrawn from <code>report.json</code> and <code>data/dj_hprc_blocks.tsv</code> by <code>python -m report.dj_figure</code>.</figcaption></figure>')
    P.h("<p>Genome by genome (the assembly's haplotypes as the k-mer median over the unit, and its core as the median over the core blocks; NGS-DOSE's level, complete copies, "
        "the copies it holds or lacks an end of, its other events and its scale; the share of 5-kb sub-blocks in which the called state equals the assembly's; the last column is what Rhie et al. 2026 report for the same assemblies):</p>")
    hdr = ["sample", "assembly (k-mer medians)", "assembly core", "partial / truncated copies in the assembly", "NGS-DOSE level", "copies", "holds (+) or lacks (−) an end", "other events", "scale", "calls = assembly", "resolved", "Rhie et al. 2026"]

    def others(t):
        v = [x for x in str(t.get("variants") or "").split(";") if x and x != "none" and x not in str(t.get("partial") or "").split(";")]
        return "; ".join(v[:6]) + (" …" if len(v) > 6 else "") or "–"
    P.table([[t["sample"], t["haplotypes"], fmt(t["assembly_core"], 0), "; ".join(x for x in (t["partial_extents"], t["truncated_extents"]) if x) or "–",
              fmt(t["reads_cn"], 2) if t["reads_cn"] is not None else "–", t.get("copies", "–"), (t.get("partial") or "–").replace("none", "–"), others(t),
              (fmt(t.get("scale_f"), 3) + (" (uncertain)" if t.get("scale_uncertain") else " (fractional)" if t.get("status") == "fractional" else "")) if t.get("scale_f") is not None else "–",
              _pct(t.get("call_concordance")) if t.get("call_concordance") is not None else "–", "yes" if t["resolved"] else "no", t["rhie"] or "–"]
             for t in sorted(A, key=lambda t: (t["assembly_mean"], t["sample"]))], hdr, numeric={2, 4, 5, 8, 9}, flagged=lambda r: r[10] == "no", wrap=True)
    P.h('<p class="small">Every block of every compared genome: <code>data/dj_hprc_blocks.tsv</code>; every junction copy the alignments found, with its contig, extent and class: '
        '<code>data/dj_hprc_copies.tsv</code>; the per-genome table: <code>data/dj_hprc.tsv</code>.</p>')

    # ------------------------------------------------------------------ what the assemblies get wrong
    odd = asm.get("oddities") or []
    P.h("<h3>What the assemblies get wrong at the junction, and what the reads show</h3>")
    trunc = [o for o in odd if o["truncated"]]
    over = [o for o in odd if o["over_blocks"] and not o["resolved"]]
    under = [o for o in odd if o["under_blocks"] and o["truncated"]]
    split = [o for o in odd if o["phasing_split"] >= 3]
    kva = [o for o in odd if o["unit_vs_blocks"]]
    unexplained = [t for t in A if t["resolved"] and (abs(t["diff_core"]) >= 0.35 or t.get("scale_uncertain") or t.get("status") == "fractional")]
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
        P.h(f'''<p>Where a resolved assembly and the reads disagree by a third of a copy or more on the core, or the call is fractional or uncertain,
{", ".join(_unexplained(t) for t in unexplained)},
neither side is confirmed: a copy the assembly does not hold, a junction lost or gained in part of the cell line (the assembly's DNA and the reads' DNA are different cultures of the same line), or a genome whose
level the cohort model misplaces are all possible. A change that one culture holds in every cell and the other in part of its cells would read as these do: a whole number in the assembly, a fraction in the reads. {"None of them has a chromosome flagged in the control regions (a chromosome 4% off its expected dosage is), so a mosaic loss or gain of a whole acrocentric is not the reason." if not any(t["flagged_chromosomes"] for t in unexplained) else "Flagged in the control regions: " + "; ".join(esc(t["sample"]) + " (" + esc(t["flagged_chromosomes"]) + ")" for t in unexplained if t["flagged_chromosomes"]) + "."}</p>''')
    n_ds = S["n_distal_start"]
    by_name = {o["name"]: o for o in offs}
    d22 = (by_name.get("distal 5-15 kb") or {}).get("states")
    d200 = (by_name.get("197-217 kb") or {}).get("states")
    P.h(f'''<p>Two structures of the junction copies themselves recur across the assemblies, and the calls find them across the cohort. <strong>Copies that begin about 22 kb into the unit</strong>
({n_ds} of {S["n_copies"]} copies, in {sum(1 for t in A if t["n_distal_start"])} of {S["n"]} genomes): the distal 22 kb is a deletion polymorphism{f"; across the cohort, at 5–15 kb, {_states(d22)}" if d22 else ""}.
And <strong>a partial copy holding the first 316 kb of the unit, in tandem with a complete copy</strong>
({S["partial_316"]} of the {S["n_partial"]} partial copies; in every gain genome with a resolved assembly): the recurrent gain is not an eleventh junction but a duplication of most of one,
which reads 11 over the first 316 kb and 10 beyond. A third is a deletion of 197–217 kb in single copies (one copy of HG00097's first haplotype lacks 193–222 kb){f": across the cohort, {_states(d200)}" if d200 else ""}.</p>''')

    # ------------------------------------------------------------------ what the method does now, and what remains
    gs = gc_slope(data.get("efficiencies"))
    cap = ((data.get("modes") or {}).get("capture") or {}).get("DJ") or {}
    P.h("<h3>What the method does now, and what remains</h3>")
    items = []
    if pinned:
        items.append(f"<li><strong>The scale.</strong> Pinned to the cohort's mode (×{fmt(sc.get('factor'), 4)})"
                     + (f"; the {k10['n']} genomes whose resolved assembly holds ten complete copies read {fmt(min(un10), 2)}–{fmt(max(un10), 2)} before it (median {fmt(float(np.median(un10)), 2)}; some were chosen for reading between steps)" if k10["n"] and un10 else "")
                     + (f". The cohort's window efficiencies fall with window GC (slope {fmt(gs[0], 2)} in log efficiency per unit GC, r = {fmt(gs[1], 2)} over {gs[2]:,} windows), as the 45S unit's do: a residual of the GC model in repeat context that the GC-rule anchors do not remove" if gs else "")
                     + ". A cohort of fewer than fifty genomes is not pinned; a genome counted alone takes the pin from a saved efficiency table of its chemistry, and without one its level rests on the GC model, a few percent low.</li>")
    items.append("<li><strong>The calls.</strong> Copies, partial copies and losses and local events are whole numbers with their breakpoints; the level remains beside them, because a genome whose level lies between two whole numbers "
                 f"throughout ({calls['scale_uncertain'] if calls else 0} here, called uncertain) is either changing in culture or misplaced by the model, and the calls cannot say which. Events shorter than about 10 kb at one copy are not called: "
                 "the small polymorphisms at 130–135, 160–165 and 225–230 kb show in the segment map and are left out of the level, not genotyped.</li>")
    items.append("<li><strong>Read counting is not the limit.</strong> " + (f"The targeted fetch reads {_pct(cap.get('median'), 2)} of the junction's reads (median over {cap.get('n', 0):,} genomes; section 3.3); " if cap.get("n") else "")
                 + "duplicate-flagged reads are counted, so the ten-fold pile-up of junction reads on GRCh38 costs nothing; a read is classified by any four of its k-mers, so sequence differences among the copies, "
                 "which lower an assembly's exact-k-mer counts by 5% per haplotype, do not lower the reads' count.</li>")
    items.append("<li><strong>Assemblies as truth, with care.</strong> Of the compared genomes the junction was fully resolved in " + f"{S['n_resolved']} of {S['n']}; for the rest the reads are the more coherent reading of the locus, "
                 "and the profile identifies which blocks the assembly broke. A screen of all 200 assembled cohort members (about four minutes a haplotype, most of it the download) would give this class a truth panel of 200.</li>")
    P.h("<ul>" + "".join(items) + "</ul>")
