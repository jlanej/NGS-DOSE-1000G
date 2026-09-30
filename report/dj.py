"""The distal junction beyond one number: its profile along the 400-kb unit in every genome, the
segments of the unit that vary between people, the integer copy states the cohort layer called
(ngsdose.segments: complete copies, partial copies, local gains and losses), and the comparison
with the HPRC release-2 assemblies of cohort members (pipeline/07_hprc_dj.sh, pipeline/hprc_dj.py).

The profiles are the cohort layer's own (`cohort_table(..., profiles=...)`): each window's estimate
divided by its efficiency, on a scale pinned to the cohort's mode, the polymorphic intervals'
efficiencies set on the cohort's comb of integers. Within a genome they read the same copy number
in every window that a junction copy holds whole; a copy that lacks part of the unit lowers the
windows it lacks, and an extra partial copy raises the windows it holds.
"""
from __future__ import annotations

import csv
import math
import re
import warnings
from collections import Counter
from pathlib import Path

import numpy as np

UNIT = 400000
SUB = 5000                       # sub-block of the segment map (matches pipeline/hprc_dj.py)
BLOCK = 20000                    # block of the profiles and of the assembly comparison
NSUB, NBLOCK = UNIT // SUB, UNIT // BLOCK
HYPERVARIABLE = 1.4              # a sub-block whose cohort SD exceeds this multiple of the median SD varies between people
MIN_WINDOWS = 2                  # windows a sub-block needs to have a value
EXPECTED = 10.0

# what Rhie et al. 2026 (bioRxiv, DJCounter) report for cohort members from the same assemblies
RHIE = {"HG00621": "DJ lost entirely (paternal), with the HSat3 distal to it and the start of the rDNA array",
        "HG01891": "DJ lost entirely (maternal), with part of the HSat3; a rare rDNA inversion at the start of the array",
        "HG01981": "9.1-9.4 copies; a partial loss of the DJ flank",
        "NA20752": "a duplicated DJ (with NA20805), the second copy partial, breakpoint inside the flank",
        "NA20805": "a duplicated DJ (with NA20752), the second copy partial, breakpoint inside the flank",
        "HG03654": "a duplicated block starting with a smaller portion of the DJ flank, down to the ACRO repeat",
        "HG00320": "partially assembled among rDNA; the partial DJ holds one arm of the palindrome, flank intact to the rDNA"}


def _blocks(starts, vals, size, n, min_windows=MIN_WINDOWS):
    idx = (starts // size).astype(int)
    out = np.full(n, np.nan)
    for i in range(n):
        m = (idx == i) & np.isfinite(vals)
        if m.sum() >= min_windows:
            out[i] = np.median(vals[m])
    return out


def profiles(lib: dict | None) -> dict | None:
    """The library's profiles of one class in blocks: per genome the 5-kb and 20-kb medians, the level on the core and
    over the whole unit, and its call."""
    if not lib or not len(lib.get("samples") or []):
        return None
    starts = np.asarray(lib["start"])
    cn = np.asarray(lib["cn"], float)
    core = np.asarray(lib["level"], bool)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        level = np.exp(np.nanmedian(np.log(cn[:, core]), axis=1))
        unit = np.exp(np.nanmedian(np.log(cn), axis=1))
    return dict(samples=list(lib["samples"]), sub=np.array([_blocks(starts, x, SUB, NSUB) for x in cn]), block=np.array([_blocks(starts, x, BLOCK, NBLOCK, 3) for x in cn]),
                level=level, unit=unit, calls=lib.get("calls") or {}, scale=lib.get("scale") or {}, offsets=lib.get("offsets") or [],
                core_windows=core, starts=starts, cn=cn)


def segments(P: dict, rules: dict | None = None) -> dict:
    """Which parts of the unit vary between people. Cohort SD of the calibrated estimate per 5-kb sub-block; the noise
    floor is the typical sub-block's SD, and a sub-block far above it holds a structural polymorphism. The core is
    what the class's rules leave in the level (`level_exclude`; a 20-kb block is core when it overlaps none of them)."""
    S = P["sub"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)          # sub-blocks without a window (no core k-mers) are all-NaN
        sd, med = np.nanstd(S, axis=0), np.nanmedian(S, axis=0)
    valid = np.isfinite(sd)
    floor = float(np.nanmedian(sd[valid])) if valid.any() else float("nan")
    hyper = valid & (sd > HYPERVARIABLE * floor) if np.isfinite(floor) else np.zeros(NSUB, bool)
    runs, i = [], 0
    while i < NSUB:
        if hyper[i]:
            j = i
            while j + 1 < NSUB and hyper[j + 1]:
                j += 1
            runs.append((i * SUB, (j + 1) * SUB))
            i = j + 1
        else:
            i += 1
    excluded = [(int(a), int(b)) for a, b in (rules or {}).get("level_exclude", [])]
    has = lambda b: bool(np.isfinite(P["block"][:, b]).any())
    core_blocks = [b for b in range(NBLOCK) if has(b) and not any(a < (b + 1) * BLOCK and e > b * BLOCK for a, e in excluded)]
    return dict(sub_kb=[i * SUB // 1000 for i in range(NSUB)], sd=[None if not np.isfinite(x) else round(float(x), 3) for x in sd],
                median=[None if not np.isfinite(x) else round(float(x), 3) for x in med], floor=floor, hypervariable=[bool(x) for x in hyper],
                hyper_intervals=runs, core_blocks=core_blocks, excluded=excluded)


def _status(c) -> str:
    """settled, fractional or uncertain (a call made before the library knew fractions is settled or uncertain)."""
    return getattr(c, "status", "uncertain" if c.uncertain else "settled")


def cohort_stats(P: dict, seg: dict) -> dict:
    """The cohort's levels against whole numbers, its calls (complete copies, partial copies, breakpoints, the genomes
    whose scale leaves their integers open) and the block states."""
    def closeness(x):
        x = x[np.isfinite(x)]
        r = np.abs(x - np.round(x))
        near = np.abs(x - np.round(x)) < 0.35
        return dict(n=int(len(x)), within_0_2=float(np.mean(r <= 0.2)), within_0_3=float(np.mean(r <= 0.3)), between=float(np.mean(r > 0.3)),
                    mode_sd=float(np.std(x[near] - np.round(x[near]), ddof=1)) if near.sum() > 2 else None)
    B = P["block"]
    d = np.diff(B, axis=1)
    out = dict(scale=P["scale"], offsets=P["offsets"], core=closeness(P["level"]), whole=closeness(P["unit"]), core_median=float(np.nanmedian(P["level"])),
               breakpoints={(i + 1) * BLOCK // 1000: int(np.sum(np.abs(d[:, i]) >= 0.6)) for i in range(NBLOCK - 1)})
    calls = P["calls"]
    if calls:
        C = [calls[s] for s in P["samples"] if s in calls]
        f = np.array([c.scale for c in C])
        bp = Counter()
        for c in C:
            for a, b in zip(c.segments[:-1], c.segments[1:]):
                bp[(a.end // SUB) * SUB // 1000] += 1
        nseg = Counter(len(c.segments) for c in C)
        tl = np.array([c.tilt for c in C])
        # an event's extent as the sample table prints it (whole kb), so that these counts are the ones a reader gets from DJ.variants
        big = lambda c: [e for e in c.events if e.end // 1000 - e.start // 1000 >= 40]
        out["calls"] = dict(n=len(C), copies={int(k): int(v) for k, v in sorted(Counter(c.copies for c in C).items())},
                            flat={int(k): int(v) for k, v in sorted(Counter(c.copies for c in C if not big(c)).items())},
                            with_large_event=int(sum(1 for c in C if big(c))), with_partial_copy=int(sum(1 for c in C if any(e.kind == "partial copy" for e in c.events))),
                            whole={int(k): int(v) for k, v in sorted(Counter(c.copies for c in C if not big(c) and _status(c) == "settled").items())},
                            status={k: int(v) for k, v in Counter(_status(c) for c in C).items()},
                            with_partial_loss=int(sum(1 for c in C if any(e.kind == "partial loss" for e in c.events))),
                            scale_sd=float(f.std(ddof=1)) if len(f) > 1 else None, scale_uncertain=int(sum(c.uncertain for c in C)),
                            off_integer=int(sum(1 for c in C if any(s.off_integer and s.end - s.start >= 40000 for s in c.segments))),
                            segments={int(k): int(v) for k, v in sorted(nseg.items())},
                            tilt_sd=float(tl.std(ddof=1)) if len(tl) > 1 else None, tilt_mad_sd=float(1.4826 * np.median(np.abs(tl - np.median(tl)))),
                            tilt_beyond_5=float(np.mean(np.abs(tl) > 0.05)), tilt_beyond_8=float(np.mean(np.abs(tl) > 0.08)),
                            breakpoints=[dict(kb=int(k), n=int(v)) for k, v in sorted(bp.items())], sigma_median=float(np.median([c.sigma for c in C])))
    return out


def _binomial_two_sided(k: int, n: int, p: float = 0.5) -> float:
    from math import comb
    pk = [comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(n + 1)]
    return float(min(1.0, sum(x for x in pk if x <= pk[k] * (1 + 1e-9))))


def inheritance(P: dict, ped: dict | None, rules: dict | None, min_trios: int = 30, n_boot: int = 2000, seed: int = 1, flank: int = 5000) -> list[dict] | None:
    """Whether the junction's common polymorphisms are passed on as a germline variant is: the test of the measurement
    that the large events cannot give. For each polymorphic interval of the class's rules, three ways.

    By regression, without a call: the child's value (the median of the interval's calibrated windows, less the
    genome's level) on the mean of its parents'. Under Mendel the slope is the value's reliability, the share of its
    variance that is not counting noise; the noise is that of the same number of core windows in the genomes whose
    call holds no event (robustly, by the median absolute deviation). By the values: a genome within 0.3 of a whole number of copies is given it, and where one
    parent lacks one copy and the other none, the child lacks one or none. By the calls: the same count on the called
    state in the middle of the interval less the state beside it (the first core position beyond the stretch the level
    leaves out, on either side; a genome whose two sides differ is not counted), in the trios whose calls are settled.
    In the last two, `new` counts the children that lack a copy where neither parent does. `child_minus_midparent` is
    the mean of the child's value less its parents' mean: zero under Mendel, and otherwise a difference between the
    libraries of the two generations in these windows."""
    if not ped or not rules or not rules.get("polymorphic") or P.get("cn") is None:
        return None
    idx = {s: i for i, s in enumerate(P["samples"])}
    trios = np.array([(idx[c], idx[q.get("father")], idx[q.get("mother")]) for c, q in ped.items()
                      if c in idx and q.get("father") in idx and q.get("mother") in idx], int)
    if len(trios) < min_trios:
        return None
    cn, starts, level, calls = P["cn"], P["starts"], P["level"], P["calls"]
    quiet = np.array([s in calls and not calls[s].events and _status(calls[s]) == "settled" and len(calls[s].segments) == 1 for s in P["samples"]])
    if quiet.sum() < 30:
        quiet = np.ones(len(P["samples"]), bool)                          # without calls, every genome; the noise is measured robustly either way
    mad_var = lambda x: float((1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))) ** 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        has = np.isfinite(cn).any(axis=0)                                  # the unit's windows that were counted (the others hold no core k-mer)
    core = np.flatnonzero(np.asarray(P["core_windows"], bool) & has)
    excluded = [(int(a), int(b)) for a, b in rules.get("level_exclude", [])]

    def value(cols):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmedian(cn[:, cols], axis=1) - level

    def count(dev):
        """Parent at -1, the other at 0: the child at -1 (passed) or at 0. Both at 0: the child at -1 (new)."""
        passed = pairs = both = new = 0
        by = dict(father=[0, 0], mother=[0, 0])
        for c, f, m in trios:
            dc, df, dm = dev[c], dev[f], dev[m]
            if dc is None or df is None or dm is None:
                continue
            for who, x, y in (("father", df, dm), ("mother", dm, df)):
                if x == -1 and y == 0 and dc in (0, -1):
                    pairs += 1
                    passed += dc == -1
                    by[who][0] += dc == -1
                    by[who][1] += 1
            if df == 0 and dm == 0:
                both += 1
                new += dc == -1
        if not pairs:
            return None
        return dict(passed=int(passed), pairs=int(pairs), p=_binomial_two_sided(int(passed), int(pairs)), both_at_zero=int(both), new=int(new),
                    father=f"{by['father'][0]} of {by['father'][1]}", mother=f"{by['mother'][0]} of {by['mother'][1]}")
    rng = np.random.default_rng(seed)
    out = []
    for iv in rules["polymorphic"]:
        a, b = iv["interval"]
        cols = np.flatnonzero((starts >= a) & (starts < b) & has)
        if len(cols) < 8:
            continue
        d = value(cols)
        noise = float(np.median([mad_var(value(core[k:k + len(cols)])[quiet]) for k in range(0, len(core) - len(cols) + 1, len(cols))]))
        ok = np.isfinite(d[trios]).all(axis=1)
        c, mid = d[trios[ok, 0]], (d[trios[ok, 1]] + d[trios[ok, 2]]) / 2
        par = np.concatenate([d[trios[ok, 1]], d[trios[ok, 2]]])
        var = float(par.var(ddof=1))
        if len(c) < min_trios or not var > noise or not mid.var() > 0:
            continue
        slope = lambda i, c=c, mid=mid: float(np.cov(mid[i], c[i])[0, 1] / mid[i].var(ddof=1))      # this interval's children and parents
        boot = [slope(rng.integers(0, len(c), len(c))) for _ in range(n_boot)]
        rel = 1 - noise / var
        full = slope(np.arange(len(c)))
        gap = c - mid                                                      # under Mendel a child is, on average, the mean of its parents
        row = dict(name=iv.get("name") or f"{a // 1000}-{b // 1000} kb", start=int(a), end=int(b), windows=int(len(cols)), n_trios=int(len(c)),
                   sd_parents=round(var ** 0.5, 4), noise_sd=round(noise ** 0.5, 4), reliability=round(rel, 4), slope=round(full, 4),
                   slope_lo=round(float(np.percentile(boot, 2.5)), 4), slope_hi=round(float(np.percentile(boot, 97.5)), 4), of_expected=round(full / rel, 4),
                   child_minus_midparent=round(float(gap.mean()), 4), child_minus_midparent_se=round(float(gap.std(ddof=1) / len(gap) ** 0.5), 4))
        row["by_value"] = count([None if not np.isfinite(x) or abs(x - round(x)) >= 0.3 else int(round(x)) for x in d])
        # by the calls: the state in the middle of the interval less the state beside it
        out_of = next(((x, y) for x, y in excluded if x <= a and b <= y), (a, b))
        sides = [q for q in (out_of[0] - flank, out_of[1] + flank) if 0 <= q < UNIT]
        dev = []
        for s in P["samples"]:
            k = calls.get(s)
            if k is None or _status(k) != "settled":
                dev.append(None)
                continue
            beside = {k.state_at(q) for q in sides}
            here = k.state_at((a + b) // 2)
            dev.append(None if len(beside) != 1 or None in beside or here is None else int(here - beside.pop()))
        row["by_calls"] = count(dev)
        out.append(row)
    return out or None


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def call_states(call, size=SUB) -> np.ndarray:
    """A call's state per sub-block: that of the segment holding most of it."""
    n = UNIT // size
    st = np.full(n, np.nan)
    for k in range(n):
        a0, b0 = k * size, (k + 1) * size
        st[k] = max(call.segments, key=lambda g: max(0, min(g.end, b0) - max(g.start, a0))).state
    return st


def assemblies(dir_, P: dict, seg: dict, by_sample: dict, log=lambda m: None) -> dict | None:
    """The assemblies against the profiles and the calls, for every cohort member whose two haplotypes were screened."""
    d = Path(dir_)
    if not (d / "haplotypes.tsv").exists():
        return None
    H = read_tsv(d / "haplotypes.tsv")
    C = read_tsv(d / "copies.tsv") if (d / "copies.tsv").exists() else []
    idx = {s: i for i, s in enumerate(P["samples"])}
    haps: dict[str, list[dict]] = {}
    for h in H:
        haps.setdefault(h["sample"], []).append(h)
    copies: dict[str, list[dict]] = {}
    for c in C:
        copies.setdefault(c["sample"], []).append(c)
    core_blocks = seg["core_blocks"]
    r4 = lambda x: None if x is None or not np.isfinite(x) else round(float(x), 4)
    per_sample, block_rows, oddities = [], [], []
    for s, hs in sorted(haps.items()):
        if len(hs) != 2 or s not in idx:
            continue
        i = idx[s]
        sub = np.array([[float(h[f"b{k}"]) if h[f"b{k}"] not in ("", None) else np.nan for k in range(NSUB)] for h in hs])   # 2 x 80
        asm_sub = sub.sum(axis=0)                                                            # copies per 5-kb sub-block
        per = BLOCK // SUB
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            asm_block = np.array([np.nanmedian(asm_sub[b * per:(b + 1) * per]) if np.isfinite(asm_sub[b * per:(b + 1) * per]).any() else np.nan for b in range(NBLOCK)])
        reads_block = P["block"][i]
        cps = copies.get(s, [])
        cls = [c["class"] for c in cps]
        n_trunc, n_partial = cls.count("truncated at contig end"), cls.count("partial")
        resolved = bool(cps) and n_trunc == 0 and n_partial <= 1 and len({c["haplotype"] for c in cps}) == 2
        by_hap = {h["haplotype"]: int(h["unit_median"]) for h in hs}
        # a haplotype's k-mer median over the whole unit against the mode of its block medians over the core: with many
        # copies, nucleotide differences among them pull the unit median below the block medians
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            hap_blocks = {h["haplotype"]: np.array([np.nanmedian(sub[k][b * per:(b + 1) * per]) if np.isfinite(sub[k][b * per:(b + 1) * per]).any() else np.nan for b in range(NBLOCK)])
                          for k, h in enumerate(hs)}

        def block_mode(v):
            v = np.round(v[core_blocks] if core_blocks else v)
            v = v[np.isfinite(v)]
            return int(np.bincount(v.astype(int)).argmax()) if len(v) else None
        unit_vs_blocks = "; ".join(f"{hp}: {by_hap[hp]} over the unit, {block_mode(hap_blocks[hp])} over most core blocks" for hp in sorted(by_hap)
                                   if block_mode(hap_blocks[hp]) is not None and by_hap[hp] < block_mode(hap_blocks[hp]))
        asm_core = float(np.nanmedian(asm_block[core_blocks])) if core_blocks else float(np.nanmedian(asm_block))
        reads_core = float(np.nanmedian(reads_block[core_blocks])) if core_blocks else float(np.nanmedian(reads_block))
        diff = reads_block - asm_block
        okb = np.isfinite(diff)
        r = by_sample.get(s, {})
        call = P["calls"].get(s)
        rec = dict(sample=s, group=r.get("dj_group", ""), haplotypes=" + ".join(f"{hp}={by_hap[hp]}" for hp in sorted(by_hap)), assembly_total=sum(by_hap.values()),
                   assembly_core=r4(asm_core), assembly_mean=r4(np.nanmean(asm_block)), reads_cn=r.get("DJ.cn"), reads_step=r.get("DJ.step"),
                   reads_core=r4(reads_core), reads_mean=r4(np.nanmean(reads_block)), diff_core=r4(reads_core - asm_core), diff_mean=r4(np.nanmean(diff[okb])) if okb.any() else None,
                   block_concordance=r4(np.mean(np.round(reads_block[okb]) == np.round(asm_block[okb]))) if okb.any() else None,
                   n_copies=len(cps), n_partial=n_partial, n_truncated=n_trunc, n_distal_start=cls.count("distal-start"), resolved=resolved,
                   phasing_split=max(by_hap.values()) - min(by_hap.values()) if by_hap else 0,
                   unit_vs_blocks=unit_vs_blocks, flagged_chromosomes=r.get("flagged_chromosomes") or "",
                   partial_extents="; ".join(f"{c['haplotype']} {int(c['dj_start']) // 1000}-{int(c['dj_end']) // 1000} kb" for c in cps if c["class"] == "partial"),
                   truncated_extents="; ".join(f"{c['haplotype']} {int(c['dj_start']) // 1000}-{int(c['dj_end']) // 1000} kb ({float(c['contig_len']) / 1e6:.2f} Mb contig)" for c in cps if c["class"] == "truncated at contig end"),
                   rhie=RHIE.get(s, ""))
        if call is not None:
            st = call_states(call)
            oks = np.isfinite(asm_sub)
            asm_partial = [(int(c["dj_start"]) // 1000, int(c["dj_end"]) // 1000) for c in cps if c["class"] == "partial"]
            from ngsdose import segments as seglib
            called_partial = [(a // 1000, b // 1000) for _, a, b in seglib.partial_copies(call, seg.get("excluded") or [])]
            # a copy that lacks an end is, in an assembly, a partial copy holding the rest
            called_partial += [(b // 1000, UNIT // 1000) if a == 0 else (0, a // 1000) for _, a, b in seglib.partial_copies(call, seg.get("excluded") or [], kind="partial loss")]
            # an end is matched within 12 kb; inside the distal 30 kb, where copies begin at 3 or at 22 kb and the calls' intervals are
            # polymorphic, an end is read as the unit's start. `off` is how far the matched ends outside that stretch lie from the assembly's
            near = lambda x, y, tol=12: abs(x - y) <= tol or (x <= 30 and y <= 30)
            pairs = [((a, b), min(((x, y) for x, y in called_partial if near(a, x) and near(b, y)), key=lambda q: abs(q[0] - a) + abs(q[1] - b), default=None)) for a, b in asm_partial]
            found = sum(1 for _, q in pairs if q is not None)
            off = [abs(u - v) for (a, b), q in pairs if q is not None for u, v in ((a, q[0]), (b, q[1])) if not (u <= 30 and v <= 30) and not (u >= UNIT // 1000 - 1 and v >= UNIT // 1000 - 1)]
            bulk = seglib._bulk(call.segments, seg.get("excluded") or [])
            rec.update(copies=call.copies, bulk=bulk, bulk_equal=bool(bulk == int(round(asm_core))), partial=r.get("DJ.partial"), variants=r.get("DJ.variants"), scale_f=round(call.scale, 3),
                       tilt=round(call.tilt, 3),
                       scale_uncertain=bool(call.uncertain), call_gap=None if not np.isfinite(call.gap) else round(call.gap, 1),
                       status=_status(call), off=r.get("DJ.off"), off_z=r.get("DJ.off_z"), fractional=r.get("DJ.fractional"),
                       call_concordance=r4(np.mean(st[oks] == asm_sub[oks])) if oks.any() else None, call_sub_blocks=int(oks.sum()),
                       call_equal=int(np.sum(st[oks] == asm_sub[oks])), assembly_partial_found=f"{found} of {len(asm_partial)}" if asm_partial else "",
                       assembly_partial_off_kb=max(off) if off else None)
        per_sample.append(rec)
        calls_block = np.round(_blocks(np.arange(NSUB) * SUB, call_states(call), BLOCK, NBLOCK, 1)) if call is not None else np.full(NBLOCK, np.nan)
        for b in range(NBLOCK):
            block_rows.append(dict(sample=s, block_kb=b * BLOCK // 1000, assembly=None if not np.isfinite(asm_block[b]) else float(asm_block[b]),
                                   reads=None if not np.isfinite(reads_block[b]) else round(float(reads_block[b]), 3),
                                   call=None if not np.isfinite(calls_block[b]) else int(calls_block[b]), core=b in core_blocks))
        # what the reads say about the assembly's breaks: blocks where a truncated or duplicated fragment moves the assembly's count away from the reads
        over = [b for b in range(NBLOCK) if okb[b] and asm_block[b] - reads_block[b] >= 0.75]
        under = [b for b in range(NBLOCK) if okb[b] and reads_block[b] - asm_block[b] >= 0.75]
        if n_trunc or over or under or rec["phasing_split"] >= 3 or rec["unit_vs_blocks"]:
            oddities.append(dict(sample=s, truncated=rec["truncated_extents"], over_blocks=[b * BLOCK // 1000 for b in over], under_blocks=[b * BLOCK // 1000 for b in under],
                                 phasing_split=rec["phasing_split"], haplotypes=rec["haplotypes"], unit_vs_blocks=rec["unit_vs_blocks"], resolved=resolved,
                                 over_max=float(np.nanmax(asm_block[okb] - reads_block[okb])) if okb.any() else None, under_max=float(np.nanmax(reads_block[okb] - asm_block[okb])) if okb.any() else None))
    if not per_sample:
        return None
    T = per_sample
    A = np.array([[t["assembly_core"], t["reads_core"], t["assembly_mean"], t["reads_mean"], t["reads_cn"] if t["reads_cn"] is not None else np.nan] for t in T], float)
    res = np.array([t["resolved"] for t in T])

    def agree(m):
        if m.sum() < 2:
            return dict(n=int(m.sum()))
        a = A[m]
        d_core, d_mean = a[:, 1] - a[:, 0], a[:, 3] - a[:, 2]
        out = dict(n=int(m.sum()), core_diff_mean=float(d_core.mean()), core_diff_sd=float(d_core.std(ddof=1)), mean_diff_mean=float(d_mean.mean()), mean_diff_sd=float(d_mean.std(ddof=1)),
                   core_within_0_5=int(np.sum(np.abs(d_core) < 0.5)))
        if m.sum() >= 3 and np.isfinite(a[:, 4]).sum() >= 3 and a[:, 2].std() > 0:
            out["r_mean_vs_cn"] = float(np.corrcoef(a[:, 2], a[:, 4])[0, 1])
        return out
    blocks_all = [b for b in block_rows if b["assembly"] is not None and b["reads"] is not None]
    dev = [b for b in blocks_all if abs(b["assembly"] - EXPECTED) >= 0.5]
    same = sum(1 for b in dev if np.sign(b["reads"] - EXPECTED) == np.sign(b["assembly"] - EXPECTED) and abs(b["reads"] - EXPECTED) >= 0.5)
    ten = [b for b in blocks_all if abs(b["assembly"] - EXPECTED) < 0.5]
    false_alarm = sum(1 for b in ten if abs(b["reads"] - EXPECTED) >= 0.5)
    ten_genomes = [t for t in T if t["resolved"] and abs(t["assembly_core"] - EXPECTED) < 0.5]
    factor = float((P["scale"] or {}).get("factor") or 1.0)
    known10 = dict(n=len(ten_genomes), reads_cn=[t["reads_cn"] for t in ten_genomes if t["reads_cn"] is not None],
                   reads_unpinned=[round(t["reads_cn"] / factor, 3) for t in ten_genomes if t["reads_cn"] is not None])
    partial_316 = sum(1 for t in T for ext in t["partial_extents"].split("; ") if ext.endswith("3-316 kb"))
    called = [t for t in T if t.get("copies") is not None]

    def calls_agree(ts):
        n, e = sum(t["call_sub_blocks"] for t in ts), sum(t["call_equal"] for t in ts)
        same = [t for t in ts if t["bulk_equal"]]                                        # the reads and the assembly hold the same state over most of the core
        ns, es = sum(t["call_sub_blocks"] for t in same), sum(t["call_equal"] for t in same)
        return dict(n_genomes=len(ts), sub_blocks=n, equal=e, concordance=(e / n if n else None), bulk_equal=len(same),
                    bulk_differs=[t["sample"] for t in ts if not t["bulk_equal"]], concordance_same_level=(es / ns if ns else None), sub_blocks_same_level=ns)
    certain = [t for t in called if t["resolved"] and t["status"] == "settled"]
    found = [t["assembly_partial_found"] for t in called if t["resolved"] and t["assembly_partial_found"]]
    stats = dict(n=len(T), all=agree(np.ones(len(T), bool)), resolved=agree(res), n_resolved=int(res.sum()),
                 blocks=dict(n=len(blocks_all), deviating=len(dev), deviating_confirmed=same, ten=len(ten), ten_false_alarm=false_alarm,
                             diff_sd=float(np.std([b["reads"] - b["assembly"] for b in blocks_all], ddof=1)) if len(blocks_all) > 2 else None,
                             diff_mad_sd=float(1.4826 * np.median(np.abs(np.array([b["reads"] - b["assembly"] for b in blocks_all]) - np.median([b["reads"] - b["assembly"] for b in blocks_all])))) if blocks_all else None,
                             concordance=float(np.mean([np.round(b["reads"]) == np.round(b["assembly"]) for b in blocks_all])) if blocks_all else None),
                 known10=known10, partial_316=partial_316, n_partial=sum(t["n_partial"] for t in T), n_truncated=sum(t["n_truncated"] for t in T),
                 n_distal_start=sum(t["n_distal_start"] for t in T), n_copies=sum(t["n_copies"] for t in T),
                 copy_classes={k: sum(1 for c in C if c["class"] == k) for k in sorted({c["class"] for c in C})}, n_haplotypes=len(H))
    if called:
        stats["calls"] = dict(resolved=calls_agree([t for t in called if t["resolved"]]), resolved_certain=calls_agree(certain),
                              fragmented=calls_agree([t for t in called if not t["resolved"]]),
                              uncertain=[t["sample"] for t in called if t["scale_uncertain"]],
                              fractional=[t["sample"] for t in called if t["status"] == "fractional"],
                              partial_found=sum(int(x.split(" of ")[0]) for x in found), partial_in_assemblies=sum(int(x.split(" of ")[1]) for x in found),
                              partial_off_kb=max((t["assembly_partial_off_kb"] for t in called if t["resolved"] and t.get("assembly_partial_off_kb") is not None), default=None),
                              partial_found_all=sum(int(t["assembly_partial_found"].split(" of ")[0]) for t in called if t["assembly_partial_found"]),
                              partial_in_all=sum(int(t["assembly_partial_found"].split(" of ")[1]) for t in called if t["assembly_partial_found"]),
                              partial_off_kb_all=max((t["assembly_partial_off_kb"] for t in called if t.get("assembly_partial_off_kb") is not None), default=None))
    log(f"[report] DJ vs HPRC assemblies: {len(T)} genomes ({int(res.sum())} with a resolved assembly), {len(blocks_all)} blocks; reads-assembly core diff "
        f"{stats['resolved'].get('core_diff_mean', float('nan')):+.3f} ± {stats['resolved'].get('core_diff_sd', float('nan')):.3f} (resolved)"
        + (f"; calls equal to the assembly in {100 * stats['calls']['resolved_certain']['concordance']:.1f}% of 5-kb sub-blocks ({stats['calls']['resolved_certain']['n_genomes']} resolved genomes with a settled call)"
           if called and stats["calls"]["resolved_certain"]["concordance"] is not None else ""))
    return dict(samples=T, blocks=block_rows, copies=C, oddities=oddities, stats=stats)


def block_table(P: dict, by: dict) -> list[dict]:
    rows = []
    for i, s in enumerate(P["samples"]):
        r = by.get(s, {})
        row = dict(sample=s, cn=r.get("DJ.cn"), cn_unit=r.get("DJ.cn_unit"), copies=r.get("DJ.copies"), partial=r.get("DJ.partial"), variants=r.get("DJ.variants"),
                   scale_f=r.get("DJ.scale_f"), tilt=r.get("DJ.tilt"), call=r.get("DJ.call"), call_gap=r.get("DJ.call_gap"),
                   off=r.get("DJ.off"), off_z=r.get("DJ.off_z"), fractional=r.get("DJ.fractional"))
        for b in range(NBLOCK):
            v = P["block"][i][b]
            row[f"b{b * BLOCK // 1000}"] = None if not np.isfinite(v) else round(float(v), 3)
        rows.append(row)
    return rows


def call_table(P: dict) -> list[dict]:
    """One row per genome and segment (kind `segment`: the whole number called, the mean on the genome's scale, and
    `raw`, the mean as the reads give it), and one per stretch that reads a fraction of a copy off the whole numbers
    (kind `fraction`: the whole number called there, the mean, and z against the cohort)."""
    cols = ("sample", "kind", "start", "end", "windows", "state", "mean", "se", "raw", "off_integer", "z", "call")
    out = []
    for s in P["samples"]:
        c = P["calls"].get(s)
        if c is None:
            continue
        status = getattr(c, "status", "uncertain" if c.uncertain else "settled")
        out += [dict(sample=s, kind="segment", **g.as_dict(), off_integer=g.off_integer, z=None, call=status) for g in c.segments]
        out += [dict(sample=s, kind="fraction", start=f.start, end=f.end, windows=f.windows, state=f.state, mean=round(f.state + f.offset, 3), se=None, raw=None,
                     off_integer=True, z=f.z, call=status) for f in getattr(c, "fractions", [])]
    return [{k: r.get(k) for k in cols} for r in out]


def parse_fractional(text, z=None) -> list[tuple[float, int, int, float | None]]:
    """'+0.53:106-400kb;-0.45:361-400kb' (and the z beside it, '+4.2;-4.1') -> [(0.53, 106, 400, 4.2), (-0.45, 361, 400, -4.1)]."""
    zs = [x for x in str(z or "").split(";") if x not in ("", "none")]
    out = []
    for k, item in enumerate(x for x in str(text or "").split(";") if x not in ("", "none")):
        m = re.fullmatch(r"([+-]\d+(?:\.\d+)?):(\d+)-(\d+)kb", item.strip())
        if m:
            try:
                zz = float(zs[k]) if k < len(zs) else None
            except ValueError:
                zz = None
            out.append((float(m.group(1)), int(m.group(2)), int(m.group(3)), zz))
    return out


def _num(r, k) -> float:
    try:
        return float(r.get(k))
    except (TypeError, ValueError):
        return float("nan")


def off_whole(rows: list[dict], ped: dict | None, eff: dict | None, compared=None, min_trios: int = 30, n_boot: int = 2000, seed: int = 1) -> dict | None:
    """The genomes that sit off whole numbers: a level some way from the whole number called (`DJ.off`, `DJ.off_z`), a
    stretch of the unit a fraction of a copy off (`DJ.fractional`), or a level between two whole numbers (uncertain).
    A change in part of the cells would leave any of them; so would a library unlike the cohort's. What can be said
    from the cohort is said here: how many there are against what the scatter of the scales gives by chance, which
    way they lie, whom they are found in (the generation, the sex, the release batch), whether they come with the
    other signs of a culture that has changed, how far the offset is shared within families and handling groups
    (which is the part of it that the libraries make), and for each genome what its relatives read."""
    called = [r for r in rows if r.get("DJ.call") not in (None, "", "NA") and np.isfinite(_num(r, "DJ.off_z"))]
    if not called:
        return None
    fr = dict((eff or {}).get("DJ", {}).get("fractions") or {})
    lz = float(fr.get("level_z", 3.0))
    by = {r["sample"]: r for r in called}
    off = {s: _num(r, "DJ.off") for s, r in by.items()}
    z = {s: _num(r, "DJ.off_z") for s, r in by.items()}
    status = Counter(r["DJ.call"] for r in called)
    clear = [r for r in called if r["DJ.call"] != "uncertain"]
    level = [r for r in clear if r["DJ.call"] == "fractional" and abs(z[r["sample"]]) >= lz - 0.005]          # by the status the library gave (the z is written to two places)
    steps = [r for r in clear if parse_fractional(r.get("DJ.fractional"))]
    normal = lambda x: 1 - math.erf(x / math.sqrt(2))                                  # two-sided tail of a normal scatter
    x = np.array([off[r["sample"]] for r in clear])
    sd = float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else float("nan")
    out = dict(n=len(called), status={k: int(v) for k, v in status.items()}, rules=fr, spread_copies=sd,
               level=dict(n=len(level), below=sum(1 for r in level if off[r["sample"]] < 0), above=sum(1 for r in level if off[r["sample"]] > 0),
                          by_chance=float(len(clear) * normal(lz)), of=len(clear),
                          beyond=[dict(copies=c, n=int(np.sum(np.abs(x - np.median(x)) > c)), below=int(np.sum(x - np.median(x) < -c)), above=int(np.sum(x - np.median(x) > c)),
                                       by_chance=float(len(x) * normal(c / sd))) for c in (0.3, 0.4, 0.5)] if sd > 0 else []),
               steps=dict(n=len(steps), events=sum(len(parse_fractional(r["DJ.fractional"])) for r in steps), looked_for=bool(fr.get("steps"))))
    bp = Counter()
    for r in steps:
        for _, a, b, _ in parse_fractional(r["DJ.fractional"]):
            for q in (a, b):
                if 0 < q < UNIT // 1000:
                    bp[int(round(q / 20.0) * 20)] += 1
    out["steps"]["breakpoints"] = [dict(kb=k, n=v) for k, v in sorted(bp.items())]
    out["uncertain_off"] = dict(n=status.get("uncertain", 0), beyond=sum(1 for r in called if r["DJ.call"] == "uncertain" and abs(z[r["sample"]]) >= lz - 0.005))
    odd = lambda r: r["DJ.call"] in ("fractional", "uncertain")
    # whom they are found in
    ped = ped or {}
    trio = {c: q for c, q in ped.items() if c in by and q.get("father") in by and q.get("mother") in by}
    fathers, mothers = {q["father"] for q in trio.values()}, {q["mother"] for q in trio.values()}
    role = lambda s: "child" if s in trio and s not in fathers | mothers else "father" if s in fathers else "mother" if s in mothers else "other"

    def share(sel, name):
        g = [r for r in called if sel(r)]
        return dict(group=name, n=len(g), fractional=sum(1 for r in g if r["DJ.call"] == "fractional"), uncertain=sum(1 for r in g if r["DJ.call"] == "uncertain"))
    out["by_role"] = [share(lambda r, k=k: role(r["sample"]) == k, k) for k in ("father", "mother", "child", "other")]
    out["by_sex"] = [share(lambda r, k=k: r.get("sex_inferred") == k, {"M": "men", "F": "women"}[k]) for k in ("M", "F")]
    batches = sorted({str(r.get("ngspca.batch")) for r in called if r.get("ngspca.batch") not in (None, "", "unknown")})
    out["by_batch"] = [share(lambda r, k=k: str(r.get("ngspca.batch")) == k, k) for k in batches]
    # the other signs of a culture that has changed: an X or a Y lost in part of the cells, a chromosome off its dosage
    def sign(r):
        x_, y_, si = _num(r, "truth.chrX"), _num(r, "truth.chrY"), r.get("sex_inferred")
        return bool((si == "F" and x_ < 1.85) or (si == "M" and y_ < 0.85) or r.get("flagged_chromosomes") not in (None, "", "NA"))
    a_, b_ = sum(1 for r in called if odd(r) and sign(r)), sum(1 for r in called if odd(r) and not sign(r))
    c_, d_ = sum(1 for r in called if not odd(r) and sign(r)), sum(1 for r in called if not odd(r) and not sign(r))
    from .report import fisher_exact
    out["other_signs"] = dict(off=a_ + b_, off_with=a_, settled=c_ + d_, settled_with=c_, p=fisher_exact(a_, b_, c_, d_) if a_ + b_ and c_ + d_ else None)
    # how far the offset is shared: within families (spouses share no genes), and between the groups the samples were handled in
    full = [(c, q["father"], q["mother"]) for c, q in trio.items() if all(by[k]["DJ.call"] != "uncertain" for k in (c, q["father"], q["mother"]))]
    if len(full) >= min_trios:
        ch = np.array([off[c] for c, _, _ in full])
        fa = np.array([off[f] for _, f, _ in full])
        mo = np.array([off[m] for _, _, m in full])
        mid = (fa + mo) / 2
        rng = np.random.default_rng(seed)
        slope = lambda i: float(np.cov(mid[i], ch[i])[0, 1] / mid[i].var(ddof=1))
        boot = [slope(rng.integers(0, len(ch), len(ch))) for _ in range(n_boot)]
        rs = float(np.corrcoef(fa, mo)[0, 1])
        bs = [float(np.corrcoef(fa[i], mo[i])[0, 1]) for i in (rng.integers(0, len(ch), len(ch)) for _ in range(n_boot))]
        out["shared"] = dict(n_trios=len(full), spouses_r=rs, spouses_lo=float(np.percentile(bs, 2.5)), spouses_hi=float(np.percentile(bs, 97.5)),
                             slope=slope(np.arange(len(ch))), slope_lo=float(np.percentile(boot, 2.5)), slope_hi=float(np.percentile(boot, 97.5)),
                             child_father_r=float(np.corrcoef(ch, fa)[0, 1]), child_mother_r=float(np.corrcoef(ch, mo)[0, 1]))
    pops = {}
    for r in clear:
        if r.get("pop"):
            pops.setdefault(r["pop"], []).append(off[r["sample"]])
    pops = {k: v for k, v in pops.items() if len(v) >= 30}
    if len(pops) >= 5:
        allv = np.concatenate([np.array(v) for v in pops.values()])
        between = sum(len(v) * (np.mean(v) - allv.mean()) ** 2 for v in pops.values()) / (len(allv) - 1)
        m = sorted(((float(np.mean(v)), k, len(v)) for k, v in pops.items()))
        sup = {}
        for r in clear:
            if r.get("superpop") and r.get("pop") in pops:
                sup.setdefault(r["superpop"], []).append(off[r["sample"]])
        bsup = sum(len(v) * (np.mean(v) - allv.mean()) ** 2 for v in sup.values()) / (len(allv) - 1) if sup else float("nan")
        out["by_population"] = dict(n=len(pops), share=float(between / allv.var(ddof=1)), share_superpop=float(bsup / allv.var(ddof=1)),
                                    lowest=[dict(pop=k, mean=v, n=n, superpop=next((r.get("superpop") for r in clear if r.get("pop") == k), "")) for v, k, n in m[:3]],
                                    highest=[dict(pop=k, mean=v, n=n, superpop=next((r.get("superpop") for r in clear if r.get("pop") == k), "")) for v, k, n in m[-3:]])
    # each genome off whole numbers, with what its relatives read
    kids = {}
    for c, q in trio.items():
        for k in (q["father"], q["mother"]):
            kids.setdefault(k, []).append(c)
    compared = compared or {}
    table = []
    for r in called:
        if not odd(r):
            continue
        s = r["sample"]
        q = trio.get(s, {})
        rel = [(w, q[w]) for w in ("father", "mother") if q.get(w) in by]
        rel += [("child", c) for c in kids.get(s, [])]
        rel += [("spouse", (trio[c]["mother"] if trio[c]["father"] == s else trio[c]["father"])) for c in kids.get(s, [])[:1]]
        ev = parse_fractional(r.get("DJ.fractional"), r.get("DJ.fractional_z"))
        same = [f"{w} {k}" for w, k in rel if abs(z[k]) >= 2 and z[k] * z[s] > 0] if abs(z[s]) >= lz - 0.005 else []
        shared_step = [f"{w} {k}" for w, k in rel for _, a, b, _ in ev for _, a2, b2, _ in parse_fractional(by[k].get("DJ.fractional"))
                       if any(abs(u - v) <= 15 for u in (a, b) if 0 < u < UNIT // 1000 for v in (a2, b2) if 0 < v < UNIT // 1000)]
        table.append(dict(sample=s, pop=r.get("pop"), sex=r.get("sex_inferred"), role=role(s), batch=r.get("ngspca.batch"), call=r["DJ.call"], level=_num(r, "DJ.cn"),
                          copies=r.get("DJ.copies"), off=off[s], off_z=z[s], fractional=r.get("DJ.fractional") if ev else "none", fractional_z=r.get("DJ.fractional_z") if ev else "none",
                          whole_numbers=r.get("DJ.variants"), scale_f=r.get("DJ.scale_f"), lean=r.get("DJ.tilt"), call_gap=r.get("DJ.call_gap"),
                          relatives="; ".join(f"{w} {k} {off[k]:+.2f}" + (f" [{by[k].get('DJ.fractional')}]" if parse_fractional(by[k].get("DJ.fractional")) else "") for w, k in rel),
                          relatives_off_the_same_way="; ".join(sorted(set(same + shared_step))), assembly=compared.get(s, ""),
                          rank=max([abs(z[s])] + [abs(e[3]) for e in ev if e[3] is not None])))
    table.sort(key=lambda t_: -t_["rank"])
    out["table"] = table
    lv = [t_ for t_ in table if t_["call"] == "fractional" and abs(t_["off_z"]) >= lz - 0.005]
    out["level"]["with_relatives"] = sum(1 for t_ in lv if t_["relatives"])
    out["level"]["relatives_same_way"] = sum(1 for t_ in lv if t_["relatives_off_the_same_way"])
    st = [t_ for t_ in table if t_["call"] == "fractional" and t_["fractional"] != "none"]
    # a step of fractional height in a parent: is it in a child, at a whole number or at a fraction?
    passed = looked = 0
    for t_ in st:
        for c in kids.get(t_["sample"], []):
            looked += 1
            mine = [q_ for _, a, b, _ in parse_fractional(t_["fractional"]) for q_ in (a, b) if 0 < q_ < UNIT // 1000]
            theirs = [q_ for _, a, b, _ in parse_fractional(by[c].get("DJ.fractional")) for q_ in (a, b) if 0 < q_ < UNIT // 1000]
            theirs += [int(q_) for it in str(by[c].get("DJ.partial") or "").split(";") for m_ in [re.fullmatch(r"[+-]\d+:(\d+)-(\d+)kb", it.strip())] if m_
                       for q_ in m_.groups() if 0 < int(q_) < UNIT // 1000]
            passed += any(abs(u - v) <= 15 for u in mine for v in theirs)
    out["steps"].update(children_looked_at=looked, children_with_the_breakpoint=passed)
    return out


def group_of(r: dict) -> str:
    """By the calls: uncertain (another whole number explains the profile nearly as well), carrier (described against
    other than ten copies), partial copy (a copy that holds or lacks an end of the unit), ten. Without calls, by the
    level: carrier, between, zero."""
    if r.get("DJ.copies") not in (None, "", "NA"):
        if r.get("DJ.call") == "uncertain":
            return "uncertain"
        if r.get("DJ.call") == "fractional":
            return "fractional"
        if int(r["DJ.copies"]) != int(EXPECTED):
            return "carrier"
        if r.get("DJ.partial") not in (None, "", "none"):
            return "partial copy"
        return "ten"
    step = r.get("DJ.step")
    if step is None or not np.isfinite(step):
        return ""
    if abs(step) >= 0.75 and abs(step - round(step)) <= 0.3:
        return "carrier"
    return "between" if abs(step - round(step)) > 0.3 else "zero"


def run(rows: list[dict], lib: dict | None, eff: dict, dj_steps: dict, assemblies_dir, log=lambda m: None, rules: dict | None = None,
        ped: dict | None = None) -> dict | None:
    """The whole DJ analysis for the report: profiles of every genome in blocks, the segment map, the cohort's levels
    and calls, and (with `assemblies_dir`) the comparison with the HPRC assemblies."""
    P = profiles(lib)
    if P is None:
        return None
    seg = segments(P, rules)
    stats = cohort_stats(P, seg)
    by = {r["sample"]: r for r in rows}
    for r in rows:
        r["dj_group"] = group_of(r)
    out = dict(n=len(P["samples"]), segments=seg, stats=stats, blocks=block_table(P, by), calls=call_table(P), rules=bool(rules),
               chain=dict((rules or {}).get("segments") or {}))
    inh = inheritance(P, ped, rules)
    if inh:
        out["inheritance"] = inh
        log("[report] DJ polymorphic intervals in the trios: "
            + "; ".join(f"{r['name']}: slope of child on midparent {r['slope']:.2f} ({r['slope_lo']:.2f}-{r['slope_hi']:.2f}) against a reliability of {r['reliability']:.2f}"
                        + "".join(f", one copy passed on in {r[k]['passed']} of {r[k]['pairs']} by the {w}" for k, w in (("by_value", "values"), ("by_calls", "calls")) if r.get(k)) for r in inh))
    out["segments"]["core_kb"] = [b * BLOCK // 1000 for b in seg["core_blocks"]]
    if assemblies_dir:
        asm = assemblies(assemblies_dir, P, seg, by, log)
        if asm:
            out["assemblies"] = asm
            for t in asm["samples"]:
                if t["sample"] in by:
                    by[t["sample"]]["DJ.hprc_core"] = round(t["assembly_core"], 3)
                    by[t["sample"]]["DJ.hprc_mean"] = round(t["assembly_mean"], 3)
    sc = stats["scale"] or {}
    log(f"[report] DJ profiles: {len(P['samples'])} genomes; scale {sc.get('rule')}" + (f" (x{sc['factor']:.4f})" if sc.get("factor") else "")
        + f"; within 0.2 of an integer: core level {100 * stats['core']['within_0_2']:.1f}%, whole unit {100 * stats['whole']['within_0_2']:.1f}%"
        + (f"; complete copies {stats['calls']['copies']}; {stats['calls']['with_partial_copy']} genomes with a partial copy" if stats.get("calls") else ""))
    return out
