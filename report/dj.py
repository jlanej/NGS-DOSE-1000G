"""The distal junction beyond one number: its profile along the 400-kb unit in every genome, the
segments of the unit that vary between people, a level set on the stable core, and the comparison
with the HPRC release-2 assemblies of cohort members (pipeline/07_hprc_dj.sh, pipeline/hprc_dj.py).

Within a genome the cohort-calibrated window estimates C_w / exp(a_w) (a_w the cohort's window
efficiency) read the same copy number in every window of a junction copy that is present whole;
a copy that lacks part of the unit lowers the windows it lacks, and an extra partial copy raises
the windows it holds. The profile in blocks therefore shows what a single median hides.
"""
from __future__ import annotations

import csv
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


def grab(estimate: dict):
    """The DJ windows of one estimate as a compact array (start, cn) of the usable windows, for the stream."""
    dj = (estimate.get("classes") or {}).get("DJ") or {}
    W = [(int(w["start"]), float(w["cn"])) for w in dj.get("windows") or [] if w.get("usable") and w.get("cn") is not None and np.isfinite(w["cn"])]
    return np.array(W, dtype=np.float64).reshape(-1, 2) if W else None


def _blocks(starts, vals, size, n, min_windows=MIN_WINDOWS):
    idx = (starts // size).astype(int)
    out = np.full(n, np.nan)
    for i in range(n):
        m = idx == i
        if m.sum() >= min_windows:
            out[i] = np.median(vals[m])
    return out


def profiles(samples: list[str], dj_windows: list, eff: dict | None) -> dict | None:
    """Per genome: calibrated window estimates -> 5-kb and 20-kb block medians and the whole-unit level.
    `eff` is the cohort's DJ efficiency table ({start, a}); without it (a single-sample run) the raw
    window estimates are used and the level is the single-sample one."""
    if not any(w is not None for w in dj_windows):
        return None
    a_by = {}
    if eff:
        a_by = {int(s): float(a) for s, a in zip(eff["start"], eff["a"]) if a is not None and np.isfinite(a)}
    subs, blocks, level, keep = [], [], [], []
    for s, W in zip(samples, dj_windows):
        if W is None:
            continue
        st, cn = W[:, 0], W[:, 1]
        if a_by:
            ok = np.array([int(x) in a_by for x in st])
            st, cn = st[ok], cn[ok] / np.exp(np.array([a_by[int(x)] for x in st[ok]]))
        if len(cn) < 20 or not (cn > 0).all():
            continue
        keep.append(s)
        subs.append(_blocks(st, cn, SUB, NSUB))
        blocks.append(_blocks(st, cn, BLOCK, NBLOCK, 3))
        level.append(float(np.exp(np.median(np.log(cn)))))
    if not keep:
        return None
    return dict(samples=keep, sub=np.array(subs), block=np.array(blocks), level=np.array(level))


def segments(P: dict) -> dict:
    """Which parts of the unit vary between people. Cohort SD of the calibrated estimate per 5-kb
    sub-block; the noise floor is the typical sub-block's SD, and a sub-block far above it holds a
    structural polymorphism. The core is the rest: the level is set there."""
    S = P["sub"]
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)          # sub-blocks without a window (no core k-mers) are all-NaN
            sd, med = np.nanstd(S, axis=0), np.nanmedian(S, axis=0)
    valid = np.isfinite(sd)
    floor = float(np.nanmedian(sd[valid])) if valid.any() else float("nan")
    hyper = valid & (sd > HYPERVARIABLE * floor) if np.isfinite(floor) else np.zeros(NSUB, bool)
    # runs of hyper-variable sub-blocks, as intervals
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
    core_sub = valid & ~hyper
    core_blocks = [b for b in range(NBLOCK) if all(core_sub[b * (BLOCK // SUB) + k] for k in range(BLOCK // SUB) if valid[b * (BLOCK // SUB) + k])
                   and any(valid[b * (BLOCK // SUB) + k] for k in range(BLOCK // SUB))]
    core_level = np.array([float(np.exp(np.nanmedian(np.log(P["sub"][i][core_sub])))) if np.isfinite(P["sub"][i][core_sub]).sum() >= 10 else np.nan for i in range(len(P["samples"]))])
    return dict(sub_kb=[i * SUB // 1000 for i in range(NSUB)], sd=[None if not np.isfinite(x) else round(float(x), 3) for x in sd],
                median=[None if not np.isfinite(x) else round(float(x), 3) for x in med], floor=floor, hypervariable=[bool(x) for x in hyper],
                hyper_intervals=runs, core_blocks=core_blocks, core_level=core_level)


def cohort_stats(P: dict, seg: dict, mode_level: float) -> dict:
    """Integer closeness of the whole-unit and core levels, the pinned scale, block-state frequencies
    and the breakpoint tally of the cohort's profiles."""
    lvl, core = P["level"], seg["core_level"]
    pin = EXPECTED / float(np.nanmedian(core)) if np.isfinite(core).any() and np.nanmedian(core) > 0 else float("nan")
    B = P["block"] * pin
    d = np.diff(B, axis=1)
    tally = [int(np.sum(np.abs(d[:, i]) >= 0.6)) for i in range(NBLOCK - 1)]
    states = {}
    for b in range(NBLOCK):
        v = np.round(B[:, b])
        v = v[np.isfinite(v)]
        states[b * BLOCK // 1000] = {int(k): int(np.sum(v == k)) for k in np.unique(v)}

    def closeness(x):
        x = x[np.isfinite(x)]
        r = np.abs(x - np.round(x))
        return dict(n=int(len(x)), within_0_2=float(np.mean(r <= 0.2)), within_0_3=float(np.mean(r <= 0.3)), between=float(np.mean(r > 0.3)),
                    mode_sd=float(np.std(x[np.abs(x - np.round(x)) < 0.35] - np.round(x[np.abs(x - np.round(x)) < 0.35]), ddof=1)) if np.sum(np.abs(x - np.round(x)) < 0.35) > 2 else None)
    return dict(pin=pin, mode_level=mode_level, core_median=float(np.nanmedian(core)), level_median=float(np.nanmedian(lvl)),
                whole=closeness(lvl * pin), core=closeness(core * pin), breakpoints={(i + 1) * BLOCK // 1000: t for i, t in enumerate(tally)}, states=states,
                block_sd=[float(np.nanstd(B[:, b])) for b in range(NBLOCK)])


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def assemblies(dir_, P: dict, seg: dict, pin: float, by_sample: dict, log=lambda m: None) -> dict | None:
    """The assemblies against the profiles, for every cohort member whose two haplotypes were screened."""
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
    per_sample, block_rows, oddities = [], [], []
    for s, hs in sorted(haps.items()):
        if len(hs) != 2 or s not in idx:
            continue
        i = idx[s]
        sub = np.array([[float(h[f"b{k}"]) if h[f"b{k}"] not in ("", None) else np.nan for k in range(NSUB)] for h in hs])   # 2 x 80
        asm_sub = sub.sum(axis=0)                                                            # copies per 5-kb sub-block
        asm_block = np.array([np.nanmedian(asm_sub[b * 4:(b + 1) * 4]) if np.isfinite(asm_sub[b * 4:(b + 1) * 4]).any() else np.nan for b in range(NBLOCK)])
        reads_block = P["block"][i] * pin
        raw_block = P["block"][i]
        cps = copies.get(s, [])
        cls = [c["class"] for c in cps]
        n_trunc, n_partial = cls.count("truncated at contig end"), cls.count("partial")
        resolved = bool(cps) and n_trunc == 0 and n_partial <= 1 and len({c["haplotype"] for c in cps}) == 2
        by_hap = {h["haplotype"]: int(h["unit_median"]) for h in hs}
        # a haplotype's k-mer median over the whole unit against the mode of its block medians over the core: with many
        # copies, nucleotide differences among them pull the unit median below the block medians
        hap_blocks = {h["haplotype"]: np.array([np.nanmedian(sub[k][b * 4:(b + 1) * 4]) if np.isfinite(sub[k][b * 4:(b + 1) * 4]).any() else np.nan for b in range(NBLOCK)]) for k, h in enumerate(hs)}
        def block_mode(v):
            v = np.round(v[core_blocks] if core_blocks else v); v = v[np.isfinite(v)]
            return int(np.bincount(v.astype(int)).argmax()) if len(v) else None
        unit_vs_blocks = "; ".join(f"{hp}: {by_hap[hp]} over the unit, {block_mode(hap_blocks[hp])} over most core blocks" for hp in sorted(by_hap)
                                   if block_mode(hap_blocks[hp]) is not None and by_hap[hp] < block_mode(hap_blocks[hp]))
        asm_core = float(np.nanmedian(asm_block[core_blocks])) if core_blocks else float(np.nanmedian(asm_block))
        reads_core = float(np.nanmedian(reads_block[core_blocks])) if core_blocks else float(np.nanmedian(reads_block))
        diff = reads_block - asm_block
        okb = np.isfinite(diff)
        r = by_sample.get(s, {})
        r4 = lambda x: None if x is None or not np.isfinite(x) else round(float(x), 4)
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
        per_sample.append(rec)
        for b in range(NBLOCK):
            block_rows.append(dict(sample=s, block_kb=b * BLOCK // 1000, assembly=None if not np.isfinite(asm_block[b]) else float(asm_block[b]),
                                   reads_pinned=None if not np.isfinite(reads_block[b]) else round(float(reads_block[b]), 3),
                                   reads=None if not np.isfinite(raw_block[b]) else round(float(raw_block[b]), 3), core=b in core_blocks))
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
    blocks_all = [b for b in block_rows if b["assembly"] is not None and b["reads_pinned"] is not None]
    dev = [b for b in blocks_all if abs(b["assembly"] - EXPECTED) >= 0.5]
    same = sum(1 for b in dev if np.sign(b["reads_pinned"] - EXPECTED) == np.sign(b["assembly"] - EXPECTED) and abs(b["reads_pinned"] - EXPECTED) >= 0.5)
    ten = [b for b in blocks_all if abs(b["assembly"] - EXPECTED) < 0.5]
    false_alarm = sum(1 for b in ten if abs(b["reads_pinned"] - EXPECTED) >= 0.5)
    # the ten-copy genomes by assembly: what the reads read there sets the scale independently of the cohort's mode
    ten_genomes = [t for t in T if t["resolved"] and abs(t["assembly_core"] - EXPECTED) < 0.5]
    known10 = dict(n=len(ten_genomes), reads_cn=[t["reads_cn"] for t in ten_genomes if t["reads_cn"] is not None],
                   reads_core_unpinned=[t["reads_core"] / pin for t in ten_genomes])
    partial_316 = sum(1 for t in T for ext in t["partial_extents"].split("; ") if ext.endswith("3-316 kb"))
    stats = dict(n=len(T), all=agree(np.ones(len(T), bool)), resolved=agree(res), n_resolved=int(res.sum()),
                 blocks=dict(n=len(blocks_all), deviating=len(dev), deviating_confirmed=same, ten=len(ten), ten_false_alarm=false_alarm,
                             diff_sd=float(np.std([b["reads_pinned"] - b["assembly"] for b in blocks_all], ddof=1)) if len(blocks_all) > 2 else None,
                             diff_mad_sd=float(1.4826 * np.median(np.abs(np.array([b["reads_pinned"] - b["assembly"] for b in blocks_all]) - np.median([b["reads_pinned"] - b["assembly"] for b in blocks_all])))) if blocks_all else None,
                             concordance=float(np.mean([np.round(b["reads_pinned"]) == np.round(b["assembly"]) for b in blocks_all])) if blocks_all else None),
                 known10=known10, partial_316=partial_316, n_partial=sum(t["n_partial"] for t in T), n_truncated=sum(t["n_truncated"] for t in T),
                 n_distal_start=sum(t["n_distal_start"] for t in T), n_copies=sum(t["n_copies"] for t in T),
                 copy_classes={k: sum(1 for c in C if c["class"] == k) for k in sorted({c["class"] for c in C})}, n_haplotypes=len(H))
    log(f"[report] DJ vs HPRC assemblies: {len(T)} genomes ({int(res.sum())} with a resolved assembly), {len(blocks_all)} blocks; reads-assembly core diff "
        f"{stats['resolved'].get('core_diff_mean', float('nan')):+.3f} ± {stats['resolved'].get('core_diff_sd', float('nan')):.3f} (resolved)")
    return dict(samples=T, blocks=block_rows, copies=C, oddities=oddities, stats=stats)


def block_table(P: dict, seg: dict, pin: float) -> list[dict]:
    rows = []
    for i, s in enumerate(P["samples"]):
        r = dict(sample=s, level=round(float(P["level"][i]), 4), core_level=None if not np.isfinite(seg["core_level"][i]) else round(float(seg["core_level"][i]), 4),
                 core_level_pinned=None if not np.isfinite(seg["core_level"][i]) else round(float(seg["core_level"][i] * pin), 4))
        for b in range(NBLOCK):
            v = P["block"][i][b] * pin
            r[f"b{b * BLOCK // 1000}"] = None if not np.isfinite(v) else round(float(v), 3)
        rows.append(r)
    return rows


def group_of(step) -> str:
    """carrier: a whole copy or more from the cohort's level; between: more than 0.3 from an integer; zero: near the level."""
    if step is None or not np.isfinite(step):
        return ""
    if abs(step) >= 0.75 and abs(step - round(step)) <= 0.3:
        return "carrier"
    return "between" if abs(step - round(step)) > 0.3 else "zero"


def run(rows: list[dict], grab: list[dict], eff: dict, dj_steps: dict, assemblies_dir, log=lambda m: None) -> dict | None:
    """The whole DJ analysis for the report: profiles of every genome, the segment map, the cohort's
    integer closeness on the whole unit and on the core, and (with `assemblies_dir`) the comparison
    with the HPRC assemblies. Adds DJ.cn_core (the core level, on the cohort's scale) to the rows."""
    samples = [g["sample"] for g in grab]
    P = profiles(samples, [g.get("dj") for g in grab], (eff or {}).get("DJ"))
    if P is None:
        return None
    seg = segments(P)
    mode_level = float(dj_steps["median"]) if dj_steps and dj_steps.get("median") is not None else float(np.nanmedian(P["level"]))
    stats = cohort_stats(P, seg, mode_level)
    by = {r["sample"]: r for r in rows}
    for i, s in enumerate(P["samples"]):
        if s in by:
            v = seg["core_level"][i]
            by[s]["DJ.cn_core"] = None if not np.isfinite(v) else round(float(v), 4)
    for r in rows:
        r["dj_group"] = group_of(r.get("DJ.step"))
    out = dict(n=len(P["samples"]), segments={k: v for k, v in seg.items() if k != "core_level"}, stats=stats, blocks=block_table(P, seg, stats["pin"]))
    out["segments"]["core_kb"] = [b * BLOCK // 1000 for b in seg["core_blocks"]]
    if assemblies_dir:
        asm = assemblies(assemblies_dir, P, seg, stats["pin"], by, log)
        if asm:
            out["assemblies"] = asm
            for t in asm["samples"]:
                if t["sample"] in by:
                    by[t["sample"]]["DJ.hprc_core"] = round(t["assembly_core"], 3)
                    by[t["sample"]]["DJ.hprc_mean"] = round(t["assembly_mean"], 3)
    log(f"[report] DJ profiles: {len(P['samples'])} genomes; hyper-variable segments {', '.join(f'{a // 1000}-{b // 1000} kb' for a, b in seg['hyper_intervals']) or 'none'}; "
        f"core level median {stats['core_median']:.3f} (pin ×{stats['pin']:.4f}); within 0.2 of an integer: whole unit {100 * stats['whole']['within_0_2']:.1f}%, core {100 * stats['core']['within_0_2']:.1f}%")
    return out
