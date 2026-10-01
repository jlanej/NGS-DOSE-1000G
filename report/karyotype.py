"""Chromosomes in copies across the cohort: what `ngsdose.karyotype` read in every genome, summed up for the page.

The readings come from the cohort layer (`cohort.cohort_table(..., karyotype=..., karyotypes=...)`): per genome, every
chromosome's level in copies with its standard error and its distance from the nearest whole number, the stretches of
a chromosome that hold another level, and the sex chromosomes on their own scales. This module counts them, sets them
against what else is known of the same genomes (the sex the pedigree gives, NGS-PCA's coverage ratios, the distal
junction's copies, the genome's relatives), and checks them from the inside: the two arms of a chromosome called off,
and the side of two copies on which nothing is expected."""
from __future__ import annotations

import math
import re
from collections import Counter

import numpy as np

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
ACROCENTRIC = ("chr13", "chr14", "chr15", "chr21", "chr22")
PLAIN = ("46,XX", "46,XY")


def _natural(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def _num(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _mad(x) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else float("nan")


def _q(x, qs=(10, 50, 90)) -> list[float]:
    x = np.asarray(x, float)
    return [float(v) for v in np.percentile(x, qs)] if len(x) else []


def _gauss_tail(z: float) -> float:
    """One-sided tail of the normal distribution."""
    return 0.5 * math.erfc(z / math.sqrt(2))


def event_rows(samples, readings) -> list[dict]:
    """Every event of every genome: a chromosome, an arm or a stretch that is not at its expected whole number."""
    out = []
    for s, rd in zip(samples, readings):
        if rd is None:
            continue
        for e in rd.events:
            out.append(dict(sample=s, chrom=e.chrom, span=e.span, start=e.start, end=e.end, regions=e.regions, copies=round(e.copies, 3), se=round(e.se, 4),
                            delta=round(e.delta, 3), z=round(e.z, 1), cells="every cell" if e.whole else "part of the cells", label=e.label(), karyotype=rd.karyotype()))
    return out


def chromosome_table(readings, model, rules: dict) -> list[dict]:
    """Per chromosome: the regions it is read from, how well its level is known in a typical genome, from what share of
    the cells a gain or a loss of one copy is called, and what was found."""
    from ngsdose import karyotype as K
    tab = K.table(model.names, model.arms)
    level_z, min_off = float(rules.get("level_z", K.LEVEL_Z)), float(rules.get("min_off", K.MIN_OFF))
    rs = [r for r in readings if r is not None]
    out = []
    for c in tab.chromosomes():
        on = (tab.chrom == c) & np.isfinite(model.sd)
        groups = [(c, lambda ch: True)] if c in AUTOSOMES else [(f"{c[3:]}, one copy", lambda ch: ch.whole == 1)] + ([(f"{c[3:]}, two copies", lambda ch: ch.whole == 2)] if c == "chrX" else [])
        for label, want in groups:
            chs = [r.chromosomes[c] for r in rs if c in r.chromosomes and r.chromosomes[c].status != "not read" and want(r.chromosomes[c])]
            if not chs:
                continue
            se = float(np.median([ch.se for ch in chs]))
            ev = [e for r in rs for e in r.events if e.chrom == c and c in r.chromosomes and want(r.chromosomes[c])]
            settled = [ch.copies for ch in chs if ch.status == "settled"]
            out.append(dict(chrom=c, label=label[3:] if label.startswith("chr") else label, regions=int(np.median([ch.regions for ch in chs])), model_regions=int(on.sum()),
                            places=int(np.median([sum(g.loci for g in ch.segments) for ch in chs])), n=len(chs), se=se, spread=_mad(settled),
                            seen_from=max(min_off, level_z * se), whole=sum(1 for e in ev if e.span == "whole" and e.whole),
                            part=sum(1 for e in ev if e.span == "whole" and not e.whole), stretches=sum(1 for e in ev if e.span != "whole"),
                            both_arms=sum(1 for ch in chs if len(ch.arms) == 2)))
    return out


def read_one(model, vec, rules: dict | None = None):
    """One genome's single-copy regions (`ngsdose.karyotype.gather`) read against a model, under the bundle's rules."""
    from ngsdose import karyotype as K
    if vec is None:
        return None
    rules = rules or {}
    return K.read(model, vec[0], vec[1], tau=float(rules.get("tau", K.TAU)), level_z=float(rules.get("level_z", K.LEVEL_Z)), min_off=float(rules.get("min_off", K.MIN_OFF)))


def _same(e, f) -> bool:
    """Two events of one genome read from two sets of regions: the same chromosome, the same direction, and a whole
    chromosome in both or spans that overlap by half of the shorter."""
    if e.chrom != f.chrom or (e.delta > 0) != (f.delta > 0):
        return False
    if e.span == "whole" or f.span == "whole":
        return e.span == f.span or min(e.end, f.end) - max(e.start, f.start) > 0.5 * max(e.end - e.start, f.end - f.start, 1)
    return min(e.end, f.end) - max(e.start, f.start) > 0.5 * min(e.end - e.start, f.end - f.start)


def windows_check(samples, readings, base: dict, windows: set) -> dict | None:
    """The genomes counted with the karyotype windows, read twice against the same model: from all their regions, and
    from the regions the earlier counts hold (`base`). What both see, and at what levels; what only the windows show."""
    both, only_w, only_b, flips = [], [], [], 0
    n = 0
    se_ratio = []
    for s, rd in zip(samples, readings):
        b = base.get(s)
        if s not in windows or rd is None or b is None:
            continue
        n += 1
        flips += rd.complement != b.complement
        for c in AUTOSOMES:
            if c in rd.chromosomes and c in b.chromosomes and rd.chromosomes[c].se and b.chromosomes[c].se and np.isfinite(rd.chromosomes[c].se) and np.isfinite(b.chromosomes[c].se):
                se_ratio.append(b.chromosomes[c].se / rd.chromosomes[c].se)
        used = set()
        for e in rd.events:
            m = next((j for j, f in enumerate(b.events) if j not in used and _same(e, f)), None)
            if m is None:
                only_w.append(dict(sample=s, label=e.label(), chrom=e.chrom, span=e.span, copies=e.copies, delta=e.delta, z=e.z, regions=e.regions,
                                   base_copies=b.chromosomes[e.chrom].copies if e.chrom in b.chromosomes else None, base_regions=b.chromosomes[e.chrom].regions if e.chrom in b.chromosomes else 0))
            else:
                used.add(m)
                both.append((e.delta, b.events[m].delta, e.se, b.events[m].se))
        only_b += [dict(sample=s, label=f.label(), chrom=f.chrom, copies=f.copies, delta=f.delta, z=f.z, with_windows=rd.chromosomes[f.chrom].copies if f.chrom in rd.chromosomes else None)
                   for j, f in enumerate(b.events) if j not in used]
    if not n:
        return None
    d = np.array(both, float).reshape(-1, 4)
    return dict(n=n, both=len(d), r=float(np.corrcoef(d[:, 0], d[:, 1])[0, 1]) if len(d) > 2 else None,
                diff_sd=_mad(d[:, 0] - d[:, 1]) if len(d) else None, diff_max=float(np.abs(d[:, 0] - d[:, 1]).max()) if len(d) else None,
                diff_z=_mad((d[:, 0] - d[:, 1]) / np.hypot(d[:, 2], d[:, 3])) if len(d) else None,
                only_windows=sorted(only_w, key=lambda x: -abs(x["z"])), only_base=sorted(only_b, key=lambda x: -abs(x["z"])), complement_differs=int(flips),
                se_ratio=_q(se_ratio, (10, 50, 90)))


def tilt_summary(samples, readings, model, windows: set) -> dict | None:
    """The fit across chromosomes (a common shift and a slope on the chromosomes' GC, each chromosome against the others):
    how the slope spreads over the genomes, how often it stands beyond its own error, and the genomes at its ends with what
    their GC-rich chromosomes would read without it."""
    gc = getattr(model, "gc", None) or {}
    rs = [(s, r) for s, r in zip(samples, readings) if r is not None and getattr(r, "tilt", None) is not None and r.gc_tilt_z is not None]
    if not gc or not rs:
        return None
    out = dict(n=len(rs), gc={c: gc[c] for c in sorted(gc, key=_natural)}, groups=[])
    for label, grp in (("with the windows", [x for x in rs if x[0] in windows]), ("without", [x for x in rs if x[0] not in windows])):
        if len(grp) >= 20:
            t = np.array([r.gc_tilt for _, r in grp])
            z = np.array([r.gc_tilt_z for _, r in grp])
            out["groups"].append(dict(label=label, n=len(grp), sd=_mad(t), q=_q(t, (1, 50, 99)), se=float(np.median(np.abs(t / z))),
                                      beyond3=int((np.abs(z) > 3).sum()), expected3=float(len(z) * 2 * _gauss_tail(3)), z_sd=_mad(z)))
    pick = [x for x in rs if x[0] in windows] or rs
    ends = sorted(pick, key=lambda x: -abs(x[1].gc_tilt_z))[:6]
    rich = [c for c in ("chr19", "chr22", "chr17", "chr16") if c in gc]
    out["ends"] = [dict(sample=s, tilt=r.gc_tilt, z=r.gc_tilt_z, karyotype=r.karyotype(),
                        chromosomes=[dict(chrom=c, copies=r.chromosomes[c].copies, without=float(r.chromosomes[c].copies * math.exp(r.chromosomes[c].tilt)),
                                          se=r.chromosomes[c].se) for c in rich if c in r.chromosomes and np.isfinite(r.chromosomes[c].copies)]) for s, r in ends]
    return out


def sex_chromosomes(rows, samples, readings, ped) -> dict:
    by = {r["sample"]: r for r in rows}
    comp = Counter()
    table, mism = [], []
    x_lost, y_lost, x_gain, y_gain, y_most = [], [], [], [], []
    for s, rd in zip(samples, readings):
        if rd is None or not rd.complement:
            continue
        rep = {"1": "M", "2": "F", 1: "M", 2: "F"}.get((ped.get(s) or {}).get("sex"), (ped.get(s) or {}).get("sex") or "")
        comp[(rd.complement, rep)] += 1
        ev = [e for e in rd.events if e.chrom in ("chrX", "chrY")]
        looks = "F" if "Y" not in rd.complement else "M"
        odd = rd.complement not in ("XX", "XY") or bool(ev) or (rep in ("M", "F") and rep != looks)
        if rep in ("M", "F") and rep != looks:
            # a man read as X with a Y left in part of the cells: a line that lost its Y in most of them, not another person
            if rep == "M" and any(e.chrom == "chrY" and not e.whole and e.delta > 0 for e in ev):
                y_most.append(dict(sample=s, y=rd.y))
            else:
                mism.append(s)
        for e in ev:
            if e.span != "whole" or e.whole:
                continue
            (x_lost if e.chrom == "chrX" and e.delta < 0 else x_gain if e.chrom == "chrX" else y_lost if e.delta < 0 else y_gain).append(abs(e.delta))
        if odd:
            r = by.get(s, {})
            table.append(dict(sample=s, pop=r.get("pop", ""), reported=rep, complement=rd.complement, x=rd.x, y=rd.y, karyotype=rd.karyotype(), status=rd.status,
                              note=rd.note, ngspca_x=_num(r.get("ngspca.chrX")), whole=rd.complement not in ("XX", "XY") and not ev))
    table.sort(key=lambda t: (t["complement"] in ("XX", "XY"), -abs((t["x"] or 0) - round(t["x"] or 0)) - abs((t["y"] or 0) - round(t["y"] or 0))))
    # against NGS-PCA's whole-chromosome coverage ratio, where the cohort table carries it
    xs = [(_num(r.get("chrX.copies")), _num(r.get("ngspca.chrX"))) for r in rows]
    xs = np.array([p for p in xs if all(np.isfinite(p))], float)
    ngs = None
    if len(xs) >= 20:
        two = xs[:, 0] > 1.5
        ngs = dict(n=len(xs), r=float(np.corrcoef(xs[:, 0], xs[:, 1])[0, 1]),
                   r_two=float(np.corrcoef(xs[two, 0], xs[two, 1])[0, 1]) if two.sum() >= 20 else None)
    return dict(complements=[dict(complement=k, reported=rep, n=v) for (k, rep), v in sorted(comp.items(), key=lambda kv: -kv[1])],
                by_complement={k: sum(v for (kk, _), v in comp.items() if kk == k) for k in sorted({k for k, _ in comp})},
                not_plain=table, mismatch=mism, y_lost_in_most=y_most,
                part=dict(x_lost=dict(n=len(x_lost), q=_q(x_lost)), y_lost=dict(n=len(y_lost), q=_q(y_lost)), x_gained=dict(n=len(x_gain), q=_q(x_gain)), y_gained=dict(n=len(y_gain), q=_q(y_gain))),
                ngspca=ngs)


def autosomes(samples, readings) -> dict:
    by_chrom = {c: dict(chrom=c, whole=0, part=0, shares=[]) for c in AUTOSOMES}
    stretches, clones, losses, gains = [], [], 0, 0
    for s, rd in zip(samples, readings):
        if rd is None:
            continue
        ev = [e for e in rd.events if e.chrom in AUTOSOMES]
        for e in ev:
            if e.span == "whole":
                d = by_chrom[e.chrom]
                d["whole" if e.whole else "part"] += 1
                d["shares"].append(e.delta)
                gains += e.delta > 0
                losses += e.delta < 0
            else:
                stretches.append(dict(sample=s, chrom=e.chrom, span=e.span, start=e.start, end=e.end, mb=(e.end - e.start) / 1e6, copies=e.copies, delta=e.delta, z=e.z,
                                      cells="every cell" if e.whole else "part of the cells", regions=e.regions))
        # several chromosomes gained in the same share of the cells: one clone that carries them all
        part = [e for e in ev if e.span == "whole" and not e.whole]
        if len(part) >= 3:
            sh = np.array([abs(e.delta) for e in part])
            if sh.max() - sh.min() <= 0.06:
                clones.append(dict(sample=s, chromosomes=[e.label() for e in part], share=float(np.median(sh)), spread=float(sh.max() - sh.min()), karyotype=rd.karyotype()))
    rows = [dict(chrom=d["chrom"], whole=d["whole"], part=d["part"], lo=min(d["shares"], key=abs) if d["shares"] else None, hi=max(d["shares"], key=abs) if d["shares"] else None,
                 gained=sum(1 for x in d["shares"] if x > 0), lost=sum(1 for x in d["shares"] if x < 0)) for d in by_chrom.values()]
    # stretches that recur: the same part of a chromosome, the same direction, in several genomes
    recur = []
    for c in AUTOSOMES:
        st = sorted((x for x in stretches if x["chrom"] == c), key=lambda x: x["start"])
        used = [False] * len(st)
        for i, a in enumerate(st):
            if used[i]:
                continue
            grp = [a]
            for j in range(i + 1, len(st)):
                b = st[j]
                if not used[j] and (b["delta"] > 0) == (a["delta"] > 0) and min(a["end"], b["end"]) - max(a["start"], b["start"]) > 0.3 * min(a["end"] - a["start"], b["end"] - b["start"]):
                    grp.append(b)
                    used[j] = True
            if len(grp) >= 3:
                recur.append(dict(chrom=c, start=min(g["start"] for g in grp), end=max(g["end"] for g in grp), n=len(grp), gain=grp[0]["delta"] > 0,
                                  samples=[g["sample"] for g in grp], common=(max(g["start"] for g in grp), min(g["end"] for g in grp))))
    stretches.sort(key=lambda x: (-abs(x["z"])))
    return dict(by_chromosome=rows, stretches=stretches, clones=clones, recurrent=sorted(recur, key=lambda r: -r["n"]), gained=int(gains), lost=int(losses),
                genomes=sum(1 for rd in readings if rd is not None and any(e.chrom in AUTOSOMES for e in rd.events)))


def checks(rows, samples, readings, ped) -> dict:
    by = {r["sample"]: r for r in rows}
    idx = {s: i for i, s in enumerate(samples)}
    rs = [r for r in readings if r is not None]
    # 1. the arms of every whole chromosome called off: a whole chromosome shows in both alike
    arm = []
    for s, rd in zip(samples, readings):
        if rd is None:
            continue
        for e in rd.events:
            ch = rd.chromosomes[e.chrom]
            if e.span == "whole" and e.chrom in AUTOSOMES and len(ch.arms) == 2:
                (p, sp, _), (q, sq, _) = ch.arms["p"], ch.arms["q"]
                arm.append(dict(sample=s, chrom=e.chrom, copies=e.copies, p=p, q=q, z=(p - q) / math.hypot(sp, sq), small=abs(e.delta) < 0.1,
                                both=abs(p - 2) > 2 * sp and abs(q - 2) > 2 * sq))
    arms = None
    if arm:
        d = np.array([(a["p"] - 2, a["q"] - 2) for a in arm])
        arms = dict(n=len(arm), same_side=int((np.sign(d[:, 0]) == np.sign(d[:, 1])).sum()), r=float(np.corrcoef(d[:, 0], d[:, 1])[0, 1]) if len(arm) > 2 else None,
                    z_sd=_mad([a["z"] for a in arm]), beyond3=sum(1 for a in arm if abs(a["z"]) > 3),
                    small=sum(1 for a in arm if a["small"]), small_same=sum(1 for a, dd in zip(arm, d) if a["small"] and np.sign(dd[0]) == np.sign(dd[1])),
                    small_both=sum(1 for a in arm if a["small"] and a["both"]))
    # 2. the quiet side: a chromosome below two copies. Losses of an autosome are rare in a culture, so the lower tail of z
    # is the measurement's own scatter, and can be set against the normal distribution
    z = np.array([ch.z for r in rs for c, ch in r.chromosomes.items() if c in AUTOSOMES and ch.z is not None and ch.whole == 2], float)
    tails = [dict(z=t, below=int((z < -t).sum()), above=int((z > t).sum()), expected=float(_gauss_tail(t) * len(z))) for t in (3, 4, 5)] if len(z) else []
    # 3. the estimator's own flags (a chromosome whose pooled depth departs by 4%): what they are here
    fl = {"n": 0, "whole": 0, "stretch": 0, "neither": 0, "stretches": []}
    for s, rd in zip(samples, readings):
        f = by.get(s, {}).get("flagged_chromosomes")
        if rd is None or not f or str(f) in ("NA", "nan", "None"):
            continue
        for c in (f if isinstance(f, (list, tuple)) else str(f).split(",")):
            c = c.strip()
            if not c:
                continue
            fl["n"] += 1
            ev = [e for e in rd.events if e.chrom == c]
            if any(e.span == "whole" for e in ev):
                fl["whole"] += 1
            elif ev:
                fl["stretch"] += 1
                fl["stretches"].append(f"{s} {ev[0].label()}")
            else:
                fl["neither"] += 1
    # 4. the distal junction: one on every acrocentric short arm, so a whole acrocentric chromosome more is a junction more
    dj = []
    for s, rd in zip(samples, readings):
        if rd is None:
            continue
        r = by.get(s, {})
        level = _num(r.get("DJ.cn"))
        whole = [e for e in rd.events if e.chrom in ACROCENTRIC and e.span == "whole"]
        other = [e for e in rd.events if e.chrom in ACROCENTRIC and e.span != "whole" and abs(e.delta) >= 0.3]
        if (whole or other) and np.isfinite(level):
            dj.append(dict(sample=s, events=", ".join(e.label() for e in whole + other), whole_chromosome=bool(whole), extra=sum(e.delta for e in whole), dj=level,
                           dj_call=r.get("DJ.call"), dj_copies=r.get("DJ.copies")))
    # 5. relatives: an event that a child shares with a parent was in the germ line, not made in the culture
    shared, n_rel = [], 0
    for s, rd in zip(samples, readings):
        p = ped.get(s) or {}
        if rd is None or not rd.events:
            continue
        for role in ("father", "mother"):
            o = p.get(role)
            if o in idx and readings[idx[o]] is not None:
                n_rel += 1
                for e in rd.events:
                    for f in readings[idx[o]].events:
                        if f.chrom == e.chrom and (f.delta > 0) == (e.delta > 0) and min(e.end, f.end) - max(e.start, f.start) > 0.5 * max(e.end - e.start, f.end - f.start, 1):
                            shared.append(dict(child=s, parent=o, role=role, child_event=e.label(), parent_event=f.label(), child_copies=e.copies, parent_copies=f.copies))
    return dict(arms=arms, tails=dict(n=int(len(z)), rows=tails), flags=fl, dj=sorted(dj, key=lambda d: -abs(d["extra"])), relatives=dict(pairs=n_rel, shared=shared))


def alleles(path, samples=None) -> dict | None:
    """The calls against the alleles of the same reads (analysis/karyotype/allele_balance.py writes the table): the events
    with enough heterozygous sites, the two shares of the cells and how they agree."""
    import csv
    try:
        with open(path) as fh:
            rows = list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))
    except OSError:
        return None
    keep = set(samples) if samples is not None else None
    # where two copies are expected; a line read as one copy with a second in part of the cells has heterozygous sites only
    # where the second is, and is listed apart (`one_copy`)
    ev = [dict(sample=r["sample"], event=r["event"], sites=int(r["sites"]), by_depth=float(r["by_depth"]), by_alleles=float(r["by_alleles"]))
          for r in rows if r["chrom"] and r["by_alleles"] not in ("", "NA") and not r.get("note") and (keep is None or r["sample"] in keep)]
    one = [dict(sample=r["sample"], event=r["event"], sites=int(r["sites"]), by_depth=float(r["by_depth"]), by_alleles=float(r["by_alleles"]))
           for r in rows if r["chrom"] and r["by_alleles"] not in ("", "NA") and r.get("note") and (keep is None or r["sample"] in keep)]
    quiet = [float(r["d"]) for r in rows if not r["chrom"] and r["d"] not in ("", "NA") and (keep is None or r["sample"] in keep)]
    few = [dict(sample=r["sample"], event=r["event"], sites=int(r["sites"]), note=r["note"]) for r in rows if r["chrom"] and r["by_alleles"] in ("", "NA") and (keep is None or r["sample"] in keep)]
    if len(ev) < 3:
        return None
    d = np.array([(e["by_depth"], e["by_alleles"]) for e in ev])
    return dict(rows=sorted(ev, key=lambda e: -e["by_depth"]), n=len(ev), genomes=len({e["sample"] for e in ev}), r=float(np.corrcoef(d.T)[0, 1]),
                diff_sd=_mad(d[:, 0] - d[:, 1]), diff_median=float(np.median(d[:, 1] - d[:, 0])), quiet_d=float(np.median(quiet)) if quiet else None, too_few=few,
                one_copy=one, near=int((np.abs(d[:, 0] - d[:, 1]) < 0.1).sum()),
                not_borne=[e for e in ev if e["by_depth"] >= 0.1 and e["by_alleles"] < 0.5 * e["by_depth"]])


def run(rows, kout: dict | None, ped: dict, rules: dict | None = None, log=None, alleles_path=None) -> dict | None:
    """The summary for the page, or None when no genome was read. `kout`: what `cohort_table` filled (`karyotypes`)."""
    if not kout or kout.get("model") is None:
        return None
    samples, readings, model, info = kout["samples"], kout["readings"], kout["model"], kout.get("info") or {}
    rs = [r for r in readings if r is not None]
    if not rs:
        return None
    from ngsdose import karyotype as K
    rules = rules or {}
    tab = K.table(model.names, model.arms)
    # the regions a typical genome is read from (a genome counted with fewer than the model knows is read from those it has)
    per = lambda pick: int(np.median([sum(ch.regions for c, ch in r.chromosomes.items() if pick(c) and ch.status != "not read") for r in rs]))
    kinds = {"A": per(lambda c: c in AUTOSOMES), "X": per(lambda c: c == "chrX"), "Y": per(lambda c: c == "chrY")}
    u = model.u[np.isfinite(model.u)]
    out = dict(n=len(rs), status=dict(Counter(r.status for r in rs)), info=info,
               rules=dict(tau=float(rules.get("tau", K.TAU)), level_z=float(rules.get("level_z", K.LEVEL_Z)), min_off=float(rules.get("min_off", K.MIN_OFF)), max_se=K.MAX_SE,
                          min_loci=K.MIN_LOCI),
               model=dict(regions=kinds["A"] + kinds["X"] + kinds["Y"], model_regions=int(np.isfinite(model.sd).sum()), autosomal=kinds["A"], x=kinds["X"], y=kinds["Y"],
                          components=model.k, n=model.n, saved=info.get("model") == "saved",
                          u=float(np.median(u)) if len(u) else None, phi={k: float(v) for k, v in model.phi.items()},
                          floor={k: float(v) for k, v in (getattr(model, "floor", None) or {}).items()}, local=float(getattr(model, "local", 0.0) or 0.0),
                          error=(model.info or {}).get("error"), single_values=sorted(str(c) for c, s0 in zip(tab.chrom, tab.start) if s0 < 0)),
               plain=sum(1 for r in rs if r.karyotype() in PLAIN), noise=_q([r.noise for r in rs], (1, 50, 99)),
               whole_number=sum(1 for r in rs if any(e.whole for e in r.events) or r.complement not in ("XX", "XY", "")),
               chromosomes=chromosome_table(readings, model, rules), sex=sex_chromosomes(rows, samples, readings, ped), autosomes=autosomes(samples, readings),
               checks=checks(rows, samples, readings, ped), karyotypes=[dict(karyotype=k, n=v) for k, v in Counter(r.karyotype() for r in rs).most_common(12)])
    win = kout.get("windows") or set()
    out["tilt"] = tilt_summary(samples, readings, model, win)
    if win:
        isw = [s in win for s in samples]
        w_r, b_r = [r for r, w in zip(readings, isw) if w], [r for r, w in zip(readings, isw) if not w]
        out["layouts"] = dict(with_windows=dict(n=sum(1 for r in w_r if r is not None), regions=int(np.median([r.regions for r in w_r if r is not None])),
                                                chromosomes=chromosome_table(w_r, model, rules), status=dict(Counter(r.status for r in w_r if r is not None))),
                              without=dict(n=sum(1 for r in b_r if r is not None), regions=int(np.median([r.regions for r in b_r if r is not None])) if any(r is not None for r in b_r) else 0,
                                           chromosomes=chromosome_table(b_r, model, rules), status=dict(Counter(r.status for r in b_r if r is not None))) if any(r is not None for r in b_r) else None)
        out["windows"] = windows_check(samples, readings, kout.get("base") or {}, win)
    if alleles_path:
        out["alleles"] = alleles(alleles_path, set(samples))
    if log:
        s = out["status"]
        log(f"[report] chromosomes: {out['n']} genomes read from {out['model']['regions']} regions; {out['plain']} are 46,XX or 46,XY throughout, "
            f"{s.get('fractional', 0)} hold a change in part of the cells, {s.get('uncertain', 0)} are uncertain")
    return out


def flags(row: dict) -> list[str]:
    """What the samples table says of a genome whose chromosomes are not plain."""
    k = row.get("karyotype")
    if not k or str(k) in ("NA", "None", "nan"):
        return []
    out = []
    if k not in PLAIN:
        out.append(f"karyotype {k}" + (" (part of the cells in brackets)" if "[" in str(k) else ""))
    if row.get("karyotype.status") == "uncertain":
        out.append("chromosomes read too coarsely to settle")
    note = row.get("karyotype.note")
    if note and str(note) not in ("NA", "None", "nan"):
        out.append(str(note))
    return out
