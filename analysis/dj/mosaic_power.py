"""What the calls make of a change in part of the cells, on the cohort's own profiles.

A change in a fraction of the cells is put into genomes that have none: forty genomes at ten copies throughout have
their windows raised or lowered by that fraction of a copy (the whole unit, or a stretch of it), are called again, and
are judged with the rest of the cohort as they would have been. The noise, the scales, the leans and whatever else the
cohort's profiles hold are therefore the real ones. For each fraction: what the chain of whole numbers calls, how
often the call is fractional or uncertain, and what height is reported.

usage (from the results repository's root; about ten minutes): python3 analysis/dj/mosaic_power.py [--cache cache/scan] [--genomes 40]
"""
import argparse, os, sys, warnings
import numpy as np
from ngsdose import cohort, resources, segments as seg
from ngsdose.tables import load_result

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="cache/scan"); ap.add_argument("--genomes", type=int, default=40); ap.add_argument("--seed", type=int, default=7)
a = ap.parse_args()
warnings.simplefilter("ignore", RuntimeWarning)
B = resources.Bundle()
rules = B.calibration()
rule = rules["DJ"]; sp = rule["segments"]; fp = rule.get("fractions") or {}
poly = [tuple(iv) for iv in rule["level_exclude"]]
files = sorted(f for f in os.listdir(a.cache) if f.endswith(".estimate.json.gz"))
prof, eff = {}, None
rows, eff, _ = cohort.cohort_table((load_result(f"{a.cache}/{f}") for f in files), B.anchors(), rules=rules, profiles=prof, n_control_pcs=0, log=lambda m: print(m, file=sys.stderr))
P = prof["DJ"]
names, starts, V, calls = list(P["samples"]), np.array(P["start"]), P["cn"].astype(float), P["calls"]
U = int(P["end"][-1]); W = int(P["end"][0] - P["start"][0])
gc = np.array([np.nan if g is None else g for g in eff["DJ"]["gc"]], float)
wsd = np.array([np.nan if g is None else g for g in eff["DJ"]["window_sd"]], float); level = np.array(eff["DJ"]["level"], bool)
ref = np.nanmedian(wsd[level]); rel = np.where(level, np.clip(wsd / ref, 0.5, 1.5), 1.0)
rng = np.random.default_rng(a.seed)
plain = [i for i, s in enumerate(names) if s in calls and calls[s].status == "settled" and not calls[s].events and len(calls[s].segments) == 1 and calls[s].copies == 10
         and abs(calls[s].off_z) < 1.5]
pick = sorted(rng.choice(plain, a.genomes, replace=False))
print(f"{len(names)} genomes; {len(plain)} at ten copies throughout, settled, within 1.5 SDs of their whole number; {len(pick)} of them take the change")


def call_again(v):
    c = seg.segment(starts, v, rel, window=W, unit_length=U, tau=float(sp["tau"]), min_windows=int(sp["min_windows"]), scale_sd=float(sp["scale_sd"]), tilt=bool(sp["tilt"]),
                    tilt_sd=float(sp["tilt_sd"]))
    c = seg.describe(c, poly, int(sp.get("min_core", seg.MIN_CORE)), expected=10)
    c.readings = []
    return c


def run(lo, hi, frac):
    Vs, Cs = V.copy(), [calls.get(s) for s in names]
    inside = (starts >= lo) & (starts < hi)
    for i in pick:
        Vs[i, inside] = V[i, inside] * (10 + frac) / 10
        Cs[i] = call_again(Vs[i])
    seg.find_fractions(Cs, starts, Vs, rel, window=W, unit_length=U, leave_out=poly, gc=gc, scale_sd=float(sp["scale_sd"]), tilt_sd=float(sp["tilt_sd"]), tau=float(sp["tau"]),
                       level_z=float(fp.get("level_z", 3)), event_z=float(fp.get("event_z", 4)), min_height=float(fp.get("min_height", 0.25)), min_windows=int(fp.get("min_windows", 100)))
    return [Cs[i] for i in pick]


def report(name, lo, hi, sign):
    print(f"\n{name}")
    print("in this share of the cells   whole number called there (share of genomes)   settled   fractional   uncertain   the change seen (fractional, uncertain or a whole copy)   height reported (median)")
    whole = lo == 0 and hi >= U
    for frac in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        C = run(lo, hi, sign * frac)
        mid = (lo + hi) // 2 if not whole else 300000
        st = {}
        for c in C:
            k = c.state_at(mid) - c.copies if not whole else c.copies - 10
            st[k] = st.get(k, 0) + 1
        status = [c.status for c in C]
        hts = []
        for c in C:
            if whole:
                hts.append(c.copies - 10 + c.off)
            else:
                near = lambda f: (lo == 0 or abs(f.start - lo) <= 15000) and (hi >= U or abs(f.end - hi) <= 15000)
                other = lambda f: (hi >= U and f.start == 0 and abs(f.end - lo) <= 15000) or (lo == 0 and f.end == U and abs(f.start - hi) <= 15000)      # the rest, described against the stretch
                for f in c.fractions:
                    if near(f):
                        hts.append(f.height(c.copies))
                    elif other(f):
                        hts.append(-f.height(c.copies))
        called = lambda c: any(e.delta == sign and (lo == 0 or abs(e.start - lo) <= 15000) and (hi >= U or abs(e.end - hi) <= 15000) for e in c.events) if not whole else c.copies == 10 + sign
        seen = sum(1 for c in C if c.status != "settled" or called(c))
        print(f"   {frac:.1f}                        {str({int(k): round(v / len(C), 2) for k, v in sorted(st.items())}):34s}   {status.count('settled') / len(C):5.2f}   "
              f"{status.count('fractional') / len(C):5.2f}       {status.count('uncertain') / len(C):5.2f}        {seen / len(C):5.2f}"
              f"                                                  {np.median(hts) if hts else float('nan'):+.2f}")


report("a whole junction lost", 0, U, -1)
report("a whole junction gained", 0, U, +1)
report("a copy of the first 316 kb", 0, 316000, +1)
report("a copy that lacks the first 110 kb (the first 110 kb lost)", 0, 110000, -1)
report("a stretch inside the unit lost, 240-320 kb", 240000, 320000, -1)
