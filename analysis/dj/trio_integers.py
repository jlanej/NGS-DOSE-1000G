"""Whole numbers of junction copies in the trios: Mendelian consistency and transmission by the parent's sex, by the level
(over the whole unit, and on the core) and by the calls (docs/data/dj_blocks.tsv: cn_unit, cn, copies).
usage: python3 analysis/dj/trio_integers.py [--tol 0.3]"""
import argparse, csv
ap = argparse.ArgumentParser(); ap.add_argument("--tol", type=float, default=0.3); a = ap.parse_args()
P = {r["sample"]: r for r in csv.DictReader(open("docs/data/dj_blocks.tsv"), delimiter="\t")}
trios = list(csv.DictReader(open("docs/data/trios.tsv"), delimiter="\t"))
def haps(s):
    return {(x, s - x) for x in range(-2, 3) if -1 <= s - x <= 2 and -1 <= x <= 2}
for label, key in (("level over the whole unit", "cn_unit"), ("level on the core", "cn"), ("calls (the copies a settled call is described against)", "copies")):
    lvl = {s: float(r[key]) for s, r in P.items() if r.get(key) not in (None, "", "NA") and (key != "copies" or r.get("call") == "settled")}
    st = {s: int(round(v - 10)) for s, v in lvl.items()}; res = {s: abs(v - round(v)) for s, v in lvl.items()}
    n = unc = ok = bad = 0; tr = {}
    for t in trios:
        c, f, m = t["child"], t["father"], t["mother"]
        if not all(x in lvl for x in (c, f, m)):
            continue
        n += 1
        if max(res[c], res[f], res[m]) > a.tol:
            unc += 1; continue
        possible = {x + y for (x, x2) in haps(st[f]) for (y, y2) in haps(st[m]) for x, y in ((x, y), (x2, y), (x, y2), (x2, y2))}
        ok += st[c] in possible; bad += st[c] not in possible
        for parent, sp, so in (("father", st[f], st[m]), ("mother", st[m], st[f])):
            if sp in (-1, 1) and so == 0:
                k = (parent, sp); tr.setdefault(k, [0, 0])
                if st[c] == sp: tr[k][0] += 1
                elif st[c] == 0: tr[k][1] += 1
    tt, nn = sum(v[0] for v in tr.values()), sum(v[0] + v[1] for v in tr.values())
    print(f"{label}: {n} trios, {unc} with a member more than {a.tol} from a whole number; {ok} consistent, {bad} Mendelian-impossible; single-carrier parents transmit {tt} of {nn}: "
          + "; ".join(f"{p} {s:+d} {v[0]} of {v[0] + v[1]}" for (p, s), v in sorted(tr.items())))
