"""The genomes whose level lies between two whole numbers: what their profiles look like (a run of blocks at another whole
number, the whole unit shifted, or noise), and what the calls make of them; from docs/data/dj_blocks.tsv.
usage: python3 analysis/dj/between_profiles.py [--unit] [--show HG00232 ...]"""
import argparse, csv, re
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--unit", action="store_true", help="select by the level over the whole unit (cn_unit) instead of the core's (cn)")
ap.add_argument("--show", nargs="*", default=[]); a = ap.parse_args()
P = {r["sample"]: r for r in csv.DictReader(open("docs/data/dj_blocks.tsv"), delimiter="\t")}
keys = [k for k in next(iter(P.values())) if re.fullmatch(r"b\d+", k)]
prof = lambda s: np.array([float(P[s][k]) if P[s][k] not in ("", "NA", "None") else np.nan for k in keys])
col = "cn_unit" if a.unit else "cn"
sel = [s for s, r in P.items() if r[col] not in ("", "NA") and abs(float(r[col]) - round(float(r[col]))) > 0.3]
kinds, by_call = {}, {}
for s in sel:
    p = prof(s); base = round(np.nanmedian(p)); dev = p - 10; big = np.abs(dev) >= 0.5
    kind = "partial (a run at another whole number)" if (np.sum(np.abs(dev) <= 0.3) >= 3 and big.sum() >= 3) else ("flat shift" if np.sum(np.abs(dev) <= 0.3) < 3 else "noisy")
    kinds[kind] = kinds.get(kind, 0) + 1
    r = P[s]
    ev = [x for x in str(r.get("variants") or "").split(";") if x and x != "none"]
    large = [x for x in ev if (lambda m: m and int(m.group(2)) - int(m.group(1)) >= 40)(re.search(r":(\d+)-(\d+)kb", x))]
    lean = r.get("tilt") not in ("", "NA", None) and abs(float(r["tilt"])) > 0.05
    call = ("an uncertain call" if r.get("call") == "uncertain" else "a copy that holds or lacks an end" if r.get("partial") not in ("", "none", "NA", None)
            else "another event of 40 kb or more" if large else "other copies than ten throughout" if r.get("copies") not in ("", "NA", None, "10")
            else "ten copies and a lean of more than 5%" if lean else "ten copies, the scale a few percent off")
    by_call[call] = by_call.get(call, 0) + 1
print(f"{len(sel)} genomes more than 0.3 from a whole number ({col}): by the look of the profile, " + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())))
print("what the calls make of them: " + ", ".join(f"{k} {v}" for k, v in sorted(by_call.items(), key=lambda kv: -kv[1])))
for s in a.show:
    print(s, P[s]["copies"], P[s]["partial"], P[s]["variants"], " ".join(f"{x:5.1f}" for x in prof(s)))
