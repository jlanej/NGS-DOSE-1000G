"""What the between-step genomes' profiles look like: a run of blocks at another integer (a partial variant), the whole
unit shifted, or noise; from docs/data/dj_blocks.tsv. usage: python3 analysis/dj/between_profiles.py [--samples FILE] [--show HG00232 ...]"""
import argparse, csv
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--samples", help="TSV with `sample` and `class` columns (default: every genome more than 0.3 from an integer)")
ap.add_argument("--show", nargs="*", default=[]); a = ap.parse_args()
P = {r["sample"]: r for r in csv.DictReader(open("docs/data/dj_blocks.tsv"), delimiter="\t")}
keys = [k for k in next(iter(P.values())) if k.startswith("b")]
def prof(s):
    return np.array([float(P[s][k]) if P[s][k] else np.nan for k in keys])
if a.samples:
    sel = [r["sample"] for r in csv.DictReader(open(a.samples), delimiter="\t") if r.get("class", "between") == "between" and r["sample"] in P]
else:
    sel = [s for s, r in P.items() if r["core_level_pinned"] and abs(float(r["core_level_pinned"]) - round(float(r["core_level_pinned"]))) > 0.3]
kinds = {}
for s in sel:
    p = prof(s); dev = p - 10; big = np.abs(dev) >= 0.5
    kind = "partial (a run at another integer)" if (np.sum(np.abs(dev) <= 0.3) >= 3 and big.sum() >= 3) else ("flat shift" if np.sum(np.abs(dev) <= 0.3) < 3 else "noisy")
    kinds[kind] = kinds.get(kind, 0) + 1
print(f"{len(sel)} between-step genomes: " + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())))
for s in a.show:
    print(s, " ".join(f"{x:5.1f}" for x in prof(s)))
