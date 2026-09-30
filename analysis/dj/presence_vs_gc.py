"""Is the junction's GC-dependent deficit sequence divergence? Exact-k-mer presence of the panel's k-mers in one assembled
haplotype (HG00097 hap1: `ngs-dose panel --report`, committed here) per 250-bp window, against window GC and against the
raw (uncalibrated) window estimate averaged over the genomes at the cohort's mode.
usage: python3 analysis/dj/presence_vs_gc.py [--report analysis/dj/HG00097_hap1.report.tsv.gz] [--cache cache/scan]"""
import argparse, csv, gzip, json, os
import numpy as np
from ngsdose import io, resources
ap = argparse.ArgumentParser(); ap.add_argument("--report", default="analysis/dj/HG00097_hap1.report.tsv.gz"); ap.add_argument("--cache", default="cache/scan")
ap.add_argument("--max-genomes", type=int, default=150); a = ap.parse_args()
core = io.load_panel(resources.Bundle(resources.default_bundle()).panel).classes["DJ"].kmer_pos
bg = {}
with gzip.open(a.report, "rt") as fh:
    for line in fh:
        if line[0] != "#":
            c, p, _, _, b = line.split("\t"); bg[int(p)] = int(b)
bgc = np.array([bg[p] for p in core])
E = json.load(open("docs/data/efficiencies.json"))["DJ"]; starts = np.array(E["start"], int); gc = np.array([np.nan if v is None else v for v in E["gc"]], float)
row = {int(s): i for i, s in enumerate(starts)}
pres = np.full(len(starts), np.nan)
for i, s in enumerate(starts):
    m = (core >= s) & (core < s + 250)
    if m.sum() >= 20:
        pres[i] = np.minimum(bgc[m], 5).mean() / 5
zero = [r["sample"] for r in csv.DictReader(open("docs/data/cohort.tsv"), delimiter="\t") if r.get("DJ.step") not in ("", "NA") and abs(float(r["DJ.step"])) < 0.2][:a.max_genomes]
acc, cnt = np.zeros(len(starts)), np.zeros(len(starts))
for s in zero:
    f = f"{a.cache}/{s}.estimate.json.gz"
    if not os.path.exists(f):
        continue
    for w in json.load(gzip.open(f, "rt"))["classes"]["DJ"]["windows"]:
        i = row.get(int(w["start"]))
        if i is not None and w.get("usable") and w.get("cn") is not None:
            acc[i] += w["cn"]; cnt[i] += 1
raw = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
ok = np.isfinite(pres) & np.isfinite(raw) & np.isfinite(gc)
print(f"windows {ok.sum()}; exact-k-mer presence per haplotype: mean {pres[ok].mean():.3f}")
print(f"corr(presence, window GC) = {np.corrcoef(pres[ok], gc[ok])[0, 1]:+.3f}; corr(raw estimate over {len(zero)} mode genomes, window GC) = {np.corrcoef(raw[ok], gc[ok])[0, 1]:+.3f}; corr(raw estimate, presence) = {np.corrcoef(raw[ok], pres[ok])[0, 1]:+.3f}")
for lo, hi in ((0, 0.35), (0.35, 0.40), (0.40, 0.45), (0.45, 0.50), (0.50, 0.60), (0.60, 1)):
    m = ok & (gc >= lo) & (gc < hi)
    if m.sum():
        print(f"  GC {lo:.2f}-{hi:.2f}: {m.sum():4d} windows, presence {pres[m].mean():.3f}, raw estimate {raw[m].mean():.2f}")
