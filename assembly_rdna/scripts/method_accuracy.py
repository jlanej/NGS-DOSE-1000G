#!/usr/bin/env python
"""Which is more accurate for rDNA copy number: the HPRC r2 assembly or NGS-DOSE short reads?

Four tests, none of which treats either method as the truth:
  1. Inheritance. A child's true dosage follows its parents'; measurement error does not. The parents are
     measured once (by NGS-DOSE, and separately by the independent Hall et al. 2021 pipeline) and each
     method's value for the child is correlated with the midparent. Paired bootstrap over trios.
     5S, which the assemblies hold whole, is the control: there the test should find no difference.
  2. A known truth. The distal junction is present once per acrocentric short arm, 10 per diploid genome.
  3. An orthogonal assay. Potapova et al. 2025 (Cell Genomics 5:101031, Table S1) measured total rDNA by
     ddPCR in 12 LCLs, 9 of them in the 1000 Genomes 30x set. NGS-DOSE counts for the 7 not yet in the
     cohort run were made here in fetch mode (potapova/run_fetch.sh) and calibrated with the cohort's saved
     window efficiencies (potapova/est/cohort.tsv).
  5. A second sequencing pipeline. Google Health's GIAB NovaSeq 6000 PCR-free 30x BAMs (Baid et al. 2020; bwa-mem
     0.7.17 on hs38DH) for HG002/HG003/HG004 (ddPCR-measured, not in 1000 Genomes) and for the CEPH trio NA12878/
     NA12891/NA12892, which NYGC also sequenced; counted in fetch mode (novaseq/run_fetch.sh) and calibrated with the
     NYGC cohort's window efficiencies. The CEPH trio shows what the pipeline change alone does.
  4. Per array. Potapova's FISH gives each array's share of the total (Table S2); the assemblies' units on
     contigs assigned to each acrocentric are compared with it.
Writes tables/method_accuracy.json, tables/potapova_comparison.tsv, tables/potapova_arrays.tsv and
figures/fig4_method_accuracy.png, figures/fig5_truths_and_assays.png.
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(BASE)
TAB, FIG, POT = (os.path.join(BASE, d) for d in ("tables", "figures", "potapova"))
rng = np.random.default_rng(20260925)

hap = pd.read_csv(f"{TAB}/haplotype_rdna.tsv", sep="\t")
hh = hap[hap.source.isin(["hprc", "hpp"]) & (hap.method == "hifiasm")]
co = pd.read_csv(f"{REPO}/docs/data/cohort.tsv", sep="\t").set_index("sample")
ped = pd.read_csv(f"{REPO}/meta/20130606_g1k_3202_samples_ped_population.txt", sep=" ").set_index("SampleID")
hall = pd.read_csv(f"{REPO}/meta/hall2021_MOESM1.txt", sep="\t").set_index("Sample")
hall.columns = [c.strip() for c in hall.columns]

asm = hh.groupby("person").agg(n=("assembly", "size"), a45=("n_18S", "sum"), a5=("n_5S_units", "sum"), aDJ=("n_DJ_copies", "sum"))
asm = asm[asm.n == 2].join(ped[["FatherID", "MotherID", "Superpopulation"]], how="inner")

# ---------------- 1. inheritance ----------------
def midp(s, kids):
    return (s.reindex(kids.FatherID).values + s.reindex(kids.MotherID).values) / 2

kids = asm[asm.FatherID.astype(str) != "0"].copy()
parents = {"NGS-DOSE": {"45S": co["rDNA45S.cn"], "5S": co["rDNA5S.cn"]},
           "Hall 2021": {"45S": hall["HC.18S.CN"] * 2, "5S": hall["HC.5S.CN"] * 2}}
child = {"45S": ("a45", co["rDNA45S.cn"]), "5S": ("a5", co["rDNA5S.cn"])}

def centre(d, c):
    return d[c] - d.groupby("Superpopulation")[c].transform("mean")

def paired(cls, pby, B=5000):
    acol, ngs = child[cls]
    d = kids.copy()
    d["mid"] = midp(parents[pby][cls], d)
    d["ngs"] = ngs.reindex(d.index).values
    d = d.dropna(subset=[acol, "ngs", "mid"]).copy()
    for c in (acol, "ngs", "mid"):
        d[c] = centre(d, c)
    x, y, z = d[acol].values, d["ngs"].values, d["mid"].values
    n = len(d)
    ra, rn = np.corrcoef(x, z)[0, 1], np.corrcoef(y, z)[0, 1]
    bs = np.array([(np.corrcoef(x[i], z[i])[0, 1], np.corrcoef(y[i], z[i])[0, 1]) for i in (rng.integers(0, n, n) for _ in range(B))])
    dd = bs[:, 1] - bs[:, 0]
    return dict(cls=cls, parents=pby, n=n, r_assembly=ra, r_assembly_lo=np.percentile(bs[:, 0], 2.5), r_assembly_hi=np.percentile(bs[:, 0], 97.5),
                r_ngsdose=rn, r_ngsdose_lo=np.percentile(bs[:, 1], 2.5), r_ngsdose_hi=np.percentile(bs[:, 1], 97.5),
                diff=rn - ra, diff_lo=np.percentile(dd, 2.5), diff_hi=np.percentile(dd, 97.5), p_ngsdose_better=float((dd > 0).mean()),
                reliability_ratio=(ra / rn) ** 2, points=d[[acol, "ngs", "mid"]].rename(columns={acol: "assembly"}))

inh = [paired(c, p) for c in ("45S", "5S") for p in ("NGS-DOSE", "Hall 2021")]
def inh_r(child_vals, pby):
    d = kids.copy(); d["mid"] = midp(parents[pby]["45S"], d); d["c"] = child_vals.reindex(d.index).values
    d["a"] = d["a45"]; d["ngs"] = co["rDNA45S.cn"].reindex(d.index).values
    d = d.dropna(subset=["a", "ngs", "mid", "c"]).copy()
    for c in ("c", "mid"):
        d[c] = centre(d, c)
    return dict(n=len(d), r=float(np.corrcoef(d.c, d.mid)[0, 1]))
inh18 = {p: inh_r(co["rDNA45S.18S.flat"], p) for p in ("NGS-DOSE", "Hall 2021")}
# assembly alone, with Hall-measured parents: the largest set (includes African families)
d = kids.copy(); d["mid"] = midp(parents["Hall 2021"]["45S"], d); d = d.dropna(subset=["a45", "mid"]).copy()
d["a45"], d["mid"] = centre(d, "a45"), centre(d, "mid")
r_big = np.corrcoef(d.a45, d.mid)[0, 1]
bsb = [np.corrcoef(d.a45.values[i], d.mid.values[i])[0, 1] for i in (rng.integers(0, len(d), len(d)) for _ in range(5000))]

# ---------------- 2. known truth: distal junction ----------------
dj = asm.join(co[["DJ.cn"]], how="inner").dropna(subset=["DJ.cn"])
dj_stats = dict(n=len(dj), assembly_mean=dj.aDJ.mean(), assembly_sd=dj.aDJ.std(), assembly_exact10=float((dj.aDJ == 10).mean()),
                ngsdose_mean=dj["DJ.cn"].mean(), ngsdose_sd=dj["DJ.cn"].std(),
                ngsdose_within_half=float((abs(dj["DJ.cn"] - co["DJ.cn"].median()) < 0.5).mean()),
                assembly_all_9_11=float(asm.aDJ.between(9, 11).mean()), assembly_all_n=int(len(asm)), cohort_DJ_median=float(co["DJ.cn"].median()))
dj["ngs_step"] = (dj["DJ.cn"] - co["DJ.cn"].median()).round().astype(int)
dj["asm_step"] = (dj.aDJ - 10).astype(int)
dj_stats.update(ngs_steps=int((dj.ngs_step != 0).sum()), ngs_steps_confirmed=int(((dj.ngs_step != 0) & (dj.ngs_step == dj.asm_step)).sum()),
                ngs_normal=int((dj.ngs_step == 0).sum()), ngs_normal_asm_10=int(((dj.ngs_step == 0) & (dj.asm_step == 0)).sum()),
                ngs_normal_asm_low=int(((dj.ngs_step == 0) & (dj.asm_step < 0)).sum()), ngs_normal_asm_high=int(((dj.ngs_step == 0) & (dj.asm_step > 0)).sum()))

# ---------------- 3. ddPCR (Potapova Table S1) ----------------
s1 = pd.read_csv(f"{POT}/potapova2025_tableS1.tsv", sep="\t").set_index("sample")
extra = pd.read_csv(f"{POT}/est/cohort.tsv", sep="\t").set_index("sample") if os.path.exists(f"{POT}/est/cohort.tsv") else pd.DataFrame()
ngs = co["rDNA45S.cn"].combine_first(extra["rDNA45S.cn"]) if len(extra) else co["rDNA45S.cn"]
flat = co["rDNA45S.18S.flat"].combine_first(extra["rDNA45S.18S.flat"]) if len(extra) else co["rDNA45S.18S.flat"]
src = pd.Series("cohort run", index=co.index)
if len(extra):
    src = src.combine_first(pd.Series("fetched here", index=extra.index))
NOV = os.path.join(BASE, "novaseq")
gg = pd.read_csv(f"{NOV}/est/google/cohort.tsv", sep="\t").set_index("sample")
gn = pd.read_csv(f"{NOV}/est/nygc/cohort.tsv", sep="\t").set_index("sample")
pc = s1.copy()
pc["ngsdose"] = ngs.reindex(pc.index)
pc["ngsdose_source"] = src.reindex(pc.index)
pc["ratio18S_flat"] = flat.reindex(pc.index)
pc["assembly_18S"] = asm.a45.reindex(pc.index)
hg002 = hap[hap.assembly.str.startswith("hg002")]
# T2T-HG002 v1.1 represents the centres of 9 of its 10 rDNA arrays as N-gaps (marbl/HG002 README; Hansen et al. 2025), so its
# 18S count is the flanking units only, not an assembly measurement: kept for display, excluded from every statistic.
pc["hg002_scaffold_18S"] = np.nan
pc.loc["HG002", "hg002_scaffold_18S"] = hg002.n_18S.sum()
pc["hall"] = (hall["HC.18S.CN"] * 2).reindex(pc.index)
pc["ngsdose_google"] = gg["rDNA45S.cn"].reindex(pc.index)
pc["ratio18S_flat_google"] = gg["rDNA45S.18S.flat"].reindex(pc.index)
pc["ngsdose_any"] = pc.ngsdose.combine_first(pc.ngsdose_google)
pc["pipeline"] = np.where(pc.ngsdose.notna(), "NYGC 1000G NovaSeq", np.where(pc.ngsdose_google.notna(), "Google GIAB NovaSeq", ""))
pc.to_csv(f"{TAB}/potapova_comparison.tsv", sep="\t")

def agree(col, ref="ddpcr"):
    x = pc.dropna(subset=[col, ref])
    lr = np.log(x[col] / x[ref])
    return dict(n=len(x), median_ratio=float(np.exp(np.median(lr))), sd_log=float(lr.std(ddof=1)) if len(x) > 1 else None,
                mean_abs_pct=float(100 * np.mean(np.abs(x[col] / x[ref] - 1))),
                r=float(np.corrcoef(x[col], x[ref])[0, 1]) if len(x) > 2 else None, samples=list(x.index))
dd_stats = {c: agree(c) for c in ["ngsdose", "conkord", "ratio18S_flat", "assembly_18S", "hall"]}
dd_stats["ngsdose_vs_conkord"] = agree("ngsdose", "conkord")
for c in ["ngsdose_google", "ratio18S_flat_google", "ngsdose_any"]:
    dd_stats[c] = agree(c)
_z = pc.dropna(subset=["ngsdose_any"])
dd_stats["conkord_same12"] = dict(n=len(_z), r=float(np.corrcoef(_z.conkord, _z.ddpcr)[0, 1]))
_x = pc.dropna(subset=["ngsdose"])
dd_stats["conkord_same9"] = dict(n=len(_x), r=float(np.corrcoef(_x.conkord, _x.ddpcr)[0, 1]),
                                 mean_abs_pct=float(100 * np.mean(abs(_x.conkord / _x.ddpcr - 1))), median_ratio=float(np.median(_x.conkord / _x.ddpcr)))
dd_stats["residual_sd_log"] = {c: float(np.log(_x[c] / _x.ddpcr).std(ddof=1)) for c in ["ngsdose", "conkord", "ratio18S_flat"]}
dd_stats["bias_pct"] = {c: float(100 * (np.exp(np.log(_x[c] / _x.ddpcr).mean()) - 1)) for c in ["ngsdose", "conkord", "ratio18S_flat"]}
_y = _x.drop(index="HG02053", errors="ignore")
dd_stats["without_HG02053"] = {c: dict(r=float(np.corrcoef(_y[c], _y.ddpcr)[0, 1]), mean_abs_pct=float(100 * np.mean(abs(_y[c] / _y.ddpcr - 1))))
                               for c in ["ngsdose", "conkord", "ratio18S_flat"]}
def rci(x, y):
    r = float(np.corrcoef(x, y)[0, 1]); n = len(x); se = 1 / np.sqrt(n - 3)
    return dict(n=n, r=r, lo=float(np.tanh(np.arctanh(r) - 1.96 * se)), hi=float(np.tanh(np.arctanh(r) + 1.96 * se)))
_a = pc.dropna(subset=["assembly_18S"])
dd_stats["r_ci"] = {"ngsdose_any": rci(_z.ngsdose_any, _z.ddpcr), "conkord": rci(pc.conkord, pc.ddpcr),
                    "ratio18S_any": rci(_z.ratio18S_flat.combine_first(_z.ratio18S_flat_google), _z.ddpcr),
                    "assembly": rci(_a.assembly_18S, _a.ddpcr), "ngsdose_same_asm": rci(_a.ngsdose_any, _a.ddpcr),
                    "conkord_same_asm": rci(_a.conkord, _a.ddpcr)}
pc["ratio18S_any"] = pc.ratio18S_flat.combine_first(pc.ratio18S_flat_google)
dd_stats["ratio18S_any"] = agree("ratio18S_any")
# ddPCR's own replicate spread, for scale
dd_stats["ddpcr_cv_median"] = float((s1.ddpcr_sd / s1.ddpcr).median())

# ---------------- 4. per array: FISH shares vs assembled units per chromosome ----------------
s2 = pd.read_csv(f"{POT}/potapova2025_tableS2.tsv", sep="\t")
s2["chrom"] = "chr" + s2.array.str[:2]
fish = s2.groupby(["sample", "chrom"]).fish_units.sum().rename("fish").reset_index()
fish["fish_share"] = fish.fish / fish.groupby("sample").fish.transform("sum")
arr = pd.read_csv(f"{TAB}/arrays_all.tsv.gz", sep="\t")
arr["person"] = arr.assembly.map(hap.set_index("assembly").person)
arr.loc[arr.assembly.str.startswith("hg002"), "person"] = np.nan   # rDNA gapped by design (see above)
pa = arr[arr.person.isin(fish["sample"])].copy()
per = pa.dropna(subset=["chrom"]).groupby(["person", "chrom"]).n_18S.sum().rename("assembled").reset_index()
tot = pa.groupby("person").n_18S.sum().rename("assembled_total")
placed = pa.dropna(subset=["chrom"]).groupby("person").n_18S.sum().rename("assembled_placed")
pm = fish.merge(per, left_on=["sample", "chrom"], right_on=["person", "chrom"], how="inner").drop(columns="person")
pm = pm.merge(tot, left_on="sample", right_index=True).merge(placed, left_on="sample", right_index=True)
pm["asm_share_of_placed"] = pm.assembled / pm.assembled_placed
pm.to_csv(f"{TAB}/potapova_arrays.tsv", sep="\t", index=False)
arr_stats = dict(samples=sorted(pm["sample"].unique()), n_chrom=len(pm),
                 r_share=float(np.corrcoef(pm.fish_share, pm.asm_share_of_placed)[0, 1]) if len(pm) > 2 else None,
                 spearman_share=float(stats.spearmanr(pm.fish_share, pm.asm_share_of_placed)[0]) if len(pm) > 2 else None,
                 placed_fraction={s: float(placed.get(s, 0) / tot.get(s, np.nan)) for s in pm["sample"].unique()})

ft = fish.groupby("sample").fish.sum().rename("fish_total").to_frame().join(pc[["conkord", "ddpcr", "ngsdose", "ngsdose_google", "assembly_18S"]])
ft.to_csv(f"{TAB}/potapova_fish_totals.tsv", sep="\t")
def _cmp(a, b):
    x = ft.dropna(subset=[a, b])
    return dict(n=len(x), r=float(np.corrcoef(x[a], x[b])[0, 1]) if len(x) > 2 else None, median_ratio=float(np.median(x[a] / x[b])))
ft["ngsdose_any"] = ft.ngsdose.combine_first(ft.ngsdose_google)
fish_tot = dict(fish_minus_conkord_max=float((ft.fish_total - ft.conkord).abs().max()),
                ngsdose=_cmp("ngsdose", "fish_total"), ngsdose_any=_cmp("ngsdose_any", "fish_total"), assembly=_cmp("assembly_18S", "fish_total"))
pm["asm_units_placed"] = pm.assembled
arr_stats["r_units"] = float(np.corrcoef(pm.fish, pm.assembled)[0, 1]) if len(pm) > 2 else None
# ---------------- 5. second NovaSeq pipeline: same people, two pipelines ----------------
BR = {"rDNA45S.cn": "45S, NGS-DOSE", "rDNA45S.18S.flat": "45S, 18S depth ratio", "rDNA5S.cn": "5S, NGS-DOSE",
      "DJ.cn": "distal junction", "truth.auto": "autosomal control", "insert_median": "insert size", "gc_rel_65": "GC-rich coverage (65%)"}
ceph = [s_ for s_ in gn.index if s_ in gg.index]
br = pd.DataFrame([dict(sample=s_, quantity=q, nygc=gn.loc[s_, q], google=gg.loc[s_, q], pct=100 * (gg.loc[s_, q] / gn.loc[s_, q] - 1))
                   for s_ in ceph for q in BR])
br.to_csv(f"{TAB}/novaseq_bridge.tsv", sep="\t", index=False)
nov = pd.concat([gg.assign(pipeline="Google GIAB NovaSeq"), gn.assign(pipeline="NYGC 1000G NovaSeq")])
nov.reset_index()[["sample", "pipeline", "depth", "insert_median", "gc_rel_65", "truth.auto", "truth.chrX", "truth.chrY", "DJ.cn",
                   "rDNA45S.cn", "rDNA45S.cn_single", "rDNA45S.18S.flat", "rDNA5S.cn"]].to_csv(f"{TAB}/novaseq_estimates.tsv", sep="\t", index=False)
bridge = {q: dict(pct=[float(v) for v in br[br.quantity == q].pct], mean_abs_pct=float(br[br.quantity == q].pct.abs().mean())) for q in BR}
bridge["samples"] = ceph
bridge["google_DJ"] = [float(v) for v in gg["DJ.cn"]]
bridge["google_gc65"] = {k: float(v) for k, v in gg.gc_rel_65.items()}
bridge["cohort_gc65_max"] = float(co.gc_rel_65.max()); bridge["cohort_insert_median"] = float(co.insert_median.median())
bridge["google_insert"] = {k: float(v) for k, v in gg.insert_median.items()}
out = dict(inheritance=[{k: v for k, v in r.items() if k != "points"} for r in inh],
           assembly_hall_parents=dict(n=len(d), r=r_big, lo=float(np.percentile(bsb, 2.5)), hi=float(np.percentile(bsb, 97.5))),
           distal_junction=dj_stats, ddpcr=dd_stats, inheritance_18S_ratio=inh18, arrays=arr_stats, fish_totals=fish_tot, novaseq=bridge)
json.dump(out, open(f"{TAB}/method_accuracy.json", "w"), indent=1, default=float)

# ---------------- figures ----------------
try:
    apply_figure_style(sizes=(8, 7, 6))
except NameError:
    plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
CA, CN, CC, CH = "#b15928", "#1f78b4", "#6a3d9a", "#999999"
PN, PG = "#404040", "#e66101"   # pipelines: NYGC 1000 Genomes NovaSeq, Google GIAB NovaSeq
GOOGLE = dict(marker="s", s=16, facecolor="white", edgecolor=CN, linewidths=1.0)   # NGS-DOSE on the Google pipeline: open square
def letter(ax, s, x=-0.22, y=1.13):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")

# Fig 4: inheritance
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.8), gridspec_kw=dict(width_ratios=[1, 1, 1.2]))
r45 = inh[0]
p = r45["points"]
lim = (min(p.min().min(), -160) * 1.05, max(p.max().max(), 190) * 1.05)
for j, (col, lab, c) in enumerate([("assembly", "the assembly", CA), ("ngs", "NGS-DOSE", CN)]):
    a = ax[j]
    a.scatter(p["mid"], p[col], s=12, color=c, lw=0)
    b = np.polyfit(p["mid"], p[col], 1); xs = np.linspace(p["mid"].min(), p["mid"].max(), 10)
    a.plot(xs, np.polyval(b, xs), color="k", lw=0.8)
    rv = r45["r_assembly" if col == "assembly" else "r_ngsdose"]
    a.set_title(f"Child measured by {lab}\nPearson r = {rv:.2f} ({r45['n']} trios)", loc="left")
    a.set_xlabel("parents' mean, NGS-DOSE 45S"); a.set_ylabel("child's 45S copies" if j == 0 else "")
    a.set_ylim(*lim); letter(a, "ab"[j], -0.25, 1.2)
a = ax[2]
ys = np.arange(len(inh))[::-1]
for y, r in zip(ys, inh):
    a.errorbar(r["r_assembly"], y + 0.13, xerr=[[r["r_assembly"] - r["r_assembly_lo"]], [r["r_assembly_hi"] - r["r_assembly"]]], fmt="o", ms=3.5, color=CA, lw=1)
    a.errorbar(r["r_ngsdose"], y - 0.13, xerr=[[r["r_ngsdose"] - r["r_ngsdose_lo"]], [r["r_ngsdose_hi"] - r["r_ngsdose"]]], fmt="s", ms=3.2, color=CN, lw=1)
a.set_yticks(ys); a.set_yticklabels([f"{r['cls']}, parents:\n{r['parents']} ({r['n']})" for r in inh])
a.axvline(0, color=CH, lw=0.6); a.set_xlim(-0.4, 1.05); a.set_ylim(-0.6, len(inh) - 0.2)
a.set_xlabel("Pearson r with parents' mean")
a.plot([], [], "o", color=CA, label="child by assembly"); a.plot([], [], "s", color=CN, label="child by NGS-DOSE")
a.legend(frameon=False, loc="upper left", bbox_to_anchor=(-0.02, 1.02), fontsize=6)
a.set_title("45S: NGS-DOSE wins; 5S: tie", loc="left"); letter(a, "c", -0.55, 1.2)
fig.tight_layout(w_pad=2.2)
fig.savefig(f"{FIG}/fig4_method_accuracy.png", dpi=300)
plt.close(fig)

# Fig 5: known truths and orthogonal assays
fig, axs = plt.subplots(2, 3, figsize=(7.4, 5.4))
ax = axs.ravel()
a = ax[0]
jit = rng.uniform(-0.12, 0.12, len(dj))
col = np.where(dj.ngs_step != 0, CC, CN)
a.scatter(dj["DJ.cn"], dj.aDJ + jit, s=10, c=col, lw=0)
a.axhline(10, color=CH, lw=0.6, ls="--"); a.axvline(10, color=CH, lw=0.6, ls="--")
for k in (-1, 1):
    a.axvline(dj_stats["cohort_DJ_median"] + k, color=CH, lw=0.5, ls=":")
a.set_xlabel("NGS-DOSE DJ copies"); a.set_ylabel("assembly DJ copies"); a.set_yticks([8, 9, 10, 11])
a.set_title(f"Distal junction (truth 10), n={dj_stats['n']}:\nsteps agree; NGS-DOSE {dj_stats['ngsdose_mean']:.2f}", loc="left")
a.text(0.03, 0.97, "purple: NGS-DOSE step", transform=a.transAxes, fontsize=5.5, va="top", color=CC)
letter(a, "a", -0.3)
a = ax[1]
x = pc.dropna(subset=["ngsdose"]); xg = pc.dropna(subset=["ngsdose_google"])
lo, hi = 400, 740
a.plot([lo, hi], [lo, hi], color=CH, lw=0.8)
a.scatter(x.ddpcr, x.ratio18S_flat, marker="o", s=9, color="#e31a1c", lw=0, label="18S depth ratio", zorder=3)
a.scatter(xg.ddpcr, xg.ratio18S_flat_google, marker="o", s=10, facecolor="white", edgecolor="#e31a1c", lw=0.8, zorder=3)
a.scatter(pc.ddpcr, pc.conkord, marker="^", s=12, color=CC, alpha=0.8, lw=0, label="CONKORD (their reads)", zorder=3)
a.errorbar(x.ddpcr, x.ngsdose, xerr=x.ddpcr_sd, fmt="s", ms=3.5, color=CN, lw=0.6, label="NGS-DOSE", zorder=4)
a.errorbar(xg.ddpcr, xg.ngsdose_google, xerr=xg.ddpcr_sd, fmt="none", ecolor=CN, lw=0.6, zorder=4)
a.scatter(xg.ddpcr, xg.ngsdose_google, zorder=5, **GOOGLE)
for s_, r_ in xg.iterrows():
    a.annotate(s_, (r_.ddpcr, r_.ngsdose_google), xytext=(4, -8), textcoords="offset points", fontsize=5, color=PG)
a.text(0.97, 0.03, "filled: NYGC 1000G NovaSeq\nopen: Google GIAB NovaSeq", transform=a.transAxes, ha="right", va="bottom", fontsize=5.3, color=PG)
a.annotate("HG02053", (713, x.loc["HG02053", "ngsdose"]), xytext=(-38, -11), textcoords="offset points", fontsize=5.5)
a.set_xlim(lo, hi); a.set_ylim(lo, hi)
a.set_xlabel("ddPCR 45S copies (±SD)"); a.set_ylabel("short-read 45S copies")
st, sc, sa12 = dd_stats["ngsdose"], dd_stats["conkord_same12"], dd_stats["ngsdose_any"]
a.set_title(f"ddPCR, n = {sa12['n']}: Pearson r\nNGS-DOSE {sa12['r']:.2f}, CONKORD {sc['r']:.2f}", loc="left")
a.legend(frameon=False, loc="upper left", fontsize=5.5); letter(a, "b", -0.3)
a = ax[2]
xa = pc.dropna(subset=["assembly_18S"])
a.plot([0, hi], [0, hi], color=CH, lw=0.8)
for s_, r_ in xa.iterrows():
    a.plot([r_.ddpcr, r_.ddpcr], [r_.assembly_18S, r_.ngsdose if pd.notna(r_.ngsdose) else r_.assembly_18S], color="#dddddd", lw=0.8, zorder=1)
a.errorbar(xa.ddpcr, xa.assembly_18S, xerr=xa.ddpcr_sd, fmt="o", ms=4, color=CA, lw=0.6, label=f"HPRC assembly ({len(xa)})", zorder=3)
xn = xa.dropna(subset=["ngsdose"]); xgn = xa.dropna(subset=["ngsdose_google"])
for s_, r_ in xgn.iterrows():
    a.plot([r_.ddpcr, r_.ddpcr], [r_.assembly_18S, r_.ngsdose_google], color="#dddddd", lw=0.8, zorder=1)
a.scatter(xn.ddpcr, xn.ngsdose, marker="s", s=14, color=CN, label=f"NGS-DOSE, NYGC ({len(xn)})", zorder=3)
h2 = pc.loc["HG002"]
a.plot([h2.ddpcr, h2.ddpcr], [h2.hg002_scaffold_18S, h2.ngsdose_google], color="#dddddd", lw=0.8, ls=":", zorder=1)
a.scatter([h2.ddpcr], [h2.hg002_scaffold_18S], marker="o", s=18, facecolor="white", edgecolor=CA, linewidths=0.9, zorder=3,
          label="HG002 T2T (rDNA gapped; excluded)")
a.scatter([h2.ddpcr], [h2.ngsdose_google], label="NGS-DOSE, Google (HG002)", zorder=4, **GOOGLE)
a.annotate("HG002: 9 of 10 arrays\nare N-gaps by design", (h2.ddpcr, h2.hg002_scaffold_18S), xytext=(-86, 10),
           textcoords="offset points", fontsize=5, color=CA)
a.set_xlim(0, hi); a.set_ylim(0, hi)
sa = dd_stats["assembly_18S"]
a.set_xlabel("ddPCR 45S copies (±SD)"); a.set_ylabel("45S copies (diploid)")
a.set_title(f"HPRC assemblies hold {100*sa['median_ratio']:.0f}%\nof ddPCR ({len(xa)} people)", loc="left")
a.legend(frameon=False, loc="upper left", fontsize=5.5); letter(a, "c", -0.3)
a = ax[3]
ftn = ft.dropna(subset=["ngsdose"]); fta = ft.dropna(subset=["assembly_18S"])
a.plot([0, hi], [0, hi], color=CH, lw=0.8)
a.scatter(ftn.fish_total, ftn.ngsdose, marker="s", s=14, color=CN, label=f"NGS-DOSE, NYGC ({len(ftn)})", zorder=3)
ftg = ft.dropna(subset=["ngsdose_google"])
a.scatter(ftg.fish_total, ftg.ngsdose_google, label=f"NGS-DOSE, Google ({len(ftg)})", zorder=4, **GOOGLE)
a.scatter(fta.fish_total, fta.assembly_18S, marker="o", s=16, color=CA, label=f"assembly ({len(fta)})", zorder=3)
a.set_xlim(0, hi); a.set_ylim(0, hi)
a.set_xlabel("FISH, summed over 10 arrays"); a.set_ylabel("45S copies (diploid)")
fn = fish_tot["ngsdose_any"]
a.set_title(f"FISH totals (= CONKORD) vs\nNGS-DOSE: Pearson r = {fn['r']:.2f}", loc="left")
a.legend(frameon=False, loc="upper left", fontsize=5.5); letter(a, "d", -0.3)
a.text(0.97, 0.03, "FISH copies = FISH share\n× CONKORD total, summed", transform=a.transAxes, ha="right", va="bottom", fontsize=5.3, color=CH)
a = ax[4]
for s_, g in pm.groupby("sample"):
    a.scatter(g.fish, g.assembled, s=12, label=s_)
m_ = max(pm.fish.max(), pm.assembled.max()) * 1.05
a.plot([0, m_], [0, m_], color=CH, lw=0.8); a.set_xlim(0, m_); a.set_ylim(0, m_)
a.set_xlabel("FISH units, both homologues"); a.set_ylabel("assembled units, same chromosome")
a.set_title(f"Per chromosome, units:\nPearson r = {arr_stats['r_units']:.2f} (n = {arr_stats['n_chrom']})", loc="left")
a.legend(frameon=False, fontsize=5, loc="upper left"); letter(a, "e", -0.3)
a.text(0.97, 0.97, "units on one chromosome\n(both homologues)", transform=a.transAxes, ha="right", va="top", fontsize=5.3, color=CH)
a = ax[5]
for s_, g in pm.groupby("sample"):
    a.scatter(g.fish_share, g.asm_share_of_placed, s=12, label=s_)
a.plot([0, 0.5], [0, 0.5], color=CH, lw=0.8); a.set_xlim(0, 0.5); a.set_ylim(0, 0.55)
a.set_xlabel("FISH share of the rDNA"); a.set_ylabel("share of placed assembled units")
a.set_title(f"Per chromosome, shares:\nPearson r = {arr_stats['r_share']:.2f}", loc="left"); letter(a, "f", -0.3)
a.text(0.97, 0.97, "share = fraction of the person's\nrDNA on that chromosome", transform=a.transAxes, ha="right", va="top", fontsize=5.3, color=CH)
fig.tight_layout(h_pad=1.6, w_pad=0.9)
fig.savefig(f"{FIG}/fig5_truths_and_assays.png", dpi=300)
plt.close(fig)

# Fig 6: a second NovaSeq pipeline
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.9), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
a = ax[0]
rows = ["rDNA45S.cn", "rDNA45S.18S.flat", "rDNA5S.cn", "DJ.cn", "truth.auto", None, "gc_rel_65", "insert_median"]
mk = dict(zip(ceph, ["o", "s", "D"]))
for y, q in enumerate(rows[::-1]):
    if q is None:
        continue
    b_ = br[br.quantity == q]
    for _, r_ in b_.iterrows():
        a.scatter(r_.pct, y, marker=mk[r_["sample"]], s=16, color=PG if q in ("gc_rel_65", "insert_median") else (CN if "NGS-DOSE" in BR[q] else PN),
                  lw=0, zorder=3, alpha=0.9)
a.axhline(2, color=CH, lw=0.5, ls=":")
a.set_yticks([y for y, q in enumerate(rows[::-1]) if q]); a.set_yticklabels([BR[q] for q in rows[::-1] if q])
a.axvline(0, color=CH, lw=0.8)
a.set_xlabel("Google minus NYGC, same person (%)")
for s_ in ceph:
    a.scatter([], [], marker=mk[s_], s=16, color=PN, label=s_)
a.legend(frameon=False, fontsize=5.5, loc="center right", bbox_to_anchor=(1.0, 0.36))
import matplotlib.transforms as mtr
a.text(0.02, 1.55, "library properties", transform=mtr.blended_transform_factory(a.transAxes, a.transData), fontsize=5.5, color=PG, va="bottom")
m45, m18 = bridge["rDNA45S.cn"]["mean_abs_pct"], bridge["rDNA45S.18S.flat"]["mean_abs_pct"]
a.set_title(f"CEPH trio, two pipelines: 45S moves\n{m45:.1f}% (NGS-DOSE), {m18:.0f}% (18S ratio)", loc="left")
letter(a, "a", -0.62, 1.2)
a = ax[1]
lo, hi = 400, 740
a.plot([lo, hi], [lo, hi], color=CH, lw=0.8)
a.errorbar(x.ddpcr, x.ngsdose, xerr=x.ddpcr_sd, fmt="s", ms=3.5, color=PN, lw=0.6, label=f"NYGC 1000G NovaSeq ({len(x)})", zorder=3)
a.errorbar(xg.ddpcr, xg.ngsdose_google, xerr=xg.ddpcr_sd, fmt="s", ms=4, mfc="white", mec=PG, mew=1.0, color=PG, lw=0.6,
           label=f"Google GIAB NovaSeq ({len(xg)})", zorder=4)
for s_, r_ in xg.iterrows():
    a.annotate(s_, (r_.ddpcr, r_.ngsdose_google), xytext=(4, -8), textcoords="offset points", fontsize=5, color=PG)
a.set_xlim(lo, hi); a.set_ylim(lo, hi)
a.set_xlabel("ddPCR 45S copies (±SD)"); a.set_ylabel("NGS-DOSE 45S copies")
sg = dd_stats["ngsdose_google"]
a.set_title(f"Against ddPCR: level {sg['median_ratio']:.2f}\n(Google), {st['median_ratio']:.2f} (NYGC)", loc="left")
a.legend(frameon=False, fontsize=5.5, loc="upper left"); letter(a, "b", -0.3, 1.2)
a = ax[2]
a.scatter(co.insert_median, co.gc_rel_65, s=4, color="#cccccc", lw=0, label=f"NYGC cohort ({len(co)})", zorder=1)
a.scatter(gn.insert_median, gn.gc_rel_65, s=18, color=PN, lw=0, label="NYGC, CEPH trio", zorder=3)
a.scatter(gg.insert_median, gg.gc_rel_65, s=18, facecolor="white", edgecolor=PG, linewidths=1.0, label="Google, all 6", zorder=4)
for s_, xy in {"NA12891": (4, -3), "NA12892": (4, -3), "HG003": (4, -6)}.items():
    a.annotate(s_, (gg.loc[s_, "insert_median"], gg.loc[s_, "gc_rel_65"]), xytext=xy, textcoords="offset points", fontsize=5, color=PG)
a.annotate("HG002, HG004,\nNA12878", (gg.loc["HG002", "insert_median"], gg.loc["HG002", "gc_rel_65"]), xytext=(-12, 9), textcoords="offset points", fontsize=5, color=PG)
a.set_xlabel("median insert size (bp)"); a.set_ylabel("coverage at 65% GC (relative)")
a.set_title("Google libraries: shorter inserts,\ntwo beyond the NYGC GC range", loc="left")
a.legend(frameon=False, fontsize=5.5, loc="lower right"); letter(a, "c", -0.3, 1.2)
fig.tight_layout(w_pad=1.4)
fig.savefig(f"{FIG}/fig6_second_pipeline.png", dpi=300)
plt.close(fig)
print(json.dumps({k: out["ddpcr"][k] for k in ["ngsdose", "ngsdose_google", "ngsdose_any", "conkord_same12", "ratio18S_flat_google"]}, default=float, indent=0)[:2500])
print(json.dumps(bridge, default=float)[:1500]); print(json.dumps(fish_tot, default=float))

# ---------------- Fig 7: the whole comparison in one self-contained figure, and its companion table ----------------
import matplotlib.gridspec as gridspec
rc = dd_stats["r_ci"]
C18 = "#e31a1c"
# bootstrap CI for the 18S ratio's inheritance r (same trios, parents by NGS-DOSE)
_d = kids.copy(); _d["mid"] = midp(parents["NGS-DOSE"]["45S"], _d); _d["c"] = co["rDNA45S.18S.flat"].reindex(_d.index).values
_d["ngs"] = co["rDNA45S.cn"].reindex(_d.index).values
_d = _d.dropna(subset=["a45", "ngs", "mid", "c"]).copy()
for _c in ("c", "mid"):
    _d[_c] = centre(_d, _c)
_bs = [np.corrcoef(_d.c.values[k], _d.mid.values[k])[0, 1] for k in (rng.integers(0, len(_d), len(_d)) for _ in range(5000))]
inh18["NGS-DOSE"].update(lo=float(np.percentile(_bs, 2.5)), hi=float(np.percentile(_bs, 97.5)))

fig = plt.figure(figsize=(7.4, 11.0))
gs = gridspec.GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 1.3], hspace=0.62, wspace=0.62,
                       left=0.09, right=0.98, top=0.945, bottom=0.01)
# a: scatter, every method against ddPCR
a = fig.add_subplot(gs[0, 0])
a.plot([0, 750], [0, 750], color=CH, lw=0.8)
a.scatter(pc.ddpcr, pc.ratio18S_any, marker="o", s=10, color=C18, lw=0, label="18S depth ratio", zorder=3)
a.scatter(pc.ddpcr, pc.conkord, marker="^", s=12, color=CC, lw=0, label="CONKORD", zorder=3)
a.scatter(pc.ddpcr, pc.ngsdose, marker="s", s=13, color=CN, lw=0, label="NGS-DOSE, NYGC", zorder=4)
a.scatter(pc.ddpcr, pc.ngsdose_google, marker="s", s=13, facecolor="white", edgecolor=CN, linewidths=0.9, label="NGS-DOSE, Google", zorder=4)
a.scatter(pc.ddpcr, pc.assembly_18S, marker="o", s=14, color=CA, lw=0, label="HPRC assembly", zorder=4)
a.scatter([pc.loc["HG002", "ddpcr"]], [pc.loc["HG002", "hg002_scaffold_18S"]], marker="o", s=14, facecolor="white", edgecolor=CA, linewidths=0.9, zorder=4)
a.annotate("HG002 T2T: rDNA\ngapped, excluded", (pc.loc["HG002", "ddpcr"], pc.loc["HG002", "hg002_scaffold_18S"]), xytext=(-62, -4), textcoords="offset points", fontsize=5.3, color=CA)
a.set_xlim(400, 740); a.set_ylim(0, 750)
a.set_xlabel("ddPCR 45S copies"); a.set_ylabel("estimate, 45S copies")
a.set_title("All methods vs ddPCR\n(grey line: equal)", loc="left")
_h, _l = a.get_legend_handles_labels()
fig.legend(_h, ["18S depth ratio", "CONKORD", "NGS-DOSE (NYGC reads)", "NGS-DOSE (Google reads)", "HPRC assembly"], loc="upper center",
           ncol=5, frameon=False, fontsize=6.3, bbox_to_anchor=(0.5, 0.998), handletextpad=0.3, columnspacing=1.4)
_tab = [("NGS-DOSE", rc["ngsdose_any"], CN), ("CONKORD", rc["conkord"], CC), ("18S ratio", rc["ratio18S_any"], C18), ("assembly", rc["assembly"], CA)]
a.text(562, 470, "Pearson r (n)", fontsize=5.5, fontweight="bold", va="top")
for k_, (lab_, d_, c_) in enumerate(_tab):
    a.text(562, 430 - 42 * k_, f"{lab_}", fontsize=5.5, color=c_, va="top")
    a.text(738, 430 - 42 * k_, f"{d_['r']:.2f} ({d_['n']})", fontsize=5.5, color=c_, va="top", ha="right")
letter(a, "a", -0.34, 1.2)
# b: level against ddPCR, every method, every person
a = fig.add_subplot(gs[0, 1:3])
rows_a = [("NGS-DOSE, NYGC", pc.ngsdose, CN, "s", True), ("NGS-DOSE, Google", pc.ngsdose_google, CN, "s", False),
          ("CONKORD", pc.conkord, CC, "^", True), ("18S depth ratio", pc.ratio18S_any, C18, "o", True),
          ("HPRC assembly", pc.assembly_18S, CA, "o", True)]
for y, (lab, v, c, m, filled) in enumerate(rows_a[::-1]):
    rr = (v / pc.ddpcr).dropna()
    yj = y + rng.uniform(-0.13, 0.13, len(rr))
    a.scatter(rr, yj, marker=m, s=16, color=c if filled else "white", edgecolor=c, linewidths=0.9, zorder=3)
    a.plot([rr.median()] * 2, [y - 0.32, y + 0.32], color="k", lw=1.2, zorder=4)
    a.text(1.33, y, f"{rr.median():.2f}  (n = {len(rr)})", va="center", fontsize=6)
a.text(1.33, len(rows_a) - 0.45, "median", fontsize=6)
a.set_yticks(range(len(rows_a))); a.set_yticklabels([r[0] for r in rows_a[::-1]])
a.axvline(1, color=CH, lw=0.8, ls="--")
a.set_xlim(0, 1.62); a.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]); a.set_ylim(-0.5, len(rows_a) - 0.1)
a.set_xlabel("estimate ÷ ddPCR, per person (1 = agrees; black tick = median)")
a.set_title(f"Level: NGS-DOSE {100*(1-dd_stats['ngsdose_any']['median_ratio']):.0f}% below ddPCR, 18S ratio "
            f"{100*(dd_stats['ratio18S_any']['median_ratio']-1):.0f}% above;\nassemblies hold {100*dd_stats['assembly_18S']['median_ratio']:.0f}%", loc="left")
letter(a, "b", -0.25, 1.2)
# c: inheritance scatter
a = fig.add_subplot(gs[1, 0])
pp = inh[0]["points"]
for col, c, m, lab in [("ngs", CN, "s", "NGS-DOSE"), ("assembly", CA, "o", "assembly")]:
    a.scatter(pp["mid"], pp[col], marker=m, s=10, color=c, lw=0, alpha=0.9, zorder=3)
    b_ = np.polyfit(pp["mid"], pp[col], 1); xs = np.linspace(pp["mid"].min(), pp["mid"].max(), 10)
    a.plot(xs, np.polyval(b_, xs), color=c, lw=0.9)
    a.text(xs[-1] + 8, np.polyval(b_, xs[-1]), f"{lab}\nPearson r\n= {np.corrcoef(pp['mid'], pp[col])[0, 1]:.2f}", color=c, fontsize=5.5, va="center")
a.axhline(0, color=CH, lw=0.5); a.axvline(0, color=CH, lw=0.5)
a.set_xlabel("parents' mean, NGS-DOSE 45S\n(minus population mean)"); a.set_ylabel("child's 45S (minus population mean)")
a.set_title(f"Child vs parents,\n{inh[0]['n']} trios", loc="left")
a.set_xlim(pp["mid"].min() - 15, pp["mid"].max() + 95)
letter(a, "c", -0.34, 1.2)
# d: inheritance r with CIs
a = fig.add_subplot(gs[1, 1])
r45, r5 = inh[0], inh[2]; i18 = inh18["NGS-DOSE"]
rows_c = [("45S  NGS-DOSE", r45["r_ngsdose"], r45["r_ngsdose_lo"], r45["r_ngsdose_hi"], CN, "s"),
          ("18S ratio", i18["r"], i18["lo"], i18["hi"], C18, "o"),
          ("assembly", r45["r_assembly"], r45["r_assembly_lo"], r45["r_assembly_hi"], CA, "o"),
          ("5S  NGS-DOSE", r5["r_ngsdose"], r5["r_ngsdose_lo"], r5["r_ngsdose_hi"], CN, "s"),
          ("assembly", r5["r_assembly"], r5["r_assembly_lo"], r5["r_assembly_hi"], CA, "o")]
ypos = [5, 4, 3, 1.5, 0.5]
for y, (lab, r_, lo_, hi_, c, m) in zip(ypos, rows_c):
    a.errorbar(r_, y, xerr=[[r_ - lo_], [hi_ - r_]], fmt=m, ms=4, color=c, lw=1)
a.set_yticks(ypos); a.set_yticklabels([r[0] for r in rows_c])
a.axhline(2.25, color=CH, lw=0.5, ls=":")
a.set_xlim(-0.2, 1.05); a.axvline(0, color=CH, lw=0.6); a.set_xlabel("Pearson r, child vs parents' mean\n(95% bootstrap CI)")
a.set_title("Inheritance: short reads\nwin for 45S; 5S tie", loc="left")
letter(a, "d", -0.72, 1.2)
# e: ranking people against ddPCR, r with 95% CI
a = fig.add_subplot(gs[1, 2])
rows_b = [(f"NGS-DOSE ({rc['ngsdose_any']['n']})", rc["ngsdose_any"], CN), (f"CONKORD ({rc['conkord']['n']})", rc["conkord"], CC),
          (f"18S ratio ({rc['ratio18S_any']['n']})", rc["ratio18S_any"], C18),
          (f"NGS-DOSE ({rc['ngsdose_same_asm']['n']})", rc["ngsdose_same_asm"], CN), (f"assembly ({rc['assembly']['n']})", rc["assembly"], CA)]
for y, (lab, d_, c) in enumerate(rows_b[::-1]):
    a.errorbar(d_["r"], y, xerr=[[d_["r"] - d_["lo"]], [d_["hi"] - d_["r"]]], fmt="o", ms=4, color=c, lw=1)
a.set_yticks(range(len(rows_b))); a.set_yticklabels([r[0] for r in rows_b[::-1]])
a.axhline(1.5, color=CH, lw=0.5, ls=":")
a.text(-0.97, 1.6, "all 12 ddPCR lines", fontsize=5.3, color=CH, va="bottom")
a.text(-0.97, 1.4, f"the {rc['assembly']['n']} with an HPRC assembly", fontsize=5.3, color=CH, va="top")
a.axvline(0, color=CH, lw=0.6); a.set_xlim(-1, 1.05); a.set_ylim(-0.5, len(rows_b) - 0.5)
a.set_xlabel("Pearson r with ddPCR\n(95% CI, Fisher z)")
a.set_title("Ranking people vs ddPCR:\nNGS-DOSE highest", loc="left")
letter(a, "e", -0.62, 1.2)
# f: 45S and 5S, assembly vs NGS-DOSE, same people
a = fig.add_subplot(gs[2, 0])
av = asm.join(co[["rDNA45S.cn", "rDNA5S.cn"]], how="inner").dropna(subset=["rDNA45S.cn", "rDNA5S.cn"])
a.plot([0, 750], [0, 750], color=CH, lw=0.8)
a.scatter(av["rDNA45S.cn"], av.a45, s=10, color=CA, lw=0, label=f"45S (Pearson r = {np.corrcoef(av['rDNA45S.cn'], av.a45)[0,1]:.2f})", zorder=3)
a.scatter(av["rDNA5S.cn"], av.a5, s=10, marker="D", facecolor="white", edgecolor=CA, linewidths=0.8,
          label=f"5S (Pearson r = {np.corrcoef(av['rDNA5S.cn'], av.a5)[0,1]:.3f})", zorder=3)
a.set_xlim(0, 720); a.set_ylim(0, 720)
a.set_xlabel("NGS-DOSE copies"); a.set_ylabel("assembled copies (both haplotypes)")
a.set_title(f"Assemblies hold 5S whole,\n45S by half ({len(av)} people)", loc="left")
a.legend(frameon=False, fontsize=5.5, loc="upper left", handletextpad=0.3)
letter(a, "f", -0.34, 1.2)
# g: known truth, distal junction
a = fig.add_subplot(gs[2, 1])
bins = np.arange(7.55, 11.6, 0.1)
dn_ = co["DJ.cn"].dropna()
a.hist(dn_.clip(7.6, 11.5), bins=bins, weights=np.ones(len(dn_)) / len(dn_), color=CN, alpha=0.85, label=f"NGS-DOSE ({len(dn_)})", zorder=2)
vc = asm.aDJ.value_counts(normalize=True).sort_index()
vin = vc[(vc.index >= 7.6) & (vc.index <= 11.5)]
a.bar(vin.index, vin.values, width=0.1, color=CA, alpha=0.9, label=f"assembly ({len(asm)})", zorder=3)
n_out = int(((asm.aDJ < 7.6) | (asm.aDJ > 11.5)).sum())
a.axvline(10, color=CH, lw=0.6, ls=":", zorder=1); a.set_xlabel("distal-junction copies (truth 10)"); a.set_ylabel("fraction of people")
a.set_xlim(7.5, 11.6)
a.set_title(f"Known truth: assemblies count\n10 exactly; NGS-DOSE {dj_stats['cohort_DJ_median']:.2f}", loc="left")
a.legend(frameon=False, fontsize=5.5, loc="upper left", handlelength=1.2)
if n_out:
    a.text(0.02, 0.8, f"{n_out} assembled people\nbelow 7.6 not shown", transform=a.transAxes, fontsize=5.3, color=CA, va="top")
letter(a, "g", -0.34, 1.2)
# h: same person, two pipelines
a = fig.add_subplot(gs[2, 2])
b45 = br[br.quantity == "rDNA45S.cn"].set_index("sample").pct; b18 = br[br.quantity == "rDNA45S.18S.flat"].set_index("sample").pct
xx = np.arange(len(ceph))
a.bar(xx - 0.18, b45.reindex(ceph), width=0.34, color=CN, label="NGS-DOSE")
a.bar(xx + 0.18, b18.reindex(ceph), width=0.34, color=C18, label="18S depth ratio")
a.axhline(0, color="k", lw=0.6)
a.set_xticks(xx); a.set_xticklabels(ceph, fontsize=6)
a.set_ylabel("45S change, Google vs NYGC (%)")
a.set_ylim(-6.5, 34); a.set_yticks([-5, 0, 5, 10, 15, 20, 25])
a.legend(frameon=False, fontsize=5.5, loc="upper left", ncol=2, columnspacing=0.8, handlelength=1.2)
a.set_title("Same person, two pipelines:\nNGS-DOSE stable", loc="left")
letter(a, "h", -0.36, 1.2)
# glossary spanning the bottom row
a = fig.add_subplot(gs[2, :]); a.axis("off")
a.set_position([0.03, 0.004, 0.95, 0.25])
gl = [("What each measure is", ""),
      ("ddPCR", "Droplet digital PCR on each cell line's DNA (Potapova et al. 2025). Counts rDNA copies without sequencing: the independent reference in a, b and e. Replicate CV about 5%."),
      ("NGS-DOSE", "This method: k-mer counts at rDNA windows, calibrated on the NYGC cohort. NYGC reads = 1000 Genomes 30× CRAMs; Google reads = a second NovaSeq pipeline (GIAB; HG002–HG004 and the CEPH trio)."),
      ("CONKORD", "Potapova et al.'s own short-read k-mer estimate, from their reads. Only available for their 12 lines."),
      ("18S depth ratio", "Read depth on the 18S gene ÷ autosomal depth: the estimator of the UK Biobank literature, computed here from the same reads as NGS-DOSE."),
      ("HPRC assembly", "18S genes found in a person's two release-2 hifiasm haplotype assemblies. T2T-HG002 v1.1 is shown open and excluded: 9 of its 10 rDNA arrays are N-gaps by design."),
      ("Pearson r, Spearman ρ", "r: Pearson's linear correlation coefficient, used for every correlation in this figure. ρ: Spearman's rank correlation, given in the overview table as a check that no single person drives r."),
      ("Inheritance", "Pearson r of the child's value with the mean of its parents' (parents by NGS-DOSE), both as deviations from the superpopulation mean. Error in the child's value lowers it; nothing can raise it."),
      ("FISH shares", "Fluorescence of each rDNA array as a fraction of the cell's total (Potapova et al.). Independent of sequencing, but a fraction, not a count."),
      ("FISH copies, totals", "Share × CONKORD total. Summed over a person's arrays they return CONKORD (within 2 copies), so FISH totals are not an independent measurement."),
      ("Per-chromosome units, shares", "Units: FISH copies on one chromosome pair vs assembled 18S genes on contigs assigned to it. Shares: each as a fraction of the person's total. "
                                       "Assembly vs FISH: Pearson r = %.2f (units), %.2f (shares), %d chromosomes in %d people (Figure 5e, f)." % (arr_stats["r_units"], arr_stats["r_share"], arr_stats["n_chrom"], len(arr_stats["samples"])))]
import textwrap
yy = 0.99
for k, v in gl:
    if not v:
        a.text(0.0, yy, k, fontsize=7, fontweight="bold", va="top", transform=a.transAxes); yy -= 0.07; continue
    lines = textwrap.wrap(v, 112)
    a.text(0.0, yy, k, fontsize=6.2, fontweight="bold", va="top", transform=a.transAxes)
    a.text(0.235, yy, "\n".join(lines), fontsize=6.2, va="top", transform=a.transAxes, linespacing=1.25)
    yy -= 0.03 + 0.042 * len(lines)
fig.savefig(f"{FIG}/fig7_at_a_glance.png", dpi=300)
plt.close(fig)

# companion table: one row per test, one column per method
SP = {c: float(stats.spearmanr(pc[c], pc.ddpcr, nan_policy="omit")[0]) for c in ["ngsdose_any", "conkord", "ratio18S_any", "assembly_18S"]}
dd_stats["spearman"] = SP
def f2(v, d=2):
    return "—" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{d}f}"
sc_rows = [
    dict(test="Level vs ddPCR (median estimate ÷ ddPCR)", truth="ddPCR", n=f"{dd_stats['ngsdose_any']['n']}; asm. {dd_stats['assembly_18S']['n']}",
         ngsdose=f2(dd_stats["ngsdose_any"]["median_ratio"]), conkord=f2(dd_stats["conkord"]["median_ratio"]),
         ratio18S=f2(dd_stats["ratio18S_any"]["median_ratio"]), assembly=f2(dd_stats["assembly_18S"]["median_ratio"]), best="CONKORD, then NGS-DOSE"),
    dict(test="Ranking vs ddPCR (Pearson r)", truth="ddPCR", n=f"{rc['ngsdose_any']['n']}; asm. {rc['assembly']['n']}",
         ngsdose=f2(rc["ngsdose_any"]["r"]), conkord=f2(rc["conkord"]["r"]), ratio18S=f2(rc["ratio18S_any"]["r"]),
         assembly=f2(rc["assembly"]["r"]), best="NGS-DOSE (CIs overlap)"),
    dict(test="Ranking vs ddPCR (Spearman ρ)", truth="ddPCR", n=f"{rc['ngsdose_any']['n']}; asm. {rc['assembly']['n']}",
         ngsdose=f2(SP["ngsdose_any"]), conkord=f2(SP["conkord"]), ratio18S=f2(SP["ratio18S_any"]), assembly=f2(SP["assembly_18S"]), best="NGS-DOSE ≈ CONKORD"),
    dict(test=f"Ranking vs ddPCR, the {rc['assembly']['n']} assembled people (Pearson r)", truth="ddPCR", n=str(rc['assembly']['n']),
         ngsdose=f2(rc["ngsdose_same_asm"]["r"]), conkord=f2(rc["conkord_same_asm"]["r"]), ratio18S="—", assembly=f2(rc["assembly"]["r"]), best="NGS-DOSE"),
    dict(test="Inheritance, 45S (Pearson r, child vs parents' mean)", truth="parents", n=str(inh[0]["n"]),
         ngsdose=f2(inh[0]["r_ngsdose"]), conkord="—", ratio18S=f2(inh18["NGS-DOSE"]["r"]), assembly=f2(inh[0]["r_assembly"]), best="NGS-DOSE = 18S ratio > assembly"),
    dict(test="Inheritance, 5S (Pearson r)", truth="parents", n=str(inh[2]["n"]), ngsdose=f2(inh[2]["r_ngsdose"]), conkord="—", ratio18S="—",
         assembly=f2(inh[2]["r_assembly"]), best="tie"),
    dict(test="Distal junction, known to be 10 (mean copies)", truth="10", n=f"{dj_stats['n']}", ngsdose=f2(dj_stats["ngsdose_mean"]), conkord="—",
         ratio18S="—", assembly=f2(dj_stats["assembly_mean"]), best="assembly (level); both see steps"),
    dict(test="Same person, second pipeline: mean |change| in 45S (%)", truth="0", n=str(len(ceph)), ngsdose=f2(bridge["rDNA45S.cn"]["mean_abs_pct"], 1),
         conkord="—", ratio18S=f2(bridge["rDNA45S.18S.flat"]["mean_abs_pct"], 1), assembly="—", best="NGS-DOSE"),
    dict(test="Per-chromosome share vs FISH (Pearson r)", truth="FISH shares", n=f"{arr_stats['n_chrom']} chrom.", ngsdose="totals only", conkord="—",
         ratio18S="—", assembly=f2(arr_stats["r_share"]), best="neither"),
]
pd.DataFrame(sc_rows).to_csv(f"{TAB}/scorecard.tsv", sep="\t", index=False)
json.dump(dict(out, scorecard=sc_rows), open(f"{TAB}/method_accuracy.json", "w"), indent=1, default=float)
