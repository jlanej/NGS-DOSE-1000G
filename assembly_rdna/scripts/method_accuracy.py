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
pc = s1.copy()
pc["ngsdose"] = ngs.reindex(pc.index)
pc["ngsdose_source"] = src.reindex(pc.index)
pc["ratio18S_flat"] = flat.reindex(pc.index)
pc["assembly_18S"] = asm.a45.reindex(pc.index)
hg002 = hap[hap.assembly.str.startswith("hg002")]
pc.loc["HG002", "assembly_18S"] = hg002.n_18S.sum()
pc["hall"] = (hall["HC.18S.CN"] * 2).reindex(pc.index)
pc.to_csv(f"{TAB}/potapova_comparison.tsv", sep="\t")

def agree(col, ref="ddpcr"):
    x = pc.dropna(subset=[col, ref])
    lr = np.log(x[col] / x[ref])
    return dict(n=len(x), median_ratio=float(np.exp(np.median(lr))), sd_log=float(lr.std(ddof=1)) if len(x) > 1 else None,
                mean_abs_pct=float(100 * np.mean(np.abs(x[col] / x[ref] - 1))),
                r=float(np.corrcoef(x[col], x[ref])[0, 1]) if len(x) > 2 else None, samples=list(x.index))
dd_stats = {c: agree(c) for c in ["ngsdose", "conkord", "ratio18S_flat", "assembly_18S", "hall"]}
dd_stats["ngsdose_vs_conkord"] = agree("ngsdose", "conkord")
_x = pc.dropna(subset=["ngsdose"])
dd_stats["conkord_same9"] = dict(n=len(_x), r=float(np.corrcoef(_x.conkord, _x.ddpcr)[0, 1]),
                                 mean_abs_pct=float(100 * np.mean(abs(_x.conkord / _x.ddpcr - 1))), median_ratio=float(np.median(_x.conkord / _x.ddpcr)))
dd_stats["residual_sd_log"] = {c: float(np.log(_x[c] / _x.ddpcr).std(ddof=1)) for c in ["ngsdose", "conkord", "ratio18S_flat"]}
dd_stats["bias_pct"] = {c: float(100 * (np.exp(np.log(_x[c] / _x.ddpcr).mean()) - 1)) for c in ["ngsdose", "conkord", "ratio18S_flat"]}
_y = _x.drop(index="HG02053", errors="ignore")
dd_stats["without_HG02053"] = {c: dict(r=float(np.corrcoef(_y[c], _y.ddpcr)[0, 1]), mean_abs_pct=float(100 * np.mean(abs(_y[c] / _y.ddpcr - 1))))
                               for c in ["ngsdose", "conkord", "ratio18S_flat"]}
# ddPCR's own replicate spread, for scale
dd_stats["ddpcr_cv_median"] = float((s1.ddpcr_sd / s1.ddpcr).median())

# ---------------- 4. per array: FISH shares vs assembled units per chromosome ----------------
s2 = pd.read_csv(f"{POT}/potapova2025_tableS2.tsv", sep="\t")
s2["chrom"] = "chr" + s2.array.str[:2]
fish = s2.groupby(["sample", "chrom"]).fish_units.sum().rename("fish").reset_index()
fish["fish_share"] = fish.fish / fish.groupby("sample").fish.transform("sum")
arr = pd.read_csv(f"{TAB}/arrays_all.tsv.gz", sep="\t")
arr["person"] = arr.assembly.map(hap.set_index("assembly").person)
arr.loc[arr.assembly.str.startswith("hg002"), "person"] = "HG002"
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

out = dict(inheritance=[{k: v for k, v in r.items() if k != "points"} for r in inh],
           assembly_hall_parents=dict(n=len(d), r=r_big, lo=float(np.percentile(bsb, 2.5)), hi=float(np.percentile(bsb, 97.5))),
           distal_junction=dj_stats, ddpcr=dd_stats, arrays=arr_stats)
json.dump(out, open(f"{TAB}/method_accuracy.json", "w"), indent=1, default=float)

# ---------------- figures ----------------
try:
    apply_figure_style(sizes=(8, 7, 6))
except NameError:
    plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
CA, CN, CC, CH = "#b15928", "#1f78b4", "#6a3d9a", "#999999"
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
    a.set_title(f"Child measured by {lab}\nr = {rv:.2f} ({r45['n']} trios)", loc="left")
    a.set_xlabel("parents' mean, NGS-DOSE 45S"); a.set_ylabel("child's 45S copies" if j == 0 else "")
    a.set_ylim(*lim); letter(a, "ab"[j], -0.25, 1.2)
a = ax[2]
ys = np.arange(len(inh))[::-1]
for y, r in zip(ys, inh):
    a.errorbar(r["r_assembly"], y + 0.13, xerr=[[r["r_assembly"] - r["r_assembly_lo"]], [r["r_assembly_hi"] - r["r_assembly"]]], fmt="o", ms=3.5, color=CA, lw=1)
    a.errorbar(r["r_ngsdose"], y - 0.13, xerr=[[r["r_ngsdose"] - r["r_ngsdose_lo"]], [r["r_ngsdose_hi"] - r["r_ngsdose"]]], fmt="s", ms=3.2, color=CN, lw=1)
a.set_yticks(ys); a.set_yticklabels([f"{r['cls']}, parents:\n{r['parents']} ({r['n']})" for r in inh])
a.axvline(0, color=CH, lw=0.6); a.set_xlim(-0.4, 1.05); a.set_ylim(-0.6, len(inh) - 0.2)
a.set_xlabel("child vs parents' mean, r")
a.plot([], [], "o", color=CA, label="child by assembly"); a.plot([], [], "s", color=CN, label="child by NGS-DOSE")
a.legend(frameon=False, loc="upper left", bbox_to_anchor=(-0.02, 1.02), fontsize=6)
a.set_title("45S: NGS-DOSE wins; 5S: tie", loc="left"); letter(a, "c", -0.55, 1.2)
fig.tight_layout(w_pad=2.2)
fig.savefig(f"{FIG}/fig4_method_accuracy.png", dpi=300)
plt.close(fig)

# Fig 5: known truths and orthogonal assays
fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
a = ax[0, 0]
jit = rng.uniform(-0.12, 0.12, len(dj))
col = np.where(dj.ngs_step != 0, CC, CN)
a.scatter(dj["DJ.cn"], dj.aDJ + jit, s=12, c=col, lw=0)
a.axhline(10, color=CH, lw=0.6, ls="--"); a.axvline(10, color=CH, lw=0.6, ls="--")
for k in (-1, 1):
    a.axvline(dj_stats["cohort_DJ_median"] + k, color=CH, lw=0.5, ls=":")
a.set_xlabel("NGS-DOSE distal-junction copies"); a.set_ylabel("assembly distal-junction copies")
a.set_yticks([8, 9, 10, 11])
a.set_title(f"DJ (truth 10): all {dj_stats['ngs_steps']} NGS-DOSE steps recur in\nthe assemblies; NGS-DOSE reads {dj_stats['ngsdose_mean']:.2f}, not 10", loc="left")
a.text(0.02, 0.95, f"purple: NGS-DOSE step\n(dotted: cohort level ±1)", transform=a.transAxes, fontsize=5.8, va="top", color="#555555")
letter(a, "a")
a = ax[0, 1]
x = pc.dropna(subset=["ngsdose"])
lo, hi = 420, 740
a.plot([lo, hi], [lo, hi], color=CH, lw=0.8)
a.scatter(x.ddpcr, x.ratio18S_flat, marker="o", s=14, facecolor="white", edgecolor="#e31a1c", label="18S depth ratio", zorder=3)
a.scatter(x.ddpcr, x.conkord, marker="^", s=16, color=CC, label="CONKORD (k-mer)", zorder=3)
a.errorbar(x.ddpcr, x.ngsdose, xerr=x.ddpcr_sd, fmt="s", ms=4, color=CN, lw=0.6, label="NGS-DOSE", zorder=4)
a.annotate("HG02053", (713, x.loc["HG02053", "ngsdose"]), xytext=(-40, -12), textcoords="offset points", fontsize=5.5)
a.set_xlim(lo, hi); a.set_ylim(lo, hi)
a.set_xlabel("ddPCR, Potapova et al. 2025 (±SD)"); a.set_ylabel("short-read 45S copies")
st, sc, sf = dd_stats["ngsdose"], dd_stats["conkord_same9"], dd_stats["ratio18S_flat"]
a.set_title(f"ddPCR (n={st['n']}): r = {st['r']:.2f} NGS-DOSE, {sc['r']:.2f} CONKORD,\n"
            f"{sf['r']:.2f} 18S ratio; NGS-DOSE {100*(st['median_ratio']-1):+.0f}%", loc="left")
a.legend(frameon=False, loc="upper left", fontsize=5.8); letter(a, "b")
a = ax[1, 0]
xa = pc.dropna(subset=["assembly_18S"])
a.plot([0, hi], [0, hi], color=CH, lw=0.8)
a.errorbar(xa.ddpcr, xa.assembly_18S, xerr=xa.ddpcr_sd, fmt="o", ms=4, color=CA, lw=0.6, label="HPRC r2 assembly")
xn = pc.loc[xa.index].dropna(subset=["ngsdose"])
a.scatter(xn.ddpcr, xn.ngsdose, marker="s", s=14, color=CN, label="NGS-DOSE, same people", zorder=3)
a.annotate("HG002 (T2T v1.1)", (xa.loc["HG002", "ddpcr"], xa.loc["HG002", "assembly_18S"]), xytext=(-20, 8), textcoords="offset points", fontsize=5.5)
a.set_xlim(0, hi); a.set_ylim(0, hi)
sa = dd_stats["assembly_18S"]
a.set_xlabel("ddPCR, Potapova et al. 2025"); a.set_ylabel("45S copies (diploid)")
a.set_title(f"Assemblies hold {100*sa['median_ratio']:.0f}% of the ddPCR copies\n(n={sa['n']}, incl. curated HG002)", loc="left")
a.legend(frameon=False, loc="upper left", fontsize=5.8); letter(a, "c")
a = ax[1, 1]
for s_, g in pm.groupby("sample"):
    a.scatter(g.fish_share, g.asm_share_of_placed, s=14, label=s_)
a.plot([0, 0.5], [0, 0.5], color=CH, lw=0.8); a.set_xlim(0, 0.5); a.set_ylim(0, 0.55)
a.set_xlabel("FISH: chromosome's share of the rDNA"); a.set_ylabel("assembly: share of placed units")
a.set_title(f"Per chromosome, assembly vs FISH:\nr = {arr_stats['r_share']:.2f} ({arr_stats['n_chrom']} chromosomes, 5 people)", loc="left")
a.legend(frameon=False, fontsize=5.5, loc="upper right", ncol=2); letter(a, "d")
fig.tight_layout(h_pad=1.5, w_pad=1.0)
fig.savefig(f"{FIG}/fig5_truths_and_assays.png", dpi=300)
plt.close(fig)
print(json.dumps({k: out[k] for k in ["ddpcr", "arrays"]}, default=float, indent=0)[:3000])
