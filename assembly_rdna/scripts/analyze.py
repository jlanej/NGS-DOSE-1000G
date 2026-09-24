#!/usr/bin/env python
"""What HPRC release-2 assemblies hold of the rDNA, and whether it tracks short-read copy number.

Inputs: work_units/*.summary.json and *.arrays.tsv (annotate_units.py), tables/assemblies.tsv, tables/acro_contigs.tsv, and from the results
repository docs/data/cohort.tsv (NGS-DOSE per-sample estimates) and the 1000 Genomes pedigree.
Everything that joins to short reads is recomputed from whatever cohort.tsv holds, so rerunning
after `regenerate.sh` has added counts updates every test.

Writes tables/haplotype_rdna.tsv, tables/person_vs_ngsdose.tsv, tables/trio_haplotypes.tsv,
tables/arrays_all.tsv.gz, tables/tests.json and figures/fig{1,2,3}_*.png
"""
import glob, json, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(BASE)
TAB, FIG = os.path.join(BASE, "tables"), os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)
rng = np.random.default_rng(20260924)
UNIT_BP = 44838


def boot_r(x, y, B=5000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 5:
        return dict(n=n, r=None, lo=None, hi=None, rho=None, p=None)
    bs = [np.corrcoef(x[i], y[i])[0, 1] for i in (rng.integers(0, n, n) for _ in range(B))]
    return dict(n=n, r=float(np.corrcoef(x, y)[0, 1]), lo=float(np.nanpercentile(bs, 2.5)), hi=float(np.nanpercentile(bs, 97.5)),
                rho=float(stats.spearmanr(x, y)[0]), p=float(stats.pearsonr(x, y)[1]))


# ---------------- per haplotype ----------------
summ = pd.DataFrame([json.load(open(f)) for f in glob.glob(f"{BASE}/work_units/*.summary.json")])
asm = pd.read_csv(f"{TAB}/assemblies.tsv", sep="\t")
summ = summ.merge(asm[["assembly_name", "sample", "haplotype", "method", "version", "phasing", "source", "censat_rdna_bp"]],
                  left_on="assembly", right_on="assembly_name", how="left").drop(columns="assembly_name").rename(columns={"sample": "person"})
summ["group"] = np.select([summ.method.eq("verkko"), summ.phasing.eq("hic"), summ.version.eq("0.19.7")],
                          ["Verkko", "hifiasm Hi-C (0.19.8-9)", "hifiasm trio (0.19.7)"], "hifiasm trio (0.19.9)")
summ.loc[~summ.source.isin(["hprc", "hpp"]), "group"] = "reference"
summ.to_csv(f"{TAB}/haplotype_rdna.tsv", sep="\t", index=False)
h = summ[summ.source.isin(["hprc", "hpp"])].copy()
hh = h[h.method == "hifiasm"]

arr = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{BASE}/work_units/*.arrays.tsv") if os.path.getsize(f) > 50], ignore_index=True)
acro = pd.read_csv(f"{TAB}/acro_contigs.tsv", sep="\t")[["assembly_name", "contig", "chrom", "level"]]
arr = arr.merge(acro, left_on=["assembly", "contig"], right_on=["assembly_name", "contig"], how="left").drop(columns="assembly_name")
arr["place"] = arr.level.fillna("unplaced").map({"chromosome": "acrocentric, chromosome-length", "random": "acrocentric, unlocalised piece",
                                                 "unplaced": "unplaced contig"})
arr["closed"] = arr.left.str.startswith("flank") & arr.right.str.startswith("flank")
arr["DJ_side"] = arr.left.str.contains("DJ") | arr.right.str.contains("DJ")
arr = arr.merge(summ[["assembly", "source", "method", "group"]], on="assembly")
arr.to_csv(f"{TAB}/arrays_all.tsv.gz", sep="\t", index=False, compression="gzip")
ah = arr[arr.source.isin(["hprc", "hpp"]) & (arr.method == "hifiasm")]
sides = pd.concat([ah.left.str.split("+").str[0], ah.right.str.split("+").str[0]]).value_counts()

# ---------------- per person, with short reads ----------------
co = pd.read_csv(f"{REPO}/docs/data/cohort.tsv", sep="\t").set_index("sample")
ped = pd.read_csv(f"{REPO}/meta/20130606_g1k_3202_samples_ped_population.txt", sep=" ").set_index("SampleID")
dip = hh.groupby("person").agg(n_hap=("assembly", "size"), a18S=("n_18S", "sum"), a28S=("n_28S", "sum"), a5S=("n_5S_units", "sum"),
                               aDJ=("DJ_copy_eq", "sum"), censat_bp=("censat_rdna_bp", "sum"), phasing=("phasing", "first"),
                               version=("version", "first"))
dip = dip[dip.n_hap == 2].join(co[["rDNA45S.cn", "rDNA45S.18S.flat", "rDNA5S.cn", "DJ.cn"]], how="inner")
dip["ratio45"] = dip.a18S / dip["rDNA45S.cn"]
dip["ratio5"] = dip.a5S / dip["rDNA5S.cn"]
dip.to_csv(f"{TAB}/person_vs_ngsdose.tsv", sep="\t")

t = hh[hh.phasing == "trio"].pivot_table(index="person", columns="haplotype", values=["n_18S", "n_5S_units"], aggfunc="first")
t.columns = [f"{a}_{b}" for a, b in t.columns]
t = t.join(ped[["FatherID", "MotherID"]], how="left")
t = t[t.FatherID.isin(co.index) & t.MotherID.isin(co.index) & t.index.isin(co.index)].copy()
for v in ["rDNA45S.cn", "rDNA5S.cn"]:
    t["father_" + v] = co.loc[t.FatherID, v].values
    t["mother_" + v] = co.loc[t.MotherID, v].values
    t["child_" + v] = co.loc[t.index, v].values
t.to_csv(f"{TAB}/trio_haplotypes.tsv", sep="\t")

# ceiling for a transmitted haplotype against the transmitting parent's diploid total, perfect measurement
sim = []
for _ in range(2000):
    X = rng.lognormal(np.log(46), 0.5, (max(len(t), 5), 5, 2))
    mom = X.sum((1, 2))
    hap = X[np.arange(len(X))[:, None], np.arange(5)[None, :], rng.integers(0, 2, (len(X), 5))].sum(1)
    sim.append(np.corrcoef(hap, mom)[0, 1])

tests = dict(
    n_haplotypes=int(len(h)), n_hifiasm=int(len(hh)), n_people_counted=int(len(dip)), n_trios=int(len(t)),
    per_haplotype={c: hh[c].describe(percentiles=[.1, .5, .9]).round(2).to_dict() for c in
                   ["n_18S", "n_28S", "n_arrays", "largest_array", "n_closed", "units_in_closed", "n_5S_units", "n_DJ_copies",
                    "n_periods", "n_periods_identical"]},
    by_group=h.groupby("group")[["n_18S", "n_arrays", "n_closed", "n_5S_units", "n_DJ_copies"]].median().to_dict(),
    frac_hap_with_18S=float((hh.n_18S > 0).mean()), frac_hap_5S_single_closed=float(((hh.n_5S_arrays == 1) & hh.largest_5S_left.eq("flank")
                                                                                   & hh.largest_5S_right.eq("flank")).mean()),
    array_sides={k: int(v) for k, v in sides.items()},
    units_by_place=ah.groupby("place").n_18S.sum().astype(int).to_dict(),
    arrays_n=int(len(ah)), arrays_closed=int(ah.closed.sum()), arrays_closed_DJ=int((ah.closed & ah.DJ_side).sum()),
    closed_units=ah[ah.closed].n_18S.describe().round(1).to_dict(),
    periods_identical=float(ah.n_identical.sum() / ah.n_periods.sum()),
    references=summ[summ.group == "reference"].set_index("assembly")[["n_18S", "n_28S", "n_arrays", "n_closed", "units_in_closed",
                                                                         "n_5S_units", "n_DJ_copies", "n_periods", "n_periods_identical"]].to_dict("index"),
    r_45S=boot_r(dip.a18S, dip["rDNA45S.cn"]), r_28S=boot_r(dip.a28S, dip["rDNA45S.cn"]), r_45S_vs_18Sratio=boot_r(dip.a18S, dip["rDNA45S.18S.flat"]),
    r_censat=boot_r(dip.censat_bp, dip["rDNA45S.cn"]), r_5S=boot_r(dip.a5S, dip["rDNA5S.cn"]), r_DJ=boot_r(dip.aDJ, dip["DJ.cn"]),
    ratio45=dip.ratio45.describe().round(3).to_dict(), ratio5=dip.ratio5.describe().round(4).to_dict(),
    ratio45_by_phasing=dip.groupby("phasing").ratio45.median().round(3).to_dict(),
    ratio45_vs_cn=float(stats.spearmanr(dip.ratio45, dip["rDNA45S.cn"])[0]) if len(dip) > 4 else None,
    hap_parent={f"{cls}_{hap}~{par}": boot_r(t[f"{col}_{hap}"], t[f"{par}_{v}"])
                for cls, col, v in [("45S", "n_18S", "rDNA45S.cn"), ("5S", "n_5S_units", "rDNA5S.cn")]
                for hap in ["mat", "pat"] for par in ["mother", "father"]},
    hap_parent_ceiling=dict(median=float(np.median(sim)), lo=float(np.percentile(sim, 2.5)), hi=float(np.percentile(sim, 97.5))),
)
json.dump(tests, open(f"{TAB}/tests.json", "w"), indent=1, default=float)

# ---------------- figures ----------------
try:
    apply_figure_style(sizes=(8, 7, 6))
except NameError:
    plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
GC = {"hifiasm trio (0.19.7)": "#1f78b4", "hifiasm trio (0.19.9)": "#a6cee3", "hifiasm Hi-C (0.19.8-9)": "#33a02c", "Verkko": "#b15928"}
ORDER = list(GC)
def letter(ax, s, x=-0.2, y=1.12):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")

# Fig 1: what the assemblies hold
fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.7), gridspec_kw=dict(width_ratios=[1.5, 1.1, 0.8, 1.0]))
a = ax[0]
for j, g in enumerate(ORDER):
    y = h[h.group == g].n_18S.values
    a.scatter(j + rng.uniform(-0.3, 0.3, len(y)), y, s=4, color=GC[g], alpha=0.6, lw=0)
    a.hlines(np.median(y), j - 0.38, j + 0.38, color="k", lw=1.3)
ref = summ.set_index("assembly")
exp_hap = float(np.median(dip["rDNA45S.cn"]) / 2) if len(dip) else None
if exp_hap:
    a.axhline(exp_hap, color="#555555", ls="--", lw=0.8)
    a.text(3.45, exp_hap + 6, f"short reads, median\nper haplotype ({exp_hap:.0f})", fontsize=6, ha="right", va="bottom", color="#555555")
a.axhline(ref.loc["chm13v2.0_maskedY_rCRS", "n_18S"], color="#999999", ls=":", lw=0.8)
a.text(-0.45, ref.loc["chm13v2.0_maskedY_rCRS", "n_18S"] + 5, "CHM13 (model arrays)", fontsize=6, color="#777777")
a.set_xticks(range(4)); a.set_xticklabels(["trio\n0.19.7", "trio\n0.19.9", "Hi-C", "Verkko"]); a.set_xlabel("hifiasm phasing / version")
a.set_ylabel("18S gene copies per haplotype"); a.set_title("Every haplotype holds rDNA", loc="left"); letter(a, "a", -0.3)
a = ax[1]
cnt = ah.n_18S.clip(upper=20).value_counts().sort_index()
a.bar(cnt.index, cnt.values, color="#1f78b4", width=0.8)
a.set_yscale("log"); a.set_xlabel("units per array (20 = 20 or more)"); a.set_ylabel("arrays (hifiasm)")
a.set_title("Arrays are 1-10 units", loc="left"); letter(a, "b")
a = ax[2]
lab = {"contig_end": "contig\nend", "flank": "flank", "gap": "gap"}
vals = [sides.get(k, 0) / sides.sum() for k in lab]
a.bar(range(3), vals, color=["#e31a1c", "#1f78b4", "#999999"])
for i, v in enumerate(vals):
    a.text(i, v + 0.02, f"{100*v:.0f}%", ha="center", fontsize=6.5)
a.set_xticks(range(3)); a.set_xticklabels(list(lab.values())); a.set_ylim(0, 1); a.set_ylabel("share of array ends")
a.set_title("Array ends", loc="left"); letter(a, "c")
a = ax[3]
pl = ah.groupby("place").n_18S.sum()
pl = pl.reindex(["acrocentric, chromosome-length", "acrocentric, unlocalised piece", "unplaced contig"]).fillna(0)
a.barh(range(3), pl.values / pl.sum(), color="#6a3d9a")
for i, v in enumerate(pl.values / pl.sum()):
    a.text(v + 0.02, i, f"{100*v:.0f}%", va="center", fontsize=6.5)
a.set_yticks(range(3)); a.set_yticklabels(["chromosome-\nlength acro.", "acrocentric\npiece", "unplaced\ncontig"]); a.set_xlim(0, 1)
a.set_xlabel("share of units"); a.set_title("Where units sit", loc="left"); letter(a, "d", -0.55)
fig.tight_layout(w_pad=0.8)
fig.savefig(f"{FIG}/fig1_what_assemblies_hold.png", dpi=300)
plt.close(fig)

# Fig 2: against short reads
fig, ax = plt.subplots(1, 3, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1, 1, 1.3]))
cc = dip.phasing.map({"trio": "#1f78b4", "hic": "#33a02c"})
a = ax[0]
a.scatter(dip["rDNA45S.cn"], dip.a18S, s=10, c=cc, lw=0)
m = max(dip["rDNA45S.cn"].max(), dip.a18S.max()) * 1.05
a.plot([0, m], [0, m], color="#bbbbbb", lw=0.8); a.plot([0, m], [0, m * dip.ratio45.median()], color="#555555", lw=0.8, ls="--")
a.set_xlim(0, m); a.set_ylim(0, m)
a.set_xlabel("NGS-DOSE 45S copies"); a.set_ylabel("assembled 18S copies (diploid)")
r = tests["r_45S"]; a.set_title(f"45S: {100*dip.ratio45.median():.0f}% assembled\nr = {r['r']:.2f}", loc="left"); letter(a, "a", -0.3, 1.2)
a.scatter([], [], s=10, c="#1f78b4", label="trio-phased"); a.scatter([], [], s=10, c="#33a02c", label="Hi-C-phased")
a.legend(frameon=False, loc="upper left", fontsize=6)
a = ax[1]
a.scatter(dip["rDNA5S.cn"], dip.a5S, s=10, c=cc, lw=0)
m5 = max(dip["rDNA5S.cn"].max(), dip.a5S.max()) * 1.05
a.plot([0, m5], [0, m5], color="#bbbbbb", lw=0.8); a.set_xlim(0, m5); a.set_ylim(0, m5)
a.set_xlabel("NGS-DOSE 5S copies"); a.set_ylabel("assembled 5S units (diploid)")
r5 = tests["r_5S"]; a.set_title(f"5S: ratio {dip.ratio5.median():.2f}\nr = {r5['r']:.3f}", loc="left"); letter(a, "b", -0.3, 1.2)
a = ax[2]
hp = tests["hap_parent"]; ceil = tests["hap_parent_ceiling"]
keys = [("45S", "mat", "mother"), ("45S", "mat", "father"), ("45S", "pat", "father"), ("45S", "pat", "mother"),
        ("5S", "mat", "mother"), ("5S", "mat", "father"), ("5S", "pat", "father"), ("5S", "pat", "mother")]
ys = np.arange(len(keys))[::-1]
a.axvspan(ceil["lo"], ceil["hi"], color="#eeeeee", zorder=0); a.axvline(ceil["median"], color="#bbbbbb", lw=0.8, zorder=0)
a.axvline(0, color="#999999", lw=0.6)
for y, (c, hap, par) in zip(ys, keys):
    v = hp[f"{c}_{hap}~{par}"]
    own = (hap == "mat") == (par == "mother")
    a.errorbar(v["r"], y, xerr=[[v["r"] - v["lo"]], [v["hi"] - v["r"]]], fmt="o", ms=3.5, lw=1,
               color=("#1f78b4" if c == "45S" else "#b15928"), mfc=("#1f78b4" if c == "45S" else "#b15928") if own else "white")
a.set_yticks(ys); a.set_yticklabels([f"{c} {'maternal' if h_=='mat' else 'paternal'} hap. vs {p}" for c, h_, p in keys])
a.set_xlim(-0.65, 1.0); a.set_xlabel("Pearson r (95% bootstrap CI)")
a.set_title(f"Each haplotype tracks its\nown parent ({tests['n_trios']} trios)", loc="left"); letter(a, "c", -0.85, 1.2)
a.set_ylim(-0.6, len(keys) - 0.1)
a.text(ceil["median"], len(keys) - 0.45, "perfect measure", fontsize=5.5, ha="center", va="center", color="#777777")
fig.tight_layout(w_pad=0.6)
fig.savefig(f"{FIG}/fig2_against_short_reads.png", dpi=300)
plt.close(fig)

# Fig 3: unit-level and control checks
periods = []
for f in glob.glob(f"{BASE}/work_units/*.arrays.tsv"):
    pass
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
a = ax[0]
grp = {"HPRC hifiasm": hh, "CHM13 (model arrays)": summ[summ.assembly.str.startswith("chm13")],
       "HG002 v1.1 (curated)": summ[summ.assembly.str.startswith("hg002")]}
vals = [(k, v.n_periods_identical.sum() / max(v.n_periods.sum(), 1)) for k, v in grp.items()]
a.bar(range(3), [v for _, v in vals], color=["#1f78b4", "#999999", "#555555"])
for i, (_, v) in enumerate(vals):
    a.text(i, v + 0.02, f"{100*v:.1f}%", ha="center", fontsize=6.5)
a.set_xticks(range(3)); a.set_xticklabels([k.replace(" (", "\n(") for k, _ in vals]); a.set_ylim(0, 1)
a.set_ylabel("unit periods identical to another"); a.set_title("Units are not copies of a template", loc="left"); letter(a, "a", -0.3)
a = ax[1]
pm = ah.period_median.dropna()
a.hist(pm / 1000, bins=np.arange(38, 52, 0.5), color="#1f78b4")
c13 = arr[arr.assembly.str.startswith("chm13")].period_median
for v in c13:
    a.axvline(v / 1000, color="#999999", lw=0.8, ls=":")
a.set_xlabel("unit period, array median (kb)"); a.set_ylabel("arrays (hifiasm)")
a.text(47.2, a.get_ylim()[1] * 0.9, "dotted: the five\nCHM13 arrays", fontsize=6, color="#777777")
a.set_title("Unit lengths are realistic", loc="left"); letter(a, "b")
a = ax[2]
dj = hh.n_DJ_copies.value_counts().sort_index()
a.bar(dj.index, dj.values / dj.sum(), color="#6a3d9a")
a.axvline(5, color="#555555", ls="--", lw=0.8)
a.set_xlabel("DJ copies per haplotype"); a.set_ylabel("share of haplotypes")
a.set_title("Distal junctions: 5 expected", loc="left"); letter(a, "c")
fig.tight_layout(w_pad=0.8)
fig.savefig(f"{FIG}/fig3_unit_checks.png", dpi=300)
plt.close(fig)
print(json.dumps({k: tests[k] for k in ["n_haplotypes", "n_people_counted", "n_trios", "r_45S", "r_5S"]}, default=float))
