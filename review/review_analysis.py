"""Independent re-analysis of the NGS-DOSE 1000 Genomes results (as of 2026-09-23, 735 genomes, 149 trios).

Reads only files committed in this repository (docs/data/cohort.tsv, docs/report.json, meta/*), writes tables and
figures to review/.  Run from the repository root:  python review/review_analysis.py
"""
import json, os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import matplotlib
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "review")
rng = np.random.default_rng(20260923)

co = pd.read_csv(f"{ROOT}/docs/data/cohort.tsv", sep="\t")
ped = pd.read_csv(f"{ROOT}/meta/20130606_g1k_3202_samples_ped_population.txt", sep=" ")
qc = pd.read_csv(f"{ROOT}/meta/ngspca_sample_qc.tsv", sep="\t")
S = (co.merge(ped[["SampleID", "FatherID", "MotherID", "Sex", "Population", "Superpopulation"]],
              left_on="sample", right_on="SampleID", how="left")
       .merge(qc[["SAMPLE_ID", "RELEASE_BATCH"]], left_on="sample", right_on="SAMPLE_ID", how="left")
       .set_index("sample").copy())
parents = set(ped.FatherID) | set(ped.MotherID)
S["gen"] = np.where(S.FatherID.astype(str) != "0", "child", np.where(S.index.isin(parents), "parent", "unrel"))
S["batch698"] = (S.RELEASE_BATCH == 698).astype(int)
S["is_child"] = (S.gen == "child").astype(int)
tr = S[(S.gen == "child") & S.FatherID.isin(S.index) & S.MotherID.isin(S.index)]

# ---------- 1. batch versus generation ----------
feat = {
    "28S / 18S (GC model)": np.log(S["rDNA45S.28S"] / S["rDNA45S.18S"]),
    "28S / 18S (no GC model)": np.log(S["rDNA45S.28S.flat"] / S["rDNA45S.18S.flat"]),
    "5'ETS / IGS (GC model)": np.log(S["rDNA45S.5ETS"] / S["rDNA45S.IGS"]),
    "18S depth ratio / calibrated 45S": np.log(S["rDNA45S.18S.flat"] / S["rDNA45S.cn"]),
    "calibrated 45S copy number": np.log(S["rDNA45S.cn"]),
    "18S depth ratio (published)": np.log(S["rDNA45S.18S.flat"]),
    "insert-size median": np.log(S["insert_median"]),
    "library rate at 65% GC": np.log(S["gc_rel_65"]),
    "held-out autosomal (truth 2)": np.log(S["truth.auto"]),
    "distal junction (truth 10)": np.log(S["DJ.cn"]),
}
rows = []
for k, v in feat.items():
    S["_y"] = v
    f = smf.ols("_y ~ batch698 + is_child + C(Population)", data=S).fit()
    rows.append(dict(quantity=k, n=int(f.nobs),
                     batch_effect_pct=100 * (np.exp(f.params["batch698"]) - 1), batch_se_pct=100 * f.bse["batch698"],
                     batch_p=f.pvalues["batch698"],
                     generation_effect_pct=100 * (np.exp(f.params["is_child"]) - 1), generation_se_pct=100 * f.bse["is_child"],
                     generation_p=f.pvalues["is_child"]))
batch_tab = pd.DataFrame(rows)
batch_tab.to_csv(f"{OUT}/batch_vs_generation.tsv", sep="\t", index=False, float_format="%.4g")

# ---------- 2. transmission reliabilities, re-derived ----------
def centred(m):
    x = S[m].astype(float)
    return x - x.groupby(S.Population).transform("mean")

def _est(c, f, mo):
    mid = (f + mo) / 2
    b = np.cov(c, mid)[0, 1] / np.var(mid, ddof=1)
    rho = np.corrcoef(f, mo)[0, 1]
    s = np.std(c, ddof=1) / np.std(np.r_[f, mo], ddof=1)
    return b - rho * (1 - b), (b / s) - rho * (1 - b / s), rho, s

def trio_stats(m, B=4000):
    x = centred(m)
    c, f, mo = x.loc[tr.index].values, x.loc[tr.FatherID].values, x.loc[tr.MotherID].values
    ok = ~(np.isnan(c) | np.isnan(f) | np.isnan(mo)); c, f, mo = c[ok], f[ok], mo[ok]
    e = _est(c, f, mo)
    bs = np.array([_est(c[i], f[i], mo[i]) for i in (rng.integers(0, len(c), len(c)) for _ in range(B))])
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    return dict(metric=m, n_trios=len(c), R=e[0], R_lo=lo[0], R_hi=hi[0], R_rescaled=e[1], R_rescaled_lo=lo[1],
                R_rescaled_hi=hi[1], spousal_r=e[2], sd_child_over_parent=e[3])

LAB = {"rDNA45S.cn": "45S, calibrated", "rDNA45S.cn_single": "45S, single-sample anchor",
       "rDNA45S.18S.flat": "45S, 18S depth ratio (published)", "rDNA5S.cn": "5S, calibrated",
       "DJ.cn": "distal junction", "HSat2.mass_Mb": "HSat2 mass", "ACRO.mass_Mb": "ACRO mass",
       "HSat3.mass_Mb": "HSat3 mass", "aSatHOR.mass_Mb": "alpha-satellite HOR mass", "TEL.mass_Mb": "telomeric repeat",
       "truth.auto": "held-out autosomal (no true variance)", "chrM.copies": "mtDNA per cell (culture)",
       "chrEBV.copies": "EBV per cell (culture)"}
rel = pd.DataFrame([trio_stats(m) for m in LAB]).assign(label=lambda d: d.metric.map(LAB))
rel.to_csv(f"{OUT}/transmission_reanalysis.tsv", sep="\t", index=False, float_format="%.4g")

def paired(a, b, B=4000):
    xa, xb = centred(a), centred(b)
    def R(x, i):
        c = x.loc[tr.index].values[i]; f = x.loc[tr.FatherID].values[i]; mo = x.loc[tr.MotherID].values[i]
        return _est(c, f, mo)[0]
    idx = np.arange(len(tr))
    dd = [R(xa, i) - R(xb, i) for i in (rng.integers(0, len(tr), len(tr)) for _ in range(B))]
    return [R(xa, idx) - R(xb, idx), *np.percentile(dd, [2.5, 97.5])]
paired_cn_flat = paired("rDNA45S.cn", "rDNA45S.18S.flat")

# sex-specific single-parent slopes, 45S calibrated
x = centred("rDNA45S.cn"); sexrows = []
for par, col in [("father", "FatherID"), ("mother", "MotherID")]:
    for sx, lab in [(1, "son"), (2, "daughter")]:
        t = tr[tr.Sex == sx]; c = x.loc[t.index].values; p = x.loc[t[col]].values
        bs = [np.cov(c[i], p[i])[0, 1] / np.var(p[i], ddof=1) for i in (rng.integers(0, len(c), len(c)) for _ in range(4000))]
        sexrows.append(dict(pair=f"{par} to {lab}", n=len(c), slope=np.cov(c, p)[0, 1] / np.var(p, ddof=1),
                            lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
sex_tab = pd.DataFrame(sexrows)
sex_tab.to_csv(f"{OUT}/transmission_by_sex_reanalysis.tsv", sep="\t", index=False, float_format="%.4g")

# ---------- 3. population structure, pilot replicate precision ----------
b25 = S[S.RELEASE_BATCH == 2504]
fpop = smf.ols('np.log(Q("rDNA45S.cn")) ~ C(Population)', data=b25).fit()
rep = pd.DataFrame(json.load(open(f"{ROOT}/docs/report.json"))["replicates"]["points"])
cal, fl = rep[rep.si == 0].reset_index(drop=True), rep[rep.si == 1].reset_index(drop=True)
def wcv(df, centre):
    lr = np.log(df.y / df.x); lr = lr - lr.mean() if centre else lr
    return np.sqrt(np.mean(lr ** 2) / 2)
bsr = np.array([wcv(fl.iloc[i], True) / wcv(cal.iloc[i], False) for i in (rng.integers(0, 12, 12) for _ in range(5000))])

summary = dict(n_genomes=len(S), n_trios=len(tr),
               batch_counts={str(k): v for k, v in pd.crosstab(S.gen, S.RELEASE_BATCH).to_dict().items()},
               paired_R_cn_minus_flat=paired_cn_flat,
               pop_R2=fpop.rsquared, pop_F_p=fpop.f_pvalue, pop_n=int(fpop.nobs),
               replicate_within_cv_cal=wcv(cal, False), replicate_within_cv_flat_centred=wcv(fl, True),
               replicate_cv_ratio_ci=list(np.percentile(bsr, [2.5, 50, 97.5])),
               r_cn_flat=float(np.corrcoef(S["rDNA45S.cn"], S["rDNA45S.18S.flat"])[0, 1]))
json.dump(summary, open(f"{OUT}/review_summary.json", "w"), indent=1, default=float)

# ---------- figures ----------
try:
    apply_figure_style(sizes=(8, 7, 6))   # figure-style helpers, when run in the analysis kernel
except NameError:
    plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})

GROUPS = [("parents, 2504 batch", (S.gen == "parent") & (S.RELEASE_BATCH == 2504), "#8c8c8c"),
          ("parents, 698 batch", (S.gen == "parent") & (S.RELEASE_BATCH == 698), "#d95f02"),
          ("children, 698 batch", (S.gen == "child") & (S.RELEASE_BATCH == 698), "#1f78b4")]
panels = [("28S / 18S (no GC model)", "28S / 18S reads\n(no GC model)"),
          ("5'ETS / IGS (GC model)", "5'ETS / IGS\n(GC model)"),
          ("18S depth ratio / calibrated 45S", "18S depth ratio /\ncalibrated 45S"),
          ("insert-size median", "insert-size median (bp)")]
fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.9))
for ax, (key, ylab), letter in zip(axes, panels, "abcd"):
    v = np.exp(feat[key])
    for j, (g, mask, col) in enumerate(GROUPS):
        y = v[mask].dropna().values
        xj = j + rng.uniform(-0.28, 0.28, len(y))
        ax.scatter(xj, y, s=4 if len(y) > 50 else 12, color=col, alpha=0.45 if len(y) > 50 else 0.95, lw=0)
        ax.hlines(np.median(y), j - 0.35, j + 0.35, color="k", lw=1.4, zorder=3)
    r = batch_tab.set_index("quantity").loc[key]
    ax.set_title(f"batch {r.batch_effect_pct:+.1f}%, p={r.batch_p:.1g}\ngeneration {r.generation_effect_pct:+.1f}%, p={r.generation_p:.2f}",
                 loc="left", fontsize=6.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(["parents\n2504", "parents\n698", "children\n698"])
    ax.set_ylabel(ylab); ax.margins(x=0.08, y=0.06)
    ax.text(-0.38, 1.30, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
fig.tight_layout(w_pad=1.0)
fig.savefig(f"{OUT}/fig_batch_signature.png", dpi=300)
plt.close(fig)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1.6, 1]))
rr = rel.iloc[::-1].reset_index(drop=True)
for i, r in rr.iterrows():
    a1.errorbar(r.R, i + 0.15, xerr=[[r.R - r.R_lo], [r.R_hi - r.R]], fmt="o", ms=3.5, color="#1f78b4", lw=1)
    a1.errorbar(r.R_rescaled, i - 0.15, xerr=[[r.R_rescaled - r.R_rescaled_lo], [r.R_rescaled_hi - r.R_rescaled]],
                fmt="s", ms=3, color="#d95f02", lw=1)
a1.axvline(0, color="#bbbbbb", lw=0.8); a1.axvline(1, color="#bbbbbb", lw=0.8, ls="--")
a1.set_yticks(range(len(rr))); a1.set_yticklabels(rr.label)
a1.set_xlabel("transmission reliability R")
a1.plot([], [], "o", color="#1f78b4", label="as measured")
a1.plot([], [], "s", color="#d95f02", label="children rescaled to\nparents' spread")
a1.legend(frameon=False, loc="lower right", markerscale=0.8, fontsize=6.5)
a1.set_title(f"Inherited: rDNA and most satellites ({len(tr)} trios)", loc="left")
a1.text(-0.75, 1.07, "a", transform=a1.transAxes, fontsize=10, fontweight="bold")
st = sex_tab.iloc[::-1].reset_index(drop=True)
for i, r in st.iterrows():
    a2.errorbar(r.slope, i, xerr=[[r.slope - r.lo], [r.hi - r.slope]], fmt="o", ms=4,
                color="#6a3d9a" if r.pair == "mother to son" else "#555555", lw=1)
    a2.text(1.14, i, f"n={r.n}", va="center", fontsize=6.5)
a2.set_yticks(range(len(st))); a2.set_yticklabels(st.pair); a2.set_xlim(0, 1.35); a2.set_ylim(-0.6, 3.6)
a2.set_xlabel("parent-offspring slope (45S)")
a2.set_title("Mother to son reads low", loc="left")
a2.text(-0.5, 1.07, "b", transform=a2.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_transmission.png", dpi=300)
plt.close(fig)
print("done")
