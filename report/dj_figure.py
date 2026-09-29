"""The distal-junction figure: the cohort's segment map, every compared genome's profile against its
HPRC assembly, and the per-genome agreement. Drawn by the report when matplotlib is installed
(docs/dj_assemblies.png); `python -m report.dj_figure --report docs/report.json --data docs/data -o FILE`
redraws it from the published tables."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

EXPECTED = 10.0


def render(dj: dict, out) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    asm = dj["assemblies"]
    samples = sorted(asm["samples"], key=lambda t: (t["assembly_mean"], t["sample"]))
    order = [t["sample"] for t in samples]
    blocks = {}
    for b in asm["blocks"]:
        blocks.setdefault(b["sample"], {})[b["block_kb"]] = b
    kbs = sorted({b["block_kb"] for b in asm["blocks"]})
    A = np.array([[blocks[s][k]["assembly"] if blocks[s][k]["assembly"] is not None else np.nan for k in kbs] for s in order], float)
    R = np.array([[blocks[s][k]["reads"] if blocks[s][k]["reads"] is not None else np.nan for k in kbs] for s in order], float)
    K = np.array([[blocks[s][k]["call"] if blocks[s][k].get("call") is not None else np.nan for k in kbs] for s in order], float)
    rows_per = 3 if np.isfinite(K).any() else 2
    seg = dj["segments"]
    n = len(order)
    fig = plt.figure(figsize=(13, 3.6 + 0.14 * rows_per * n))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.0, 0.14 * rows_per * n], width_ratios=[1.35, 1], hspace=0.32, wspace=0.28)

    # A: the cohort's segment map
    ax = fig.add_subplot(gs[0, 0])
    kb = np.array(seg["sub_kb"], float)
    sd = np.array([np.nan if v is None else v for v in seg["sd"]], float)
    for a, b in seg["hyper_intervals"]:
        ax.axvspan(a / 1000, b / 1000, color="#f4c7a1", alpha=0.7, lw=0)
    ax.plot(kb + 2.5, sd, color="#1f4e79", lw=1.6)
    ax.axhline(seg["floor"], color="#888", lw=1, ls="--")
    ax.set_xlim(0, 400)
    ax.set_xlabel("position in the distal-junction unit (kb)")
    ax.set_ylabel("cohort SD of the calibrated estimate\nper 5-kb sub-block (copies)")
    ax.set_title(f"A. Where the junction varies between people ({dj['n']:,} genomes)", loc="left", fontsize=10.5, fontweight="bold")
    ax.text(398, seg["floor"], " noise floor", va="bottom", ha="right", fontsize=8, color="#666")
    y0, y1 = ax.get_ylim()
    for i, (a, b) in enumerate(seg["hyper_intervals"]):
        ax.annotate(f"{a // 1000}–{b // 1000} kb", ((a + b) / 2000, y1), xytext=(0, -4 - 11 * (i % 2)), textcoords="offset points",
                    ha="left" if a < 20000 else "center", va="top", fontsize=7.5, color="#8a4b08", annotation_clip=False)

    # B: per-genome agreement
    ax = fig.add_subplot(gs[0, 1])
    res = np.array([t["resolved"] for t in samples])
    x = np.array([t["assembly_mean"] for t in samples]); y = np.array([t["reads_mean"] for t in samples])
    lo, hi = min(x.min(), y.min()) - 0.3, max(x.max(), y.max()) + 0.3
    ax.plot([lo, hi], [lo, hi], color="#bbb", lw=1, zorder=1)
    ax.scatter(x[res], y[res], s=34, color="#1f4e79", label="assembly resolved", zorder=3)
    ax.scatter(x[~res], y[~res], s=34, facecolor="white", edgecolor="#c0504d", lw=1.4, label="assembly fragmented at the junction", zorder=3)
    labelled = sorted(((xi, yi, t["sample"]) for t, xi, yi in zip(samples, x, y) if abs(yi - xi) > 0.35 or t["group"] == "carrier"), key=lambda v: (v[1], v[0]))
    placed = []                                            # push a label up when it would sit on the previous one
    for xi, yi, name in labelled:
        yl = yi
        for px, py in placed:
            if abs(xi - px) < 0.45 and abs(yl - py) < 0.11:
                yl = py + 0.11
        placed.append((xi, yl))
        ax.annotate(name, (xi, yi), xytext=(xi + 0.06, yl), textcoords="data", fontsize=7, color="#333", va="center",
                    arrowprops=dict(arrowstyle="-", color="#bbb", lw=0.6) if abs(yl - yi) > 0.02 else None)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("assembly: copies per 20-kb block, mean over the unit")
    ax.set_ylabel("NGS-DOSE, on the ten-copy scale")
    st = asm["stats"]["resolved"]
    ax.set_title(f"B. Genome by genome ({asm['stats']['n']} compared)", loc="left", fontsize=10.5, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    if st.get("n", 0) >= 2:
        ax.text(0.98, 0.04, f"resolved (n = {st['n']}): reads − assembly {st['mean_diff_mean']:+.2f} ± {st['mean_diff_sd']:.2f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#333")

    # C: the profiles, assembly and reads side by side
    ax = fig.add_subplot(gs[1, :])
    M = np.full((rows_per * n, len(kbs)), np.nan)
    M[0::rows_per], M[1::rows_per] = A, R
    if rows_per == 3:
        M[2::3] = K
    norm = TwoSlopeNorm(vmin=7, vcenter=EXPECTED, vmax=13)
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", norm=norm, interpolation="nearest")
    ax.set_yticks(np.arange(0, rows_per * n, rows_per) + (rows_per - 1) / 2)
    ax.set_yticklabels([f"{t['sample']}{'' if t['resolved'] else ' ⚠'}  {t['assembly_mean']:.1f} | {t['reads_mean']:.1f}" for t in samples], fontsize=7.5)
    for i in range(n):
        ax.axhline(rows_per * i + rows_per - 0.5, color="white", lw=1.6)
    ax.set_xticks(np.arange(len(kbs))[::2])
    ax.set_xticklabels([f"{k}" for k in kbs[::2]], fontsize=8)
    ax.set_xlabel("20-kb block of the unit (kb); each genome: assembly (upper row), NGS-DOSE's profile (middle)" + (", its integer call (lower)" if rows_per == 3 else ""))
    core = set(seg.get("core_kb", []))
    for j, k in enumerate(kbs):
        if k not in core:
            ax.axvspan(j - 0.5, j + 0.5, ymin=0, ymax=1, facecolor="none", edgecolor="#8a4b08", hatch="//", lw=0, alpha=0.25)
    cb = fig.colorbar(im, ax=ax, fraction=0.012, pad=0.01)
    cb.set_label("copies (10 = one per acrocentric short arm)", fontsize=8)
    ax.set_title("C. Along the unit: HPRC assembly (k-mer block medians, both haplotypes), NGS-DOSE's calibrated profile, and the integer states called from it\n"
                 "⚠ assembly fragmented at the junction; hatched blocks are the segments left out of the level", loc="left", fontsize=9.5, fontweight="bold")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return Path(out).name


def from_files(report_json, data_dir) -> dict:
    d = json.load(open(report_json))["dj"]
    with open(Path(data_dir) / "dj_hprc_blocks.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["block_kb"] = int(r["block_kb"])
        r["assembly"] = float(r["assembly"]) if r["assembly"] not in ("", "None", "NA") else None
        r["reads"] = float(r["reads"]) if r["reads"] not in ("", "None", "NA") else None
        r["call"] = float(r["call"]) if r.get("call") not in ("", "None", "NA", None) else None
    d["assemblies"]["blocks"] = rows
    return d


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(argv)
    print(render(from_files(a.report, a.data), a.out))


if __name__ == "__main__":
    main()
