#!/usr/bin/env python3
"""Evaluate the 1000 Genomes pilot: four trios (NYGC 30x, NovaSeq 2x150) and an independent, older
library of each of the twelve samples (HGSVC 2015: HiSeq 2500 2x126, ~75x, nine samples; Illumina
Platinum pedigree: HiSeq 2000 2x101, ~50x, the CEU trio).

    python pilot/evaluate_pilot.py            # writes pilot_report.md + tables here

Inputs are the engine's counts files in counts_nygc/ and counts_replicates/ next to this script; they
were produced by `run_pilot.sh` straight from the public CRAMs (no download of whole files). Also read:
hall2021_MOESM1.txt here (or the repository's meta/ copy) and, for the cohort figures the report quotes,
the repository's docs/report.json when present. Every input is checked before anything is estimated, and
every output is written only once all of them have been computed, so a failed run leaves no partial set.
"""
import csv
import json
import sys
from io import StringIO as io_text
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
from ngsdose import cohort, estimate, io, resources  # noqa: E402

TRIOS = [("NA12878", "NA12891", "NA12892", "CEU"), ("NA19240", "NA19239", "NA19238", "YRI"),
         ("HG00514", "HG00512", "HG00513", "CHS"), ("HG00733", "HG00731", "HG00732", "PUR")]
SAMPLES = [s for t in TRIOS for s in t[:3]]          # each counted from NYGC and from an older library
DESIGN = "[NGS-DOSE docs/DESIGN.md](https://github.com/jlanej/NGS-DOSE/blob/main/docs/DESIGN.md)"
SEX = dict(NA12878="F", NA12891="M", NA12892="F", NA19240="F", NA19239="M", NA19238="F",
           HG00514="F", HG00512="M", HG00513="F", HG00733="F", HG00731="M", HG00732="F")


def run(dirname, bundle, panel, units, feats, anchors=None):
    out, tables = {}, {}
    for f in sorted((HERE / dirname).glob("*.json.gz")):
        c = io.load_counts(f)
        L = estimate.nearest_table(c, None)["l"]
        if L not in tables:
            tables[L] = bundle.region_tables(L)
        out[c["sample"]] = estimate.estimate_sample(c, panel, units, feats, region_tables=tables[L], regions=bundle.regions(), anchors=anchors)
    return out


def check_counts():
    """Every pilot sample has both counts files; a missing one would silently shrink every n below."""
    missing = [f"{d}/{s}.json.gz" for d in ("counts_nygc", "counts_replicates") for s in SAMPLES if not (HERE / d / f"{s}.json.gz").exists()]
    if missing:
        sys.exit(f"evaluate_pilot: {len(missing)} counts file(s) missing in {HERE}: {', '.join(missing)} (run_pilot.sh counts them)")


def hall_table():
    """Hall et al. 2021 Supplementary Data 1, from here or the repository's meta/ copy; None if neither exists."""
    for path in (HERE / "hall2021_MOESM1.txt", HERE.parent / "meta" / "hall2021_MOESM1.txt"):
        if path.exists():
            with open(path, newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            lack = [c for c in ("Sample", "HC.18S.CN") if not rows or c not in rows[0]]
            if lack:
                sys.exit(f"evaluate_pilot: {path} has no column {' or '.join(lack)}: it is not Hall et al.'s Supplementary Data 1 "
                         "(an error page saved in its place?); copy the repository's meta/hall2021_MOESM1.txt over it")
            return {r["Sample"]: r for r in rows}
    return None


def cohort_report():
    """The cohort page's numbers (docs/report.json), for the figures this report quotes from it; None without it."""
    path = HERE.parent / "docs" / "report.json"
    return json.loads(path.read_text()) if path.exists() else None


def release_batches():
    """sample -> release batch ('2504' or '698') from NGS-PCA's QC table, when the repository has it."""
    path = HERE.parent / "meta" / "ngspca_sample_qc.tsv"
    if not path.exists():
        return {}
    with open(path, newline="") as fh:
        return {r["SAMPLE_ID"].strip(): r.get("RELEASE_BATCH", "").strip() for r in csv.DictReader(fh, delimiter="\t")}


def window_differences(ny, rep, common, cls="rDNA45S"):
    """d[pair, window] = log(C_replicate / C_nygc), each library under its own GC model."""
    w0 = ny[common[0]]["classes"][cls]["windows"]
    start = np.array([w["start"] for w in w0]); end = np.array([w["end"] for w in w0]); gc = np.array([w["gc"] for w in w0])
    D = np.full((len(common), len(w0)), np.nan)
    for i, s in enumerate(common):
        for j, (x, y) in enumerate(zip(ny[s]["classes"][cls]["windows"], rep[s]["classes"][cls]["windows"])):
            if x["cn"] and y["cn"]:
                D[i, j] = np.log(y["cn"] / x["cn"])
    return start, end, gc, D


def consensus_windows(D, gc, max_abs=0.06, max_spread=0.15, min_pairs=3):
    """Moderate-GC windows on which two library chemistries agree, reproducibly across individuals."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(D, axis=0)
        spread = np.nanpercentile(D, 90, axis=0) - np.nanpercentile(D, 10, axis=0)
    n = np.sum(~np.isnan(D), axis=0)
    return (gc >= 0.40) & (gc <= 0.60) & (np.abs(med) < max_abs) & (spread < max_spread) & (n >= min_pairs)


def main():
    write_anchors = "--write-anchors" in sys.argv
    check_counts()
    H, rep_json = hall_table(), cohort_report()
    outputs, anchors_out = {}, None                  # file name -> text, and anchors.json's text, written at the end
    bundle = resources.Bundle()                      # NGSDOSE_RESOURCES, or the installed checkout's resources/GRCh38
    panel, units, feats = io.load_panel(bundle.panel), bundle.units(), bundle.features()
    ny = run("counts_nygc", bundle, panel, units, feats)
    hg = run("counts_replicates", bundle, panel, units, feats)
    lines = ["# 1000 Genomes pilot: results", "",
             f"{len(ny)} NYGC samples (4 trios), {len(hg)} independent-library replicates (HGSVC and Illumina Platinum). "
             "Generated by `evaluate_pilot.py`; every number below is recomputed from the committed counts files"
             + (", except the cohort figures quoted in sections 1 and 4, which come from the cohort page's docs/report.json." if rep_json else "."), ""]

    # ---------------- per-sample table
    cal_ny = cohort.calibrate(list(ny.values()), "rDNA45S")
    cn_ny = dict(zip(cal_ny.samples, np.exp(cal_ny.c)))
    rows = []
    for s, r in ny.items():
        k = r["classes"]
        rows.append(dict(sample=s, sex=SEX.get(s, "?"), depth=r["depth_equiv"], insert=r["insert_median"], dup=r["ctrl_dup_frac"],
                         gc65=r["gc_rel"]["65"], auto=r["truth_regions"]["auto"]["cn"], chrX=r["truth_regions"]["chrX"]["cn"],
                         chrY=(r["truth_regions"].get("chrY") or {}).get("cn", float("nan")), chrM=(r["truth_regions"].get("chrM") or {}).get("cn", float("nan")),
                         chrEBV=(r["truth_regions"].get("chrEBV") or {}).get("cn", float("nan")),
                         DJ=k["DJ"]["cn"], rDNA5S=k["rDNA5S"]["cn"], rDNA45S_anchor=k["rDNA45S"]["cn_anchor"], rDNA45S_cal=cn_ny[s],
                         r18S_flat=k["rDNA45S"]["features"]["18S"]["cn_flat"], r28S_18S=k["rDNA45S"]["features"]["28S"]["cn"] / k["rDNA45S"]["features"]["18S"]["cn"]))
    fh = io_text()
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    outputs["pilot_samples.tsv"] = fh.getvalue()

    lines += ["## 1. Known-truth controls (NYGC)", "",
              "| sample | sex | depth | dup-flag | GC rate @65% | held-out autosomal (truth 2) | chrX (truth 1 / 2) | DJ (truth 10) |",
              "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        lines.append(f"| {r['sample']} | {r['sex']} | {r['depth']:.1f} | {r['dup']:.3f} | {r['gc65']:.2f} | {r['auto']:.3f} | {r['chrX']:.3f} | {r['DJ']:.2f} |")
    a = np.array([r["auto"] for r in rows]); d = np.array([r["DJ"] for r in rows])
    xm = np.array([r["chrX"] for r in rows if r["sex"] == "M"]); xf = np.array([r["chrX"] for r in rows if r["sex"] == "F"])
    lines += ["", f"Held-out autosomal: mean {a.mean():.4f}, SD {a.std(ddof=1):.4f}. chrX: males {xm.mean():.3f} (SD {xm.std(ddof=1):.3f}), "
                  f"females {xf.mean():.3f} (SD {xf.std(ddof=1):.3f}). DJ: mean {d.mean():.2f}, SD {d.std(ddof=1):.2f}, range {d.min():.2f}-{d.max():.2f}.", ""]
    if hg:
        ha = np.array([r["truth_regions"]["auto"]["cn"] for r in hg.values()]); hd = np.array([r["classes"]["DJ"]["cn"] for r in hg.values()])
        hxm = np.array([r["truth_regions"]["chrX"]["cn"] for s, r in hg.items() if SEX.get(s) == "M"])
        hxf = np.array([r["truth_regions"]["chrX"]["cn"] for s, r in hg.items() if SEX.get(s) == "F"])
        # women whose NYGC culture has lost an X in part of its cells (below 1.85 copies) are set aside for the female-X range
        lost = sorted(r["sample"] for r in rows if r["sex"] == "F" and r["chrX"] < 1.85)
        xf_int = np.array([r["chrX"] for r in rows if r["sex"] == "F" and r["sample"] not in lost])
        lost_txt = "; ".join(f"{s} has lost an X chromosome in part of its culture - {next(r['chrX'] for r in rows if r['sample'] == s):.2f} copies in "
                             f"the NYGC DNA, {hg[s]['truth_regions']['chrX']['cn']:.2f} in the older DNA" for s in lost if s in hg)
        sphase = ""
        dj_x = ((rep_json or {}).get("biology") or {}).get("DJ_vs_chrX_female")
        if dj_x and dj_x.get("n"):
            meta = rep_json.get("meta", {})
            sphase = (f" The cohort report has checked the first part of that prediction. On {meta.get('as_of', '?')}, with {meta.get('n', 0):,} genomes "
                      f"counted, it compared {dj_x['n']} women whose cell line has kept both X chromosomes (chrX 1.85-2.15 copies) and found "
                      f"r = {dj_x['r']:.2f} (95% CI {dj_x['r_lo']:.2f} to {dj_x['r_hi']:.2f}) between DJ and chrX. "
                      + (("The two deficits do not track each other, which weakens a shared S-phase cause without ruling out a small one; "
                         "the control-PC part of the prediction is untested.")
                         if dj_x["r_lo"] < 0 < dj_x["r_hi"] else
                         "They move together, as the hypothesis predicts." if dj_x["r_lo"] > 0 else
                         "They move in opposite directions, against the hypothesis."))
        lines += [f"The same controls in the {len(hg)} older libraries: autosomal {ha.mean():.4f} (SD {ha.std(ddof=1):.4f}); chrX males {hxm.mean():.3f}, "
                  f"females {hxf.mean():.3f} (SD {hxf.std(ddof=1):.3f}); DJ {hd.mean():.2f} (SD {hd.std(ddof=1):.2f}, range {hd.min():.2f}-{hd.max():.2f}).", "",
                  "Two things to read off these numbers. " + (f"{lost_txt} - so the two DNA batches of a cell line are not the same biological sample, and "
                  "disagreement between libraries bounds measurement error from above. And " if lost_txt else "")
                  + f"the NYGC libraries read both the female X ({xf_int.min():.2f}-{xf_int.max():.2f}"
                  + (f", {', '.join(lost)} aside" if lost else "") + f") and the DJ ({d.min():.2f}-{d.max():.2f}) a few percent low while the older libraries "
                  "read them near 2 and 10: the inactive X and the acrocentric short arms replicate late, and DNA from a culture with more cells in S phase "
                  "under-represents late-replicating sequence. That is a hypothesis, not a result. It predicts that DJ, female X and the leading control PC "
                  f"move together across a cohort ({DESIGN}, section 13, \"The S-phase hypothesis\")." + sphase, ""]

    # ---------------- trios
    lines += ["## 2. Transmission (4 trios; descriptive only at this n)", "",
              "| class | child | child | father | mother | midparent | child - midparent | within [0, F+M] |", "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    val = {"rDNA45S": cn_ny, "rDNA5S": {s: r["classes"]["rDNA5S"]["cn"] for s, r in ny.items()}, "DJ": {s: r["classes"]["DJ"]["cn"] for s, r in ny.items()}}
    for cls, v in val.items():
        for c, f, m, pop in TRIOS:
            if all(x in v for x in (c, f, m)):
                mid = (v[f] + v[m]) / 2
                lines.append(f"| {cls} | {c} ({pop}) | {v[c]:.1f} | {v[f]:.1f} | {v[m]:.1f} | {mid:.1f} | {v[c] - mid:+.1f} ({100 * (v[c] / mid - 1):+.1f}%) | {'yes' if v[c] <= v[f] + v[m] else 'NO'} |")
    lines.append("")

    # ---------------- cross-library replicates
    common = [s for s in hg if s in ny]
    cn_hg_x = {}
    lines += ["## 3. Independent libraries of the same cell line", "",
              "NYGC (TruSeq PCR-free, NovaSeq 2x150, 2019) against an older library of the same sample: HGSVC (HiSeq 2500 2x126, "
              "2015, bwakit + postalt) or, for the CEU trio, Illumina Platinum (HiSeq 2000 2x101, 2012-13). Same resources and the "
              "same sinks, learned from NYGC bwa-mem scans; all counts are fetch mode, so their capture on the bwakit realignments of "
              "the older libraries is assumed, not measured. Read length, insert size, GC behaviour, depth and alignment pipeline all "
              "differ (the aligner is bwa-mem, ALT-aware, against the same GRCh38 analysis set in both, in another version and "
              "post-processing). The cell-line cultures differ too, so disagreement bounds measurement error from above.", ""]
    if len(common) >= 3:
        cal_hg = cohort.calibrate([hg[s] for s in common], "rDNA45S")
        cn_hg_own = dict(zip(cal_hg.samples, np.exp(cal_hg.c)))
        a_ny = cal_ny.a
        cal_x = cohort.calibrate([hg[s] for s in common], "rDNA45S", a_fixed=a_ny)
        cn_hg_x = dict(zip(cal_x.samples, np.exp(cal_x.c)))
        est = {
            "18S, no GC model (published practice)": (lambda r: r["classes"]["rDNA45S"]["features"]["18S"]["cn_flat"],) * 2,
            "28S, no GC model": (lambda r: r["classes"]["rDNA45S"]["features"]["28S"]["cn_flat"],) * 2,
            "all windows, no GC model": (lambda r: r["classes"]["rDNA45S"]["cn_all_flat"],) * 2,
            "18S, fragment-GC model": (lambda r: r["classes"]["rDNA45S"]["features"]["18S"]["cn"],) * 2,
            "28S, fragment-GC model": (lambda r: r["classes"]["rDNA45S"]["features"]["28S"]["cn"],) * 2,
            "moderate-GC (40-60%) windows, fragment-GC model (headline rule before anchors.json; `--gc-rule-anchors`)": (lambda r: r["classes"]["rDNA45S"]["cn_anchor"],) * 2,
        }
        table = []
        for name, (fa, fb) in est.items():
            x = np.array([fa(ny[s]) for s in common]); y = np.array([fb(hg[s]) for s in common])
            table.append((name, x, y))
        table.append(("calibrated, efficiencies learned per library", np.array([cn_ny[s] for s in common]), np.array([cn_hg_own[s] for s in common])))
        table.append(("calibrated, NYGC efficiencies applied to the replicates", np.array([cn_ny[s] for s in common]), np.array([cn_hg_x[s] for s in common])))
        table.append(("rDNA5S, fragment-GC model", np.array([ny[s]["classes"]["rDNA5S"]["cn"] for s in common]), np.array([hg[s]["classes"]["rDNA5S"]["cn"] for s in common])))
        table.append(("rDNA5S, no GC model", np.array([ny[s]["classes"]["rDNA5S"]["cn_all_flat"] for s in common]), np.array([hg[s]["classes"]["rDNA5S"]["cn_all_flat"] for s in common])))
        table.append(("DJ, fragment-GC model", np.array([ny[s]["classes"]["DJ"]["cn"] for s in common]), np.array([hg[s]["classes"]["DJ"]["cn"] for s in common])))
        lines += ["| estimator | mean log ratio replicate/NYGC | SD of log ratio (pair noise) | Pearson r |", "| --- | --- | --- | --- |"]
        rep = {}
        for name, x, y in table:
            lr = np.log(y / x)
            rep[name] = dict(bias=float(lr.mean()), sd=float(lr.std(ddof=1)), r=float(np.corrcoef(x, y)[0, 1]))
            lines.append(f"| {name} | {lr.mean():+.4f} | {lr.std(ddof=1):.4f} | {np.corrcoef(x, y)[0, 1]:.4f} |")
        # the single-sample headline as `ngsdose estimate` reports it today, with the bundle's anchors.json - which was learned from
        # these same pairs (section 3, window by window), so this figure is in sample
        anchors = bundle.anchors()
        if anchors.get("rDNA45S"):
            ny_a = run("counts_nygc", bundle, panel, units, feats, anchors=anchors)
            hg_a = run("counts_replicates", bundle, panel, units, feats, anchors=anchors)
            xa = np.array([ny_a[s]["classes"]["rDNA45S"]["cn"] for s in common]); ya = np.array([hg_a[s]["classes"]["rDNA45S"]["cn"] for s in common])
            la = np.log(ya / xa)
            lines += ["", "This script estimates with the fragment-GC 40-60% anchor rule throughout, because the shipped anchors.json "
                          f"(the {sum((e - s) // 250 for s, e in anchors['rDNA45S'])} windows below) was learned from these same {len(common)} pairs. "
                          f"With anchors.json the single-sample headline reads {la.mean():+.4f} (SD {la.std(ddof=1):.4f}, r = {np.corrcoef(xa, ya)[0, 1]:.3f}), "
                          "but that is in sample; the out-of-sample figure is the held-out single-sample result below."]
        lines += ["", "Library properties: " + "; ".join(f"{s}: NYGC {ny[s]['depth_equiv']:.0f}x ins {ny[s]['insert_median']} GC65 {ny[s]['gc_rel']['65']:.2f} / "
                                                          f"replicate {hg[s]['depth_equiv']:.0f}x read {hg[s]['read_length']} ins {hg[s]['insert_median']} GC65 {hg[s]['gc_rel']['65']:.2f}" for s in common), ""]
        outputs["pilot_replicates.tsv"] = ("sample\t" + "\t".join(f"{n} [NYGC]\t{n} [replicate]" for n, _, _ in table) + "\n"
                                           + "".join(s + "\t" + "\t".join(f"{x[i]:.3f}\t{y[i]:.3f}" for _, x, y in table) + "\n"
                                                     for i, s in enumerate(common)))
        outputs["pilot_replicates.json"] = json.dumps(rep, indent=1)

    # ---------------- where the libraries agree, window by window
    if len(common) >= 6:
        start, end, gc, D = window_differences(ny, hg, common)
        fam = {s: next(t[3] for t in TRIOS if s in t[:3]) for s in common}
        mod = (gc >= 0.40) & (gc <= 0.60)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(D, axis=0)
            spread = np.nanpercentile(D, 90, axis=0) - np.nanpercentile(D, 10, axis=0)
        lines += ["### Window by window", "",
                  f"d = log(C_replicate / C_NYGC) per 250-bp window of the 45S unit, each library under its own GC model, {len(common)} pairs. "
                  f"Over the {int(mod.sum())} moderate-GC (40-60%) windows the median d ranges from {np.nanmin(med[mod]):+.2f} to {np.nanmax(med[mod]):+.2f}, "
                  f"yet for a given window it is nearly the same in every individual (median 10-90% spread across pairs {np.nanmedian(spread[mod]):.2f}): "
                  "the disagreement is a deterministic property of sequence x chemistry - dropout zones one fragment wide around particular "
                  "elements, strand-asymmetric at fine scale - not noise, and not something a genome-wide GC curve can carry.", ""]
        def intervals(mask):
            iv = []
            for s0, e0 in zip(start[mask], end[mask]):
                if iv and iv[-1][1] == s0:
                    iv[-1][1] = int(e0)
                else:
                    iv.append([int(s0), int(e0)])
            return iv

        held = []
        for f in sorted(set(fam.values())):
            train = [i for i, s in enumerate(common) if fam[s] != f]
            test = [s for s in common if fam[s] == f]
            if len(train) < 3 or not test:
                continue
            iv_f = intervals(consensus_windows(D[train], gc))
            # each library calibrated on its own (window efficiencies are chemistry-specific); only the
            # level is shared, through anchors that were chosen without this family
            c_ny = cohort.calibrate(list(ny.values()), "rDNA45S", anchors=iv_f)
            v_ny, v_rp = dict(zip(c_ny.samples, np.exp(c_ny.c))), {}
            for rl in sorted({hg[s]["read_length"] for s in common}):
                grp = [hg[s] for s in common if hg[s]["read_length"] == rl]
                c_rp = cohort.calibrate(grp, "rDNA45S", anchors=iv_f)
                v_rp.update(zip(c_rp.samples, np.exp(c_rp.c)))
            for s in test:
                # the same anchors without any cohort: ratio of sums over the anchor windows of one sample
                def single(r, iv=iv_f):
                    w = [x for x in r["classes"]["rDNA45S"]["windows"] if x["cn"] and any(a0 <= x["start"] and x["end"] <= b0 for a0, b0 in iv)]
                    return 2 * sum(x["obs"] for x in w) / sum(x["exp"] for x in w)
                held.append((s, f, int(c_ny.anchor.sum()), float(v_ny[s]), float(v_rp[s]), single(ny[s]), single(hg[s])))
        if held:
            lr = np.log([h[4] / h[3] for h in held])
            lr1 = np.log([h[6] / h[5] for h in held])
            lines += ["**Consensus windows** - moderate GC, |median d| < 0.06, 10-90% spread < 0.15 - are windows on which the two "
                      "chemistries agree. They were learned with one family left out; each library was then calibrated separately "
                      "(its own window efficiencies) with those windows as the anchor, and the left-out family's pairs compared. "
                      "The cross-library level is therefore out of sample:", "",
                      "| sample | family held out | anchor windows | NYGC | replicate | log ratio | single-sample NYGC | single-sample replicate | log ratio |",
                      "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
            for (s, f, n, a, h, a1, h1), x, x1 in zip(held, lr, lr1):
                lines.append(f"| {s} | {f} | {n} | {a:.1f} | {h:.1f} | {x:+.3f} | {a1:.1f} | {h1:.1f} | {x1:+.3f} |")
            r_cal = np.corrcoef([h[3] for h in held], [h[4] for h in held])[0, 1]
            r_one = np.corrcoef([h[5] for h in held], [h[6] for h in held])[0, 1]
            lines += ["", f"Held-out agreement, calibrated per library: mean log ratio {lr.mean():+.4f}, SD {lr.std(ddof=1):.4f}, r = {r_cal:.4f} "
                          f"(n = {len(held)}). Single sample, no cohort (ratio of sums over the same anchor windows): mean {lr1.mean():+.4f}, "
                          f"SD {lr1.std(ddof=1):.4f}, r = {r_one:.4f}.", ""]
            outputs["pilot_heldout.tsv"] = ("sample\tfamily_held_out\tn_anchor_windows\tnygc\treplicate\tnygc_18S_flat\treplicate_18S_flat\n"
                                            + "".join(f"{s}\t{f}\t{n}\t{a:.3f}\t{h:.3f}\t{ny[s]['classes']['rDNA45S']['features']['18S']['cn_flat']:.3f}\t"
                                                      f"{hg[s]['classes']['rDNA45S']['features']['18S']['cn_flat']:.3f}\n" for s, f, n, a, h, _a1, _h1 in held))
            rep["calibrated per library, consensus anchors learned with the family held out"] = dict(
                bias=float(lr.mean()), sd=float(lr.std(ddof=1)), r=float(r_cal))
            rep["single sample, consensus anchors learned with the family held out"] = dict(
                bias=float(lr1.mean()), sd=float(lr1.std(ddof=1)), r=float(r_one))
            outputs["pilot_replicates.json"] = json.dumps(rep, indent=1)
        clean = consensus_windows(D, gc)
        iv = intervals(clean)
        lines += [f"All pairs together give {int(clean.sum())} consensus windows ({sum(e - s for s, e in iv):,} bp in {len(iv)} intervals): "
                  + ", ".join(f"{s / 1000:.2f}-{e / 1000:.2f} kb" for s, e in iv) + "."
                  + (" These intervals are the rDNA45S anchors shipped in NGS-DOSE resources/GRCh38/anchors.json (written by "
                     "`python pilot/evaluate_pilot.py --write-anchors`). Unless `--gc-rule-anchors` is given, `ngsdose estimate` and "
                     "`ngsdose cohort` use them to set the absolute level of 45S copy number, together with the fragment-GC 40-60% rule. "
                     "The estimates in this report do not use them: they use the GC rule alone (sections 1-2) or anchors learned with the "
                     "family held out (the table above)." if [list(x) for x in bundle.anchors().get("rDNA45S", [])] == iv else ""), ""]
        if write_anchors:
            anchors_out = json.dumps({"rDNA45S": {
                "intervals": iv, "window": 250,
                "rule": "fragment GC 40-60%, |median d| < 0.06 and 10-90% spread < 0.15, d = log(C_replicate/C_NYGC) per window",
                "derived_from": f"{len(common)} sample pairs: NYGC NovaSeq 2x150 vs HGSVC HiSeq 2500 2x126 / Platinum HiSeq 2000 2x101 "
                                "(NGS-DOSE-1000G: python pilot/evaluate_pilot.py --write-anchors)"}}, indent=1)

    # ---------------- Hall et al. 2021
    if H is not None:
        ss = [s for s in ny if s in H]
        if len(ss) >= 3:
            mine = np.array([ny[s]["classes"]["rDNA45S"]["features"]["18S"]["cn_flat"] / 2 for s in ss])
            theirs = np.array([float(H[s]["HC.18S.CN"]) for s in ss])
            dup_c = np.array([ny[s]["ctrl_dup_frac"] for s in ss]); dup_r = np.array([ny[s]["classes"]["rDNA45S"]["dup_flag_frac"] for s in ss])
            pred = mine * (1 - dup_r) / (1 - dup_c)
            pred_den = mine / (1 - dup_c)
            batch = release_batches()
            in_batch = {b: sum(batch.get(s) == b for s in H) for b in ("2504", "698")}
            table_txt = f"Their table holds {len(H):,} genomes"
            if batch:
                table_txt += (", all from the 2,504 of the first release batch and none of the 698 added later" if in_batch["2504"] == len(H)
                              else f", {in_batch['2504']:,} of them from the first release batch of 2,504")
            absent = [s for s in ny if s not in H and batch.get(s) == "2504"]
            coh = ((rep_json or {}).get("hall") or {})
            coh_txt = ""
            if coh.get("n") and coh.get("flat_ratio"):
                r_c = (coh.get("flat") or {}).get("r")
                coh_txt = (f" The cohort page repeats this comparison on {coh['n']:,} genomes ({(rep_json.get('meta') or {}).get('as_of', '?')}): "
                           + (f"r = {r_c:.3f}, " if r_c is not None else "") + f"mean ratio {coh['flat_ratio']:.3f}"
                           + (f", {coh['dup_denominator_ratio']:.3f} with the duplicate-free denominator" if coh.get("dup_denominator_ratio") else "")
                           + (f", {coh['dup_corrected_ratio']:.3f} with duplicates dropped on both sides" if coh.get("dup_corrected_ratio") else "") + ".")
            lines += ["## 4. Against Hall, Turner & Queitsch (Sci Rep 2021), the same CRAMs", "",
                      f"{len(ss)} shared samples. Their 18S value (mean `samtools depth` over a 145-bp segment of the 18S, reads re-aligned from "
                      "FASTQ, divided by their chromosome-1 depth, which matches our duplicate-excluded depth (an inference from their published values, "
                      "not their stated method)) against a comparable depth ratio from our "
                      f"counts (the callable 50-bp bins of the 18S, no GC model), halved to their per-haploid scale: r = {np.corrcoef(mine, theirs)[0, 1]:.4f}, mean ratio "
                      f"theirs/ours {np.mean(theirs / mine):.3f}. Dividing ours by the unflagged fraction of control reads, 1/(1 - dup_control), the "
                      f"denominator their values imply, brings the mean ratio to {np.mean(theirs / pred_den):.3f} (SD {np.std(theirs / pred_den, ddof=1):.3f}); "
                      f"dropping duplicate-flagged reads on both sides, (1 - dup_rDNA)/(1 - dup_control), brings it to {np.mean(theirs / pred):.3f} "
                      f"(SD {np.std(theirs / pred, ddof=1):.3f}). Most of the offset between the two pipelines fits the duplicate flag (Finding 1 in "
                      f"{DESIGN}, section 2). ({table_txt}. That leaves {len(ss)} shared pilot samples: {', '.join(ss)}"
                      + (f"; {', '.join(absent)}, also of that batch, {'is' if len(absent) == 1 else 'are'} missing from their table" if absent else "")
                      + ".)" + coh_txt, "",
                      "| sample | control dup-flag | rDNA dup-flag | Hall 18S | ours/2 | ours/2 with dup exclusion | ours/2, duplicate-free denominator |",
                      "| --- | --- | --- | --- | --- | --- | --- |"]
            for i, s in enumerate(ss):
                lines.append(f"| {s} | {dup_c[i]:.3f} | {dup_r[i]:.3f} | {theirs[i]:.1f} | {mine[i]:.1f} | {pred[i]:.1f} | {pred_den[i]:.1f} |")
            lines.append("")
    # ---------------- chrY, and the state of the culture
    t = lambda r, k: (r["truth_regions"].get(k) or {}).get("cn", float("nan"))
    if all("chrY" in r["truth_regions"] for r in ny.values()):
        ym = np.array([t(r, "chrY") for s_, r in ny.items() if SEX.get(s_) == "M"]); yf = np.array([t(r, "chrY") for s_, r in ny.items() if SEX.get(s_) == "F"])
        lines += ["## 5. chrY, and the state of the culture", "",
                  f"chrY ({next(iter(ny.values()))['truth_regions']['chrY'].get('n_regions', 40)} X-degenerate regions; truth 1 or 0), NYGC: males {ym.mean():.3f} (SD {ym.std(ddof=1):.3f}, range {ym.min():.3f}-{ym.max():.3f}), "
                  f"females {yf.mean():.4f} (max {yf.max():.4f})." + ("" if not hg else
                  " Older libraries: males {:.3f}, females {:.4f}.".format(np.mean([t(r, "chrY") for s_, r in hg.items() if SEX.get(s_) == "M"]),
                                                                          np.mean([t(r, "chrY") for s_, r in hg.items() if SEX.get(s_) == "F"]))), "",
                  "Mitochondrial genomes and EBV episomes per cell are not truths but covariates: they describe the culture the DNA came from. "
                  "The two libraries of a cell line were made years apart from different cultures, so their disagreement here is biology, and it is "
                  "large next to anything the known-truth controls show:", "",
                  "| sample | chrM copies, NYGC | older library | EBV copies, NYGC | older library | 45S NYGC / older (calibrated) |", "| --- | --- | --- | --- | --- | --- |"]
        ratio45, dm, de = [], [], []
        for s_ in ny:
            r2 = hg.get(s_)
            c45 = (cn_ny[s_] / cn_hg_x[s_]) if (r2 is not None and s_ in cn_hg_x) else float("nan")
            lines.append(f"| {s_} | {t(ny[s_], 'chrM'):.0f} | {t(r2, 'chrM') if r2 else float('nan'):.0f} | {t(ny[s_], 'chrEBV'):.1f} | "
                         f"{t(r2, 'chrEBV') if r2 else float('nan'):.1f} | {c45:.3f} |")
            if r2 is not None and np.isfinite(c45):
                ratio45.append(np.log(c45)); dm.append(np.log(t(ny[s_], "chrM") / t(r2, "chrM"))); de.append(np.log(t(ny[s_], "chrEBV") / t(r2, "chrEBV")))
        if len(ratio45) >= 6:
            lines += ["", f"Between the two cultures of the same line, mitochondrial content differs by a factor of {np.exp(np.std(dm, ddof=1)):.2f} (SD of the log ratio) "
                          f"and EBV load by {np.exp(np.std(de, ddof=1)):.2f}, against {np.exp(np.std(ratio45, ddof=1)):.3f} for calibrated 45S. Correlation of the 45S log ratio "
                          f"with the chrM log ratio: {np.corrcoef(ratio45, dm)[0, 1]:+.2f}; with the EBV log ratio: {np.corrcoef(ratio45, de)[0, 1]:+.2f} "
                          f"(n = {len(ratio45)}" + ("; |r| > 0.58 would be nominally significant" if len(ratio45) == 12 else "") + "). A look, not a test: the cohort has the trios to ask whether "
                          "a child's departure from the midparent tracks the state of its culture.", ""]
    outputs["pilot_report.md"] = "\n".join(lines)
    for name, text in outputs.items():
        (HERE / name).write_text(text)
    if write_anchors and anchors_out:
        out = resources.default_bundle() / "anchors.json"
        out.write_text(anchors_out)
        print("wrote", out, file=sys.stderr)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
