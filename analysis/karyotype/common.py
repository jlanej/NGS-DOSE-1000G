"""Shared by the scripts of analysis/karyotype: every genome's single-copy regions from the cached estimates, read as the
page reads them, and the cohort read against its own model. Run the scripts from the repository's root, after
`regenerate.sh` has filled `cache/` (it estimates counts_scan into cache/scan and counts_karyotype into cache/karyotype)."""
import os
import pickle
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from ngsdose import karyotype as K
from ngsdose import resources
from ngsdose.tables import load_result

warnings.simplefilter("ignore", RuntimeWarning)
BUNDLE = resources.Bundle()
KAR = BUNDLE.karyotype()
ARMS, GC = KAR["arms"], KAR.get("gc")


def _gather(directory):
    files = sorted(f for f in os.listdir(directory) if f.endswith(".estimate.json.gz")) if os.path.isdir(directory) else []

    def one(f):
        r = load_result(os.path.join(directory, f))
        return r["sample"], K.gather(r, KAR["control_names"])
    with ThreadPoolExecutor(8) as ex:
        return dict(ex.map(one, files))


def vectors(cache="cache/scan", windows="cache/karyotype"):
    """(samples, names, matrix, with_windows): each genome's regions as the page reads them, from its estimate with the
    karyotype windows where it has one and from the scan's where it has not; NaN where a genome lacks a region."""
    scan, win = _gather(cache), _gather(windows)
    samples = sorted(scan)
    names, Y = K.assemble([win.get(s) or scan[s] for s in samples])
    return samples, names, Y.astype(float), {s for s in samples if win.get(s)}


def cohort(cache="cache/scan", windows="cache/karyotype", keep="cache/karyotype.analysis.pkl", log=lambda m: print(m, file=sys.stderr)):
    """The cohort read against its own model: dict(samples, names, Y, with_windows, readings, model, info). Kept in `keep`
    and used again while it is newer than both estimate directories."""
    newest = max(os.path.getmtime(d) for d in (cache, windows) if os.path.isdir(d))
    if keep and os.path.exists(keep) and os.path.getmtime(keep) > newest:
        return pickle.load(open(keep, "rb"))
    samples, names, Y, win = vectors(cache, windows)
    readings, model, info = K.cohort([(names, y) for y in Y], ARMS, rules=KAR.get("rules"), log=log, gc=GC)
    out = dict(samples=samples, names=names, Y=Y, with_windows=win, readings=readings, model=model, info=info)
    if keep:
        pickle.dump(out, open(keep, "wb"))
    return out


def aligned(C, rows=None):
    """The genomes' regions in the model's order, the model's table, and the clean genomes (no event, settled), of `rows`
    (indexes into C["samples"]; all by default)."""
    model = C["model"]
    pos = {n: i for i, n in enumerate(C["names"])}
    col = np.array([pos.get(n, -1) for n in model.names])
    Ym = np.where(col[None, :] >= 0, C["Y"][:, np.maximum(col, 0)], np.nan)
    rows = range(len(C["samples"])) if rows is None else rows
    clean = [i for i in rows if C["readings"][i] is not None and not C["readings"][i].events and C["readings"][i].status == "settled"]
    return Ym, K.table(model.names, model.arms), clean


def shuffled_places(tab, model, present, rng):
    """The model's table with each chromosome's places (a window's pieces, kept together) in random order along it: a
    stretch of real change is scattered, so what the chain still calls is chance."""
    start, end = tab.start.copy(), tab.end.copy()
    for c in tab.chromosomes():
        idx = np.flatnonzero(present & (tab.chrom == c) & np.isfinite(model.sd) & (tab.start >= 0))
        if len(idx) < 2:
            continue
        idx = idx[np.argsort(tab.mid[idx], kind="stable")]
        first, _, _, nreg = K.loci(tab.start[idx], tab.end[idx], np.zeros(len(idx)), model.sd[idx])
        for slot, p in enumerate(rng.permutation(len(first))):
            j = idx[first[p]:first[p] + nreg[p]]
            base = 1_000_000 * (slot + 1)
            start[j], end[j] = base + (tab.start[j] - tab.start[j[0]]), base + (tab.end[j] - tab.start[j[0]])
    return K.Table(tab.names, tab.chrom, start, end, tab.arm, tab.kind)


mad = lambda x: float(1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x))))
