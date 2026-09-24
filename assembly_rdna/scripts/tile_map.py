#!/usr/bin/env python
"""Shared rDNA-content measure (the ONE implementation; used for whole assemblies in
validate_full.py and for the extracted sequence of all haplotypes in the next stage).

Measure
  1. Cut sequence into non-overlapping 5-kb tiles on a grid anchored at contig position 0
     (tile k = contig[k*5000, (k+1)*5000)). For extracted records named
     "contig:start-end" (1-based inclusive) the record is placed back on contig
     coordinates, so tiles inside an extracted region are identical to the tiles of the
     whole assembly; tiles cut by a region edge are partial (tile_len < 5000).
     Tiles shorter than MIN_TILE (500 bp) or with > 50% N are skipped.
  2. Map tiles with minimap2 -x asm20 -c --secondary=no against
       target A: 45S unit KY962518.1 concatenated with itself (circular unit), and
       target B: 5S unit X12811.1 doubled + distal-junction core (CHM13 chr21:2,708,299-3,108,298),
     run separately.
  3. Per tile and class (45S / 5S / DJ): aln_frac = union of query intervals of
     alignments with identity >= 0.90 (identity = PAF matches / alignment block length)
     divided by tile length. A tile is class sequence if aln_frac >= 0.50.
  4. Secondary, run-based summary: passing tiles of a class are joined into runs along the
     contig (gaps of <= 1 tile allowed); {cls}_bp_runs_ge50kb counts passing-tile bp in runs
     >= 50 kb. Validation showed the DJ class is unspecific per tile (~44-47 Mb of
     interspersed-repeat hits genome-wide) while DJ runs >= 50 kb are confined to acrocentric
     short arms / unplaced p-arm contigs, so use DJ_bp_runs_ge50kb for DJ.

Outputs for a run with prefix P:
  P.hits.tsv.gz   one row per (tile, class) with any alignment: seq, contig, tile_start,
                  tile_end (0-based half-open contig coords), tile_len, n_frac, cls,
                  aln_frac, best_identity, n_aln, pass
  P.summary.json  tiles considered / skipped, bp tiled, passing bp per class

Usage:
  python scripts/tile_map.py --fasta in.fa.gz [--fasta ...] --out P [--threads 8]
"""
import argparse
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile

import pandas as pd
import pysam

TILE = 5000
MIN_TILE = 500
MAX_N_FRAC = 0.5
MIN_IDENT = 0.90
MIN_ALN_FRAC = 0.50
MM2_ARGS = ["-x", "asm20", "-c", "--secondary=no"]
UNIT_45S = "/Users/Kitty/git/NGS-DOSE/work/ref/KY962518.1.fa"
UNIT_5S = "/Users/Kitty/git/NGS-DOSE/work/ref/X12811.1.fa"
DJ_CORE = "/Users/Kitty/git/NGS-DOSE/work/ref/DJ.chm13_chr21.fa"
HERE = os.path.dirname(os.path.abspath(__file__))
TARGET_DIR = os.path.join(os.path.dirname(HERE), "data", "targets")
REGION_RE = re.compile(r"^(.+):(\d+)-(\d+)$")


def _seq(path):
    with pysam.FastxFile(path) as fh:
        r = next(iter(fh))
        return r.sequence.upper()


def build_targets(outdir=TARGET_DIR):
    os.makedirs(outdir, exist_ok=True)
    a = os.path.join(outdir, "target_45S_2x.fa")
    b = os.path.join(outdir, "target_5S2x_DJ.fa")
    if not os.path.exists(a):
        s = _seq(UNIT_45S)
        with open(a, "w") as fh:
            fh.write(f">45S_KY962518.1_2x\n{s + s}\n")
    if not os.path.exists(b):
        s5, dj = _seq(UNIT_5S), _seq(DJ_CORE)
        with open(b, "w") as fh:
            fh.write(f">5S_X12811.1_2x\n{s5 + s5}\n>DJ_chm13_chr21_2708299_3108298\n{dj}\n")
    return {"45S": a, "5S_DJ": b}


def target_class(tname):
    return "45S" if tname.startswith("45S") else ("5S" if tname.startswith("5S") else "DJ")


def placement(name):
    """Record name -> (contig, 0-based contig offset of record start)."""
    m = REGION_RE.match(name)
    if m:
        return m.group(1), int(m.group(2)) - 1
    return name, 0


def write_tiles(fastas, tile_fa):
    """Write grid tiles to tile_fa. Returns (tile table, stats)."""
    meta = []
    st = {"records": 0, "bp_input": 0, "tiles_total": 0, "tiles_skipped_N": 0,
          "tiles_skipped_short": 0, "bp_tiled": 0}
    with open(tile_fa, "w") as out:
        for fa in fastas:
            with pysam.FastxFile(fa) as fh:
                for rec in fh:
                    contig, off = placement(rec.name)
                    seq = rec.sequence
                    L = len(seq)
                    st["records"] += 1
                    st["bp_input"] += L
                    g0 = off
                    g1 = off + L
                    k = (g0 // TILE) * TILE
                    while k < g1:
                        s, e = max(k, g0), min(k + TILE, g1)
                        k += TILE
                        tl = e - s
                        st["tiles_total"] += 1
                        if tl < MIN_TILE:
                            st["tiles_skipped_short"] += 1
                            continue
                        t = seq[s - off:e - off]
                        nf = (t.count("N") + t.count("n")) / tl
                        if nf > MAX_N_FRAC:
                            st["tiles_skipped_N"] += 1
                            continue
                        tid = f"{contig}|{s}|{e}"
                        out.write(f">{tid}\n{t}\n")
                        meta.append((tid, rec.name, contig, s, e, tl, round(nf, 4)))
                        st["bp_tiled"] += tl
    df = pd.DataFrame(meta, columns=["tile", "seq", "contig", "tile_start", "tile_end",
                                     "tile_len", "n_frac"])
    return df, st


def run_minimap2(target, tile_fa, paf, threads):
    cmd = ["minimap2", *MM2_ARGS, "-t", str(threads), target, tile_fa]
    with open(paf, "w") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)


def _union(iv):
    iv = sorted(iv)
    tot, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if cs is not None:
        tot += ce - cs
    return tot


def score_paf(paf):
    """PAF -> per (tile, class): aln_frac, best_identity, n_aln."""
    agg = {}
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            q, qlen, qs, qe, tn = f[0], int(f[1]), int(f[2]), int(f[3]), f[5]
            ident = int(f[9]) / int(f[10])
            key = (q, target_class(tn))
            a = agg.setdefault(key, {"qlen": qlen, "iv": [], "best": 0.0, "n": 0})
            a["n"] += 1
            a["best"] = max(a["best"], ident)
            if ident >= MIN_IDENT:
                a["iv"].append((qs, qe))
    rows = []
    for (q, c), a in agg.items():
        af = _union(a["iv"]) / a["qlen"]
        rows.append((q, c, round(af, 4), round(a["best"], 4), a["n"], af >= MIN_ALN_FRAC))
    return pd.DataFrame(rows, columns=["tile", "cls", "aln_frac", "best_identity", "n_aln", "pass"])


RUN_MIN = 50_000


def class_runs(hits, cls, min_len=RUN_MIN, max_gap=TILE):
    """Runs of passing tiles of one class: list of (contig, start, end, n_tiles, bp)."""
    x = hits[(hits["cls"] == cls) & (hits["pass"].astype(bool))]
    out = []
    for c, g in x.groupby("contig"):
        g = g.sort_values("tile_start")
        cur = None
        for s, e, l in zip(g.tile_start, g.tile_end, g.tile_len):
            if cur is None or s > cur[1] + max_gap:
                if cur is not None:
                    out.append((c, *cur))
                cur = [s, e, 1, l]
            else:
                cur[1] = max(cur[1], e)
                cur[2] += 1
                cur[3] += l
        if cur is not None:
            out.append((c, *cur))
    return [r for r in out if r[2] - r[1] >= min_len]


def run(fastas, out_prefix, threads=8, keep_tiles=False):
    if isinstance(fastas, str):
        fastas = [fastas]
    targets = build_targets()
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(out_prefix))) as td:
        tile_fa = os.path.join(td, "tiles.fa")
        meta, st = write_tiles(fastas, tile_fa)
        hits = []
        for tk, tp in targets.items():
            paf = os.path.join(td, f"{tk}.paf")
            run_minimap2(tp, tile_fa, paf, threads)
            hits.append(score_paf(paf))
        if keep_tiles:
            os.replace(tile_fa, out_prefix + ".tiles.fa")
    hits = pd.concat(hits, ignore_index=True)
    hits = meta.merge(hits, on="tile", how="inner").drop(columns="tile")
    hits = hits.sort_values(["contig", "tile_start", "cls"])
    hits.to_csv(out_prefix + ".hits.tsv.gz", sep="\t", index=False, compression="gzip")
    hits["pass"] = hits["pass"].astype(bool)
    passing = hits.loc[hits["pass"]]
    st["params"] = {"tile": TILE, "run_min": RUN_MIN, "min_tile": MIN_TILE, "max_n_frac": MAX_N_FRAC,
                    "min_identity": MIN_IDENT, "min_aln_frac": MIN_ALN_FRAC,
                    "minimap2": " ".join(MM2_ARGS), "targets": targets}
    st["inputs"] = fastas
    for c in ("45S", "5S", "DJ"):
        p = passing[passing.cls == c]
        st[f"{c}_tiles"] = int(len(p))
        st[f"{c}_bp"] = int(p.tile_len.sum())
        rr = class_runs(hits, c)
        st[f"{c}_n_runs_ge50kb"] = len(rr)
        st[f"{c}_bp_runs_ge50kb"] = int(sum(r[4] for r in rr))
    with open(out_prefix + ".summary.json", "w") as fh:
        json.dump(st, fh, indent=1)
    return hits, st


def load_hits(out_prefix, passing_only=True):
    h = pd.read_csv(out_prefix + ".hits.tsv.gz", sep="\t")
    h["pass"] = h["pass"].astype(bool)
    return h.loc[h["pass"]] if passing_only else h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", action="append", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    _, st = run(a.fasta, a.out, a.threads)
    print(json.dumps({k: v for k, v in st.items() if k not in ("params", "inputs")}))


if __name__ == "__main__":
    main()
