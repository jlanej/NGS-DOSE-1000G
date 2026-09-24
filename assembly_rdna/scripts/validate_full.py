#!/usr/bin/env python
"""Step 4: validate the targeted extraction on whole haplotype assemblies.

For each validation assembly:
  1. download the whole bgzipped FASTA to data/full/ (deleted afterwards unless --keep);
  2. apply the shared measure (tile_map.run) to the whole assembly and, separately, to the
     extracted sequence data/seq/<name>.fa.gz  -> data/validation/<name>.{full,extracted}.*
  3. compare: class bp genome-wide, share inside the extracted regions (a tile counts as
     inside when >= 50% of it lies in data/regions/<name>.bed), every passing tile outside
     the regions (with CenSat label and chromosome assignment at that position), and 45S
     tiles inside vs outside CenSat rDNA intervals.

Writes tables/validation.tsv and tables/validation_outside_hits.tsv.

Usage: python scripts/validate_full.py [--names N1 N2 ...] [--threads 8] [--keep] [--redo-extracted]
"""
import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, META, TABLES, http_get, is_rdna_label, read_censat, s3_url  # noqa
import tile_map  # noqa

DEFAULT = [  # method / version / phasing
    "HG00642_mat_hprc_r2_v1.0.1",   # hifiasm 0.19.7 trio
    "HG01433_pat_hprc_r2_v1.0.1",   # hifiasm 0.19.7 trio
    "HG00235_hap2_hprc_r2_v1.0.1",  # hifiasm 0.19.9 hic
    "HG00097_hap1_hprc_r2_v1.0.1",  # hifiasm 0.19.9 hic
    "NA20762_hap1_hprc_r2_v1.0.1",  # verkko 2.2.1 hic
    "hg002v1.1.pat",                # verkko 1.1 trio (Q100 curated HG002)
    "HG01358_mat_hprc_r2_v1.0.1",   # hifiasm 0.19.9 trio
    "HG00140_hap2_hprc_r2_v1.0.1",  # hifiasm 0.19.9 hic
]
FULL = os.path.join(DATA, "full")
VAL = os.path.join(DATA, "validation")
ARRAY_MIN_RUN = 50_000   # >= 50 kb of consecutive 45S tiles (gap <= 1 tile) = array-like


def intervals_by_contig(rows):
    d = {}
    for c, s, e in rows:
        d.setdefault(c, []).append((s, e))
    return {c: (np.array(sorted(v))) for c, v in d.items()}


def overlap_bp(iv, s, e):
    if iv is None or len(iv) == 0:
        return 0
    m = (iv[:, 0] < e) & (iv[:, 1] > s)
    if not m.any():
        return 0
    x = iv[m]
    return int((np.minimum(x[:, 1], e) - np.maximum(x[:, 0], s)).sum())


def labels_at(cs_by, c, s, e):
    rows = cs_by.get(c, [])
    labs = sorted({lab for (a, b, lab) in rows if a < e and b > s})
    return ",".join(labs) if labs else "none"


def download(name, url):
    dest = os.path.join(FULL, f"{name}.fa.gz")
    if not os.path.exists(dest):
        t0 = time.time()
        n = http_get(url, stream_to=dest, timeout=600)
        print(f"downloaded {name}: {n / 1e9:.2f} GB in {time.time() - t0:.0f}s", flush=True)
    return dest


def validate(name, man, threads, redo_extracted=False):
    m = man.loc[name]
    full = os.path.join(FULL, f"{name}.fa.gz")
    pf, pe = os.path.join(VAL, f"{name}.full"), os.path.join(VAL, f"{name}.extracted")
    if not os.path.exists(pf + ".hits.tsv.gz"):
        tile_map.run([full], pf, threads)
    if redo_extracted or not os.path.exists(pe + ".hits.tsv.gz"):
        tile_map.run([os.path.join(DATA, "seq", f"{name}.fa.gz")], pe, threads)
    hf, he = tile_map.load_hits(pf), tile_map.load_hits(pe)
    reg = pd.read_csv(os.path.join(DATA, "regions", f"{name}.bed"), sep="\t", header=None,
                      names=["contig", "start", "end", "reason"])
    reg_by = intervals_by_contig(zip(reg.contig, reg.start, reg.end))
    cs = read_censat(m["censat"]) if isinstance(m["censat"], str) else []
    cs_by = {}
    for c, s, e, lab in cs:
        cs_by.setdefault(c, []).append((s, e, lab))
    rd_by = intervals_by_contig([(c, s, e) for c, s, e, lab in cs if is_rdna_label(lab)])
    alias = {}
    if isinstance(m["alias"], str):
        for line in open(m["alias"]):
            if not line.startswith("#"):
                f = line.rstrip("\n").split("\t")
                alias[f[0]] = f[1]
    hf = hf.copy()
    hf["in_region_bp"] = [overlap_bp(reg_by.get(c), s, e) for c, s, e in zip(hf.contig, hf.tile_start, hf.tile_end)]
    hf["inside"] = hf.in_region_bp >= 0.5 * hf.tile_len
    hf["in_censat_rdna_bp"] = [overlap_bp(rd_by.get(c), s, e) for c, s, e in zip(hf.contig, hf.tile_start, hf.tile_end)]
    row = {"assembly_name": name, "censat_available": bool(cs)}
    for cls in ("45S", "5S", "DJ"):
        x = hf[hf.cls == cls]
        y = he[he.cls == cls]
        row[f"{cls}_bp_full"] = int(x.tile_len.sum())
        row[f"{cls}_bp_full_inside_regions"] = int(x[x.inside].tile_len.sum())
        row[f"{cls}_share_inside"] = round(x[x.inside].tile_len.sum() / x.tile_len.sum(), 5) if len(x) else np.nan
        row[f"{cls}_bp_extracted_measure"] = int(y.tile_len.sum())
        row[f"{cls}_extracted_over_full"] = round(y.tile_len.sum() / x.tile_len.sum(), 5) if len(x) else np.nan
        row[f"{cls}_n_tiles_outside"] = int((~x.inside).sum())
        row[f"{cls}_bp_outside"] = int(x[~x.inside].tile_len.sum())
    x = hf[hf.cls == "45S"]
    inc = x.in_censat_rdna_bp > 0
    row["45S_bp_in_censat_rdna"] = int(x[inc].tile_len.sum())
    row["45S_bp_outside_censat_rdna"] = int(x[~inc].tile_len.sum())
    row["45S_n_tiles_outside_censat_rdna"] = int((~inc).sum())
    allr = tile_map.class_runs(hf, "45S", min_len=0)
    arr = [r for r in allr if r[2] - r[1] >= ARRAY_MIN_RUN]
    row["45S_bp_array_like_runs"] = int(sum(r[4] for r in arr))
    row["45S_n_array_like_runs"] = len(arr)
    row["45S_n_isolated_runs"] = len(allr) - len(arr)
    labs = [labels_at(cs_by, c, s, e) for c, s, e in zip(x.contig, x.tile_start, x.tile_end)]
    row["45S_bp_no_censat_label"] = int(x.tile_len[[l == "none" for l in labs]].sum())
    # DJ: per-tile class is unspecific; runs >= 50 kb are the DJ measure
    for tag, hh in (("full", hf), ("extracted", he.assign(**{"pass": True}))):
        rr = tile_map.class_runs(hh, "DJ")
        row[f"DJ_n_runs50_{tag}"] = len(rr)
        row[f"DJ_bp_runs50_{tag}"] = int(sum(r[4] for r in rr))
        if tag == "full":
            ins = [r for r in rr if overlap_bp(reg_by.get(r[0]), r[1], r[2]) >= 0.5 * (r[2] - r[1])]
            row["DJ_bp_runs50_full_inside_regions"] = int(sum(r[4] for r in ins))
            row["DJ_runs50_outside"] = ";".join(f"{r[0]}:{r[1]}-{r[2]}({alias.get(r[0], 'NA')})"
                                                for r in rr if r not in ins)
            djrun = intervals_by_contig([(r[0], r[1], r[2]) for r in rr])
    row["DJ_runs50_extracted_over_full"] = round(row["DJ_bp_runs50_extracted"] / row["DJ_bp_runs50_full"], 5) \
        if row["DJ_bp_runs50_full"] else np.nan
    # CenSat rDNA label vs 45S tiles
    rdi = [(c, s, e) for c, s, e, lab in cs if is_rdna_label(lab)]
    tiles45 = intervals_by_contig(zip(x.contig, x.tile_start, x.tile_end))
    cov = [overlap_bp(tiles45.get(c), s, e) for c, s, e in rdi]
    L = np.array([e - s for c, s, e in rdi]) if rdi else np.array([])
    row["censat_rdna_n"] = len(rdi)
    row["censat_rdna_bp"] = int(L.sum())
    small = L < 10_000
    row["censat_rdna_n_lt10kb"] = int(small.sum())
    row["censat_rdna_bp_lt10kb"] = int(L[small].sum())
    row["censat_rdna_bp_ge10kb"] = int(L[~small].sum())
    covd = np.array(cov) if cov else np.array([])
    row["censat_rdna_bp_covered_by_45S_tiles"] = int(covd.sum()) if len(covd) else 0
    row["censat_rdna_lt10kb_n_with_45S_tile"] = int((covd[small] > 0).sum()) if len(covd) else 0
    row["censat_rdna_ge10kb_share_covered"] = round(covd[~small].sum() / L[~small].sum(), 4) if (~small).any() else np.nan
    # outside hits
    out = hf[~hf.inside].copy()
    out["assembly_name"] = name
    out["censat_label"] = [labels_at(cs_by, c, s, e) for c, s, e in zip(out.contig, out.tile_start, out.tile_end)]
    out["chrom_assignment"] = out.contig.map(alias).fillna("NA")
    out["in_censat_rdna"] = out.in_censat_rdna_bp > 0
    out["in_class_run_ge50kb"] = [overlap_bp(djrun.get(c), s, e) > 0 if k == "DJ" else None
                                  for k, c, s, e in zip(out.cls, out.contig, out.tile_start, out.tile_end)]
    out = out[["assembly_name", "cls", "contig", "chrom_assignment", "tile_start", "tile_end", "tile_len",
               "n_frac", "aln_frac", "best_identity", "censat_label", "in_censat_rdna", "in_class_run_ge50kb", "in_region_bp"]]
    hf.to_csv(pf + ".hits.annotated.tsv.gz", sep="\t", index=False, compression="gzip")
    return row, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="*", default=DEFAULT)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--redo-extracted", action="store_true",
                    help="re-measure the extracted sequence (after re-running extract_regions.py)")
    a = ap.parse_args()
    os.makedirs(FULL, exist_ok=True)
    os.makedirs(VAL, exist_ok=True)
    idx = pd.read_csv(os.path.join(DATA, "index.csv")).set_index("assembly_name")
    man = pd.read_csv(os.path.join(META, "manifest.tsv"), sep="\t").set_index("assembly_name")
    asm = pd.read_csv(os.path.join(TABLES, "assemblies.tsv"), sep="\t").set_index("assembly_name")
    todo = [n for n in a.names if not os.path.exists(os.path.join(VAL, f"{n}.full.hits.tsv.gz"))]
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda n: download(n, s3_url(idx.loc[n, "assembly"])), todo))
    rows, outs = [], []
    for n in a.names:
        t0 = time.time()
        r, o = validate(n, man, a.threads, a.redo_extracted)
        for k in ("sample", "haplotype", "method", "version", "phasing", "total_length", "n_contigs"):
            r[k] = asm.loc[n, k]
        r["full_fasta_bytes"] = os.path.getsize(os.path.join(FULL, f"{n}.fa.gz")) \
            if os.path.exists(os.path.join(FULL, f"{n}.fa.gz")) else np.nan
        rows.append(r)
        outs.append(o)
        print(f"validated {n} in {time.time() - t0:.0f}s", flush=True)
    v = pd.DataFrame(rows)
    front = ["assembly_name", "sample", "haplotype", "method", "version", "phasing", "total_length", "n_contigs"]
    v = v[front + [c for c in v.columns if c not in front]]
    v.to_csv(os.path.join(TABLES, "validation.tsv"), sep="\t", index=False)
    pd.concat(outs).to_csv(os.path.join(TABLES, "validation_outside_hits.tsv"), sep="\t", index=False)
    if not a.keep:
        for n in a.names:
            p = os.path.join(FULL, f"{n}.fa.gz")
            if os.path.exists(p):
                os.remove(p)
    print(v.T.to_string())


if __name__ == "__main__":
    main()
