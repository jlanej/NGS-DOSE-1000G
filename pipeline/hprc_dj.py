"""The distal junction in HPRC release-2 assemblies: from the per-haplotype screens to two small tables.

A screen of one haplotype FASTA is two files, both made by 07_hprc_dj.sh:

  <sample>_<hap>.report.tsv   `ngs-dose panel --report`: for every k-mer of the DJ unit, how often it
                              occurs in the assembly (the bundle's DJ panel is the once-per-junction
                              "core" k-mers of the five CHM13 junctions; a complete copy contributes 1)
  <sample>_<hap>.masked.paf   minimap2 (asm20) of the assembly against the unit with everything outside
                              a core k-mer masked to N, so that the repeat elements the unit shares with
                              the rest of the genome seed nothing: each junction copy is a run of
                              alignments on one contig

`python3 pipeline/hprc_dj.py screen DIR -o meta/dj_hprc` writes

  haplotypes.tsv   one row per haplotype: the median core-k-mer count over the unit and per 5-kb
                   sub-block (a haplotype with five complete copies reads 5 in every sub-block; a
                   deletion in one copy reads 4 across it; SNVs only lower single k-mers)
  copies.tsv       one row per junction copy found by alignment: contig, extent of the unit it holds,
                   how far it sits from the contig's ends, and a class (complete; complete with internal
                   gaps; distal-start, for the common copies that begin 22 kb in; truncated at a contig
                   end, which is the assembly's break, not the genome's; partial, an internal breakpoint)

`python3 pipeline/hprc_dj.py masked-unit -o DJ.core_masked.fa` writes the masked unit the screen aligns to.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np

SUB = 5000                     # sub-block of the k-mer table
BLOCK = 20000                  # block of a copy's coverage string
MIN_KMERS = 20                 # a sub-block with fewer core k-mers has no value


def bundle(resources_dir):
    from ngsdose import io, resources
    B = resources.Bundle(resources_dir or resources.default_bundle())
    pc = io.load_panel(B.panel).classes["DJ"]
    unit = B.units()["DJ"]
    return pc.kmer_pos, len(unit), unit, io.load_panel(B.panel).k


def opener(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def read_report(path, core) -> np.ndarray:
    """Background count of every core k-mer, in the panel's order."""
    bg = {}
    with opener(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            c, p, _, _, b = line.rstrip("\n").split("\t")
            if c == "DJ":
                bg[int(p)] = int(b)
    missing = [p for p in core if p not in bg]
    if missing:
        raise SystemExit(f"{path}: {len(missing)} core k-mers are not in the report (is it a report of the bundle's DJ panel?)")
    return np.array([bg[p] for p in core], dtype=np.int32)


def haplotype_row(sample, hap, bgc, core, unit_len) -> dict:
    med = int(np.median(bgc))
    row = dict(sample=sample, haplotype=hap, n_core=len(core), unit_median=med, frac_at_median=round(float(np.mean(bgc == med)), 4),
               frac_above_median=round(float(np.mean(bgc > med)), 4), mean=round(float(bgc.mean()), 3))
    for i, s in enumerate(range(0, unit_len, SUB)):
        m = (core >= s) & (core < s + SUB)
        row[f"b{i}"] = int(np.median(bgc[m])) if m.sum() >= MIN_KMERS else ""
    return row


def copies_from_paf(path, unit_len, cover_ref=None) -> list[dict]:
    """Cluster the masked-unit alignments (>= 2 kb, identity >= 0.9) of each contig into junction copies:
    hits within 300 kb of each other on the contig, without the unit starting over, are one copy; a copy
    holds at least 30 kb of core."""
    hits: dict[str, list] = {}
    for line in open(path):
        p = line.split("\t")
        if int(p[10]) < 2000 or int(p[9]) / int(p[10]) < 0.9:
            continue
        hits.setdefault(p[0], []).append((int(p[2]), int(p[3]), int(p[7]), int(p[8]), int(p[1]), p[4]))
    out = []
    for ctg, h in hits.items():
        h.sort()
        groups, cur = [], [h[0]]
        for x in h[1:]:
            # a new copy: far along the contig, or the unit starting over. Along the contig a forward copy's hits
            # climb the unit and a reverse copy's descend it; a hit on the copy's dominant strand that jumps back
            # by 200 kb or more starts another copy (an inverted palindrome arm within a copy lies on the other strand)
            plus = sum(c[1] - c[0] for c in cur if c[5] == "+"); minus = sum(c[1] - c[0] for c in cur if c[5] == "-")
            dominant = "+" if plus >= minus else "-"
            top, bottom = max(c[3] for c in cur), min(c[2] for c in cur)
            restart = x[5] == dominant and ((dominant == "+" and x[2] < top - 200000) or (dominant == "-" and x[3] > bottom + 200000))
            if x[0] - cur[-1][1] > 300000 or restart:
                groups.append(cur)
                cur = [x]
            else:
                cur.append(x)
        groups.append(cur)
        for g in groups:
            cov = np.zeros(unit_len, bool)
            for qs, qe, ts, te, ql, st in g:
                cov[ts:te] = True
            if cov.sum() < 30000:
                continue
            qs, qe, ql = min(x[0] for x in g), max(x[1] for x in g), g[0][4]
            lo, hi = int(np.argmax(cov)), int(unit_len - np.argmax(cov[::-1]))
            blocks = [float(cov[s:s + BLOCK].mean()) for s in range(0, unit_len, BLOCK)]
            out.append(dict(contig=ctg, contig_len=ql, q_start=qs, q_end=qe, dj_start=lo, dj_end=hi, covered_core_bp=int(cov.sum()),
                            to_contig_start=qs, to_contig_end=ql - qe, strands="".join(sorted({x[5] for x in g})), blocks=blocks))
    return out


def classify(cp: dict, unit_len: int, complete_bp: int) -> str:
    full = cp["dj_start"] <= 10000 and cp["dj_end"] >= unit_len - 10000
    trunc = (cp["to_contig_start"] < 5000 and cp["dj_start"] > 10000) or (cp["to_contig_end"] < 5000 and cp["dj_end"] < unit_len - 10000)
    if full and cp["covered_core_bp"] >= 0.9 * complete_bp:
        return "complete"
    if full:
        return "complete, internal gaps"
    if 15000 <= cp["dj_start"] <= 30000 and cp["dj_end"] >= unit_len - 10000:
        return "distal-start"
    if trunc:
        return "truncated at contig end"
    return "partial"


def write_tsv(rows: list[dict], path: Path):
    cols = list(rows[0].keys()) if rows else []
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else (",".join(f"{x:.2f}" for x in v) if isinstance(v, list) else str(v)) for v in (r[c] for c in cols)) + "\n")


def screen(a):
    core, unit_len, _, _ = bundle(a.resources)
    d = Path(a.dir)
    haps = sorted({f[:-len(".report.tsv")] for f in os.listdir(d) if f.endswith(".report.tsv")})
    if not haps:
        raise SystemExit(f"{d}: no <sample>_<hap>.report.tsv")
    hrows, crows = [], []
    for h in haps:
        sample, hap = h.rsplit("_", 1)
        bgc = read_report(d / f"{h}.report.tsv", core)
        hrows.append(haplotype_row(sample, hap, bgc, core, unit_len))
        paf = d / f"{h}.masked.paf"
        if paf.exists() and paf.stat().st_size:
            cps = copies_from_paf(paf, unit_len)
            complete_bp = max((c["covered_core_bp"] for c in cps), default=0)
            for cp in cps:
                cp["class"] = classify(cp, unit_len, max(complete_bp, 250000))
                crows.append(dict(sample=sample, haplotype=hap, **cp))
        print(f"[hprc_dj] {h}: unit median {hrows[-1]['unit_median']}, {sum(1 for c in crows if c['sample'] == sample and c['haplotype'] == hap)} copies by alignment", file=sys.stderr)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(hrows, out / "haplotypes.tsv")
    write_tsv(crows, out / "copies.tsv")
    print(f"[hprc_dj] {len(hrows)} haplotypes, {len(crows)} copies -> {out}", file=sys.stderr)


def masked_unit(a):
    core, unit_len, unit, k = bundle(a.resources)
    keep = np.zeros(unit_len, bool)
    for p in core:
        keep[p:p + k] = True
    seq = "".join(c if kp else "N" for c, kp in zip(unit, keep))
    with open(a.out, "w") as fh:
        fh.write(">DJ\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
    print(f"[hprc_dj] {a.out}: {unit_len} bp, {int(keep.sum())} bp inside core k-mers kept", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen", help="reports and PAFs of a directory -> haplotypes.tsv, copies.tsv")
    s.add_argument("dir")
    s.add_argument("-o", "--out", required=True)
    s.add_argument("-r", "--resources")
    s.set_defaults(fn=screen)
    m = sub.add_parser("masked-unit", help="the DJ unit with everything outside a core k-mer masked")
    m.add_argument("-o", "--out", required=True)
    m.add_argument("-r", "--resources")
    m.set_defaults(fn=masked_unit)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
