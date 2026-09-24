#!/usr/bin/env python
"""Unit-level annotation of rDNA in the extracted sequence of one assembly.

For data/seq/<name>.fa.gz (records named contig:start-end, 1-based inclusive) this aligns the
45S unit (KY962518.1) and its 18S and 28S genes, the 5S unit (X12811.1) and the distal-junction
core (CHM13 chr21:2,708,299-3,108,298) as queries against the extracted sequence, reporting every
copy (minimap2 -x asm20, all secondary chains kept, repetitive-minimizer filter disabled).

Definitions (all on contig coordinates)
  gene copy   an 18S alignment covering >= 95% of the 18S at >= 95% identity (1,869 bp)
  28S copy    the same for the 28S (5,051 bp)
  full unit   a 45S-unit alignment covering >= 80% of KY962518.1 at >= 95% identity
  5S unit     an X12811.1 alignment covering >= 90% of the 2,231-bp unit at >= 95% identity
  DJ copy-eq  bp of DJ-core alignments at >= 90% identity / 400 kb; DJ copy = alignment >= 200 kb
  array       18S copies on one contig whose successive starts are <= 120 kb apart (a lone copy
              is an array of one)
  array side  'contig_end' if the contig ends within 100 kb beyond the outermost copy, 'gap' if an
              assembly gap (gaps.bed, or any N-run >= 10 bp in the sequence) starts within 100 kb, else 'flank' (non-rDNA sequence follows);
              'DJ' is added when a DJ-core alignment >= 50 kb lies within 1 Mb on that side
  closed      both sides 'flank'
  period      sequence from one 18S start to the next in an array, same strand; its length and
              whether it is byte-identical to another period in the same haplotype

Usage: python scripts/annotate_units.py NAME [NAME ...] [--threads 2]
Writes work_units/<NAME>.paf.gz, work_units/<NAME>.copies.tsv.gz, work_units/<NAME>.arrays.tsv,
work_units/<NAME>.summary.json
"""
import argparse, collections, gzip, hashlib, json, os, re, subprocess, sys
import numpy as np, pandas as pd, pysam

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
SEQ = os.path.join(BASE, "data", "seq")
META = os.path.join(BASE, "data", "meta")
WORK = os.path.join(BASE, "work_units")
UNIT45 = "/Users/Kitty/git/NGS-DOSE/work/ref/KY962518.1.fa"
UNIT5 = "/Users/Kitty/git/NGS-DOSE/work/ref/X12811.1.fa"
DJ = "/Users/Kitty/git/NGS-DOSE/work/ref/DJ.chm13_chr21.fa"
FEAT = {"18S": (3657, 5526), "28S": (7920, 12971)}
MM2 = ["minimap2", "-x", "asm20", "-c", "-N", "200000", "-p", "0.01", "--secondary=yes", "-f", "1000000"]
MAX_STEP, END_WIN, DJ_WIN = 120_000, 100_000, 1_000_000
REGION_RE = re.compile(r"^(.+):(\d+)-(\d+)$")


def one(path):
    with pysam.FastxFile(path) as fh:
        return next(iter(fh)).sequence.upper()


def query_fa():
    q = os.path.join(WORK, "queries.fa")
    if not os.path.exists(q):
        u = one(UNIT45)
        with open(q, "w") as fh:
            fh.write(f">unit45S\n{u}\n")
            for k, (a, b) in FEAT.items():
                fh.write(f">{k}\n{u[a:b]}\n")
            fh.write(f">unit5S\n{one(UNIT5)}\n>DJcore\n{one(DJ)}\n")
    return q


def placement(name):
    m = REGION_RE.match(name)
    return (m.group(1), int(m.group(2)) - 1) if m else (name, 0)


def contig_lengths(name):
    for p in (os.path.join(META, "fai", f"{name}.fa.gz.fai"), os.path.join(META, "fai", f"{name}.fa.fai")):
        if os.path.exists(p):
            return pd.read_csv(p, sep="\t", header=None, usecols=[0, 1], names=["c", "L"]).set_index("c").L.to_dict()
    return {}


def gaps(name):
    d = collections.defaultdict(list)
    ca = os.path.join(META, "chrom_assignment")
    for p in (f"{name}.gaps.bed", f"{name}.gaps.derived.bed"):
        p = os.path.join(ca, p)
        if os.path.exists(p):
            for line in open(p):
                f = line.split("\t")
                if len(f) >= 3:
                    d[f[0]].append((int(f[1]), int(f[2])))
    return d


def run_minimap(name, threads):
    paf = os.path.join(WORK, f"{name}.paf.gz")
    if not os.path.exists(paf):
        fa = os.path.join(SEQ, f"{name}.fa.gz")
        p = subprocess.run(MM2 + ["-t", str(threads), fa, query_fa()], capture_output=True, check=True)
        with gzip.open(paf + ".tmp", "wb") as fh:
            fh.write(p.stdout)
        os.replace(paf + ".tmp", paf)
    cols = ["q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te", "nmatch", "alen", "mapq"]
    rows = [l.split("\t")[:12] for l in gzip.open(paf, "rt") if l.strip()]
    df = pd.DataFrame(rows, columns=cols)
    for c in cols:
        if c not in ("q", "strand", "t"):
            df[c] = df[c].astype(int)
    df["ident"] = df.nmatch / df.alen
    df["qcov"] = (df.qe - df.qs) / df.qlen
    pl = df.t.map(placement)
    df["contig"] = [p[0] for p in pl]
    off = np.array([p[1] for p in pl], dtype=np.int64)
    df["cs"] = df.ts + off
    df["ce"] = df.te + off
    return df


def side(contig, pos, direction, L, gp):
    """What lies beyond an array edge at `pos`, looking in `direction` (-1 left, +1 right)."""
    if direction < 0 and pos <= END_WIN:
        return "contig_end"
    if direction > 0 and L is not None and L - pos <= END_WIN:
        return "contig_end"
    for a, b in gp.get(contig, []):
        if (direction < 0 and pos - END_WIN <= b <= pos + 1) or (direction > 0 and pos - 1 <= a <= pos + END_WIN):
            return "gap"
    return "flank"


def annotate(name, threads=2):
    df = run_minimap(name, threads)
    L = contig_lengths(name)
    gp = gaps(name)
    copies = df[((df.q == "18S") | (df.q == "28S")) & (df.qcov >= 0.95) & (df.ident >= 0.95)].copy()
    copies = pd.concat([copies,
                        df[(df.q == "unit45S") & (df.qcov >= 0.80) & (df.ident >= 0.95)],
                        df[(df.q == "unit5S") & (df.qcov >= 0.90) & (df.ident >= 0.95)]])
    copies = copies.drop_duplicates(["q", "contig", "cs", "ce"]).sort_values(["q", "contig", "cs"])
    dj = df[(df.q == "DJcore") & (df.ident >= 0.90)].drop_duplicates(["contig", "cs", "ce"])
    s18 = copies[copies.q == "18S"]
    arrays = []
    fa = pysam.FastaFile(os.path.join(SEQ, f"{name}.fa.gz"))
    rec_of = collections.defaultdict(list)
    for r in fa.references:
        c, o = placement(r)
        rec_of[c].append((o, o + fa.get_reference_length(r), r))
    def fetch(contig, a, b):
        for o, e, r in rec_of[contig]:
            if o <= a and b <= e:
                return fa.fetch(r, a - o, b - o).upper()
        return None
    for r in fa.references:                                   # N-runs >= 10 bp in the sequence itself
        c, o = placement(r)
        for m in re.finditer(r"[Nn]{10,}", fa.fetch(r)):
            gp[c].append((o + m.start(), o + m.end()))
    periods = []
    for contig, g in s18.groupby("contig"):
        g = g.sort_values("cs")
        starts = g.cs.values
        brk = np.r_[0, np.where(np.diff(starts) > MAX_STEP)[0] + 1, len(starts)]
        for i in range(len(brk) - 1):
            sub = g.iloc[brk[i]:brk[i + 1]]
            a0, a1 = int(sub.cs.min()), int(sub.ce.max())
            Lc = L.get(contig)
            left, right = side(contig, a0, -1, Lc, gp), side(contig, a1, +1, Lc, gp)
            djl = dj[(dj.contig == contig) & (dj.ce <= a0) & (dj.ce >= a0 - DJ_WIN) & (dj.ce - dj.cs >= 50_000)]
            djr = dj[(dj.contig == contig) & (dj.cs >= a1) & (dj.cs <= a1 + DJ_WIN) & (dj.ce - dj.cs >= 50_000)]
            strand = sub.strand.mode().iat[0]
            arr_id = f"{contig}:{a0}"
            ss = sub[sub.strand == strand].cs.values if strand == "+" else sub[sub.strand == strand].ce.values
            ss = np.sort(ss)
            for x, y in zip(ss[:-1], ss[1:]):
                seq = fetch(contig, int(x), int(y))
                periods.append(dict(array=arr_id, length=int(y - x),
                                    md5=hashlib.md5(seq.encode()).hexdigest() if seq and "N" not in seq else None))
            arrays.append(dict(assembly=name, array=arr_id, contig=contig, start=a0, end=a1, span=a1 - a0,
                               n_18S=len(sub), n_28S=int(((copies.q == "28S") & (copies.contig == contig) &
                                                          (copies.cs >= a0 - 15000) & (copies.ce <= a1 + 15000)).sum()),
                               strand=strand, mixed_strand=int((sub.strand != strand).sum()),
                               left=left + ("+DJ" if len(djl) else ""), right=right + ("+DJ" if len(djr) else ""),
                               contig_len=Lc))
    arr = pd.DataFrame(arrays)
    per = pd.DataFrame(periods)
    if len(per):
        dup = per.dropna(subset=["md5"]).md5.duplicated(keep=False)
        per["identical"] = False
        per.loc[dup.index, "identical"] = dup
        if len(arr):
            agg = per.groupby("array").agg(n_periods=("length", "size"), period_median=("length", "median"),
                                           period_min=("length", "min"), period_max=("length", "max"),
                                           n_identical=("identical", "sum"))
            arr = arr.merge(agg, left_on="array", right_index=True, how="left")
    # 5S arrays (same chaining, 5 kb step)
    s5 = copies[copies.q == "unit5S"].sort_values(["contig", "cs"])
    n5_arrays, n5_max, n5_sides = 0, 0, ("none", "none")
    for contig, g in s5.groupby("contig"):
        st, en = g.cs.values, g.ce.values
        brk = np.r_[0, np.where(np.diff(st) > 5000)[0] + 1, len(st)]
        n5_arrays += len(brk) - 1
        for i in range(len(brk) - 1):
            n = int(brk[i + 1] - brk[i])
            if n > n5_max:
                a0, a1 = int(st[brk[i]]), int(en[brk[i + 1] - 1])
                n5_max = n
                n5_sides = (side(contig, a0, -1, L.get(contig), gp), side(contig, a1, +1, L.get(contig), gp))
    dj_eq = float((dj.ce - dj.cs).sum() / 400_000)
    closed = arr[(arr.left.str.startswith("flank")) & (arr.right.str.startswith("flank"))] if len(arr) else arr
    summ = dict(assembly=name, n_18S=int(len(s18)), n_28S=int((copies.q == "28S").sum()),
                n_full_units=int((copies.q == "unit45S").sum()), n_5S_units=int(len(s5)),
                n_5S_arrays=n5_arrays, largest_5S_array=n5_max, largest_5S_left=n5_sides[0], largest_5S_right=n5_sides[1], DJ_copy_eq=round(dj_eq, 3),
                n_DJ_copies=int(((dj.ce - dj.cs) >= 200_000).sum()),
                n_arrays=int(len(arr)), n_arrays_ge2=int((arr.n_18S >= 2).sum()) if len(arr) else 0,
                n_arrays_ge10=int((arr.n_18S >= 10).sum()) if len(arr) else 0,
                largest_array=int(arr.n_18S.max()) if len(arr) else 0,
                n_closed=int(len(closed)), units_in_closed=int(closed.n_18S.sum()) if len(closed) else 0,
                n_closed_DJ=int(closed.left.str.contains("DJ").sum() + closed.right.str.contains("DJ").sum()) if len(closed) else 0,
                n_periods=int(len(per)), n_periods_identical=int(per.identical.sum()) if len(per) else 0,
                period_median=float(per.length.median()) if len(per) else None,
                n_periods_outside_40_50kb=int(((per.length < 40000) | (per.length > 50000)).sum()) if len(per) else 0)
    for s in ("contig_end", "gap", "flank"):
        summ[f"sides_{s}"] = int((arr.left.str.startswith(s)).sum() + (arr.right.str.startswith(s)).sum()) if len(arr) else 0
    copies[["q", "contig", "cs", "ce", "strand", "ident", "qcov"]].to_csv(os.path.join(WORK, f"{name}.copies.tsv.gz"),
                                                                        sep="\t", index=False, compression="gzip")
    arr.to_csv(os.path.join(WORK, f"{name}.arrays.tsv"), sep="\t", index=False)
    json.dump(summ, open(os.path.join(WORK, f"{name}.summary.json"), "w"), indent=1)
    return summ


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    for n in a.names:
        try:
            s = annotate(n, a.threads)
            print(json.dumps(s))
        except Exception as e:
            print(json.dumps(dict(assembly=n, error=repr(e))), file=sys.stderr)
