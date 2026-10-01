"""The calls set against the alleles of the same reads.

Depth says how many copies a chromosome is held in. The alleles say it again, independently: two copies carry the two
alleles of a heterozygous site in equal shares; a copy gained in a share f of the cells makes them (1 + f) : 1, a copy
lost in a share f, 1 : (1 - f). The distance d of the allele fraction from one half is estimated from the spread of the
sites' allele fractions less what each site's depth gives by chance, d^2 = mean (AF - 1/2)^2 - mean 1 / (4 depth), less
the same quantity on the chromosomes of the same genome that depth reads at two copies throughout (the caller's and the
aligner's share). Where two copies are expected (an autosome; an X in a woman), the share of the cells is then
f = 4d / (1 - 2d) for a gain and 4d / (1 + 2d) for a loss. Where one is (a line read as one X with a second in part
of the cells), heterozygous sites exist only in the cells that hold the second, the alleles stand f : 1, and
f = (1 - 2d) / (1 + 2d); the caller misses the sites whose second allele is rarest, so this reads high where f is
small. The agreement is reported over the events where two copies are expected.

The variant calls are made on the reads fetched for the karyotype windows (allele_calls.sh: bcftools mpileup and call on
the regions of the bundle's controls file), one VCF per genome. The events are the page's (docs/data/karyotype_events.tsv).

usage (from the repository's root): python3 analysis/karyotype/allele_balance.py --vcfs DIR [--events docs/data/karyotype_events.tsv]
           [-o meta/karyotype_allele_balance.tsv] [--min-share 0.1]
"""
import argparse
import csv
import glob
import gzip
import os
import sys

import numpy as np

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
MIN_SITES = 30


def sites(vcf):
    """Heterozygous single-nucleotide sites: (chromosome, position, fraction of the reads with the other allele, depth)."""
    out = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p[3]) != 1 or len(p[4]) != 1 or float(p[5]) < 30:
                continue
            f = dict(zip(p[8].split(":"), p[9].split(":")))
            ad = f.get("AD", "").split(",")
            if len(ad) != 2 or f.get("GT") not in ("0/1", "0|1", "1|0"):
                continue
            r, a = int(ad[0]), int(ad[1])
            if 20 <= r + a <= 120:
                out.setdefault(p[0], []).append((int(p[1]), a / (r + a), r + a))
    return {c: np.array(v, float) for c, v in out.items()}


def excess(af, dp, boot=200, seed=1):
    """The sites' spread beyond chance (d^2 before the baseline), and its bootstrap SE."""
    one = lambda a, d: float(np.mean((a - 0.5) ** 2) - np.mean(0.25 / d))
    rng = np.random.default_rng(seed)
    bs = [one(af[k], dp[k]) for k in (rng.integers(0, len(af), len(af)) for _ in range(boot))]
    return one(af, dp), float(np.std(bs))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcfs", required=True, help="directory of SAMPLE.vcf.gz")
    ap.add_argument("--events", default="docs/data/karyotype_events.tsv")
    ap.add_argument("-o", "--out", default="meta/karyotype_allele_balance.tsv")
    ap.add_argument("--min-share", type=float, default=0.1, help="events of at least this share of the cells are set against the alleles")
    a = ap.parse_args()
    events = {}
    with open(a.events) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            events.setdefault(r["sample"], []).append(r)
    rows = []
    for vcf in sorted(glob.glob(os.path.join(a.vcfs, "*.vcf.gz"))):
        s = os.path.basename(vcf)[:-len(".vcf.gz")]
        by = sites(vcf)
        ev = events.get(s, [])
        off = {e["chrom"] for e in ev}
        quiet = [excess(by[c][:, 1], by[c][:, 2])[0] for c in AUTOSOMES if c in by and c not in off and len(by[c]) >= MIN_SITES]
        if len(quiet) < 5:
            print(f"{s}: too few chromosomes at two copies with heterozygous sites; left out", file=sys.stderr)
            continue
        b0 = float(np.median(quiet))
        n_x = len(by.get("chrX", []))
        rows.append(dict(sample=s, event="chromosomes at two copies", chrom="", start="", end="", sites=int(sum(len(by[c]) for c in AUTOSOMES if c in by and c not in off)),
                         d=round(float(np.sqrt(max(b0, 0))), 4), d_se="", by_depth=0.0, by_alleles="", note=f"the genome's baseline; {n_x} heterozygous sites on chrX"))
        for e in ev:
            delta = float(e["delta"])
            c = e["chrom"]
            n0 = int(round(float(e["copies"]) - delta))           # the whole number expected
            if abs(delta) < a.min_share or c == "chrY" or n0 not in (1, 2) or (n0 == 1 and delta < 0):
                continue
            lo, hi = int(e["start"]), int(e["end"])
            v = by.get(c)
            v = v[(v[:, 0] >= lo) & (v[:, 0] < hi)] if v is not None else np.zeros((0, 3))
            row = dict(sample=s, event=e["label"], chrom=c, start=lo, end=hi, sites=len(v), d="", d_se="", by_depth=round(abs(delta), 3), by_alleles="",
                       note="" if n0 == 2 else "one copy expected: a second in part of the cells")
            if len(v) < MIN_SITES:
                row["note"] = "too few heterozygous sites" + (" (one copy carries none)" if float(e["copies"]) < 1.3 else "")
            else:
                x, se = excess(v[:, 1], v[:, 2])
                d = float(np.sqrt(max(x - b0, 0.0)))
                d = min(d, 0.49)
                share = (1 - 2 * d) / (1 + 2 * d) if n0 == 1 else 4 * d / (1 - 2 * d) if delta > 0 else 4 * d / (1 + 2 * d)
                row.update(d=round(d, 4), d_se=round(se / max(2 * d, 1e-3), 4), by_alleles=round(share, 3))
            rows.append(row)
    cols = ["sample", "event", "chrom", "start", "end", "sites", "d", "d_se", "by_depth", "by_alleles", "note"]
    with open(a.out, "w") as fh:
        fh.write("# the calls against the alleles of the same reads (analysis/karyotype/allele_balance.py): by_depth is the share of the cells the level gives, by_alleles the share the allele fractions give\n")
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    both = np.array([(r["by_depth"], r["by_alleles"]) for r in rows if r["chrom"] and r["by_alleles"] != "" and not r["note"]], float)
    print(f"{len({r['sample'] for r in rows})} genomes; {len(both)} events with enough sites: r = {np.corrcoef(both.T)[0, 1]:.3f}, "
          f"robust SD of the difference {1.4826 * np.median(np.abs((both[:, 0] - both[:, 1]) - np.median(both[:, 0] - both[:, 1]))):.3f}; written to {a.out}" if len(both) > 2 else f"{len(rows)} rows")


if __name__ == "__main__":
    main()
