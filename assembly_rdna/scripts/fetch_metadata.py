#!/usr/bin/env python
"""Step 2: fetch per-assembly metadata for every row of the HPRC r2 index and build
tables/assemblies.tsv.

For each of the 466 index rows (464 HPRC/HPP haplotypes + GRCh38 + CHM13 reference rows):
  data/meta/fai/<name>.fa.gz.fai, data/meta/gzi/<name>.fa.gz.gzi
  data/meta/chrom_assignment/<name>.{chromAlias.txt,gaps.bed,t2t_chromosomes.tsv}
  data/meta/chains/<name>_vs_GRCh38.chain.gz (+ <name>.chainsum.tsv summary)
  CenSat: reused from hprc_censat/ (read-only) when present, else data/censat/<file>.gz
  data/meta/manifest.tsv : chosen S3 key / local path / status per file type

Usage: python scripts/fetch_metadata.py [--threads 8]
"""
import argparse
import collections
import gzip
import json
import os
import re
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (ACRO, CHM13_CENSAT, DATA, HG002_CENSAT, LOCAL_CENSAT_DIR, META,  # noqa
                    TABLES, hap_label, http_get, is_rdna_label, open_text, parse_chain,
                    read_censat, read_fai, s3_list, s3_url, session)

GRCH38_LOCAL = "/Users/Kitty/git/NGS-DOSE/work/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
SUFFIX = {"censat": ".cenSat.bed", "alias": ".chromAlias.txt", "gaps": ".gaps.bed",
          "t2t": ".t2t_chromosomes.tsv", "chain": "_vs_GRCh38.chain.gz"}
DIRS = {k: os.path.join(META, d) for k, d in
        [("fai", "fai"), ("gzi", "gzi"), ("alias", "chrom_assignment"),
         ("gaps", "chrom_assignment"), ("t2t", "chrom_assignment"), ("chain", "chains")]}
BYTES = collections.Counter()


def key_of(uri):
    return uri.replace("s3://human-pangenomics/", "")


def stems(name):
    s = [name]
    short = re.sub(r"(_v\d+)(\.\d+)+$", r"\1", name)
    if short != name:
        s.append(short)
    return s


def pick(listing, name, suffix):
    for st in stems(name):
        c = [k for k, _ in listing if os.path.basename(k).startswith(st)
             and os.path.basename(k).endswith(suffix)
             and os.path.basename(k)[len(st):len(st) + 1] in (".", "_")]
        if c:
            return sorted(c)[-1]
    return None


def fetch(url, dest, gz=False):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "cached"
    data = http_get(url)
    BYTES["meta"] += len(data)
    tmp = dest + ".part"
    if gz:
        with gzip.open(tmp, "wb") as fh:
            fh.write(data)
    else:
        with open(tmp, "wb") as fh:
            fh.write(data)
    os.replace(tmp, dest)
    return "fetched"


def chain_summary(chain_path, prefix, out):
    """Per (contig, GRCh38 chr, strand): aligned bp and GRCh38 (+strand) span."""
    agg = {}
    for ch in parse_chain(chain_path, keep_q=set()):  # no blocks kept, only totals
        if ch["qStrand"] == "+":
            qs, qe = ch["qStart"], ch["qEnd"]
        else:
            qs, qe = ch["qSize"] - ch["qEnd"], ch["qSize"] - ch["qStart"]
        k = (prefix + ch["tName"], ch["qName"], ch["qStrand"])
        a = agg.setdefault(k, [0, qs, qe, ch["tStart"], ch["tEnd"]])
        a[0] += ch.get("aligned", 0)
        a[1], a[2] = min(a[1], qs), max(a[2], qe)
        a[3], a[4] = min(a[3], ch["tStart"]), max(a[4], ch["tEnd"])
    rows = [(*k, *v) for k, v in agg.items()]
    df = pd.DataFrame(rows, columns=["contig", "grch38_chr", "strand", "aligned_bp",
                                     "grch38_min", "grch38_max", "contig_min", "contig_max"])
    df.sort_values(["contig", "aligned_bp"], ascending=[True, False]).to_csv(out, sep="\t", index=False)


def ncbi_chrom_titles(accs):
    """Chromosome names for GenBank CM accessions via NCBI esummary."""
    out = {}
    r = session().get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                      params={"db": "nuccore", "id": ",".join(accs), "retmode": "json"}, timeout=60)
    r.raise_for_status()
    res = r.json()["result"]
    for uid in res["uids"]:
        e = res[uid]
        m = re.search(r"chromosome (\w+)", e["title"])
        out[e["accessionversion"]] = ("chr" + m.group(1)) if m else "chrUn"
    return out


def grch38_gaps(names, dest):
    """N-runs (>=1 bp) on the given GRCh38 chromosomes, from the local GRCh38 FASTA."""
    import pysam
    fa = pysam.FastaFile(GRCH38_LOCAL)
    with open(dest, "w") as fh:
        for pn in names:
            chrom = pn.split("#")[-1]
            if chrom not in fa.references:
                continue
            seq = fa.fetch(chrom).upper()
            for m in re.finditer(r"N+", seq):
                fh.write(f"{pn}\t{m.start()}\t{m.end()}\n")


def process(row, listings):
    name = row.assembly_name
    rec = {"assembly_name": name, "notes": []}
    asm_key = key_of(row.assembly)
    listing = listings[os.path.dirname(asm_key)]
    ann = [(k, s) for k, s in listing if "/annotation/" in k]
    # fai / gzi
    for kind, col in (("fai", "assembly_fai"), ("gzi", "assembly_gzi")):
        dest = os.path.join(DIRS[kind], f"{name}.fa.gz.{kind}")
        try:
            fetch(s3_url(getattr(row, col)), dest)
            rec[kind] = dest
        except Exception as e:  # noqa
            rec[kind] = None
            rec["notes"].append(f"{kind} missing ({e.__class__.__name__})")
    fai = read_fai(rec["fai"]) if rec.get("fai") else {}
    prefix = ""
    if fai:
        c0 = next(iter(fai))
        prefix = c0.rsplit("#", 1)[0] + "#" if "#" in c0 else ""
    rec["pansn_prefix"] = prefix
    # annotation files
    for kind in ("alias", "gaps", "t2t", "chain", "censat"):
        k = pick(ann, name, SUFFIX[kind])
        rec[kind + "_key"] = k
        if kind == "censat":
            continue
        if k is None:
            rec[kind] = None
            continue
        dest = os.path.join(DIRS[kind], os.path.basename(k))
        try:
            fetch(s3_url(k), dest)
            rec[kind] = dest
        except Exception as e:  # noqa
            rec[kind] = None
            rec["notes"].append(f"{kind} fetch failed ({e})")
    # CenSat: local reuse, S3 fetch, or documented substitute for the extramural rows
    rec["censat"], rec["censat_source"] = None, None
    for st in stems(name):
        p = os.path.join(LOCAL_CENSAT_DIR, f"{st}.cenSat.bed.gz")
        if os.path.exists(p):
            rec["censat"], rec["censat_source"] = p, "local:hprc_censat"
            break
    if rec["censat"] is None and rec["censat_key"]:
        dest = os.path.join(DATA, "censat", os.path.basename(rec["censat_key"]) + ".gz")
        fetch(s3_url(rec["censat_key"]), dest, gz=True)
        rec["censat"], rec["censat_source"] = dest, "s3:" + rec["censat_key"]
    if rec["censat"] is None and fai:
        dest = os.path.join(DATA, "censat", f"{name}.cenSat.substitute.bed.gz")
        src = None
        if name in HG002_CENSAT:
            src = "s3:" + HG002_CENSAT[name]
            txt = http_get(s3_url(HG002_CENSAT[name])).decode()
            BYTES["meta"] += len(txt)
            ren = lambda c: prefix + re.sub(r"_(PATERNAL|MATERNAL)$", "", c)  # noqa
        elif name.startswith("chm13"):
            src = "local:" + CHM13_CENSAT
            txt = open(CHM13_CENSAT).read()
            ren = lambda c: prefix + c  # noqa
        if src:
            with gzip.open(dest, "wt") as fh:
                for line in txt.splitlines():
                    if line.startswith("track") or not line.strip():
                        continue
                    f = line.split("\t")
                    f[0] = ren(f[0])
                    fh.write("\t".join(f) + "\n")
            rec["censat"], rec["censat_source"] = dest, "substitute:" + src
            rec["notes"].append(f"CenSat not in HPRC annotation; substituted {src}")
    if rec["censat"] is None:
        rec["notes"].append("CenSat missing")
    # chromosome assignment fallback for rows without chromAlias
    if rec.get("alias") is None and fai:
        dest = os.path.join(DIRS["alias"], f"{name}.chromAlias.derived.txt")
        if all(re.search(r"#chr[0-9XYM]+$", c) or re.search(r"#chr", c) for c in fai):
            amap = {c: c.split("#")[-1] for c in fai}
            how = "derived from PanSN contig names"
        else:
            accs = [c.split("#")[-1] for c in fai]
            t = ncbi_chrom_titles(accs)
            amap = {c: t.get(c.split("#")[-1], "chrUn") for c in fai}
            how = "derived from NCBI nuccore titles"
        with open(dest, "w") as fh:
            fh.write("# assembly\tucsc\tsource\n")
            for c, u in amap.items():
                fh.write(f"{c}\t{u}\t{how}\n")
        rec["alias"] = dest
        rec["notes"].append(f"chromAlias missing; {how}")
    if rec.get("gaps") is None and fai and name.startswith("GCA_000001405"):
        dest = os.path.join(DIRS["gaps"], f"{name}.gaps.derived.bed")
        if not os.path.exists(dest):
            grch38_gaps([c for c in fai if c.split("#")[-1] in ACRO], dest)
        rec["gaps"] = dest
        rec["notes"].append("gaps.bed missing; N-runs on acrocentrics computed from local GRCh38 FASTA")
    elif rec.get("gaps") is None:
        rec["notes"].append("gaps.bed missing")
    if rec.get("t2t") is None:
        rec["notes"].append("t2t_chromosomes.tsv missing")
    if rec.get("chain"):
        out = os.path.join(DIRS["chain"], f"{name}.chainsum.tsv")
        if not os.path.exists(out):
            chain_summary(rec["chain"], prefix, out)
        rec["chainsum"] = out
    else:
        rec["chainsum"] = None
        rec["notes"].append("GRCh38 chain missing")
    return rec


def summarise(row, rec):
    name = row.assembly_name
    out = dict(sample=row.sample_id, haplotype=hap_label(name), hap_index=row.haplotype,
               assembly_name=name, method=row.assembly_method,
               version=row.assembly_method_version, phasing=row.phasing, source=row.source)
    fai = read_fai(rec["fai"]) if rec.get("fai") else {}
    out["total_length"] = sum(fai.values()) if fai else None
    out["n_contigs"] = len(fai) if fai else None
    alias = {}
    if rec.get("alias"):
        for line in open(rec["alias"]):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            alias[f[0]] = f[1]
    t2t = {}
    if rec.get("t2t"):
        tt = pd.read_csv(rec["t2t"], sep="\t")
        t2t = dict(zip(tt.sequence_id, tt.level))
    gaps = collections.Counter()
    if rec.get("gaps"):
        for line in open(rec["gaps"]):
            f = line.split("\t")
            if len(f) >= 3:
                gaps[f[0]] += 1
    cs = read_censat(rec["censat"]) if rec.get("censat") else []
    rd = [r for r in cs if is_rdna_label(r[3])]
    out["censat_rdna_n"] = len(rd) if rec.get("censat") else None
    out["censat_rdna_bp"] = sum(e - s for _, s, e, _ in rd) if rec.get("censat") else None
    for chrom in ACRO:
        prim = [c for c, u in alias.items() if u == chrom]
        rand = [c for c, u in alias.items() if u.startswith(chrom + "_") and u.endswith("_random")]
        out[f"{chrom}_contigs"] = ",".join(prim) if prim else ""
        out[f"{chrom}_t2t"] = ",".join(str(t2t.get(c, "no")) for c in prim) if rec.get("t2t") else "NA"
        out[f"{chrom}_gaps"] = ",".join(str(gaps.get(c, 0)) for c in prim) if rec.get("gaps") else "NA"
        out[f"{chrom}_n_random"] = len(rand)
        out[f"{chrom}_rdna_bp"] = sum(e - s for c, s, e, _ in rd if c in set(prim) | set(rand))
    unassigned = {c for c, u in alias.items() if u.startswith("chrUn")}
    out["unplaced_rdna_bp"] = sum(e - s for c, s, e, _ in rd if c in unassigned)
    out["censat_source"] = rec.get("censat_source")
    missing = [k for k in ("fai", "gzi", "alias_key", "gaps_key", "t2t_key", "chain", "censat_key")
               if not rec.get(k)]
    out["missing_files"] = ",".join(m.replace("_key", "") for m in missing)
    out["notes"] = "; ".join(rec["notes"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    for d in set(DIRS.values()) | {os.path.join(DATA, "censat"), TABLES}:
        os.makedirs(d, exist_ok=True)
    idx_path = os.path.join(DATA, "index.csv")
    idx = pd.read_csv(idx_path)
    t0 = time.time()
    dirs = sorted({os.path.dirname(key_of(u)) for u in idx.assembly})
    listings = {}
    with ThreadPoolExecutor(a.threads) as ex:
        futs = {ex.submit(s3_list, d + "/"): d for d in dirs}
        for f in as_completed(futs):
            listings[futs[f]] = f.result()
    json.dump(listings, open(os.path.join(META, "s3_listing.json"), "w"))
    recs, fails = {}, []
    with ThreadPoolExecutor(a.threads) as ex:
        futs = {ex.submit(process, r, listings): r.assembly_name for r in idx.itertuples()}
        for f in as_completed(futs):
            n = futs[f]
            try:
                recs[n] = f.result()
            except Exception as e:  # noqa
                fails.append(n)
                recs[n] = {"assembly_name": n, "notes": [f"FAILED: {e}"]}
                print("FAIL", n, e, file=sys.stderr)
    pd.DataFrame([{k: (v if not isinstance(v, list) else "; ".join(v)) for k, v in r.items()}
                  for r in recs.values()]).to_csv(os.path.join(META, "manifest.tsv"), sep="\t", index=False)
    rows = [summarise(r, recs[r.assembly_name]) for r in idx.itertuples()]
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(TABLES, "assemblies.tsv"), sep="\t", index=False)
    print(json.dumps({"n_rows": len(tab), "failures": fails, "meta_bytes_fetched": BYTES["meta"],
                      "seconds": round(time.time() - t0, 1)}))


if __name__ == "__main__":
    main()
