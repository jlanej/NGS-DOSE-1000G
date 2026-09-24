#!/usr/bin/env python
"""Step 3: targeted extraction of rDNA-bearing / acrocentric short-arm / 5S sequence for
every assembly in data/meta/manifest.tsv (written by fetch_metadata.py).

Region classes (4th BED column; ';'-joined after merging overlapping regions)
  rdna_censat   every CenSat interval whose label contains "rDNA" (pure "GAP" excluded),
                merged within 1 kb, padded 200 kb each side, clipped to the contig.
  acro_parm     every contig assigned to chr13/14/15/21/22 (chromosome-level "chrN" or
                unlocalised "chrN_<acc>_random" in chromAlias). The short-arm end is taken
                from the GRCh38 chain strand (majority aligned bp to chrN), else from the
                side nearest the CenSat rDNA, else the contig start. From that end the
                segment runs max(12 Mb, D + 1 Mb) into the contig, D = distance to the far
                edge of the farthest CenSat rDNA interval or gap lying within 25 Mb of that
                end; clipped to the contig. Unlocalised contigs are only taken when they
                can hold short-arm sequence: CenSat rDNA on them, or a chain alignment to
                chrN reaching GRCh38 < centromere end + 2 Mb, or no chain alignment to chrN.
  unplaced_rdna whole contig, for contigs not assigned to a chromosome (chrUn / absent from
                chromAlias) that carry CenSat rDNA, plus every other unassigned contig that
                the shared tile measure (tile_map.py) finds to carry a passing 45S or 5S tile
                or a DJ run >= 50 kb (added after validation showed whole-rDNA chrUn contigs
                with no CenSat label). Assemblies without CenSat: all unassigned contigs are
                screened with tile_map.py (45S tiles).
  5S_locus      GRCh38 flanks chr1:228,400,000-228,600,000 and 228,650,000-228,850,000 are
                lifted through <name>_vs_GRCh38.chain.gz (chain with most aligned bp in the
                flank). Region = inner span between the two lifted flanks + 300 kb each side.
                No chain: all chr1-assigned contigs are tile-scanned (tile_map.py, 5S class)
                and the largest 5S tile cluster +/- 300 kb is taken (noted).

Coordinates: BED 0-based half-open; FASTA record names are samtools region strings
"contig:start-end" (1-based inclusive), i.e. name = f"{contig}:{bed_start+1}-{bed_end}".

Usage: python scripts/extract_regions.py [--threads 8] [--only NAME ...] [--force]
"""
import argparse
import bisect
import json
import os
import struct
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (ACRO, DATA, GRCH38_CEN_END, META, TABLES, is_rdna_label,  # noqa
                    merge_intervals, parse_chain, read_censat, region_name, s3_url)
import tile_map  # noqa

MB = 1_000_000
PAD_RDNA = 200_000
MERGE_RDNA = 1_000
ACRO_MIN = 12 * MB
ACRO_WINDOW = 25 * MB
ACRO_EXTRA = 1 * MB
FLANK_L = (228_400_000, 228_600_000)
FLANK_R = (228_650_000, 228_850_000)
PAD_5S = 300_000
SEQ_DIR = os.path.join(DATA, "seq")
REG_DIR = os.path.join(DATA, "regions")


def read_fai_full(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        d[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return d


def read_gzi(path):
    with open(path, "rb") as fh:
        n = struct.unpack("<Q", fh.read(8))[0]
        arr = struct.unpack(f"<{2 * n}Q", fh.read(16 * n))
    c = [0] + list(arr[0::2])
    u = [0] + list(arr[1::2])
    return c, u


def est_compressed_bytes(regions, fai, gzi, file_size=None):
    """Estimated bytes of bgzf blocks covering the regions (what remote faidx must read)."""
    c, u = gzi
    spans = []
    for contig, s, e, _ in regions:
        L, off, lb, lw = fai[contig]
        u0 = off + (s // lb) * lw + s % lb
        u1 = off + ((e - 1) // lb) * lw + (e - 1) % lb
        i0 = bisect.bisect_right(u, u0) - 1
        i1 = bisect.bisect_right(u, u1)
        c0 = c[i0]
        c1 = c[i1] if i1 < len(c) else (file_size or c[-1] + 65536)
        spans.append((c0, c1))
    spans.sort()
    tot, cs, ce = 0, None, None
    for a, b in spans:
        if cs is None or a > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = a, b
        else:
            ce = max(ce, b)
    if cs is not None:
        tot += ce - cs
    return tot


def load_alias(path):
    al = {}
    if not isinstance(path, str) or not os.path.exists(path):
        return al
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        al[f[0]] = f[1]
    return al


def load_gaps(path):
    g = {}
    if not isinstance(path, str) or not os.path.exists(path):
        return g
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            g.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return g


def acro_of(ucsc):
    """'chr13' -> ('chr13','chromosome'); 'chr13_X_random' -> ('chr13','random')."""
    for c in ACRO:
        if ucsc == c:
            return c, "chromosome"
        if ucsc.startswith(c + "_") and ucsc.endswith("_random"):
            return c, "random"
    return None, None


def lift_flank(chains, flank, prefix):
    """Lift a GRCh38 chr1 interval: per chain, aligned bp and assembly span inside flank."""
    best = None
    for ch in chains:
        if ch["qName"] != "chr1" or not ch["blocks"]:
            continue
        a, b = flank
        if ch["qStrand"] == "-":
            a, b = ch["qSize"] - flank[1], ch["qSize"] - flank[0]
        bp, tmin, tmax = 0, None, None
        for tp, qp, size in ch["blocks"]:
            s, e = max(qp, a), min(qp + size, b)
            if s >= e:
                continue
            bp += e - s
            t0, t1 = tp + (s - qp), tp + (e - qp)
            tmin = t0 if tmin is None else min(tmin, t0)
            tmax = t1 if tmax is None else max(tmax, t1)
        if bp and (best is None or bp > best["bp"]):
            best = dict(contig=prefix + ch["tName"], bp=bp, tmin=tmin, tmax=tmax,
                        strand=ch["qStrand"], chain_id=ch["id"])
    return best


def faidx_to_file(url, fai, gzi, regions, out_fa):
    """Fetch regions (list of (contig,s,e)) from a remote bgzipped FASTA to out_fa (plain)."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as rf:
        for contig, s, e in regions:
            rf.write(region_name(contig, s, e) + "\n")
        rpath = rf.name
    try:
        with open(out_fa, "w") as fh:
            subprocess.run(["samtools", "faidx", "--fai-idx", fai, "--gzi-idx", gzi, url,
                            "-r", rpath], stdout=fh, stderr=subprocess.PIPE, check=True)
    finally:
        os.unlink(rpath)


def plan(m, tmpdir):
    """Return (regions [(contig,s,e,reason)], acro_rows, info)."""
    name = m["assembly_name"]
    fai = read_fai_full(m["fai"])
    L = {k: v[0] for k, v in fai.items()}
    alias = load_alias(m.get("alias"))
    gaps = load_gaps(m.get("gaps"))
    prefix = m.get("pansn_prefix") if isinstance(m.get("pansn_prefix"), str) else ""
    info = {"notes": [], "fallback_bytes": 0}
    regs = []
    censat = read_censat(m["censat"]) if isinstance(m.get("censat"), str) else None
    rd = {}
    if censat is not None:
        for c, s, e, lab in censat:
            if is_rdna_label(lab):
                rd.setdefault(c, []).append((s, e, "rdna_censat"))
    # (a) CenSat rDNA intervals
    for c, iv in rd.items():
        for s, e, _ in merge_intervals(iv, MERGE_RDNA):
            regs.append((c, max(0, s - PAD_RDNA), min(L[c], e + PAD_RDNA), "rdna_censat"))
    # chain summary for orientation
    cs = None
    if isinstance(m.get("chainsum"), str) and os.path.exists(m["chainsum"]):
        cs = pd.read_csv(m["chainsum"], sep="\t")
    # (b) acrocentric short arms
    acro_rows = []
    for c, u in alias.items():
        chrom, level = acro_of(u)
        if chrom is None or c not in L:
            continue
        n = L[c]
        rec = dict(assembly_name=name, contig=c, chrom=chrom, level=level, length=n)
        strand, qmin, qmax, abp = None, None, None, 0
        if cs is not None:
            sub = cs[(cs.contig == c) & (cs.grch38_chr == chrom)]
            if len(sub):
                bys = sub.groupby("strand").aligned_bp.sum()
                strand = bys.idxmax()
                abp = int(bys.max())
                qmin, qmax = int(sub.grch38_min.min()), int(sub.grch38_max.max())
        if name.startswith("GCA_000001405"):
            strand, qmin, abp = "+", 0, n  # GRCh38 itself
        riv = rd.get(c, [])
        rmid = None
        if riv:
            w = sum(e - s for s, e, _ in riv)
            rmid = sum((s + e) / 2 * (e - s) for s, e, _ in riv) / w
        rdna_end = None if rmid is None else ("start" if rmid < n / 2 else "end")
        if strand is not None:
            pend, how = ("start" if strand == "+" else "end"), "chain"
        elif rdna_end is not None:
            pend, how = rdna_end, "censat_rdna"
        else:
            pend, how = "start", "default"
        rec.update(chain_strand=strand, chain_aligned_bp=abp, grch38_min=qmin, grch38_max=qmax,
                   rdna_bp=sum(e - s for s, e, _ in riv), rdna_side=rdna_end,
                   parm_end=pend, orientation_source=how,
                   orientation_conflict=bool(rdna_end and strand and rdna_end != pend))
        take = level == "chromosome" or bool(riv) or strand is None or \
            (qmin is not None and qmin < GRCH38_CEN_END[chrom] + 2 * MB)
        rec["taken"] = take
        if take:
            feats = [(s, e) for s, e, _ in riv] + gaps.get(c, [])
            if pend == "start":
                far = [e for s, e in feats if s < ACRO_WINDOW]
                D = max(far) if far else 0
                seg = (0, min(n, max(ACRO_MIN, D + ACRO_EXTRA)))
            else:
                far = [n - s for s, e in feats if e > n - ACRO_WINDOW]
                D = max(far) if far else 0
                seg = (max(0, n - max(ACRO_MIN, D + ACRO_EXTRA)), n)
            regs.append((c, seg[0], seg[1], "acro_parm"))
            rec.update(seg_start=seg[0], seg_end=seg[1])
        acro_rows.append(rec)
    # (c) unplaced contigs with rDNA
    unassigned = [c for c in L if c not in alias or alias[c].startswith("chrUn")]
    info["unplaced_tilescan_added"] = ""
    if censat is not None:
        for c in unassigned:
            if c in rd:
                regs.append((c, 0, L[c], "unplaced_rdna"))
        # Validation (8 whole assemblies) found whole-rDNA chrUn contigs with NO CenSat label.
        # Screen every remaining unassigned contig with the shared tile measure.
        rest = [c for c in unassigned if c not in rd]
        if rest:
            fa = os.path.join(tmpdir, "unplaced_rest.fa")
            faidx_to_file(s3_url(m["assembly"]), m["fai"], m["gzi"], [(c, 0, L[c]) for c in rest], fa)
            info["fallback_bytes"] += est_compressed_bytes([(c, 0, L[c], "") for c in rest], fai,
                                                           read_gzi(m["gzi"]))
            scan = os.path.join(REG_DIR, f"{name}.unplaced_scan")
            h, _ = tile_map.run([fa], scan, threads=2)
            hit = set(h[h["pass"] & h.cls.isin(["45S", "5S"])].contig)
            hit |= {r[0] for r in tile_map.class_runs(h, "DJ")}
            for c in sorted(hit, key=lambda x: list(L).index(x)):
                regs.append((c, 0, L[c], "unplaced_rdna"))
            info["unplaced_tilescan_added"] = ",".join(sorted(hit))
            if hit:
                info["notes"].append(f"{len(hit)} unplaced contig(s) without CenSat rDNA label added by "
                                     f"tile screen ({len(rest)} screened, {sum(L[c] for c in rest)} bp)")
    else:
        tot = sum(L[c] for c in unassigned)
        if unassigned and tot <= 100 * MB:
            fa = os.path.join(tmpdir, "unassigned.fa")
            faidx_to_file(s3_url(m["assembly"]), m["fai"], m["gzi"],
                          [(c, 0, L[c]) for c in unassigned], fa)
            info["fallback_bytes"] += est_compressed_bytes(
                [(c, 0, L[c], "") for c in unassigned], fai, read_gzi(m["gzi"]))
            h, _ = tile_map.run([fa], os.path.join(tmpdir, "unassigned_scan"), threads=2)
            hit = sorted(set(h[(h["pass"]) & (h.cls == "45S")].contig))
            for c in hit:
                regs.append((c, 0, L[c], "unplaced_rdna"))
            info["notes"].append(f"no CenSat: {len(unassigned)} unassigned contigs "
                                 f"({tot} bp) tile-scanned, {len(hit)} with 45S tiles")
        elif unassigned:
            info["notes"].append(f"no CenSat and {tot} bp unassigned: unplaced rDNA not screened")
    # (d) 5S locus on chr1q42
    info["5S_method"] = None
    if name.startswith("GCA_000001405"):
        c = prefix + "chr1"
        core = (FLANK_L[1], FLANK_R[0])
        regs.append((c, core[0] - PAD_5S, core[1] + PAD_5S, "5S_locus"))
        info.update({"5S_method": "identity (GRCh38)", "5S_contig": c, "5S_core": f"{core[0]}-{core[1]}"})
    elif isinstance(m.get("chain"), str) and os.path.exists(m["chain"]):
        chains = list(parse_chain(m["chain"], keep_q={"chr1"}))
        lf, rf = lift_flank(chains, FLANK_L, prefix), lift_flank(chains, FLANK_R, prefix)
        info["5S_left"] = None if lf is None else f"{lf['contig']}:{lf['tmin']}-{lf['tmax']}({lf['strand']},{lf['bp']}bp)"
        info["5S_right"] = None if rf is None else f"{rf['contig']}:{rf['tmin']}-{rf['tmax']}({rf['strand']},{rf['bp']}bp)"
        if lf and rf and lf["contig"] == rf["contig"]:
            a, b = sorted([(lf["tmin"], lf["tmax"]), (rf["tmin"], rf["tmax"])])
            core = (a[1], b[0]) if a[1] <= b[0] else (b[0], a[1])
            if core[1] - core[0] > 5 * MB:
                info["notes"].append(f"5S lift span implausible ({core[1]-core[0]} bp)")
            c = lf["contig"]
            regs.append((c, max(0, core[0] - PAD_5S), min(L[c], core[1] + PAD_5S), "5S_locus"))
            info.update({"5S_method": "chain_both_flanks", "5S_contig": c, "5S_core": f"{core[0]}-{core[1]}"})
        elif lf or rf:
            for side, x in (("left", lf), ("right", rf)):
                if x is None:
                    continue
                c = x["contig"]
                # inner edge of a single lifted flank: the end facing the array
                inner = x["tmax"] if (side == "left") == (x["strand"] == "+") else x["tmin"]
                regs.append((c, max(0, inner - PAD_5S), min(L[c], inner + PAD_5S), "5S_locus"))
            info["5S_method"] = "chain_one_flank" if not (lf and rf) else "chain_flanks_split_contigs"
            info["notes"].append(f"5S lift partial: {info['5S_method']}")
        else:
            info["5S_method"] = "lift_failed"
            info["notes"].append("5S lift failed: no chain block in either flank")
    if info["5S_method"] in (None, "lift_failed"):
        chr1 = [c for c, u in alias.items() if u == "chr1" or (u.startswith("chr1_") and u.endswith("_random"))]
        if chr1:
            fa = os.path.join(tmpdir, "chr1.fa")
            faidx_to_file(s3_url(m["assembly"]), m["fai"], m["gzi"], [(c, 0, L[c]) for c in chr1], fa)
            info["fallback_bytes"] += est_compressed_bytes([(c, 0, L[c], "") for c in chr1], fai,
                                                           read_gzi(m["gzi"]))
            h, _ = tile_map.run([fa], os.path.join(tmpdir, "chr1_scan"), threads=2)
            h5 = h[(h["pass"]) & (h.cls == "5S")]
            if len(h5):
                cl = []
                for c, g in h5.groupby("contig"):
                    cl += [(c, s, e) for s, e, _ in merge_intervals(
                        [(r.tile_start, r.tile_end, "") for r in g.itertuples()], 50_000)]
                c, s, e = max(cl, key=lambda x: x[2] - x[1])
                regs.append((c, max(0, s - PAD_5S), min(L[c], e + PAD_5S), "5S_locus"))
                was = info["5S_method"]
                info.update({"5S_method": "tile_scan_chr1" + ("" if was is None else "_after_lift_failure"),
                             "5S_contig": c, "5S_core": f"{s}-{e}"})
                info["notes"].append(f"5S located by tile scan of {len(chr1)} chr1 contigs (no usable chain)")
            else:
                info["notes"].append("5S not found by chr1 tile scan")
        else:
            info["notes"].append("5S: no chain and no chr1-assigned contig")
    # merge per contig
    merged = []
    byc = {}
    for c, s, e, r in regs:
        byc.setdefault(c, []).append((s, e, r))
    for c in byc:
        for s, e, r in merge_intervals(byc[c], 0):
            merged.append((c, s, e, r))
    order = {c: i for i, c in enumerate(fai)}
    merged.sort(key=lambda x: (order[x[0]], x[1]))
    return merged, acro_rows, info


def verify(out_gz, regions):
    fai = out_gz + ".fai"
    got = {}
    for line in open(fai):
        f = line.split("\t")
        got[f[0]] = int(f[1])
    want = {region_name(c, s, e): e - s for c, s, e, _ in regions}
    return got == want


def process(m, force=False):
    name = m["assembly_name"]
    out_gz = os.path.join(SEQ_DIR, f"{name}.fa.gz")
    done = os.path.join(SEQ_DIR, f"{name}.done.json")
    if not force and os.path.exists(done):
        return json.load(open(done))
    t0 = time.time()
    with tempfile.TemporaryDirectory(dir=os.path.join(DATA, "tmp")) as td:
        regions, acro_rows, info = plan(m, td)
        bed = os.path.join(REG_DIR, f"{name}.bed")
        with open(bed, "w") as fh:
            for c, s, e, r in regions:
                fh.write(f"{c}\t{s}\t{e}\t{r}\n")
        pd.DataFrame(acro_rows).to_csv(os.path.join(REG_DIR, f"{name}.acro.tsv"), sep="\t", index=False)
        ok = False
        info["sequence_refetched"] = True
        if os.path.exists(out_gz + ".fai"):
            try:
                ok = verify(out_gz, regions)   # unchanged region set: keep existing sequence
                info["sequence_refetched"] = not ok
            except Exception:  # noqa
                ok = False
        for attempt in range(0 if ok else 4):
            try:
                raw = os.path.join(td, "x.fa")
                faidx_to_file(s3_url(m["assembly"]), m["fai"], m["gzi"],
                              [(c, s, e) for c, s, e, _ in regions], raw)
                subprocess.run(["bgzip", "-f", "-@", "2", raw], check=True)
                os.replace(raw + ".gz", out_gz)
                subprocess.run(["samtools", "faidx", out_gz], check=True)
                if verify(out_gz, regions):
                    ok = True
                    break
                info["notes"].append(f"attempt {attempt + 1}: length mismatch")
            except Exception as e:  # noqa
                info["notes"].append(f"attempt {attempt + 1}: {str(e)[:200]}")
                time.sleep(5 * (attempt + 1))
    fai = read_fai_full(m["fai"])
    est = est_compressed_bytes(regions, fai, read_gzi(m["gzi"]))
    res = {"assembly_name": name, "status": "ok" if ok else "FAILED", "n_regions": len(regions),
           "bp_extracted": sum(e - s for _, s, e, _ in regions),
           "est_bytes_fetched": est + info.pop("fallback_bytes"),
           "seconds": round(time.time() - t0, 1)}
    for r in ("rdna_censat", "acro_parm", "unplaced_rdna", "5S_locus"):
        res[f"bp_{r}"] = sum(e - s for _, s, e, t in regions if r in t.split(";"))
    res["n_acro_contigs_taken"] = sum(1 for a in acro_rows if a.get("taken"))
    res["n_acro_contigs_skipped"] = sum(1 for a in acro_rows if not a.get("taken"))
    res["n_orientation_reversed"] = sum(1 for a in acro_rows if a.get("taken") and a["parm_end"] == "end")
    res["n_orientation_conflict"] = sum(1 for a in acro_rows if a.get("orientation_conflict"))
    res["sequence_refetched"] = info.get("sequence_refetched", True)
    res["unplaced_tilescan_added"] = info.get("unplaced_tilescan_added", "")
    res["n_unplaced_tilescan_added"] = len([x for x in res["unplaced_tilescan_added"].split(",") if x])
    for k in ("5S_method", "5S_contig", "5S_core", "5S_left", "5S_right"):
        res[k] = info.get(k)
    res["notes"] = "; ".join(info["notes"])
    if ok:
        json.dump(res, open(done, "w"))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    for d in (SEQ_DIR, REG_DIR, os.path.join(DATA, "tmp"), TABLES):
        os.makedirs(d, exist_ok=True)
    man = pd.read_csv(os.path.join(META, "manifest.tsv"), sep="\t")
    idx = pd.read_csv(os.path.join(DATA, "index.csv"))[["assembly_name", "assembly"]]
    man = man.merge(idx, on="assembly_name")
    if a.only:
        man = man[man.assembly_name.isin(a.only)]
    tile_map.build_targets()
    t0 = time.time()
    results = []
    with ThreadPoolExecutor(a.threads) as ex:
        futs = {ex.submit(process, r.to_dict(), a.force): r["assembly_name"] for _, r in man.iterrows()}
        for i, f in enumerate(as_completed(futs)):
            n = futs[f]
            try:
                results.append(f.result())
            except Exception as e:  # noqa
                results.append({"assembly_name": n, "status": "FAILED", "notes": f"exception: {e}"})
            if (i + 1) % 25 == 0:
                print(f"{i + 1}/{len(futs)} done, {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(results)
    out = os.path.join(TABLES, "extraction.tsv")
    if a.only and os.path.exists(out):
        old = pd.read_csv(out, sep="\t")
        res = pd.concat([old[~old.assembly_name.isin(res.assembly_name)], res])
    order = {n: i for i, n in enumerate(pd.read_csv(os.path.join(DATA, "index.csv")).assembly_name)}
    res = res.sort_values("assembly_name", key=lambda s: s.map(order))
    res.to_csv(out, sep="\t", index=False)
    acro = [pd.read_csv(os.path.join(REG_DIR, f"{n}.acro.tsv"), sep="\t")
            for n in res.assembly_name if os.path.exists(os.path.join(REG_DIR, f"{n}.acro.tsv"))
            and os.path.getsize(os.path.join(REG_DIR, f"{n}.acro.tsv")) > 1]
    if acro:
        pd.concat(acro).to_csv(os.path.join(TABLES, "acro_contigs.tsv"), sep="\t", index=False)
    print(json.dumps({"n": len(res), "ok": int((res.status == "ok").sum()),
                      "failed": res.loc[res.status != "ok", "assembly_name"].tolist(),
                      "est_bytes_fetched": float(res.est_bytes_fetched.sum()),
                      "wall_seconds": round(time.time() - t0, 1)}))


if __name__ == "__main__":
    main()
