#!/usr/bin/env python3
"""Fetch exactly the CRAM containers that hold a set of regions, by byte range, and write them as a small CRAM.

htslib reads a remote file with range requests that have no end: every query starts a stream that is thrown away at the
next seek, and what was in flight is lost (measured on the karyotype windows' candidates: 1,517 MB moved for 410 MB of
slices). The CRAM index says where each container lies, so the containers can be asked for exactly: the file definition
and header, the containers that hold the regions, the end-of-file container. The result is a valid CRAM of those
containers (index it with `samtools index`), which `ngs-dose count -m fetch` counts as it would the remote file: the
same reads in the same regions.

usage: slice_fetch.py URL CRAI BED OUT.cram --fai REFERENCE.fa.fai [--pad 600] [--threads 8]
  URL   the CRAM over http(s); a signed URL keeps its query
  CRAI  its index, local (the CRAM's own: an index of another file of the same sample gives other containers)
  BED   the regions (the bundle's controls.bed, or a control set's BED)
  --fai the reference's .fai, or any file whose lines start with the contig names in the CRAM header's order: the index
        numbers contigs by that order
  --reuse CRAM BED: a small CRAM made earlier by this script from the same URL and index for the regions of BED: the byte
        ranges it holds are copied from it and only the others fetched (to add regions to a fetch without paying twice)
"""
import argparse
import gzip
import http.client
import os
import sys
import threading
import time
import urllib.parse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

EOF_LEN = 38                      # the CRAM 3 end-of-file container
local = threading.local()


def connection(u):
    return (http.client.HTTPSConnection if u.scheme == "https" else http.client.HTTPConnection)(u.netloc, timeout=60)


def target(u):
    return u.path + (f"?{u.query}" if u.query else "")


def get(url, a, b, tries=6):
    """Bytes a..b inclusive, on a connection kept per thread."""
    u = urllib.parse.urlsplit(url)
    for t in range(tries):
        try:
            c = getattr(local, "conn", None)
            if c is None:
                c = local.conn = connection(u)
            c.request("GET", target(u), headers={"Range": f"bytes={a}-{b}", "Connection": "keep-alive"})
            r = c.getresponse()
            data = r.read()
            if r.status == 206 and len(data) == b - a + 1:
                return data
            raise OSError(f"status {r.status}, {len(data)} bytes of {b - a + 1}")
        except (OSError, http.client.HTTPException) as e:
            if getattr(local, "conn", None) is not None:
                local.conn.close()
            local.conn = None
            if t == tries - 1:
                raise OSError(f"bytes {a}-{b}: {e}") from None
            time.sleep(1.5 * (t + 1))


def size_of(url):
    """The file's length, from the answer to a request for its first byte (a URL signed for GET is not signed for HEAD)."""
    u = urllib.parse.urlsplit(url)
    c = connection(u)
    c.request("GET", target(u), headers={"Range": "bytes=0-0"})
    r = c.getresponse()
    r.read()
    c.close()
    total = (r.getheader("Content-Range") or "").rpartition("/")[2]
    if r.status != 206 or not total.isdigit():
        raise OSError(f"the server did not answer a range request (status {r.status}): a fetch needs range requests")
    return int(total)


def plan(crai, bed, pad, ids):
    """The regions of `bed` by contig number, each container's slices by contig, and every container's offset."""
    by = defaultdict(list)
    for line in open(bed):
        p = line.split("\t")
        if not line.startswith("#") and len(p) >= 3 and p[0] in ids:
            by[ids[p[0]]].append((int(p[1]) - pad, int(p[2]) + pad))
    slices = defaultdict(list)
    offs = set()
    with gzip.open(crai, "rt") as fh:
        for line in fh:
            s, st, span, coff, _, _ = map(int, line.split("\t"))
            offs.add(coff)
            if s in by:
                slices[s].append((st, st + span, coff))
    return by, slices, sorted(offs)


def ranges_of(by, slices, offs, total, join):
    """The containers that hold the regions, and the byte ranges that hold them (each container whole; ranges less than
    `join` bytes apart asked for as one, the containers between them included)."""
    nxt = {o: n for o, n in zip(offs, offs[1:] + [total - EOF_LEN])}
    want = set()
    for s, regs in by.items():
        sl = sorted(slices[s])
        for rs, re_ in regs:
            want.update(coff for st, en, coff in sl if st <= re_ and en >= rs)
    ranges = []
    for o in sorted(want):
        if ranges and o - ranges[-1][1] <= join:
            ranges[-1][1] = nxt[o]
        else:
            ranges.append([o, nxt[o]])
    return sorted(want), ranges


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("url")
    ap.add_argument("crai")
    ap.add_argument("bed")
    ap.add_argument("out")
    ap.add_argument("--fai", required=True, help="contig names in the CRAM header's order (the reference's .fai)")
    ap.add_argument("--pad", type=int, default=600, help="bp added either side of a region, as the engine reads them")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--join", type=int, default=65536, help="ranges closer than this many bytes are asked for as one")
    ap.add_argument("--reuse", nargs=2, metavar=("CRAM", "BED"))
    a = ap.parse_args()
    with open(a.fai) as fh:
        ids = {line.split("\t")[0].strip(): i for i, line in enumerate(fh) if line.strip()}
    by, slices, offs = plan(a.crai, a.bed, a.pad, ids)
    total = size_of(a.url)
    want, ranges = ranges_of(by, slices, offs, total, a.join)
    have = {}
    if a.reuse:
        by0, slices0, _ = plan(a.crai, a.reuse[1], a.pad, ids)
        _, old = ranges_of(by0, slices0, offs, total, a.join)
        at = offs[0]                      # the earlier file: the header, its ranges in order, the end-of-file container
        for x, y in old:
            have[(x, y)] = at
            at += y - x
        if os.path.getsize(a.reuse[0]) != at + EOF_LEN:
            raise SystemExit(f"{a.reuse[0]} is not the CRAM of {a.reuse[1]} with this index: {os.path.getsize(a.reuse[0])} bytes, {at + EOF_LEN} expected")

    def piece(r):
        x, y = r
        for (ox, oy), pos in have.items():
            if ox <= x and y <= oy:
                with open(a.reuse[0], "rb") as fh:
                    fh.seek(pos + x - ox)
                    return fh.read(y - x), 0
        return get(a.url, x, y - 1), y - x
    parts = [(0, offs[0])] + [(x, y) for x, y in ranges] + [(total - EOF_LEN, total)]
    t0 = time.time()
    with ThreadPoolExecutor(a.threads) as ex:
        got = list(ex.map(piece, parts))
    with open(a.out + ".part", "wb") as fh:
        for b, _ in got:
            fh.write(b)
    os.replace(a.out + ".part", a.out)
    size, moved = sum(len(b) for b, _ in got), sum(n for _, n in got)
    print(f"{a.out}: {len(want)} containers in {len(ranges)} ranges, {size / 1e6:.1f} MB of {total / 1e9:.2f} GB ({100 * size / total:.2f}%), "
          f"{moved / 1e6:.1f} MB fetched, {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
