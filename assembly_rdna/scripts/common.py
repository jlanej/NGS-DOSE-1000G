"""Shared helpers for the HPRC r2 assembly rDNA extraction (assembly_rdna/).

Coordinate conventions used everywhere in this directory
  * BED files: 0-based half-open [start, end).
  * Sequence names / samtools faidx region strings: 1-based inclusive,
    i.e. BED (start, end) <-> "contig:{start+1}-{end}".
"""
import gzip
import io
import os
import re
import time
import xml.etree.ElementTree as ET

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
META = os.path.join(DATA, "meta")
TABLES = os.path.join(ROOT, "tables")
S3_HTTP = "https://s3-us-west-2.amazonaws.com/human-pangenomics"
LOCAL_CENSAT_DIR = "/Users/Kitty/git/NGS-DOSE-1000G/hprc_censat"   # read-only reuse
CHM13_CENSAT = "/Users/Kitty/git/NGS-DOSE/work/ref/chm13v2.0_censat_v2.1.bed"
HG002_CENSAT = {  # T2T consortium CenSat v2.0 for the Q100 hg002v1.1 assembly
    "hg002v1.1.pat": "T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.PAT.bed",
    "hg002v1.1.mat_MT": "T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.MAT.bed",
}
ACRO = ["chr13", "chr14", "chr15", "chr21", "chr22"]
# GRCh38 centromere ends (approx., UCSC centromeres track) - used to decide whether an
# unlocalised (_random) acrocentric contig reaches the short-arm/centromere side.
GRCH38_CEN_END = {"chr13": 18_051_248, "chr14": 18_173_523, "chr15": 19_725_254,
                  "chr21": 12_915_808, "chr22": 15_054_318}

_session = None


def session():
    global _session
    if _session is None:
        _session = requests.Session()
        a = requests.adapters.HTTPAdapter(pool_connections=16, pool_maxsize=16)
        _session.mount("https://", a)
    return _session


def s3_url(key_or_uri):
    k = key_or_uri
    if k.startswith("s3://human-pangenomics/"):
        k = k[len("s3://human-pangenomics/"):]
    return f"{S3_HTTP}/{k}"


def http_get(url, retries=5, stream_to=None, timeout=120):
    """GET with retries. Returns bytes (or writes to stream_to path and returns size)."""
    last = None
    for i in range(retries):
        try:
            r = session().get(url, timeout=timeout, stream=stream_to is not None)
            if r.status_code == 404:
                raise FileNotFoundError(url)
            r.raise_for_status()
            if stream_to is None:
                return r.content
            tmp = stream_to + ".part"
            n = 0
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
                    n += len(chunk)
            os.replace(tmp, stream_to)
            return n
        except FileNotFoundError:
            raise
        except Exception as e:  # noqa
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed after {retries} tries: {url}: {last}")


def s3_list(prefix):
    """List all keys (with sizes) under an S3 prefix (paginated ListObjectsV2)."""
    out, token = [], None
    ns = "{http://s3.amazonaws.com/doc/2006-03-01/}"
    while True:
        params = {"list-type": "2", "prefix": prefix}
        if token:
            params["continuation-token"] = token
        for i in range(5):
            try:
                r = session().get(S3_HTTP, params=params, timeout=60)
                r.raise_for_status()
                break
            except Exception:
                time.sleep(2 * (i + 1))
        else:
            raise RuntimeError(f"listing failed: {prefix}")
        root = ET.fromstring(r.content)
        for c in root.findall(f"{ns}Contents"):
            out.append((c.find(f"{ns}Key").text, int(c.find(f"{ns}Size").text)))
        if root.findtext(f"{ns}IsTruncated") == "true":
            token = root.findtext(f"{ns}NextContinuationToken")
        else:
            return out


def hap_label(name):
    m = re.search(r"(?:^|[_.])(mat|pat|hap1|hap2)(?:[_.]|$)", name)
    return m.group(1) if m else "ref"


def open_text(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path), encoding="utf-8")
    return open(path)


def read_fai(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[0]] = int(f[1])
    return d


def read_censat(path):
    """Return list of (contig, start, end, label) from a CenSat BED (track lines skipped)."""
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(("track", "#", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), f[3]))
    return rows


def is_rdna_label(label):
    return "rdna" in label.lower()


def merge_intervals(iv, gap=0):
    """Merge [(s,e,tag)] (same contig) when separated by <= gap. Tags are ';'-joined sets."""
    iv = sorted(iv)
    out = []
    for s, e, t in iv:
        if out and s <= out[-1][1] + gap:
            ps, pe, pt = out[-1]
            tags = set(pt.split(";")) | set(t.split(";"))
            out[-1] = (ps, max(pe, e), ";".join(sorted(tags)))
        else:
            out.append((s, e, t))
    return out


def region_name(contig, s, e):
    """0-based half-open (s,e) -> samtools 1-based inclusive region string."""
    return f"{contig}:{s + 1}-{e}"


def parse_chain(path, keep_q=None):
    """Parse a UCSC chain file (gz). For these HPRC files t = assembly contig (no PanSN
    prefix), q = GRCh38 chromosome. Yields dicts with header fields and aligned blocks
    as (t_start, q_start, size) in + strand coordinates of t and strand-specific q
    coordinates (q coords on '-' strand are relative to the reverse complement)."""
    with gzip.open(path, "rt") as fh:
        cur = None
        for line in fh:
            if line.startswith("chain"):
                if cur is not None:
                    yield cur
                f = line.split()
                cur = dict(score=int(f[1]), tName=f[2], tSize=int(f[3]), tStrand=f[4],
                           tStart=int(f[5]), tEnd=int(f[6]), qName=f[7], qSize=int(f[8]),
                           qStrand=f[9], qStart=int(f[10]), qEnd=int(f[11]), id=f[12],
                           blocks=[] if (keep_q is None or f[7] in keep_q) else None)
                tp, qp = cur["tStart"], cur["qStart"]
                continue
            f = line.split()
            if not f or cur is None:
                continue
            size = int(f[0])
            if cur["blocks"] is not None:
                cur["blocks"].append((tp, qp, size))
            cur.setdefault("aligned", 0)
            cur["aligned"] += size
            if len(f) == 3:
                tp += size + int(f[1])
                qp += size + int(f[2])
        if cur is not None:
            yield cur
