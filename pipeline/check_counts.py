#!/usr/bin/env python3
"""Checks the cohort scripts make before they count, and before a CRAM is removed. Standard library only:
they run inside whichever image the cohort uses, which may predate anything in the ngsdose package.

    check_counts.py fetch-panels SINKS_BED PANEL...   every class of the panels has an interval in the sinks BED
    check_counts.py engine-build COUNTS               print the engine commit the counts file was written by
    check_counts.py same-reads SCAN FETCH [SINKS [CONTROLS]]
                                                      the fetch saw the scan's control and region reads, read for read;
                                                      with the sinks BED the fetch used, also every class read the scan
                                                      placed inside those sinks (and no more reads than the scan has,
                                                      unless the scan loaded panels the fetch did not: then a note).
                                                      With the controls FASTA the fetch was to use (a fetch plan's
                                                      controls.txt), the fetch must have used it; when it is not the
                                                      scan's (controls.lite200), the fetch's regions must all be regions
                                                      of the scan, and are compared one by one, with ctrl_reads against
                                                      the scan's reads in the fetch's control regions

Exit status 0: passed; 1: failed, with the reason on stderr; 2: a file could not be read, or the two were made with
region sets that cannot be compared (a fetch region the scan lacks).
"""
import bisect
import gzip
import hashlib
import json
import sys
import zlib


def _open(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def panel_classes(path):
    """The class names of a panel, from its ##class header lines."""
    names = []
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##class"):
                fields = dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:] if "=" in f)
                names.append(fields.get("name", ""))
    if not names:
        raise ValueError(f"{path}: no ##class header lines")
    return names


def sink_classes(path):
    with _open(path) as fh:
        return {p[3] for p in (line.split() for line in fh if line.strip() and not line.startswith(("#", "track", "browser"))) if len(p) > 3}


def counts(path):
    with _open(path) as fh:
        return json.load(fh)


def reads(c):
    return c["ctrl_reads"], sorted((r["name"], r["obs"]) for r in c["regions"])


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Uncomparable(ValueError):
    """The two files cannot be compared (exit status 2), as opposed to disagreeing (1)."""


def subset_reads(scan, fetch):
    """The scan's control and region reads restricted to the regions the fetch holds, for a fetch made with a subset
    of the scan's control regions (a fetch plan's controls.lite200): (ctrl_reads, sorted (name, obs)). ctrl_reads is
    the sum of obs over the regions of role control, as the engine counts it. A fetch region that the scan lacks, or
    holds with another role, means the two were made with different region sets, not a subset: Uncomparable."""
    by = {r["name"]: r for r in scan["regions"]}
    missing = [r["name"] for r in fetch["regions"] if r["name"] not in by]
    if missing:
        raise Uncomparable(f"the fetch holds {len(missing)} region(s) the scan does not (e.g. {missing[0]}): its controls "
                           f"({fetch.get('controls')}) are not a subset of the scan's ({scan.get('controls')})")
    role = [r["name"] for r in fetch["regions"] if r.get("role", "control") != by[r["name"]].get("role", "control")]
    if role:
        raise Uncomparable(f"{len(role)} region(s) have another role in the fetch than in the scan (e.g. {role[0]})")
    if not any(r.get("role", "control") == "control" for r in fetch["regions"]):
        raise Uncomparable("the fetch holds no control region")
    names = {r["name"] for r in fetch["regions"]}
    ctrl = sum(r["obs"] for n, r in by.items() if n in names and r.get("role", "control") == "control")
    return ctrl, sorted((n, r["obs"]) for n, r in by.items() if n in names)


def read_sinks(path):
    """{(class, contig): merged [start, end] intervals} of a sinks BED."""
    iv = {}
    with _open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 4 or line.startswith(("#", "track", "browser")):
                continue
            iv.setdefault((p[3], p[0]), []).append((int(p[1]), int(p[2])))
    out = {}
    for key, rows in iv.items():
        merged = []
        for s, e in sorted(rows):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[key] = ([s for s, _ in merged], merged)
    return out


def in_sinks(scan, sinks):
    """Per class, the scan's reads placed in a bin that lies wholly inside the class's sinks (the bin clipped at
    its contig's end): reads a fetch through those sinks reads without fail."""
    length = {c["name"]: c["len"] for c in scan["contigs"]}
    kind = {c["name"]: c["kind"] for c in scan["classes"]}
    n = {}
    for p in scan["placements"]:
        cls, contig = p["class"], p["contig"]
        if (cls, contig) not in sinks or contig not in length:
            continue
        width = scan["placement_bin"] if kind.get(cls) == "positional" else scan["placement_bin_compositional"]
        start, end = p["start"], min(p["start"] + width, length[contig])
        starts, merged = sinks[(cls, contig)]
        i = bisect.bisect_right(starts, start) - 1
        if i >= 0 and merged[i][0] <= start and end <= merged[i][1]:
            n[cls] = n.get(cls, 0) + p["reads"]
    return n


def class_reads(scan, fetch, sinks_path):
    """Problems with the fetch's class reads: fewer than the scan placed inside the sinks (an index that lacks a
    contig), or more than the scan has. [] when the sinks file is not the one the fetch was made with."""
    with open(sinks_path, "rb") as fh:
        if hashlib.sha256(fh.read()).hexdigest() != fetch.get("sinks_sha256"):
            print(f"note: {sinks_path} is not the sinks file the fetch was made with; its class reads are not checked",
                  file=sys.stderr)
            return []
    inside = in_sinks(scan, read_sinks(sinks_path))
    total = {c["name"]: c["reads"] for c in scan["classes"]}
    # A scan that loaded panels the fetch did not (the satellites) gave some reads to those classes, or left them
    # ambiguous, that the fetch can give to one of its own: there, more reads than the scan is noted, not failed.
    co_loaded = sorted(set(total) - {c["name"] for c in fetch["classes"]})
    bad = []
    for c in fetch["classes"]:
        name, got = c["name"], c["reads"]
        if name not in total:
            continue
        if got < inside.get(name, 0):
            bad.append(f"{name}: fetch {got:,} reads, fewer than the {inside[name]:,} the scan placed inside its sinks")
        elif got > total[name]:
            msg = f"{name}: fetch {got:,} reads, more than the scan's {total[name]:,}"
            if co_loaded:
                print(f"note: {msg} (the scan also loaded {', '.join(co_loaded[:3])}"
                      f"{' ...' if len(co_loaded) > 3 else ''}, which the fetch did not)", file=sys.stderr)
            else:
                bad.append(msg)
    return bad


def main(argv):
    cmd, args = (argv[0], argv[1:]) if argv else ("", [])
    try:
        if cmd == "fetch-panels" and len(args) >= 2:
            have = sink_classes(args[0])
            missing = [(n, p) for p in args[1:] for n in panel_classes(p) if n not in have]
            for n, p in missing:
                print(f"class {n} (panel {p}) has no interval in {args[0]}", file=sys.stderr)
            return 1 if missing else 0
        if cmd == "engine-build" and len(args) == 1:
            print(counts(args[0]).get("engine_build", ""))
            return 0
        if cmd == "same-reads" and len(args) in (2, 3, 4):
            scan, fetch = counts(args[0]), counts(args[1])
            (sc, sr), (fc, fr) = reads(scan), reads(fetch)
            if len(args) == 4:
                want = sha256_of(args[3])
                if fetch.get("controls_sha256") != want:
                    print(f"fetch and scan disagree:\n  the fetch was made with the controls {fetch.get('controls')} "
                          f"(sha256 {str(fetch.get('controls_sha256'))[:12]}), not with {args[3]} (sha256 {want[:12]})", file=sys.stderr)
                    return 1
                if fetch.get("controls_sha256") != scan.get("controls_sha256"):
                    try:
                        sc, sr = subset_reads(scan, fetch)
                    except Uncomparable as e:
                        print(f"check_counts.py same-reads: {e}", file=sys.stderr)
                        return 2
                    print(f"note: the fetch counted {len(fr)} of the scan's {len(scan['regions'])} regions (controls "
                          f"{fetch.get('controls')}, the scan's {scan.get('controls')}); compared on those", file=sys.stderr)
            bad = class_reads(scan, fetch, args[2]) if len(args) >= 3 else []
            if (sc, sr) == (fc, fr) and not bad:
                return 0
            s, f = dict(sr), dict(fr)
            diff = sorted(n for n in set(s) | set(f) if s.get(n) != f.get(n))
            print("fetch and scan disagree:", file=sys.stderr)
            if (sc, sr) != (fc, fr):
                print(f"  ctrl_reads {fc} against {sc}; {len(diff)} of {len(s)} regions differ"
                      + "".join(f"\n  {n}: fetch {f.get(n)}, scan {s.get(n)}" for n in diff[:5]), file=sys.stderr)
            for b in bad:
                print(f"  {b}", file=sys.stderr)
            return 1
    except (OSError, EOFError, zlib.error, ValueError, KeyError, TypeError) as e:
        print(f"check_counts.py {cmd}: {e}", file=sys.stderr)
        return 2
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
