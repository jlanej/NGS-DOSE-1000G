#!/usr/bin/env python3
"""Is the fathers' deficit in passing on a lost junction a property of reading a male genome?

A father's loss of a whole junction passes to few of his children, a mother's to half (docs/DJ.md, section 4.3). If the
calls read a man's genome differently from a woman's, nine-copy calls would differ in frequency between the sexes, or
the level of the genomes at ten would. This counts both, the carriers by their place in the pedigree, and the
transmissions by the child's sex. Reads docs/data/cohort.tsv and the pedigree; run from the repository root.
"""
import csv
import re
import statistics as st
from collections import Counter

PED = "meta/20130606_g1k_3202_samples_ped_population.txt"
rows = {r["sample"]: r for r in csv.DictReader(open("docs/data/cohort.tsv"), delimiter="\t")}
ped = {p["SampleID"]: p for p in csv.DictReader(open(PED), delimiter=" ")}
sex = {s: ("M" if p["Sex"] == "1" else "F") for s, p in ped.items()}


def large(r):
    """An event of 40 kb or more, as the sample table prints it."""
    for item in r["DJ.variants"].split(";"):
        m = re.fullmatch(r"([+-]\d+):(\d+)-(\d+)kb", item.strip())
        if m and int(m.group(3)) - int(m.group(2)) >= 40:
            return True
    return False


settled = {s: r for s, r in rows.items() if r["DJ.call"] == "settled"}
whole = {s: int(r["DJ.copies"]) for s, r in settled.items() if not large(r)}           # one state throughout

print("by sex (settled calls; one state throughout):")
for sx in "MF":
    S = [s for s in settled if sex[s] == sx]
    c = Counter(whole[s] for s in S if s in whole)
    ten = [float(rows[s]["DJ.cn"]) for s in S if whole.get(s) == 10]
    print(f"  {sx}: {len(S)} genomes; " + ", ".join(f"{k}: {v}" for k, v in sorted(c.items()))
          + f"; level of the genomes at ten {st.mean(ten):.3f} (SD {st.stdev(ten):.3f}, n {len(ten)})")

trio = {s: p for s, p in ped.items() if p["FatherID"] in rows and p["MotherID"] in rows and s in rows}
fathers, mothers = {p["FatherID"] for p in trio.values()}, {p["MotherID"] for p in trio.values()}
others = set(rows) - set(trio) - fathers - mothers
print("nine copies throughout, by place in the pedigree:")
for name, S in (("fathers", fathers), ("mothers", mothers), ("children", set(trio)), ("others, men", {s for s in others if sex[s] == "M"}),
                ("others, women", {s for s in others if sex[s] == "F"})):
    S = [s for s in S if s in settled]
    n = sum(1 for s in S if whole.get(s) == 9)
    print(f"  {name:14s} {n:3d} of {len(S):4d} ({100 * n / len(S):.1f}%)")

tr = Counter()
for s, p in trio.items():
    f, m = p["FatherID"], p["MotherID"]
    if s not in settled or f not in settled or m not in settled:
        continue
    for who, par, other in (("father", f, m), ("mother", m, f)):
        if whole.get(par) == 9 and whole.get(other) == 10 and whole.get(s) in (9, 10):
            tr[(who, "son" if sex[s] == "M" else "daughter", whole[s] == 9)] += 1
print("a parent at nine throughout, the other at ten, the child at nine or ten: passed on, by the child's sex")
for who in ("father", "mother"):
    for kid in ("son", "daughter"):
        y, n = tr[(who, kid, True)], tr[(who, kid, False)]
        print(f"  {who}'s loss to a {kid}: {y} of {y + n}")
