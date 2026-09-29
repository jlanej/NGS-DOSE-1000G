#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_sinks
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/sinks_%j.out
# Learn where the reads of named classes land - new ones scanned through CANDIDATE_PANELS (config.sh), or any
# other class the scans carry - from the whole-file scans, and check it on scans the learning did not see:
#
#   bash 06_learn_sinks.sh NAME CLASS [CLASS...]          # or: sbatch 06_learn_sinks.sh NAME CLASS [CLASS...]
#
# The scans in $WORK_DIR/counts_scan that carry every CLASS are shuffled (seed SINKS_SEED, default 7). The first
# SINKS_TRAIN of them (default 30, and at most half) learn the sinks (`ngsdose sinks --classes ... --stats`): the
# 10-kb bins that hold >= 1e-5 of a class's reads and >= 25 reads in any training scan, merged and padded by 1 kb.
# All the others are held out and measure what those sinks capture (`ngsdose sinks --evaluate ... --stats`). With
# fewer than 30 scans carrying the classes it stops: a candidate is scanned until they are there. It writes, in
# $WORK_DIR/sinks:
#   NAME.bed              the sink intervals of the CLASSes (the positional classes the learner always adds are left out)
#   NAME.train.stats.tsv  per interval, its share of the class and the capture curve, in the training scans (in-sample)
#   NAME.stats.tsv        the same in the held-out scans: what `ngsdose fetchplan --capture` trims by
#   NAME.heldout.tsv      the capture of every held-out scan and class
#   NAME.scans.txt        which scans trained and which were held out
#   NAME.menu.tsv         the image's fetch menu, its paths made absolute, with a row per CLASS (status experimental,
#                         tier D, these sinks and statistics) and a preset NAME; fetch them with
#                         FETCH_MENU=$WORK_DIR/sinks/NAME.menu.tsv FETCH_PRESET="core NAME" (config.sh)
#   NAME.<panel>          a copy of each candidate panel defining a CLASS, which the menu names (a panel of the image
#                         is named where it is)
# and prints each class's held-out capture: median, 10th percentile and minimum over the held-out scans. Needs an image
# whose ngsdose has `sinks --stats` (and, for the menu, resources/fetch_menu.tsv); the fae1124 image has neither.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
NAME="${1:-}"; shift || true
[ -n "$NAME" ] && [ "$#" -gt 0 ] || { echo "usage: 06_learn_sinks.sh NAME CLASS [CLASS...]" >&2; exit 2; }
case "$NAME" in *[!A-Za-z0-9._-]*|.*) echo "ERROR: NAME $NAME: letters, digits, '.', '_' and '-' only" >&2; exit 2 ;; esac
CLASSES="$*"
SINKS_TRAIN="${SINKS_TRAIN:-30}"
SINKS_SEED="${SINKS_SEED:-7}"
MIN_SCANS=30                                           # scans that carry the classes, training and held-out together
OUT="$WORK_DIR/sinks"; mkdir -p "$OUT"
TMP="$OUT/.$NAME.tmp.$$"; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

SINKS_HELP="$(ngsdose_py ngsdose sinks --help 2>/dev/null)" || true
case "$SINKS_HELP" in *--stats*) ;; *)
  echo "ERROR: this image's ngsdose has no \`sinks --stats\` (the fae1124 image predates it): learn sinks with a newer image" >&2
  echo "  (SIF); the scans need not be made again" >&2; exit 1 ;;
esac
# the held-out statistics say so in their file ('# held-out: yes'), and fetchplan says whether its expected capture is held out
HELD_OUT=(); case "$SINKS_HELP" in *--held-out*) HELD_OUT=(--held-out) ;; esac

# the scans that carry every class, shuffled and split
ngsdose_py python3 - "$WORK_DIR/counts_scan" "$SINKS_SEED" "$SINKS_TRAIN" "$MIN_SCANS" "$TMP" $CLASSES <<'PY'
import glob, gzip, json, os, random, sys, zlib
d, seed, n_train, least, tmp, want = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6:]
files, kinds, bad, other = [], {}, 0, 0
for f in sorted(glob.glob(os.path.join(d, "*.json.gz"))):
    try:
        with gzip.open(f, "rt") as fh:
            c = json.load(fh)
    except (OSError, EOFError, zlib.error, ValueError):
        bad += 1
        continue
    have = {x["name"]: x["kind"] for x in c.get("classes", ())}
    if c.get("mode") != "scan" or any(n not in have for n in want):
        other += 1
        continue
    kinds.update((n, have[n]) for n in want)
    files.append(f)
print(f"{len(files)} of {len(files) + other + bad} scans in {d} carry {' '.join(want)}" + (f"; {bad} could not be read" if bad else ""), file=sys.stderr)
if len(files) < least:
    print(f"ERROR: at least {least} scans must carry the classes before their sinks are learned: keep them in the scans "
          "(CANDIDATE_PANELS) until there are", file=sys.stderr)
    sys.exit(1)
random.Random(seed).shuffle(files)
k = min(n_train, len(files) // 2)
with open(os.path.join(tmp, "scans.txt"), "w") as fh:
    fh.write(f"# seed {seed}: {k} training scans, {len(files) - k} held out\n")
    fh.writelines(f"{'train' if i < k else 'heldout'}\t{f}\n" for i, f in enumerate(files))
with open(os.path.join(tmp, "kinds.txt"), "w") as fh:
    fh.writelines(f"{n}\t{kinds[n]}\n" for n in want)
PY
train=(); heldout=()
while IFS=$'\t' read -r role f; do
  case "$role" in train) train+=("$f") ;; heldout) heldout+=("$f") ;; esac
done < "$TMP/scans.txt"
echo "learning the sinks of $CLASSES from ${#train[@]} scans; ${#heldout[@]} held out"

# shellcheck disable=SC2086
ngsdose_py ngsdose sinks "${train[@]}" --classes $CLASSES -o "$TMP/all.bed" --stats "$TMP/all.train.stats.tsv"
keep=" $CLASSES "
awk -F'\t' -v keep="$keep" 'index(keep, " " $4 " ")' "$TMP/all.bed" > "$TMP/$NAME.bed"
awk -F'\t' -v keep="$keep" '/^#/ || $1 == "class" || index(keep, " " $1 " ")' "$TMP/all.train.stats.tsv" > "$TMP/$NAME.train.stats.tsv"
for c in $CLASSES; do
  awk -F'\t' -v c="$c" '$4 == c {n = 1} END {exit !n}' "$TMP/$NAME.bed" \
    || { echo "ERROR: no bin of $c reached the rule in the training scans: $c has no sinks (too few reads?)" >&2; exit 1; }
done
ngsdose_py ngsdose sinks "${heldout[@]}" --evaluate "$TMP/$NAME.bed" --stats "$TMP/$NAME.stats.tsv" ${HELD_OUT[@]+"${HELD_OUT[@]}"} > "$TMP/$NAME.heldout.tsv"

# the capture of each class in the held-out scans, and the menu that makes the classes fetchable
ngsdose_py python3 - "$TMP" "$OUT" "$NAME" "$FETCH_MENU" "${#train[@]}" "$CANDIDATE_PANEL_FILES" "$BUNDLE/panel.k31.tsv.gz" $EXTRA_PANELS <<'PY'
import gzip, os, shutil, sys
tmp, out, name, base, n_train, cands, panels = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6].split(), sys.argv[7:]
kinds = dict(line.rstrip("\n").split("\t") for line in open(os.path.join(tmp, "kinds.txt")))
frac, unm, bp, nint = {c: [] for c in kinds}, {c: [] for c in kinds}, dict.fromkeys(kinds, 0), dict.fromkeys(kinds, 0)
with open(os.path.join(tmp, f"{name}.heldout.tsv")) as fh:
    next(fh)
    for line in fh:
        s, c, tot, _, f, u = line.rstrip("\n").split("\t")
        if c in frac and int(tot):
            frac[c].append(float(f))
            unm[c].append(int(u) / int(tot))
for line in open(os.path.join(tmp, f"{name}.bed")):
    p = line.split("\t")
    bp[p[3].strip()] += int(p[2]) - int(p[1])
    nint[p[3].strip()] += 1

def q(v, x):
    v = sorted(v)
    return v[min(len(v) - 1, int(x * (len(v) - 1) + 0.5))] if v else float("nan")

print("class\tkind\tintervals\tbp\theldout_scans\tcapture_median\tcapture_p10\tcapture_min\tunmapped_median")
summary = {}
for c in kinds:
    summary[c] = (q(frac[c], 0.5), q(frac[c], 0.1), len(frac[c]))
    print(f"{c}\t{kinds[c]}\t{nint[c]}\t{bp[c]}\t{len(frac[c])}\t{q(frac[c], 0.5):.5f}\t{q(frac[c], 0.1):.5f}\t{min(frac[c], default=float('nan')):.5f}"
          f"\t{q(unm[c], 0.5):.5f}")
if not os.path.isfile(base):
    print(f"note: no fetch menu at {base} (an image older than the menu): {name}.menu.tsv is not written", file=sys.stderr)
    sys.exit(0)

def opened(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)

owner = {}                                               # class -> the panel that defines it
for p in panels:
    if not os.path.exists(p):
        continue
    with opened(p) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##class\t"):
                n = dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:] if "=" in f).get("name", "")
                owner.setdefault(n, p)
lost = [c for c in kinds if c not in owner]
if lost:
    print(f"note: no panel named here (the bundle panel, EXTRA_PANELS, CANDIDATE_PANELS) defines {', '.join(lost)}, so {name}.menu.tsv "
          "is not written: run again with CANDIDATE_PANELS set to the panels the scans were made with", file=sys.stderr)
    sys.exit(0)
where = {}
for c in kinds:
    p = owner[c]
    if p in cands and p not in where:                    # a candidate on the host: copied beside its sinks
        where[p] = os.path.join(out, f"{name}.{os.path.basename(p)}")
        shutil.copyfile(p, os.path.join(tmp, os.path.basename(where[p])))
lines, head, rows = [], None, []
bdir = os.path.dirname(os.path.abspath(base))
for line in open(base):
    line = line.rstrip("\n")
    if line.startswith("#") or not line.strip():
        lines.append(line)
    elif head is None:
        head = line.split("\t")
    else:
        rows.append(dict(zip(head, line.split("\t"))))
if "stats" not in head:
    head.append("stats")
for r in rows:
    for col in ("panel", "sinks", "stats"):
        v = r.get(col, "-").strip() or "-"
        r[col] = v if v in ("-", ".") or os.path.isabs(v) else os.path.join(bdir, v)
note = f"sinks by 06_learn_sinks.sh {name} from {n_train} scans"
by = {r["class"]: r for r in rows}
for c in kinds:
    med, p10, n = summary[c]
    r = by.get(c)
    if r is None:
        r = {"class": c, "group": "candidate", "kind": kinds[c], "tier": "D", "presets": "-", "truth": "-",
             "measures": f"{c} ({os.path.basename(owner[c])})", "notes": "-"}
        rows.append(r)
    r.update(status=f"experimental: {note}, no fetch compared yet", panel=where.get(owner[c], owner[c]),
             sinks=os.path.join(out, f"{name}.bed"), stats=os.path.join(out, f"{name}.stats.tsv"),
             presets=",".join([p for p in r["presets"].split(",") if p.strip() not in ("", "-", name)] + [name]))
    r["notes"] = (r["notes"] if r["notes"].strip() not in ("", "-") else "") + \
        f"{'; ' if r['notes'].strip() not in ('', '-') else ''}held-out capture median {med:.5f}, p10 {p10:.5f} over {n} scans"
with open(os.path.join(tmp, f"{name}.menu.tsv"), "w") as fh:
    fh.write(f"# written by 06_learn_sinks.sh {name} from {base}, its paths made absolute; rows {', '.join(kinds)} and the preset {name} added\n")
    fh.write("\n".join(lines) + "\n")
    fh.write(f"##preset\t{name}\t{', '.join(kinds)}: {note}, checked on held-out scans\n")
    fh.write("\t".join(head) + "\n")
    for r in rows:
        fh.write("\t".join(r.get(col, "-") for col in head) + "\n")
PY
mv -f "$TMP/scans.txt" "$OUT/$NAME.scans.txt"
rm -f "$OUT/$NAME.menu.tsv"                             # a menu of an earlier run would name the old statistics' classes
for f in "$TMP/$NAME".*; do mv -f "$f" "$OUT/"; done
echo "wrote $OUT/$NAME.bed, .stats.tsv (held out), .train.stats.tsv, .heldout.tsv, .scans.txt$( [ -s "$OUT/$NAME.menu.tsv" ] && echo ", .menu.tsv")"
if [ -s "$OUT/$NAME.menu.tsv" ]; then
  echo "To fetch these classes as well: FETCH_MENU=$OUT/$NAME.menu.tsv FETCH_PRESET=\"core $NAME\" (config.sh; costs:"
  echo "  ngsdose fetchplan --menu $OUT/$NAME.menu.tsv --preset core $NAME --crai <indexes> --contigs $REF_FASTA.fai)"
fi
