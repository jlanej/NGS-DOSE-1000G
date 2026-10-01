# The cohort run

The expanded 1000 Genomes cohort (3,202 samples, 602 trios; NYGC, TruSeq PCR-free, NovaSeq
2×150, GRCh38) is the validation cohort for NGS-DOSE: it is public, it has trios, many samples
have HPRC assemblies, and NGS-PCA has already been run on it (`../meta/ngspca/`).

| path | what |
| --- | --- |
| `../pilot/` | four trios, each sample also as an independent older library; run locally; counts files, the evaluation and plotting scripts, the report and the figure are there |
| `00_setup.sh`, `01_stage_and_dose.sh`, `01a_dispatch_staged.sh`, `01b_dose_sample.sh`, `01_count.sh`, `02_cohort.sh`, `03_compare_modes.sh` (+ `compare_modes.py`), `04_hprc_satellites.sh`, `05_report.sh`, `06_learn_sinks.sh`, `config.sh` (+ `check_counts.py`) | the full-cohort run, for SLURM or a plain loop |
| `07_hprc_dj.sh`, `hprc_dj.py` | the distal junction in the HPRC release-2 assemblies: each haplotype FASTA fetched, screened for the panel's k-mers (`ngs-dose panel --report`) and aligned to the core-masked unit (minimap2), then the two tables the report reads (`../meta/dj_hprc/`; `--dj-assemblies`). `SAMPLES="HG00097 HG01891"` restricts it; about four minutes per haplotype, dominated by the download |
| `ngs-dose.def` | fallback Apptainer definition of the container image |
| `../meta/hprc_r2_censat.keys.txt`, `hprc_satellites.py` | the HPRC release-2 CenSat annotations (S3 keys; 205 samples, 200 of them in this cohort) and the comparison of satellite estimates against them |
| `../meta/ngspca/` | NGS-PCA output for this cohort (200 coverage PCs, `AUTO_HQ_median`), produced by [NGS-PCA's 1000G example](https://github.com/jlanej/NGS-PCA/tree/master/example/1000G_highcov) |

**Upgrading a cohort that is already running** (the 1000 Genomes run, on the `sha-fae1124` image):
stop the running stager first (`scancel` its job, or end its tmux session), then start this
version's. It resumes the partial downloads the old one left, and leaves alone any modified in the
last 15 minutes. Keep the image the cohort was counted with and say so:
`export NGSDOSE_IMAGE=docker://ghcr.io/jlanej/ngs-dose:sha-fae1124 EXPECTED_ENGINE_BUILD=fae1124`.
`FETCH_PANELS` needs nothing: left unset it fetches the telomeric repeat only when the image's bundle
has its sinks, so on `sha-fae1124` it fetches 45S, 5S and DJ as before, with a note. Jobs already
queued keep the script they were submitted with (`sbatch` copies it at submission) and source the
`config.sh` of the directory they were submitted from (for this cohort, NGS-DOSE's
`example/1000G` at fae1124): they run the old checks, with no `gzip -t` on the index, no check that
the fetch saw the scan's reads and no build check. Let them finish, then check the pairs they wrote
(the image's `sinks.bed` is the one those fetches used; a failing sample is printed):

```bash
source config.sh
for f in "$WORK_DIR"/counts_fetch/*.json.gz; do s=$(basename "$f" .json.gz)
  ngsdose_py python3 "$EX_DIR/check_counts.py" same-reads "$WORK_DIR/counts_scan/$s.json.gz" "$f" \
    "$BUNDLE/sinks.bed" || echo "$s"
done
```

or `scancel` them and let the new stager submit them again. The 1,748 pairs counted before the
upgrade all pass this check.

## The pilot

```bash
REF=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa bash ../pilot/run_pilot.sh
```

counts every sample in fetch mode straight from the public CRAMs (AWS Open Data mirror for the
NYGC data, EBI for the HGSVC replicates) and writes [`pilot/pilot_report.md`](../pilot/pilot_report.md).
`python pilot/evaluate_pilot.py` (from the repository root) alone regenerates the report from the committed counts files.

The replicates are what make the pilot informative. HG00512/3/4, HG00731/2/3 and NA19238/39/40
were sequenced years earlier for the HGSVC (HiSeq 2500, 2×126, ~78×, bwakit + postalt), and the
CEU trio for Illumina's Platinum pedigree (HiSeq 2000, 2×101, ~54×): different library
preparation, chemistry, read length, insert size, depth and alignment pipeline, and a GC
response that is the reverse of NYGC's. Agreement between the two is reproducibility of the
*measurement*, not of the file — and an upper bound on its error, because the two DNA batches
come from different cultures of the cell line.

`python pilot/evaluate_pilot.py --write-anchors` additionally rewrites the method's
`resources/GRCh38/anchors.json` (the installed bundle) from the replicate pairs; `python pilot/plot_pilot.py` draws
`pilot_figure.png` (needs matplotlib).

## The cohort

Full accuracy for everyone: every CRAM is staged, scanned whole, fetched as well (about five
seconds more on a local file), verified, and removed.

Nothing is installed on the cluster: the image holds the engine, the `ngsdose` package, the
resource bundle and `aria2c`; these scripts come with this repository; the host needs Apptainer,
SLURM, `curl` and the coreutils. (From a checkout of NGS-DOSE instead: leave `SIF` unset, point
`NGSDOSE_SRC` at it, `cargo build --release`, `pip install .`.)

```bash
apptainer pull ngs-dose.sif docker://ghcr.io/jlanej/ngs-dose:sha-<commit>   # pinned: see "The image" below
git clone --depth 1 https://github.com/jlanej/NGS-DOSE-1000G   # the scripts submit jobs, so they live on the host
cd NGS-DOSE-1000G/pipeline                              # submit from here: config.sh and logs/ are relative to it
export SIF=/path/to/ngs-dose.sif WORK_DIR=/scratch/$USER/ngs_dose_1000G EXPECTED_ENGINE_BUILD=<commit>

bash 00_setup.sh                                        # reference, pedigree, manifest (sample, URL, MD5)
LIMIT=5 bash 01_stage_and_dose.sh                       # smoke test: five samples staged, counted, removed
sbatch 01_stage_and_dose.sh                             # the cohort: a few transfers at a time, one 01b_dose_sample.sh job per CRAM

MODE=scan sbatch 02_cohort.sh                           # estimate, calibrate, adjust, transmission - on the scans
MODE=fetch sbatch 02_cohort.sh                          # the same on the targeted fetches, for the comparison
bash 03_compare_modes.sh                                # sink capture per sample, sinks re-learned, fetch / scan per sample
sbatch 04_hprc_satellites.sh                            # satellite array mass against the HPRC assemblies of the same people
sbatch 07_hprc_dj.sh                                    # the distal junction copy by copy in those assemblies (SAMPLES="..." for a few)
python3 hprc_dj.py screen $WORK_DIR/hprc_dj -o ../meta/dj_hprc   # -> the tables regenerate.sh passes to the report
sbatch 05_report.sh                                     # the page: known truths, fetch vs scan, trios, cell line, satellites - from whatever exists
```

**Publishing as the run proceeds.** The counts files are small (~240 kB per scan, ~70 kB per fetch;
under 1 GB for the cohort), public-data derivatives with no reads or genotypes in them, and the
primary product of the run: from them everything else is a minute's computation. They go, as they
land, into this repository: `counts_scan/`, `counts_fetch/`, and `docs/` for the page that
`05_report.sh` (or `regenerate.sh` at this repository's root) builds from them, served by GitHub
Pages. `rsync` the two counts directories from `$WORK_DIR` into it, regenerate, commit.
The page is honest about how much of the cohort it rests on, and each section appears when the
data for it exist (three trios for inheritance, twenty for its intervals, sixty samples for the
PC sweep, assemblies for the satellites).

`01_stage_and_dose.sh` keeps a few multi-connection `aria2c` transfers going (the image's
`aria2c`; each file goes to `<file>.part` and takes its name only once whole: the CRAM checked
against the MD5 of the sequence index, its index by `gzip -t`), holds at most `MAX_LOCAL_CRAMS`
files on disk, submits one job per landed CRAM and can be re-run at any time: finished samples and
queued jobs are recognised, and each run first verifies and re-submits the CRAMs already on disk
(a sample whose job ended without counts is submitted again, up to `MAX_TRIES` jobs, 3 by default)
before it downloads anything. A job counts as ended only when `squeue` answers without it; while
`squeue` does not answer (a busy controller times out), every job is taken as alive, so nothing is
submitted twice. The stager waits for room only before a download, and stops, naming the CRAMs,
when `MAX_KEPT_CRAMS` (10) of them have been kept by failed jobs, or when the disk is full and no
job on it could free room, each on five checks in a row a minute apart: a failure that repeats (an
image without a panel, an unbound reference) is fixed at its cause, in `logs/dose_<job>.out`, and
the stager run again. One stager runs at a time (`$WORK_DIR/stager.lock`, where `flock` exists).
Compute nodes need no internet unless the manager itself runs as a job; on a cluster where they
have none, run it on a login or transfer node (`bash 01_stage_and_dose.sh` in tmux).

`01b_dose_sample.sh SAMPLE LINE [MD5]` has the calling convention of NGS-PCA's
`01b_mosdepth_sample.sh` (verify the MD5 if given, process, check the outputs, delete the CRAM),
so the aria2 download manager that staged this cohort for NGS-PCA can drive it: it needs its
per-sample job script and its "already done" test to be settable (`$WORK_DIR/counts_scan/<sample>.json.gz`
instead of the mosdepth output), nothing else. `01a_dispatch_staged.sh` is for files staged by
other means (Globus, rsync): it neither downloads nor deletes, skips files aria2 is still
writing and indexes that are not yet whole, and submits each landed CRAM once, again only when
its job has ended without counts (up to `MAX_TRIES` jobs; a sample with a job is left for the next
run while `squeue` does not answer); the dispatch mark is written only when `sbatch` accepted the
job, and a refused submission ends the run with a non-zero status. About 5
minutes on 8 CPUs (median 304 s over 1,748 scans) and 1.3 GB of memory per sample; the counts files
are ~230 kB (scan) and ~70 kB (fetch), so the whole cohort is under 1 GB.

Without staging, `01_count.sh` counts straight from the public bucket, a block of the manifest
per array task: `MODE=fetch` moves ~0.5 GB per sample (1.5 TB for the cohort, about a minute
each) and is how a biobank would run; `MODE=scan` streams the whole file (~15 GB per sample over
one connection). The array has to cover the manifest at `SAMPLES_PER_TASK` samples per task (a task
of an array from 0 that falls short stops with an error):

```bash
N=$(wc -l < $WORK_DIR/manifest.tsv)
MODE=fetch sbatch --array=0-$(( (N - 1) / ${SAMPLES_PER_TASK:-10} ))%25 01_count.sh
```

**What the full scan buys, and why this cohort should have it.** A scan is placement-independent:
it counts every class read wherever the aligner put it, and it is how sinks are learned. In this
pipeline it is the only mode that measures the ten satellite families, three of them particular to
the acrocentric short arms and pericentromeres (`EXTRA_PANELS` in `config.sh`; 200 samples of the
cohort have HPRC assemblies to hold the satellites against). `FETCH_PANELS` fetches only the
telomeric repeat, which the aligner concentrates at the chromosome ends, through the bundle's sinks.

That is a matter of what ships, not of what a fetch can do. In these CRAMs, sinks learned from 30
scans hold at least 99.8% of the reads of nine of the families in every one of 200 other genomes
(99.85% in two of three random draws; aSatHOR needs 59.6 Mb of intervals, the others 0.2-3.5 Mb
each); HSat1B reaches only 96.6%, part of it wholly unmapped in the 698-sample batch. No satellite
sinks ship in the bundle, so `FETCH_PANELS` stays with the telomeric repeat; a fetch plan can read
the families' experimental sinks, learned from 100 of these scans ("New classes, and what a fetch
reads" below), though no such fetch has yet been compared with its scan. Sinks belong to an aligner
and a reference, and are learned again for another pipeline. A class loaded in these scans can be fetched later from the public CRAMs once its
sinks are learned from them, and so can a new class once a handful of new scans has said where its
reads land. What a later fetch cannot recover without staging 48 TB again is the reads each genome
has outside the sinks, and a class whose reads do not settle in a stable set of intervals.

The scan also records, on the 1-kb grid of `mosdepth --by 1000` (10 kb for the satellite and
telomere classes), where every class read was aligned *and every other read in those bins* - which
is the data that decides how far the measurement can be simplified for a biobank:

1. **Targeted fetch** (exists): with both modes run on every sample, `03_compare_modes.sh` gives
   fetch / scan per sample across 26 populations and both sexes, and re-learns the sinks at
   n = 3,202 on the finer grid (tighter sinks, fewer requests).
2. **A depth proxy from bins a cohort already has.** NGS-PCA's mosdepth run over this cohort
   (1-kb bins, all contigs, no MAPQ filter, duplicates excluded) is kept. In NA12878, 99.8% of
   the aligned 45S reads sit in 273 of those bins. Whether a fixed set of such bins, normalised
   to the autosomal median, tracks the full estimate - and at what loss of transmission
   reliability, which `ngsdose trios --compare-to` measures directly - needs exactly three
   things: the scans, the kept mosdepth outputs, and the per-bin census that says how much of
   each bin's depth is the class and how its duplicate-flag rate differs from the genome's
   (5.5% against 10.8% in NA12878: a depth tool that drops flagged reads over-reads rDNA by the
   difference, by a different amount in every sample). All three exist after this run.
3. Anything else - fewer or smaller sinks, fewer controls, lower depth - can be tried on the
   counts files alone, because they hold positions, not summaries.

Operational notes:

- `01b_dose_sample.sh` is idempotent under a requeue and removes a CRAM (`KEEP_CRAMS=1` keeps
  them) only once both counts files exist and are whole, were written by `EXPECTED_ENGINE_BUILD`
  when that is set, and agree: the fetch's `ctrl_reads` and every region's reads must equal the
  scan's, and each class must have at least the reads the scan placed inside the sinks and no
  more than the scan's total, as they do by construction (in all 1,748 pairs so far). A fetch that
  differs - the mark of a cut index, which htslib reads without an error and which loses every
  contig past the cut - is removed and the CRAM kept; the index is `gzip -t`-tested before
  anything is counted. The CRAM is kept whenever anything failed. Counts are removed only while
  the CRAM and a whole index are there to make them again: a job that finds them without the CRAM
  (a requeue after the CRAM was removed) moves a fetch that disagrees aside to
  `counts_fetch/<sample>.json.gz.rejected`, so that the stager stages the CRAM again, and leaves a
  file of another engine build where it is. A second job for the same sample waits for the first
  (`$WORK_DIR/dispatched/<sample>.lock`, where `flock` exists). `compare_modes.py` lists any sample whose known-truth or
  dosage ratio between the modes is not 1.
  `01_count.sh` is idempotent too (finished samples are skipped), tries each sample three times,
  and exits non-zero if any sample still failed; re-submit to retry those. The engine retries
  failed intervals itself, gives up with exit status 75 when a connection has gone silent for
  five minutes (a dead HTTPS connection waits for ever rather than failing; `timeout` is a
  backstop where it exists, and `--stall-timeout 0` turns the watchdog off, e.g. for files that
  have to be recalled from tape), writes its output under a temporary name and renames it, and
  refuses a BAM/CRAM without an end-of-file marker (a truncated copy decodes without error and
  loses exactly the contigs the rDNA is on).
  `%25` caps concurrent array tasks: the bottleneck is the network, and hundreds of parallel
  readers against one S3 prefix earn throttling rather than speed.
- The manifest's second column may be a local path instead of a URL (e.g. CRAMs staged by
  NGS-PCA's download stage); fetch mode then looks for `<cram>.crai` next to it.
- With `export SIF=/path/to/ngs-dose.sif` both the engine and the Python steps run through
  Apptainer. `00_setup.sh` pulls the image if the file is not there yet (`NGSDOSE_IMAGE`,
  default `docker://ghcr.io/jlanej/ngs-dose:latest`, published by NGS-DOSE's
  `.github/workflows/container.yml` on every push to its `main`, with `:sha-<commit>`, and as
  `:X.Y.Z` on version tags; a package that has not been made public needs `apptainer remote
  login` first). Where nothing has been published, `ngs-dose.def` builds the same image from the
  root of an NGS-DOSE checkout (`apptainer build --fakeroot --build-arg NGSDOSE_BUILD=$(git
  rev-parse HEAD) ngs-dose.sif /path/to/NGS-DOSE-1000G/pipeline/ngs-dose.def`). The scripts bind
  `$WORK_DIR`, the repository, the reference directory and `$CRAM_DIR` (`APPTAINER_BINDS` in
  `config.sh`). Without a container: in the NGS-DOSE checkout, `cargo build --release` (needs
  libclang, e.g. `module load llvm`) and `pip install .`.
- **The image.** A cohort should keep one engine, and the scripts have to match the image's
  bundle. Pin `NGSDOSE_IMAGE` to `:sha-<commit>` (or `:X.Y.Z`, or a digest) rather than `:latest`,
  and set `EXPECTED_ENGINE_BUILD` to that commit (a prefix is enough): a counts file from any other
  build is then removed and its job fails with the CRAM kept, and the stager notes which build the
  counts so far carry. Before anything is staged or counted, every script that fetches checks
  each `FETCH_PANELS` class against the image's `sinks.bed`. Left unset, `FETCH_PANELS` loads the
  telomere panel only when the bundle has its sinks and otherwise leaves it out with a note; set
  explicitly, a panel whose classes have no sinks stops the script. The cohort's first 1,748
  genomes were counted with `ngs-dose:sha-fae1124`, whose bundle predates the telomere sinks: on
  that image the default fetches 45S, 5S and DJ, as those genomes were fetched (set
  `EXPECTED_ENGINE_BUILD=fae1124` as well). A newer image adds the telomeric repeat to the fetch,
  but the cohort then mixes engine builds.
- Every download (`00_setup.sh`, the stager, `01_count.sh`, `04_hprc_satellites.sh`) goes to
  `<file>.part`, is resumed after an interruption and checked (an MD5, `gzip -t`, or for the CenSat
  BED files at least four fields on every record), and only then takes its name.
- Every counts file records the SHA-256 of the panel, controls and sinks it was made with and
  the lengths of the contigs it used; `ngsdose estimate` warns if a cohort mixes resource sets
  and refuses a file aligned to another reference build.
- QC columns to look at first in `cohort.tsv`: `truth.auto` (2), `truth.chrX` and `truth.chrY`
  (sex; 1.6 and the like is mosaic loss, common in LCLs), `DJ.cn` (10), `flagged_chromosomes`
  (aneuploidy), `eof_marker`, `gc_curve_max_se`.

### New classes, and what a fetch reads

**Scan first.** A new class is never fetched from where the reference says it is. Its reads are wherever
the aligner put them: on paralogs, decoys, unplaced contigs, or nowhere (the unmapped bin), and that
differs between aligners and references (with DRAGEN and an alt-masked reference most 45S and DJ reads of
the one genome checked were unmapped). So a new class is first a k-mer panel loaded in whole-file scans;
where its reads land (its sinks) is learned from those scans and checked on other scans; only then can a
fetch read it. On this cohort:

1. **Scan the candidates.** The candidate panels are built and final: `resources/experimental/candidates/`
   of NGS-DOSE holds seven files (83 classes, 5,076,728 k-mers; `candidates.tsv` there says what each class
   measures, its truth and its tier, and `resources/build/build_candidate_panels.sh` rebuilds them byte for
   byte). The fae1124 image does not contain them, so copy that directory, or the files wanted, to the
   cluster, and check that they arrived whole: each counts file records the sha256 of the panels it was
   counted with (`panel_sha256`), and `ngsdose estimate` finds a candidate's unit by it.

   | file | k-mers | sha256 |
   | --- | --- | --- |
   | `coding-vntrs.k31.panel.tsv.gz` | 86,162 | `46bf3ddd0b48031d2d968bf433fe2a36d2cd7bdde8c71c4c1ec2acae39917a29` |
   | `macrosatellites.k31.panel.tsv.gz` | 45,362 | `385bf0d44373a563770e3e4cb6cfc9bb272d202bc582c9a87d3b06be5c5eaf6f` |
   | `multicopy-genes.k31.panel.tsv.gz` | 544,098 | `a7dafd323ae5c23dc825bd9f8375be9123d7883f906c2e20208fca4b5cfd2aa8` |
   | `nonhuman-myco.k31.panel.tsv.gz` | 3,731,523 | `ad14926aec70efba19cf87e60f795699b6a8487b58fba76df8d19db0a644450f` |
   | `nonhuman.k31.panel.tsv.gz` | 395,211 | `8c4a653669fc45979785203d27880e5d4ec89b58491e11d91912b9aff8070f4f` |
   | `rna-arrays.k31.panel.tsv.gz` | 148,026 | `c33c4a131c8c5bee0b4d260619f94a265226f49af12193401c31d89ff1b9d1d3` |
   | `sex-chromosome-arrays.k31.panel.tsv.gz` | 126,346 | `d103151637f200c147dedf90aa6a2b5640b65e316f542bbc2a554745d1ff1aa8` |

   MYCO has a file of its own (`nonhuman-myco.k31.panel.tsv.gz`, 3.7 M k-mers), which is most of the added
   memory; leave it out if culture mycoplasma is not wanted. Then, before (re)starting the stager:

   ```bash
   export CANDIDATE_PANELS=/path/to/candidates     # a directory's *panel.tsv.gz files, or files; space-separated
   ( cd /path/to/candidates && sha256sum *.k31.panel.tsv.gz )      # compare with the table above
   ```

   The candidates join `EXTRA_PANELS` in every scan from then on, and their directories join
   `APPTAINER_BINDS`. The fetch never loads them. Every script that counts first checks that no candidate
   shares a k-mer or a class name with the bundle panel, `EXTRA_PANELS` or `FETCH_PANELS`, or with another
   candidate, and stops if one does: the engine drops a k-mer that two loaded panels claim from both, which
   would change the counts of the classes the cohort already measures. The check is standard-library Python
   in the image and needs no data. On the seven final candidate files (5,076,728 k-mers, checked against the
   bundle, satellite and telomere panels) it took 4.2-4.3 s and 0.67 GB in three runs, and a pass is
   remembered in `$WORK_DIR/candidate_panels.checked` until a panel file or the image changes. The engine
   and the image stay those of the cohort. With the seven final files loaded, fae1124 scanned NA12878 (30×
   NYGC) in 119.5 s with a peak of 1.65 GB, against 106.1 s and 1.43 GB without them (one run on a shared
   machine; runs over the review rounds gave 109.6-119.5 s and 1.49-1.81 GB with them), and every shipped
   class had the same reads and placements as without them, in fae1124 and in the NGS-DOSE branch engine
   alike; so did the test fixture, with three engine builds. That follows from how the engine classifies a
   read. Each positional class (61 of the candidates) is scored on its own, so a positional candidate never
   changes a shipped class. The compositional classes compete: a read goes to the one with most hits, or
   is counted as ambiguous when the runner-up has more than a fifth of the best's hits. So only the 22
   compositional candidates (the 19 coding VNTRs, DYZ19, MYCO and EBV2) could take a read from a shipped
   satellite or TEL class; none did in NA12878 (ambiguous reads 3,214 with and without them). The check of
   each fetch against its scan would not catch it: when the scan loaded panels the fetch did not, a fetch
   class with more reads than the scan gets a note, not a failure. The genomes counted before the candidates
   were added do not carry them, and are not scanned again.

2. **Learn their sinks** once at least 30 scans carry them (it refuses fewer):

   ```bash
   bash 06_learn_sinks.sh ychr TSPY RBMY DAZ       # NAME, then the classes; or sbatch it
   ```

   The scans that carry every class are shuffled (seed 7). Thirty of them (at most half) learn the sinks
   (`ngsdose sinks --classes ... --stats`): a 10-kb neighbourhood becomes a sink where it holds at least
   1e-5 of a class's reads and at least 25 reads in any training scan, and the sink keeps the placement
   bins inside it that hold the class's reads (1 kb for a positional class, 10 kb for a compositional
   one), merged and padded by 1 kb. All the other scans are held out and say what those sinks capture
   (`ngsdose sinks --evaluate ... --stats`, with `--held-out` when the image has it). It prints, per class, the intervals, the bp, and
   the held-out capture (median, 10th percentile, minimum, and the median share left unmapped, which only
   the unmapped bin can recover). It writes `$WORK_DIR/sinks/NAME.bed`, `NAME.stats.tsv` (held out, and
   marked so: what `FETCH_CAPTURE` trims by), `NAME.train.stats.tsv` (in-sample), `NAME.heldout.tsv` (per
   scan), `NAME.scans.txt`, and `NAME.menu.tsv`: the image's fetch menu, with its paths made absolute, a
   row per class (status experimental; tier D for a class new to the menu, while a row the menu already
   has keeps its tier) and a preset NAME. A candidate panel is copied beside its sinks, so that the menu
   does not depend on `CANDIDATE_PANELS` any more. Learning needs an image whose `ngsdose` has
   `sinks --stats` (not fae1124); it reads counts files only, so it can use a newer image than the one
   counting the cohort. A class whose held-out capture is low, or whose reads are largely unmapped, stays
   scan-only. These sinks are those of the NYGC bwa-mem alignments to the GRCh38 analysis set.

3. **Fetch them** on other CRAMs of the same pipeline (NYGC bwa-mem, GRCh38 analysis set):
   `FETCH_MENU=$WORK_DIR/sinks/NAME.menu.tsv FETCH_PRESET="core NAME"` (below). With an image newer than
   fae1124, `01b_dose_sample.sh` fetches every staged CRAM of a cohort still being scanned by that plan and
   checks the fetch against the scan through the plan's sinks: that comparison is what an experimental row
   still lacks. No such fetch has been run yet.

The same script re-learns the sinks of a class the scans already carry (the satellite families, TEL) for
this cohort: `bash 06_learn_sinks.sh sat SST1 CER`. On 80 of the cohort's scans (30 to learn, 50 held out)
it gave SST1 34 intervals (607 kb, held-out capture median 0.99972, 10th percentile 0.99949) and CER 49
(1.16 Mb, 0.99976 and 0.99970).

**Customising the fetch.** Left unset, the fetch is what it has always been: the control regions and
`$BUNDLE/sinks.bed`, with the bundle panel and `FETCH_PANELS`. Set any of the following and the fetch reads
what a plan says instead. `ngsdose fetchplan` (NGS-DOSE) builds the plan once from the fetch menu
(`resources/fetch_menu.tsv` in the image: one row per class or sub-option with its panel, sinks, status,
tier and presets) into `FETCH_PLAN_DIR`, and every script that fetches uses its sinks BED, controls, panels
and flags. The image must have `ngsdose fetchplan`. The fae1124 image has not, so on it these variables
stop every script with a message. Export them, as `CANDIDATE_PANELS`: the jobs the stager submits read them
from the environment they inherit, and a job without them would fetch as before.

| variable | what |
| --- | --- |
| `FETCH_CLASSES` | options by name, e.g. `"rDNA45S rDNA5S DJ TEL SST1"`, `DXZ1` for a sub-option, or `unmapped` for the reads without a coordinate |
| `FETCH_PRESET` | presets of the menu: `core` (45S, 5S, DJ: the fetch of the fae1124 image), `core_tel` (+ TEL), `truths`, `satellites`, `biobank_lite`, `xy_arrays`, and a NAME from `06_learn_sinks.sh` |
| `FETCH_BUDGET_MB` | add options in tier order (A to D), cheapest first within a tier, until the plan would read more than this many MB of a CRAM (median over the costing indexes; the controls count); from the options above if given, otherwise from the whole menu |
| `FETCH_CAPTURE` | per class, keep only the sink intervals of highest yield until the expected capture (10th percentile of the statistics' scans) reaches this fraction, e.g. `0.995` |
| `FETCH_PLAN_ARGS` | any other `fetchplan` flags: `--fill` (a budget skips what does not fit and goes on), `--capture-class TEL=0.99`, `--capture-stat median`, `--stats FILE`, `--sinks BED ...`, `--controls $BUNDLE/controls.lite200.bed` |
| `FETCH_MENU` | the menu (default: the image's); a copy with other tiers or presets, or one written by `06_learn_sinks.sh` |
| `FETCH_CRAI`, `FETCH_COST_SAMPLES` | the indexes the plan is costed on; by default the first 5 manifest samples' `.crai` (downloaded to `$WORK_DIR/crai`), with `$REF_FASTA.fai`. A budget needs them, and with them a capture target ranks intervals per byte |
| `FETCH_PLAN_DIR` | where the plan is kept (default `$WORK_DIR/fetchplan`) |
| `FETCH_SINKS_KNOWN` | sinks BEDs of earlier fetches (space-separated), for `ngsdose estimate` (below) |

To see what a choice costs before any of it is fetched, print the plan without writing one:

```bash
source config.sh
ngsdose_py ngsdose fetchplan --menu "$FETCH_MENU" --preset biobank_lite --capture 0.995 \
  --crai "$WORK_DIR"/crai/*.crai --contigs "$REF_FASTA.fai"
```

The plan is built by the first script that needs it (`00_setup.sh`, the stager, `01a_dispatch_staged.sh`,
`01_count.sh`), under a lock, and kept with the settings and the SHA-256 of every file it came from
(`settings.txt`). A later run whose settings or files differ stops and shows the difference: a cohort's
fetches should share one plan (each counts file records the SHA-256 of the sinks BED it used). To change the
plan, remove `FETCH_PLAN_DIR` or name a new one; the fetches made so far keep theirs. Build the plan on a
node with internet (the stager's), since it downloads the costing indexes. `FETCH_PANELS` and a plan
cannot both be set. What the scripts take from the plan:

- **Controls.** The fetch reads the control FASTA the plan names (`plan.controls.txt`, `FETCH_CONTROLS`):
  `$BUNDLE/controls.fa.gz`, or `controls.lite200.fa.gz` for a plan made with
  `--controls $BUNDLE/controls.lite200.bed`. A controls set the bundle does not name is refused, since
  `ngsdose estimate` would refuse every counts file made with it. The scan always reads all 800 controls.
  Before a CRAM is deleted, `01b_dose_sample.sh` checks the fetch against the scan
  (`check_counts.py same-reads SCAN FETCH SINKS CONTROLS`): the fetch must have used the plan's controls,
  and with the lighter set each of its regions (382 with lite200) is compared with the scan's count of the
  same region. A fetch made with other controls than the plan's is removed and made again while the CRAM
  is still there; a region set that is not a subset of the scan's cannot be compared, and then both counts
  files and the CRAM are kept.
- **Count flags.** `plan.count_flags.txt` holds `--unmapped` when the unmapped bin is selected. When the
  plan loads a panel only part of whose classes it selects (SST1 alone out of the satellite panel),
  `fetchplan` asks the engine that will fetch (`count --help`) what it takes. An engine with `--classes`
  (NGS-DOSE's current one) gets `--classes=` and the selected classes: the fetch counts and reads only
  those. One with only `--allow-missing-sinks` (NGS-DOSE 7772e32, main since 645ae55) gets that flag: the
  other classes are counted only where their reads fall inside the plan, the counts list them in
  `sinks_missing_classes`, and `ngsdose estimate` marks them `no_sinks_in_fetch`. The cohort's fae1124
  takes neither, and `fetchplan` refuses such a plan for it (`--unmarked-companions` would write it anyway,
  leaving the other classes undercounted and unmarked unless `estimate` gets `--fetch-sinks`). With
  fae1124, select whole panels: `core` and `core_tel` load no class they do not select and need no flag.
  The scripts pass `--engine $NGSDOSE_BIN` to `fetchplan` without an image, and refuse a plan whose flags
  hold `--classes` when the engine that will fetch does not take it (fae1124).
- **Sub-options.** A plan that fetches `DXZ1`, `DYZ3`, `DYZ1` or `DYZ2` without the rest of its family
  counts the family only there; `ngsdose estimate` then reports `DXZ1.mass_Mb` and so on and marks the
  family `subset_only`. It knows a fetch's sinks BED only by the SHA-256 the counts file records, so
  `02_cohort.sh` and `03_compare_modes.sh` pass `FETCH_SINKS_KNOWN`, or else the plan's `plan.sinks.bed`,
  as `--fetch-sinks` when the image's `ngsdose estimate` takes it. Without it, for a fetch through a BED
  `estimate` does not know, the sub-options are `unverified` (reads, no mass), and so are the three
  families that have sub-options (`aSatHOR`, `HSat1B`, `HSat3`): NaN wherever the fetch counted them,
  even in a plan that selected no sub-option. The report (`05_report.sh`) does not pass it yet, so for
  plan fetches its page loses those three families as well as the sub-options.

What the choices cost, per 30× NYGC CRAM, is in NGS-DOSE's [`docs/fetch_examples.md`](https://github.com/jlanej/NGS-DOSE/blob/main/docs/fetch_examples.md): its
`resources/build/fetch_examples.sh` runs `fetchplan` for a set of worked examples (the presets, capture
targets, the lite controls, the X and Y arrays, budgets; the `biobank_lite --capture 0.995` plan above is
its example 14) on the `.crai` of 13 genomes of this cohort and writes each plan's table, count flags and total (MB = 1e6 bytes of CRAM slices and their containers'
compression headers, median over the 13). There, `core_tel` reads 522.7 MB, 3.09% of the file. These
costs, like the sinks, are those of this pipeline's files (NYGC bwa-mem, GRCh38 analysis set): another
pipeline's are read off its own indexes.

Capture targets cut most where a class's reads are concentrated. Every class of the menu is trimmed by
held-out statistics that ship with the image: the bundle's (`sinks.stats.tsv`, 1,375 scans of this cohort
not used to learn them) and the satellites' (`sinks.satellites.stats.tsv`, 1,648), both ranked per byte of
their CRAM slices when the plan has indexes. The sub-options are always fetched whole, and at the default 10th-percentile
statistic `HSat1B` never reaches 0.99 (97.17% with all its sinks), so a higher target keeps all of it; with
`--capture-stat median` its full capture is 0.990, so a target up to 0.990 trims it. An
interval the plan reads for one option is kept for every class that has it, at no extra cost. The controls
are the floor: the lighter `controls.lite200` (200 of the 800 control regions and all the truth and dosage
regions) reads 93.1 MB instead of 231.4. Its effect on the estimates comes from GC tables rebuilt from the
scans (45S +0.21%, SD 0.30%; 12 of the 32 aneuploidy flags missed), not yet from fetches of whole CRAMs,
and a cohort should not mix it with the full set.

At biobank scale the order is the same, with the biobank's own aligner and reference: scan a small share of
the genomes with every panel of interest (0.1-1% of a biobank is still thousands), learn and check the sinks
there (`ngsdose sinks`, as `06_learn_sinks.sh` does), and fetch the rest by a plan whose cost is read off
the biobank's own indexes. The sinks and costs above serve that pipeline's CRAMs only.

### Chromosomes in copies: the windows for every genome

NGS-DOSE 0.3.0 reads every chromosome in copies from the single-copy regions (its DESIGN.md, section 8),
and its bundle adds the *karyotype windows*: 2,265 pieces of clean single-copy sequence along every arm,
which read every autosome to about 0.006 copies (standard error at 30x) where the controls alone leave
chromosome 19 at 0.028 and chromosome 22 at 0.035. The cohort was counted before they existed. The page reads 186 genomes
with them (`counts_karyotype/`: the containers of the public CRAMs that hold the regions, fetched by exact
byte range with `analysis/karyotype/slice_fetch.py` and counted locally with `ngs-dose count -m fetch`)
and the other 3,016 from the 980 regions their counts hold. To read every genome with the windows, fetch
them, with an image of 0.3.0 or later and the cheapest class beside them (a plan needs one):

```bash
export NGSDOSE_IMAGE=docker://ghcr.io/jlanej/ngs-dose:sha-<commit>    # 0.3.0 or later, pinned
export SIF=$WORK_DIR/ngs-dose.karyotype.sif EXPECTED_ENGINE_BUILD=<commit>
bash 00_setup.sh
N=$(wc -l < $WORK_DIR/manifest.tsv)
MODE=fetch COUNTS_DIR=$WORK_DIR/counts_karyotype FETCH_CLASSES=rDNA5S FETCH_PLAN_ARGS="--controls all" \
  FETCH_PLAN_DIR=$WORK_DIR/fetchplan_karyotype sbatch --array=0-$(( (N - 1) / ${SAMPLES_PER_TASK:-10} ))%25 01_count.sh
```

`--controls all` reads every region of the bundle, the 800 controls with the windows, as a scan counts
them: 404 MB of CRAM slices per genome as the engine decodes them (median of 13 NYGC indexes; 1.3 TB for
the cohort); `--controls karyotype` (200 controls and every window, 277 MB) reads the chromosomes almost as
well and misses a few short stretches that the 800 controls define. Over HTTPS htslib moves several times
what it decodes (each query opens a range request without an end, and what is in flight at the next seek is
thrown away: 1,517 MB moved for 410 MB of slices, measured), so on a network-bound cluster stage the
containers instead: `analysis/karyotype/slice_fetch.py URL CRAI REGIONS.bed OUT.cram` writes the CRAM of
exactly the containers that hold the regions (about 450 MB per genome for every region with 600 bp of
padding), which `ngs-dose count -m fetch -i OUT.cram` counts like the remote file. Then copy the counts
into this repository's `counts_karyotype/` and run `regenerate.sh`: a genome with a file there is read from
it. With every genome counted so, NGS-DOSE's `resources/build/karyotype_model.py` can learn the bundle's
model of the windows on all 3,202 instead of 186.

### What the cohort run is for

`02_cohort.sh` ends with the table the design is waiting on — transmission reliability of every
candidate estimator (`rDNA45S.18S.flat`, the estimator in the literature; the single-sample
anchor estimate; the calibrated `cn`; 5S; DJ) and of three traits whose answer is known in
advance (`truth.auto`: no variance but error; `chrEBV.copies` and `chrM.copies`: large variance,
none of it transmitted through the nuclear genome - if these come out "reliable", families
share batches and every other reliability in the table is inflated by as much), unadjusted, adjusted on NGS-PCA's coverage PCs and
adjusted on the internal control PCs, each with a family-bootstrap confidence interval, the
spousal correlation, a permuted-family null, and the *paired* bootstrap of every estimator
against the 18S depth ratio (separate intervals are about ±0.09 at 602 trios; the paired
difference is much sharper). The number of PCs is the Marchenko–Pastur default, separately for the two PC
sets (`N_PC=20 N_CTRL_PC=5 sbatch 02_cohort.sh` overrides them), and `pcsweep.*.tsv` with its `.summary.txt` is the evidence for or
against that default: for every number of PCs, the cross-validated error of the known truths
and the reliability of the classes.

The results page (`python -m report`, run by `05_report.sh`) carries most of the checks that need
the cohort: `truth.auto`, and `truth.chrX` / `truth.chrY` by sex, with the `flagged_chromosomes`
aneuploidy screen as sample flags (3.1); `DJ.cn` and its whole-copy steps, carriers, transmission
and de novo steps (3.2); the sample-by-sample comparison with Hall et al. 2021 and the duplicate-flag
correction of the ratio between the pipelines (3.7); what the coverage PCs remove, reliability after
adjustment, and the cross-validated PC sweep (3.9); and 45S-5S, 45S-mtDNA (Gibbons et al. 2014,
2015) and DJ against the female X, the S-phase test (4). Not scripted yet:

- whether a child's departure from the midparent tracks the state of the culture it was
  sequenced from (`chrEBV.copies`, `chrM.copies`, the leading control PCs): the part of the
  non-transmitted variance that is biology of the cell line rather than measurement error;
- whether the leading control PC moves with DJ and the female X.

## Data sources and attribution

- NYGC 30× CRAMs: Byrska-Bishop et al., *Cell* 185:3426 (2022); AWS Open Data mirror `s3://1000genomes/1000G_2504_high_coverage/`.
- HGSVC and Illumina Platinum-pedigree CRAMs: IGSR data collections `hgsv_sv_discovery` and `illumina_platinum_pedigree` (EBI).
- `pilot/hall2021_MOESM1.txt`: Supplementary Data 1 of Hall, Turner & Queitsch, *Sci Rep* 11:449 (2021),
  doi:10.1038/s41598-020-80049-y, distributed under CC BY 4.0; unmodified.
