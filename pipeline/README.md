# The cohort run

The expanded 1000 Genomes cohort (3,202 samples, 602 trios; NYGC, TruSeq PCR-free, NovaSeq
2×150, GRCh38) is the validation cohort for NGS-DOSE: it is public, it has trios, many samples
have HPRC assemblies, and NGS-PCA has already been run on it (`../meta/ngspca/`).

| path | what |
| --- | --- |
| `../pilot/` | four trios, each sample also as an independent older library; run locally; counts files, the evaluation and plotting scripts, the report and the figure are there |
| `00_setup.sh`, `01_stage_and_dose.sh`, `01a_dispatch_staged.sh`, `01b_dose_sample.sh`, `01_count.sh`, `02_cohort.sh`, `03_compare_modes.sh` (+ `compare_modes.py`), `04_hprc_satellites.sh`, `05_report.sh`, `config.sh` (+ `check_counts.py`) | the full-cohort run, for SLURM or a plain loop |
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
sinks ship yet, so `FETCH_PANELS` stays with the telomeric repeat. Sinks belong to an aligner and a reference, and are learned again for
another pipeline. A class loaded in these scans can be fetched later from the public CRAMs once its
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
