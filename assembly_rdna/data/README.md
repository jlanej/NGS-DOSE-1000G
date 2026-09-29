# assembly_rdna/data — HPRC release-2 targeted rDNA extraction

Built 2026-09-24 with `scripts/fetch_metadata.py`, `scripts/extract_regions.py`,
`scripts/validate_full.py` and the shared measure `scripts/tile_map.py` (conda env `asm`:
samtools/htslib 1.24, minimap2 2.31-r1302, pysam 0.24.1, pandas 3.0.6, numpy 2.5.3, Python 3.12).

## Inputs
* `index.csv` — HPRC r2 index `assemblies_release2_v1.0.index.csv`
  (github human-pangenomics/hprc_intermediate_assembly, main branch, fetched 2026-09-24).
  466 rows = 464 haplotype assemblies of 232 people + 2 reference rows (GRCh38 no-alt
  analysis set, CHM13 v2.0 maskedY rCRS). All 466 rows were processed.
* S3 `s3://human-pangenomics/...` read over `https://s3-us-west-2.amazonaws.com/human-pangenomics/...`.
* Unit sequences: 45S KY962518.1, 5S X12811.1, DJ core CHM13 chr21:2,708,299-3,108,298
  (all from `/Users/Kitty/git/NGS-DOSE/work/ref/`). Tile-mapping targets are written to `targets/`.

## Layout
| path | content |
|---|---|
| `meta/fai/`, `meta/gzi/` | `.fai` / `.gzi` of every assembly (used locally by remote `samtools faidx`) |
| `meta/chrom_assignment/` | `chromAlias.txt`, `gaps.bed`, `t2t_chromosomes.tsv` as published (exact S3 file versions in `meta/manifest.tsv`); `*.derived.*` = our substitutes for rows without them |
| `meta/chains/` | `<name>_vs_GRCh38.chain.gz` (minigraph-cactus) + `<name>.chainsum.tsv` (aligned bp / strand / GRCh38 span per contig x GRCh38 chromosome) |
| `meta/manifest.tsv` | per assembly: chosen S3 key and local path for every file type, notes |
| `meta/s3_listing.json` | S3 listing used to resolve file names |
| `censat/` | CenSat BEDs not already in `../../hprc_censat/` (60 fetched, gzipped) and 3 documented substitutes (`*.cenSat.substitute.bed.gz`). The 400 files in `hprc_censat/` are used in place (read-only). |
| `regions/<name>.bed` | extracted regions, BED 0-based half-open, col 4 = reason(s) `rdna_censat` / `acro_parm` / `unplaced_rdna` / `5S_locus`, `;`-joined where merged regions overlap |
| `regions/<name>.acro.tsv` | every acrocentric-assigned contig: level (chromosome / random), chain strand, short-arm end, orientation source, taken or skipped, segment |
| `seq/<name>.fa.gz` (+ `.fai`, `.gzi`) | bgzipped extracted sequence; record names `contig:start-end` = samtools region strings, **1-based inclusive** (`name = contig:{bed_start+1}-{bed_end}`); verified: record set and lengths equal the BED |
| `seq/<name>.done.json` | per-assembly extraction record (also collected in `../tables/extraction.tsv`) |
| `validation/` | tile hits of the 8 whole assemblies (`*.full.*`) and of the extracted sequence (`*.extracted.*`) of those 8 plus GRCh38, CHM13 and HG06807 pat/mat, each with `summary.json` |
| `full/` | whole-assembly FASTAs for validation (removed after validation) |

## Extraction parameters (extract_regions.py)
* **rdna_censat**: CenSat intervals whose label contains `rDNA` (pure `GAP` excluded), merged within 1 kb, padded 200 kb each side, clipped to the contig.
* **acro_parm**: all contigs whose chromAlias name is `chr13/14/15/21/22` or `chrN_<acc>_random`
  (in HPRC r2 only T2T chromosomes get a chromosome-level name; other chromosome pieces are
  `_random`). Short-arm end = GRCh38 chain strand (majority aligned bp to chrN), else side of the
  CenSat rDNA, else contig start (`orientation_source`). Segment from that end of length
  max(12 Mb, D + 1 Mb), D = far edge of the farthest CenSat rDNA interval or gap within 25 Mb of
  that end, clipped. `_random` contigs are skipped only when they carry no CenSat rDNA and their
  chain alignment to chrN starts beyond the GRCh38 centromere end + 2 Mb (q-arm-only pieces).
* **unplaced_rdna**: whole contig for chrUn / unaliased contigs with CenSat rDNA, **plus** every other
  unassigned contig on which the shared tile measure finds a passing 45S or 5S tile or a DJ run
  >= 50 kb (scan hits in `regions/<name>.unplaced_scan.*`). This screen was added after validation
  pass 1 showed whole-rDNA chrUn contigs with no CenSat label at all; it added 505 contigs in 282
  assemblies (427 via 45S tiles, 54.1 Mb of 45S-tile bp; 77 via DJ runs; 1 via 5S). Without CenSat
  (GRCh38, HG06807), all unassigned contigs are screened (45S tiles) if they total <= 100 Mb; above that
  the screen is skipped and noted. Here GRCh38 had 127 unassigned contigs (4.49 Mb, 1 with 45S),
  HG06807_mat had 1 (16.6 kb, no hit) and HG06807_pat had none.
* **5S_locus**: GRCh38 flanks chr1:228,400,000-228,600,000 and 228,650,000-228,850,000 lifted through the
  GRCh38 chain (chain with most aligned bp in each flank); region = span between the two lifted
  flanks + 300 kb each side. If the flanks land on different contigs, each flank's inner edge
  ±300 kb is taken (`chain_flanks_split_contigs`). Without a chain (HG00272 hap1/hap2, HG06807
  pat/mat, CHM13) all chr1-assigned contigs are tile-scanned and the largest 5S tile cluster
  ±300 kb is taken (`tile_scan_chr1`). GRCh38: identity (228.60-228.65 Mb).
* Overlapping regions are merged per contig. Remote access: `samtools faidx --fai-idx --gzi-idx <https-url> -r regions.txt`,
  8 assemblies concurrently, up to 4 attempts each.

## Run log (2026-09-24)
* fetch_metadata.py: 466/466 rows, 1.51 GB metadata (1.1 GB of it chains), 248 s.
* extract_regions.py pass 1: 466/466 (HG06807_mat_v1 failed once on an empty-hit bug, fixed and rerun),
  ~26 min. Pass 2 (`--force`, with the unplaced tile screen): 466/466 verified, 281 assemblies whose
  region set changed were re-fetched, 1238 s. Estimated bgzf bytes needed for the final extraction
  of all 466: 9.56 GB (sum of `est_bytes_fetched` in `../tables/extraction.tsv`; pass 1 alone 9.36 GB);
  40.0 Gb of sequence, 7.9 GB (7.4 GiB) on disk in `seq/`.
* validate_full.py: 8 whole assemblies (7.02 GB downloaded in 27 min, 4 concurrent), each tile-mapped in
  ~100 s (9 threads). Whole FASTAs deleted afterwards (moved to Trash).

## Shared rDNA measure (tile_map.py)
5-kb tiles on a grid anchored at contig position 0 (extracted records are placed back on contig
coordinates, so interior tiles are identical between whole and extracted sequence); tiles < 500 bp
or > 50% N skipped. minimap2 `-x asm20 -c --secondary=no` against (A) KY962518.1 doubled and (B)
X12811.1 doubled + DJ core, run separately. A tile is class sequence (45S / 5S / DJ) if alignments
at >= 90% identity (PAF matches / block length) cover >= 50% of the tile.
Secondary run summary: `{cls}_bp_runs_ge50kb` = passing-tile bp in runs >= 50 kb (gaps <= 1 tile).
**Use `DJ_bp_runs_ge50kb`, not `DJ_bp`, for DJ**: per tile the DJ class hits 44-47 Mb per genome
(interspersed repeats; only 6-8% inside the extracted regions), whereas DJ runs >= 50 kb are 1.4-2.4 Mb
per haplotype, 96-100% inside the extracted regions (the only runs outside are a 50-80 kb DJ-like
segment at chr11 ~93.0-93.8 Mb in 4 of 8 haplotypes).

## Substitutions and missing files (see `../tables/assemblies.tsv` `missing_files`, `notes`)
* HG06807 pat/mat (verkko 2, HPRC_PLUS GenBank build): no CenSat, chromAlias, gaps, T2T table or chains.
  Chromosome names from NCBI nuccore titles; short-arm end = contig start (default, no evidence).
* hg002v1.1.pat / .mat_MT (Q100): CenSat from T2T consortium `hg002v1.1.cenSatv2.0.{PAT,MAT}.bed`
  (contigs renamed to PanSN); chromosome names from PanSN; no gaps/T2T tables (T2T by construction).
* CHM13: CenSat v2.1 from local copy (renamed); no GRCh38 chain.
* GRCh38: no CenSat; acrocentric N-runs computed from local GRCh38 FASTA as gaps.
* HG00272 hap1/hap2: no GRCh38 chain in S3 (all other files present).
