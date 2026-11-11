# downstream/ - Fig. 1 (LongBench) and Fig. 2 (Dong) on the oarfish runs

These are the downstream analyses of the original (Bambu-only) runs, adapted to the pulposeq
runs in `config-runs/` (`run_oarfish: true`, `oarfish_num_bootstraps: 20`,
`min_total_counts: 5`, pulposeq `transcriptome/` output layout). The original scripts were in
`before-oarfish-integration/` (removed from this branch; see the git history).

## Where the scripts look

| What | Path (set in `longbench/00_paths.R`, `dong/00_paths.R`) |
|---|---|
| project folder (`project`) | `C:/Users/mbarb/OneDrive/pulposeq-oarfish/` |
| LongBench run folders (pulposeq `--outdir`) | `<project>/ont_cdna`, `<project>/ont_drs`, `<project>/pb` |
| Dong run folder | `<project>/dong` |
| this repository (clone) | `<project>/pulposeq-benchmarking-2026/` |
| reference files | `<repo>/references/` (e.g. the Lexogen SIRV design .xlsx) |
| sample sheets | `<repo>/downstream/<dataset>/data/` |
| GENCODE v50 GTF (Dong `04_isoformic.R` only; too large for the repo) | `ref_gtf` in `dong/00_paths.R` |
| results and figures | `<repo>/downstream/<dataset>/<analysis>/{results,figures}/` (git-ignored) |

To move the project: edit `project` in both `00_paths.R`, `BASE` in `run_all.R` and the
`source(".../00_paths.R")` line at the top of each script.
Each analysis keeps its old script number and logic; what changed is listed below.

## What changed compared with the original (Bambu-only) scripts

**Counts that are modelled come from oarfish, with the bootstrap replicates used as edgeR intends.**
`helpers/oarfish_counts.R` (`read_oarfish()`) is the single reader:

| Level | edgeR function | What edgeR receives |
|---|---|---|
| Transcript | `catchOarfish()` | counts / read-to-transcript ambiguity (RTA) overdispersion |
| Gene | `catchOarfishGene()` | counts / gene-level overdispersion (bootstraps summed per gene) |

- The overdispersion is estimated from each sample's `.infreps.pq` bootstraps, over every
  model oarfish quantified, then restricted to pulposeq's curated set. It is >= 1
  (1 = no ambiguity); every DE table reports it in an `overdispersion` column.
- Divided counts are not rounded: they "behave much like negative binomial counts" and go into
  the unchanged edgeR quasi-likelihood test (glmQLFit robust -> glmTreat |log2FC| > 0.5, FDR 0.05)
  (edgeR reference manual; Baldoni et al. 2024a NAR 52:e13; 2024b NARGAB lqae151; 2025 NAR gkaf1305).
- Expression filters (>= 10 reads in >= n samples) are applied to the **undivided read counts**;
  edgeR then tests the divided counts.
- The sequin scripts run both levels on **divided** (primary, figures) and **undivided** counts
  (`count_mode`), so the effect of the RTA division is visible. Their DE filter follows
  `01_de_edger.R` (`filter_counts = "undivided"`); change both together.
- SIRV/ERCC accuracy uses undivided oarfish counts (nothing is tested) and reports each
  isoform's overdispersion, with the Spearman correlation between abundance error and
  overdispersion (assignment errors come with high overdispersion; capture errors do not).
- **Quantifier comparison.** Every spike-in script also runs on Bambu's counts of the same runs
  (`Rscript <script> bambu`, or `@bambu` in `run_all.R`; results in `<folder>_bambu/`), and
  `longbench/spike-ins/08_quantifier_comparison.R` sets oarfish and Bambu side by side
  (same alignments, annotation and libraries) for the Supplement.
- **Spike-in definitions** (identical for both quantifiers and both datasets): detected = >= 3
  reads in >= 2 libraries (ERCC: >= 1 read in >= 2 libraries), counts rounded per library;
  SIRV accuracy = median across libraries of log2(observed / expected), expected = 1/84;
  sequin outliers, capture efficiency vs length and robust fold-change agreement are written
  next to the DE tables.
- Read support (full-length / unique reads) is Bambu's, as before. The curated set is now chosen
  on oarfish counts, so transcripts kept on oarfish counts alone have no Bambu reads and fall in
  "neither".

**Paths.** Every pulposeq file is read from the `transcriptome/` layout; `00_paths.R` of each
dataset holds all of them (`run_dir`, `pp`). Results go to `<analysis>/results/` and
figures to `<analysis>/figures/`.

**Sample sheets.** `longbench/data/` holds the run sample sheets from `config-runs/`
(`*_ONTc`, `*_ONTd`, `*_PB`). `dong/data/samplesheet.csv` keeps the group labels the scripts use
(H1975 = TKI-resistant, HCC827 = TKI-sensitive) plus a `cell_line` column.

## Requirements

- **edgeR >= 4.99.12** (Bioconductor devel, released as edgeR 5.0.0), for
  `catchOarfish()` and `catchOarfishGene()`:

  ```r
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(version = "devel")
  BiocManager::install("edgeR")
  ```
- `readr`, `nanoparquet`, `jsonlite` (read by `catchOarfish`), plus the packages each script
  loads (tidyverse, Anaquin, ggrepel, pheatmap, patchwork, circlize, showtext, svglite,
  isoformic, plotgardener, rtracklayer, ...), as before.
- The run folders must contain `transcriptome/oarfish/` with the per-sample `.quant`,
  `.meta_info.json` and `.infreps.pq` files (published by pulposeq when `run_oarfish` is true).

## Layout and run order

Edit `base` (and `run_dir`) in `longbench/00_paths.R` and `dong/00_paths.R`, and `BASE` in
`run_all.R`. Then `Rscript run_all.R` runs everything in order, each script in a fresh R process,
with logs in `logs/<date_time>/`.

| Dataset | Analysis folder | Script | Needs |
|---|---|---|---|
| longbench | `spike-ins/` | `01_sequin_edger.R` (sequin DE + identification) | - |
| longbench | `spike-ins/` | `02_sirv_ercc.R` (SIRV E0, length bias, ERCC) | - |
| longbench | `spike-ins/` | `01` and `02` again with `bambu` (quantifier comparison) | - |
| longbench | `unique_intron_chains/` | `03_ujc_report.R`, `04_ujc_length.R`, `06_supp_tables_ujc.R` | 03 before 04 and 06 |
| longbench | `read-support/` | `05_read_support.R` | - |
| longbench | `figure-panels/` | `07_fig1_panels.R` (Fig. 1) | 01-05 |
| longbench | `nanocomp/` | NanoComp merges (shell) | - |
| dong | `diff-exp/` | `01_de_edger.R` (gene + transcript DE) | - |
| dong | `sequins/` | `02_sequin_edger.R` (sequin DE; again with `bambu`) | - |
| dong (runs last) | `../longbench/spike-ins/` | `08_quantifier_comparison.R` (oarfish vs Bambu, both datasets) | LongBench 01-02 and Dong 02, each with both quantifiers |
| dong | `correlation-cis-1Mb/` | `03_cis_circos.R` | 01 |
| dong | `isoformic/` | `04_isoformic.R` | 01, 03 |
| dong | `read-support/` | `05_read_support.R` | - |
| dong | `figure-panels/` | `06_fig2_panels.R` (Fig. 2), `07_circos_coverage_biotypes.R` | 06: 01, 02, 05 |

`read_oarfish()` caches what it reads (`spike-ins/results/oarfish_cache/`,
`diff-exp/results/oarfish_cache/`) and re-reads when the oarfish files change; delete the cache
folder to force a fresh read.
