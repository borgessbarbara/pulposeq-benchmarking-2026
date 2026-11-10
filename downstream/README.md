# downstream/ - Fig. 1 (LongBench) and Fig. 2 (Dong) on the oarfish runs

These are the downstream analyses of `before-oarfish-integration/`, adapted to the pulposeq
runs in `config-runs/` (`run_oarfish: true`, `oarfish_num_bootstraps: 50`,
`min_total_counts: 5`, outdir `results_oarfish/`, pulposeq `transcriptome/` output layout).
Each analysis keeps its old script number and logic; what changed is listed below.

## What changed compared with before-oarfish-integration/

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
- The sequin scripts also run gene level on **undivided** summed oarfish counts
  (`count_mode = "undivided"`), as a sensitivity check of the gene-level division.
  Figures use `count_mode = "divided"` only.
- SIRV/ERCC accuracy uses undivided oarfish counts (nothing is tested) and reports each
  isoform's overdispersion.
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
| longbench | `unique_intron_chains/` | `03_ujc_report.R`, `04_ujc_length.R`, `06_supp_tables_ujc.R` | 03 before 04 and 06 |
| longbench | `read-support/` | `05_read_support.R` | - |
| longbench | `figure-panels/` | `07_fig1_panels.R` (Fig. 1) | 01-05 |
| longbench | `nanocomp/` | NanoComp merges (shell) | - |
| dong | `diff-exp/` | `01_de_edger.R` (gene + transcript DE) | - |
| dong | `sequins/` | `02_sequin_edger.R` (sequin DE) | - |
| dong | `correlation-cis-1Mb/` | `03_cis_circos.R` | 01 |
| dong | `isoformic/` | `04_isoformic.R` | 01, 03 |
| dong | `read-support/` | `05_read_support.R` | - |
| dong | `figure-panels/` | `06_fig2_panels.R` (Fig. 2), `07_circos_coverage_biotypes.R` | 06: 01, 02, 05 |

`read_oarfish()` caches what it reads (`spike-ins/results/oarfish_cache/`,
`diff-exp/results/oarfish_cache/`) and re-reads when the oarfish files change; delete the cache
folder to force a fresh read.
