# downstream/ — final scripts for Fig. 1 (LongBench) and Fig. 2 (Dong)

Consolidated from `organize_scripts/` (copies of the scripts in `analysis_cross_protocol/` and `dong/de/scripts/`).
Nothing outside `downstream/` is modified: the scripts **read** the pulposeq run folders
(`ont_cdna/`, `ont_drs/`, `pb/`, `dong/`) and `refs/`, and **write** only to `downstream/*/results/` and `downstream/*/figures/`.

Every path is set once, in `00_paths.R` of each folder. To move the project, edit `base` there (and the one
`source(".../00_paths.R")` line at the top of each script).

## Folder layout

```
downstream/
├── README.md
├── run_all.R                       runs every script below in order (logs in logs/)
├── longbench/                      Fig. 1
│   ├── 00_paths.R
│   ├── 01_sequin_edger.R           sequin DE ground truth (edgeR) + sequin identification
│   ├── 02_sirv_ercc.R              SIRV / ERCC benchmark (length bias, E0, ERCC dose response)
│   ├── 03_ujc_report.R             novel transcripts grouped by unique intron chain
│   ├── 04_ujc_length.R             UIC sharing by transcript length
│   ├── 05_read_support.R           read-level support per biotype (3 LongBench runs)
│   ├── 06_supp_tables_ujc.R        UIC supplementary tables (draft ST2–ST7)
│   ├── 07_fig1_panels.R            Fig. 1 b–i, one file per panel
│   ├── helpers/  plot_theme.R, figure_helpers.R
│   ├── data/     samplesheet_ont.csv, samplesheet_ont_drs.csv, samplesheet_pb.csv,
│   │             SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx   (copied, unchanged)
│   ├── results/  (created by the scripts)
│   └── figures/fig1_panels/        (created by 07)
└── dong/                           Fig. 2
    ├── 00_paths.R
    ├── 01_de_edger.R               DE, gene + isoform level (complete, every plot kept);
    │                               also the sample Spearman matrix of Fig. 2d
    ├── 02_sequin_edger.R           sequin DE ground truth (edgeR)
    ├── 03_cis_circos.R             cis co-expression + circos (complete, every plot kept)
    ├── 04_isoformic.R              isoformic: events, cis pairs, journal panels (complete)
    ├── 05_read_support.R           read-level support per biotype (Dong run)
    ├── 06_fig2_panels.R            Fig. 2 a–f, one file per panel
    ├── 07_circos_coverage_biotypes.R  genome-wide circos: RNA coverage, annotated + novel
    │                               transcripts per strand, legend in the centre
    ├── helpers/  plot_theme.R, id_helpers.R, figure_helpers.R
    ├── results/  (created by the scripts)
    └── figures/  fig2_panels/ (06), isoformic_panels/ (04), circos_coverage/ (07)
```

## Run order

**Everything at once:** `Rscript run_all.R` (or Source it in RStudio). It runs every script below in order,
each in its own fresh R process, prints OK / FAILED / SKIPPED per script and writes the full output to
`logs/<date_time>/` (one `.log` per script + `summary.csv`). Settings at its top: `SKIP` (reuse old results),
`ONLY` (run just some scripts), `STOP_ON_ERROR` (default FALSE: a failure only skips the scripts that need it).

Run each folder top to bottom (`Rscript 01_….R`, or `source()` in RStudio). Within a folder, a script only
needs the ones with a lower number that it lists under "Needs".

| Folder | Script | Needs | Main outputs |
|---|---|---|---|
| longbench | 01_sequin_edger.R | – | `results/spikeins/sequins/` |
| longbench | 02_sirv_ercc.R | – | `results/spikeins/sirv_ercc/` |
| longbench | 03_ujc_report.R | – | `results/ujc/` |
| longbench | 04_ujc_length.R | 03 | `results/ujc/length/` |
| longbench | 05_read_support.R | – | `results/read_support/` |
| longbench | 06_supp_tables_ujc.R | 03 | `results/supp_tables/` |
| longbench | 07_fig1_panels.R | 01–05 | `figures/fig1_panels/`, `results/figure_tables/` |
| dong | 01_de_edger.R | – | `results/de_edger/` |
| dong | 02_sequin_edger.R | – | `results/sequins_edger/` |
| dong | 03_cis_circos.R | 01 | `results/cis_circos/` (incl. `paper_cis_circos.*`) |
| dong | 04_isoformic.R | 01, 03 | `results/isoformic/`, `figures/isoformic_panels/` |
| dong | 05_read_support.R | – | `results/read_support/` |
| dong | 06_fig2_panels.R | 01, 02, 05 | `figures/fig2_panels/`, `results/figure_tables/` |
| dong | 07_circos_coverage_biotypes.R | – (reads `dong/ref_transcripts`, `dong/novel_transcripts`, `dong/coverage/*.bw`) | `figures/circos_coverage/`, `results/circos_coverage/` |

`03_cis_circos.R` must run before `04_isoformic.R`: isoformic reads the cis-pair table
`results/cis_circos/circos_cis_de_pairs.csv` (the old `de/scripts/cis_top_pairs.csv` was a copy of it).

## Panel map

### Fig. 1 (LongBench)

| Panel | Content | Drawn by | File |
|---|---|---|---|
| a | pulposeq workflow schematic | vector editor | – |
| b | SIRV length bias + identification | 07_fig1_panels.R (`pB`) | `fig1b_length_bias` |
| c | ERCC dose response | 07 (`pD`) | `fig1c_ercc` |
| d | sequin ground truth, isoform level | 07 (`pC`) | `fig1d_sequin_isoform` |
| e | read-level support per protocol | 07 (`pE`) | `fig1e_read_support` |
| f | novel transcript length | 07 (`pF`) | `fig1f_novel_length` |
| g | UIC sharing | 07 (`pG`) | `fig1g_uic_sharing` |
| h | UIC sharing by length | 07 (`pH`) | `fig1h_uic_sharing_by_length` |
| i | transcript ends within an identical UIC | 07 (`pI`) | `fig1i_transcript_ends` |

The old Fig. 1e (novel transcripts by biotype × class code) is a table now:
`results/figure_tables/supp_novel_biotype_classcode_longbench.csv`.

### Fig. 2 (Dong)

| Panel | Content | Drawn by | File |
|---|---|---|---|
| a | read-level support across biotypes | 06_fig2_panels.R (`pA`) | `fig2a_read_support` |
| b | novel transcripts by biotype × class code | 06 (`pB`) | `fig2b_novel_biotype` |
| c | sequin ground truth, isoform level | 06 (`pC`) | `fig2c_sequin_isoform` |
| d | sample correlation (Spearman on edgeR logCPM; matrix from 01) | 06 (`pD`) | `fig2d_sample_correlation` |
| e | volcano, gene level | 06 (`pE`) | `fig2e_volcano` |
| f | DE isoforms of the marker genes (vertical heatmap) | 06 (`pF`) | `fig2f_isoform_heatmap` |
| g | EGFR isoform profile | 04_isoformic.R, section 7e | `isoformic_panels/fig2_panels/fig2g_profile_EGFR` |
| h | EGFR locus | 04_isoformic.R | `isoformic_panels/locus_*/EGFR.*` (right part) — see "Not yet as printed" |
| i, j, k | UBE2I / UBE2I-AS1, ELOC, NDUFA6 profiles | 04_isoformic.R, section 7e | `isoformic_panels/fig2_panels/fig2{i,j,k}_profile_*` |
| l | cis-acting pairs, circos | 03_cis_circos.R | `results/cis_circos/paper_cis_circos.*` — see "Not yet as printed" |
| m | cis-pair heatmap | 04_isoformic.R, section 7e | `isoformic_panels/fig2_panels/fig2m_cis_pairs_heatmap` |

**Lineage plasticity markers** (not in the printed Fig. 2) = `isoformic_panels/fig2_panels/marker_profile_{POU2F3,HMGA2}`
(04, section 7e; same style and size as Fig. 2k, every tested isoform shown), plus a page each in
`results/isoformic/plots/02_profiles.pdf` and `03_log2fc.pdf`. ASCL1 is not expressed (not tested) and is skipped.
The Extended Data grid is unchanged (`cs$panel_exclude`).

**Genome-wide circos** (coverage + biotypes) = `figures/circos_coverage/circos_coverage_biotypes.{pdf,png,svg}`
(07; 8 × 7.49 in; ring labels in the 12 o'clock gap, legend in the centre, no poster version).

**Extended Data Fig. 1** = `figures/isoformic_panels/ed_profiles_events_cis.*` (every event + cis-pair profile).

Each panel is written as `.pdf`, `.svg` (editable text; needs the font "Source Sans 3" installed) and a 600 dpi `.png`.

## What changed compared with the old scripts

Analyses are unchanged (same inputs, filters, tests, scoring). Changes:

1. **Paths**: from `00_paths.R`; outputs go to `downstream/` only. The LongBench novel-transcript metadata is read
   from the run folders (`<run>/novel_transcripts/…validated_metadata.csv`), which are byte-identical (md5) to the
   copies in `analysis_cross_protocol/data/`.
2. **Only the plots used in the figures** are drawn. Removed from the analysis scripts (tables kept):
   Anaquin linearity / fold-change / ROC / LOD plots (LongBench and Dong sequins, SIRV / ERCC), the 04b edgeR
   overview and the DESeq2-vs-edgeR comparison, the 06 report figures (UpSet, corroboration, disagreement rates,
   ECDF, agreement vs ends), the 07 plots, `make_fig_dong2.R` assembly summary (`pE`) and one-page figure.
   **Kept complete, with every plot**: `01_de_edger.R`, `03_cis_circos.R`, `04_isoformic.R`.
3. **04_isoformic.R** = `isoformic_dong2.R` (engine) + `isoformic_consolidated2.R` in one file, with the settings
   consolidated2 used applied directly (`modifyList`, no `isoformic_overrides`). The engine sections that this run had
   switched off (biotype summary, default plot sets, fgsea enrichment, lncRNA–known-gene pairs) and their helper
   functions are removed; everything consolidated2 wrote is still written.
4. **No DESeq2**: Fig. 2d uses the Spearman matrix `01_de_edger.R` already writes
   (`results/de_edger/transcript/sample_spearman.csv`: logCPM, prior.count 2, the isoforms it tested); the panel
   subtitle says "logCPM" instead of "VST".
5. **03_ujc_report.R / 06_supp_tables_ujc.R**: the command-line options (`optparse`) are a config list now.
6. **Panels as printed.** The final figures (`manuscript/Figure1_2.pdf`, `manuscript/Figure2.pdf`) are A3 pages
   (297 × 420 mm). Each panel is now saved at the size it has on that page, with the text sizes, legend positions,
   key sizes, line spacing and panel-border widths measured from those PDFs, and with the text and orientation as
   printed. In short:
   - **Fig. 1**: titles 10 pt, subtitles 9 pt, axis / legend text 8 pt, strips 9.5 pt (b, d) or 9 pt;
     b, d corner text 7 pt (ONT DRS text bottom-left, as printed); c, f labels 8 pt; e bar labels 5.3 pt;
     g bar and row labels 7 pt; i agreement note 8 pt under the legend.
   - **Fig. 2**: titles and subtitles 10 pt; a, c, e axis titles 10 pt; d sample names 10 pt; e totals 10 pt and
     gene labels 7.5 pt; f gene names 8.5 pt, isoforms as rows, columns H1975_2, _1, _3 | HCC827_2, _1, _3 | DE;
     g, i–k single profiles (titles / subtitles 10 pt, isoform labels 6 pt, the printed subtitles, e.g.
     "Driver locus"); m heatmap at 1.34 × the Extended Data version (names 8 pt, ρ / location 7 pt,
     colour bar "z-score (logCPM)", no biotype key).
   - Geometry (point sizes, line widths, bar widths) is what the original scripts drew.
   - One deliberate difference: Fig. 2d subtitle says "logCPM", not "VST" (no DESeq2).