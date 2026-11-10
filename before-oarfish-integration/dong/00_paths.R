## ===========================================================================
## 00_paths.R - every path used by the Dong (Fig. 2) scripts
## Sourced by each script in downstream/dong/; edit here only.
##
## Inputs  (read only): the pulposeq run folder <base>/dong/ (bambu_validated/,
##                      ref_transcripts/, novel_transcripts/, transcriptome_report/,
##                      coverage/, samplesheet.csv) and the reference GTF in <base>/refs/
## Outputs (written)  : downstream/dong/results/ and downstream/dong/figures/
## Nothing outside downstream/ is written.
## ===========================================================================

base       <- "C:/Users/mbarb/OneDrive/paper"
data_dir   <- file.path(base, "dong")                    # pulposeq --outdir of the Dong run
sheet      <- file.path(data_dir, "samplesheet.csv")     # columns: sample, group, fastq
ref_gtf    <- file.path(base, "refs/gencode.v50.primary_assembly.annotation.gtf")  # pulposeq --annotation

dg_dir     <- file.path(base, "downstream/dong")
helper_dir <- file.path(dg_dir, "helpers")
res_dir    <- file.path(dg_dir, "results")
fig_dir    <- file.path(dg_dir, "figures")

## results of each step
de_dir     <- file.path(res_dir, "de_edger")             # 01: <level>/de_*.csv, logcpm.rds, sample_spearman.csv, all_results.csv
seq_dir    <- file.path(res_dir, "sequins_edger")        # 02: sequin_de_results.csv, summary.csv
cis_dir    <- file.path(res_dir, "cis_circos")           # 03: circos_cis_de_pairs.csv, paper_cis_circos.*
iso_dir    <- file.path(res_dir, "isoformic")            # 04: tables/, plots/
iso_paper  <- file.path(fig_dir, "isoformic_panels")     # 04: journal panels (Fig. 2g-k, m; ED Fig. 1)
rs_dir     <- file.path(res_dir, "read_support")         # 05
tab_dir    <- file.path(res_dir, "figure_tables")        # 06: numbers behind the panels
panel_dir  <- file.path(fig_dir, "fig2_panels")          # 06: one file per panel
