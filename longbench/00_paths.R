## ===========================================================================
## 00_paths.R - every path used by the LongBench (Fig. 1) scripts
## Sourced by each script in downstream/longbench/; edit here only.
##
## Inputs  (read only): the pulposeq run folders <base>/{ont_cdna,ont_drs,pb}/
##                      and downstream/longbench/data/ (sample sheets + Lexogen SIRV design)
## Outputs (written)  : downstream/longbench/results/ and downstream/longbench/figures/
## Nothing outside downstream/ is written.
## ===========================================================================

base     <- "C:/Users/mbarb/OneDrive/paper"              # pulposeq run folders live here
lb_dir   <- file.path(base, "downstream/longbench")

helper_dir <- file.path(lb_dir, "helpers")
data_dir   <- file.path(lb_dir, "data")                  # samplesheet_*.csv, SIRV design .xlsx
res_dir    <- file.path(lb_dir, "results")
fig_dir    <- file.path(lb_dir, "figures")

## pulposeq run folders (one per LongBench protocol)
platforms <- c("ont_cdna", "ont_drs", "pb")
run_dir   <- setNames(file.path(base, platforms), platforms)

## results of each step
spk_dir   <- file.path(res_dir, "spikeins")                  # 01, 02
seq_dir   <- file.path(spk_dir, "sequins")                   # 01
sirv_dir  <- file.path(spk_dir, "sirv_ercc")                 # 02
ujc_dir   <- file.path(res_dir, "ujc")                       # 03, 04 (length/)
rs_dir    <- file.path(res_dir, "read_support")              # 05
supp_dir  <- file.path(res_dir, "supp_tables")               # 06
tab_dir   <- file.path(res_dir, "figure_tables")             # 07: numbers behind the panels
panel_dir <- file.path(fig_dir, "fig1_panels")               # 07: one file per panel
