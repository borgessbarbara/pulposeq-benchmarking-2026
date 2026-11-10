#!/usr/bin/env Rscript
# =============================================================================
# pulposeq | Computational resource usage across runs
# -----------------------------------------------------------------------------
# Input : pulposeq_supplementary_tables.xlsx
#           ST1-ST4 = per-process Nextflow trace metrics (one sheet per run)
#           ST5     = run-level totals
# Output: <outdir>/resources_panel.{pdf,png,svg}        assembled panel
#         <outdir>/resources_{totals,steps,profile}.{pdf,png,svg}   each part alone
#         <outdir>/source_data_resources_*.csv          source data
#
# Usage : Rscript plot_resources_panel.R <supp_tables.xlsx> [outdir] [input_stats.tsv]
#    or : source("plot_resources_panel.R")   # xlsx in the same folder as the script
#
#   input_stats.tsv is optional: tab-separated, columns `run` and `gigabases`.
#   `run` must match the ST5 column names (e.g. "LongBench ONT cDNA").
#   When you supply it, the top row gets a "CPU time per Gb input" facet.
#
# Parts (choose with PULPO_PARTS; default c("totals", "steps", "profile", "io"));
# in the assembled panel they are tagged A, B, C, ... in this order
#   totals : one bar per run - elapsed time, peak RAM, peak CPU cores,
#            data written.
#   steps  : wall time, CPU time and peak RAM of the main steps, one bar per
#            run; every step x metric cell on its own scale (needs ggh4x).
#   profile: every Nextflow process, step real time vs peak RAM (log-log),
#            point area = CPU cores of the busiest task, colour = Fig. 2a stage.
#   io     : main steps, data read vs data written (equal log axes), with the
#            written = read diagonal.
#
# Dependencies: readxl, dplyr, tidyr, ggplot2 (>= 3.5), patchwork, scales, ggh4x
# Optional    : ggrepel (labels), ragg (PNG), svglite (SVG), systemfonts (font lookup), rstudioapi
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
if (!requireNamespace("ggh4x", quietly = TRUE)) {
  stop("Package 'ggh4x' is required (independent x axes per facet): install.packages('ggh4x')")
}

# ---- 0. Arguments --------------------------------------------------------------
# Inputs are resolved in this order:
#   1. command-line arguments (Rscript plot_resources_panel.R <xlsx> [outdir] [tsv])
#   2. variables set before source(): PULPO_XLSX, PULPO_OUTDIR, PULPO_INPUT_STATS
#   3. defaults next to this script: ./pulposeq_supplementary_tables.xlsx and
#      ./figures/resources/
script_dir <- function() {
  # Rscript: the path comes in as --file=...
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  # source(): the calling frame stores the file path as `ofile`
  for (f in rev(sys.frames())) {
    if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  }
  # RStudio "Run" without source(): ask the editor for the open file
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
pick_arg <- function(i, var, default) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= i) return(args[i])
  if (exists(var, envir = globalenv())) return(get(var, envir = globalenv()))
  default
}

SCRIPT_DIR  <- script_dir()
xlsx        <- pick_arg(1, "PULPO_XLSX",
                        file.path(SCRIPT_DIR, "pulposeq_supplementary_tables.xlsx"))
outdir      <- pick_arg(2, "PULPO_OUTDIR", file.path(SCRIPT_DIR, "figures", "resources"))
input_stats <- pick_arg(3, "PULPO_INPUT_STATS", NA_character_)
message("Input : ", xlsx, "\nOutput: ", outdir)

if (!file.exists(xlsx)) {
  stop("Supplementary table not found: ", xlsx,
       "\nPut the xlsx next to the script, or set PULPO_XLSX <- \"<path>\" before source().")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Configuration ----------------------------------------------------------
# Runs: sheet with per-process metrics -> display label -> ST5 column name.
RUNS <- tribble(
  ~sheet, ~run,                    ~st5_col,
  "ST1",  "LongBench · ONT cDNA",  "LongBench ONT cDNA",
  "ST2",  "LongBench · ONT DRS",   "LongBench ONT DRS",
  "ST3",  "LongBench · PacBio",    "LongBench PacBio",
  "ST4",  "Dong · ONT cDNA",       "Dong ONT cDNA"
)

# Protocol colours as in Fig. 2 (cDNA navy, DRS green, PacBio magenta).
# Dong is the same protocol as LongBench cDNA, so it gets a lighter blue.
RUN_COLS <- setNames(c("#2B4B9B", "#1B8A4B", "#D6297F", "#6FA8DC"), RUNS$run)

# Steps shown in part (2), top to bottom in pipeline order. Every Nextflow
# process not listed here is pooled into "All other steps".
STEP_GROUPS <- tribble(
  ~step,                    ~group,
  "NANOCOMP_RAW",           "NanoComp QC (Dong only)",
  "NANOCOMP_REST",          "NanoComp QC (Dong only)",
  "NANOCOMP_MINIMAP2",      "NanoComp QC (Dong only)",
  "RESTRANDER",             "Restrander (ONT cDNA only)",
  "MINIMAP2_ALIGN",         "minimap2 alignment",
  "BAM_COVERAGE",           "Mapping coverage",
  "BAMBU",                  "Bambu assembly + quantification",
  "VALIDATE_BAMBU_OUTPUTS", "Validate Bambu outputs"
)
STEP_ORDER <- unique(STEP_GROUPS$group)
OTHER      <- "All other steps"

# Stages of the Fig. 2a subway map (colour of the points in the profile part).
# Processes missing here stop the script, so new steps are never misfiled.
STAGE_MAP <- list(
  "1 · QC and mapping" = c("NANOCOMP_RAW", "NANOCOMP_REST", "NANOCOMP_MINIMAP2",
                           "RESTRANDER", "MINIMAP2_ALIGN", "BAM_COVERAGE"),
  "2 · Assembly, quantification, coding potential" =
                         c("BAMBU", "FILTER_BAMBU_COUNTS", "GFFREAD",
                           "GFFCOMPARE", "CPC2"),
  "3 · Annotation refinement, lncRNA discovery" =
                         c("NOVEL_TRANSCRIPTS", "VALIDATE_BAMBU_COUNTS",
                           "VALIDATE_BAMBU_OUTPUTS", "VALIDATE_BAMBU_GTF",
                           "REF_TRANSCRIPTS", "ANNOTATION"),
  "4 · Visualization and reporting" =
                         c("GENOMIC_CONTEXT", "RENDER_REPORT", "MULTIQC")
)
STAGE_TBL  <- tibble(stage = rep(names(STAGE_MAP), lengths(STAGE_MAP)),
                     step  = unlist(STAGE_MAP, use.names = FALSE))
# Ordinal ramp in the pulposeq magenta family (dark -> light = stage 1 -> 4)
STAGE_COLS <- setNames(c("#6E0B47", "#B8327F", "#E08AB9", "#F5D3E6"), names(STAGE_MAP))

# Steps labelled in the profile part
STEP_LABELS <- c(
  NANOCOMP_MINIMAP2      = "NanoComp (BAM)",
  RESTRANDER             = "Restrander",
  MINIMAP2_ALIGN         = "minimap2",
  BAM_COVERAGE           = "Coverage",
  BAMBU                  = "Bambu",
  VALIDATE_BAMBU_OUTPUTS = "Validate outputs"
)

# Steps shown in the disk I/O part (one shape each; at most 5 fillable shapes)
IO_STEPS <- c("RESTRANDER", "MINIMAP2_ALIGN", "BAM_COVERAGE",
                   "BAMBU", "VALIDATE_BAMBU_OUTPUTS")

# Which parts to draw, top to bottom (set PULPO_PARTS before source() to change).
# Available: "totals", "steps", "profile", "io"
PARTS <- if (exists("PULPO_PARTS", envir = globalenv())) {
  get("PULPO_PARTS", envir = globalenv())
} else c("totals", "steps", "profile", "io")
stopifnot(all(PARTS %in% c("totals", "steps", "profile", "io")))
# Height of each part in the assembled panel (mm)
PART_H <- c(totals = 52, steps = 125, profile = 100, io = 85)

# Metrics in the steps part: column in `proc`, facet title, how to pool processes
STEP_METRICS <- tribble(
  ~id,           ~title,                              ~pool,
  "real_h",      "Realtime",                          "sum",
  "peak_rss_gb", "Peak RAM, largest task (GB)",       "max"
)

# ---- 2. Theme (matches Fig. 2 / Fig. 3 styling) --------------------------------
pick_font <- function(candidates = c("Source Sans 3", "Source Sans Pro",
                                     "Arial", "Helvetica")) {
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    fams <- unique(systemfonts::system_fonts()$family)
    hit  <- candidates[candidates %in% fams]
    if (length(hit)) return(hit[1])
  }
  "sans"
}
FONT      <- pick_font()
BASE_SIZE <- 7                      # pt; ED figures are printed at ~180 mm wide
TXT_SIZE  <- (BASE_SIZE - 1) / .pt  # geom_text size (mm) for 6-pt labels

theme_pulpo <- function(base_size = BASE_SIZE, base_family = FONT) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.4),
      strip.text       = element_text(colour = "black", margin = margin(2, 2, 2, 2)),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      axis.text        = element_text(colour = "black"),
      plot.title       = element_text(size = rel(1.25), hjust = 0, margin = margin(b = 2)),
      plot.subtitle    = element_text(colour = "grey35", hjust = 0, margin = margin(b = 4)),
      plot.title.position = "plot",
      legend.key.size  = unit(3, "mm"),
      legend.background = element_blank(),
      legend.title     = element_blank(),
      legend.text      = element_text(size = rel(0.95)),
      plot.margin      = margin(6, 8, 10, 6)   # t, r, b, l (pt): room between parts
    )
}

fmt_num <- function(x) {
  ifelse(is.na(x), NA_character_,
  ifelse(x >= 100,  sprintf("%.0f", x),
  ifelse(x >= 10,   sprintf("%.1f", x),
  ifelse(x >= 0.01, sprintf("%.2f", x), "<0.01"))))
}

# Hours -> "8h 37m" / "37m" / "45s" (as in the per-task columns of ST1-ST4).
# short = TRUE drops ":00" minutes ("4h"), for axis labels.
fmt_hm <- function(h, short = FALSE) {
  m   <- round(h * 60)
  out <- ifelse(m >= 60,
                ifelse(short & m %% 60 == 0,
                       sprintf("%dh", m %/% 60),
                       sprintf("%dh %02dm", m %/% 60, m %% 60)),
                ifelse(m >= 1, sprintf("%dm", m), sprintf("%ds", round(h * 3600))))
  out[!is.na(h) & h == 0] <- "0"
  out[is.na(h)] <- NA_character_
  out
}
# Axis breaks at round clock steps (5 min ... 24 h), at most n intervals
breaks_hm <- function(n = 4) {
  function(lim) {
    steps <- c(5, 10, 15, 30, 60, 120, 180, 240, 360, 720, 1440) / 60
    st <- steps[which(diff(range(lim)) / steps <= n)[1]]
    if (is.na(st)) st <- tail(steps, 1)
    seq(0, lim[2], by = st)
  }
}
time_axis <- function() scale_x_continuous(breaks = breaks_hm(), labels = function(x) fmt_hm(x, short = TRUE),
                                expand = expansion(mult = c(0, 0.35)))
num_axis  <- function() scale_x_continuous(breaks = breaks_extended(n = 3), labels = label_number(),
                                expand = expansion(mult = c(0, 0.35)))

# ---- 3. Read supplementary tables ----------------------------------------------
# Each ST sheet has a title, a description and a blank row above the header.
# Find the header row by its first cell rather than hard-coding `skip`.
read_st <- function(path, sheet, header_key) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE, col_types = "text",
                    .name_repair = "minimal")
  hdr <- which(raw[[1]] == header_key)[1]
  if (is.na(hdr)) stop(sprintf("Header '%s' not found in sheet %s", header_key, sheet))
  read_excel(path, sheet = sheet, skip = hdr - 1, .name_repair = "minimal")
}

## 3a. Per-process metrics (ST1-ST4)
read_process_sheet <- function(sheet, run) {
  d <- read_st(xlsx, sheet, "#")
  if (!"Realtime (h)" %in% names(d)) {
    stop("Sheet ", sheet, " has no 'Realtime (h)' column: use the updated ",
         "pulposeq_supplementary_tables.xlsx (step time incl. queue + realtime columns).")
  }
  d <- d[!is.na(d[["Nextflow process"]]) & d[["Nextflow process"]] != "Total / peak", ]
  tibble(
    run         = run,
    process     = d[["Nextflow process"]],
    tasks       = as.numeric(d[["Tasks"]]),
    # realtime: time with >= 1 task executing, queue waiting excluded; parallel
    # tasks counted once (ST1 legend). The plots use this, not "incl. queue".
    real_h      = as.numeric(d[["Realtime (h)"]]),
    queue_h     = as.numeric(d[["Step time incl. queue (h)"]]),
    cpu_h       = as.numeric(d[["CPU time (CPU-h)"]]),
    mean_cores  = as.numeric(d[["Mean CPU cores"]]),
    max_cores   = as.numeric(d[["Max CPU cores (one task)"]]),
    peak_rss_gb = as.numeric(d[["Peak RSS, max task (GB)"]]),
    read_gb     = as.numeric(d[["Data read (GB)"]]),
    written_gb  = as.numeric(d[["Data written (GB)"]])
  )
}

proc <- bind_rows(Map(read_process_sheet, RUNS$sheet, RUNS$run)) |>
  mutate(step = sub("^.*:", "", process)) |>
  left_join(STEP_GROUPS, by = "step") |>
  left_join(STAGE_TBL, by = "step") |>
  mutate(group = coalesce(group, OTHER),
         run   = factor(run, levels = RUNS$run))

unmapped <- unique(proc$step[is.na(proc$stage)])
if (length(unmapped)) {
  stop("Processes without a stage in STAGE_MAP: ", paste(unmapped, collapse = ", "))
}
proc <- proc |> mutate(stage = factor(stage, levels = names(STAGE_MAP)))

## 3b. Run-level totals (ST5)
st5 <- read_st(xlsx, "ST5", "Metric")
get_metric <- function(metric) {
  row <- st5[st5$Metric == metric, , drop = FALSE]
  if (nrow(row) != 1) stop("ST5 metric not found (or duplicated): ", metric)
  unlist(row[RUNS$st5_col], use.names = FALSE)
}

runs <- tibble(
  run           = factor(RUNS$run, levels = RUNS$run),
  st5_col       = RUNS$st5_col,
  samples       = as.numeric(get_metric("Samples")),
  elapsed_h     = as.numeric(get_metric("Total realtime")),
  cpu_h         = as.numeric(get_metric("Total CPU time")),
  peak_ram_gb   = as.numeric(get_metric("Peak RAM, single task")),
  peak_ram_step = as.character(get_metric("Step with peak RAM")),
  peak_cores    = as.numeric(get_metric("Peak CPU cores, single task")),
  peak_cpu_step = as.character(get_metric("Step with peak CPU")),
  written_tb    = as.numeric(get_metric("Total data written"))
)

## 3c. Optional: input size, for per-Gb efficiency
if (!is.na(input_stats)) {
  ins <- read.delim(input_stats, check.names = FALSE)
  stopifnot(all(c("run", "gigabases") %in% names(ins)))
  runs <- runs |>
    left_join(ins |> select(st5_col = run, gigabases), by = "st5_col") |>
    mutate(cpu_h_per_gb = cpu_h / gigabases)
}

## 3d. Consistency check: per-process sums must match the ST5 totals
chk <- proc |>
  group_by(run) |>
  summarise(cpu_sum = sum(cpu_h), .groups = "drop") |>
  left_join(runs |> select(run, cpu_h), by = "run")
bad <- with(chk, abs(cpu_sum - cpu_h) > 0.01)
if (any(bad)) {
  warning("Per-process sums differ from ST5 totals for: ",
          paste(chk$run[bad], collapse = ", "))
}

# ---- 4. Part (1): one number per run -------------------------------------------
pretty_step <- function(x) {
  lab <- unname(STEP_LABELS[toupper(x)])
  ifelse(is.na(lab), x, lab)
}
# "Bambu" if one step holds the peak in every run, else "run: step; ..."
peak_who <- function(steps) {
  if (length(unique(steps)) == 1) {
    paste(pretty_step(steps[1]), "in every run")
  } else {
    paste(runs$run, pretty_step(steps), sep = ": ", collapse = "; ")
  }
}
peak_note <- sprintf("Peak RAM: %s; peak CPU cores: %s",
                     peak_who(runs$peak_ram_step), peak_who(runs$peak_cpu_step))

TOTAL_TIME <- "Realtime, whole run"
totals_long <- bind_rows(
  tibble(run = runs$run, metric = TOTAL_TIME,                   value = runs$elapsed_h),
  if ("cpu_h_per_gb" %in% names(runs))
    tibble(run = runs$run, metric = "CPU time per Gb input (CPU-h)", value = runs$cpu_h_per_gb),
  tibble(run = runs$run, metric = "Peak RAM, largest task (GB)",  value = runs$peak_ram_gb),
  tibble(run = runs$run, metric = "Peak CPU cores, largest task", value = runs$peak_cores),
  tibble(run = runs$run, metric = "Data written, whole run (TB)", value = runs$written_tb)
) |>
  mutate(metric = factor(metric, levels = unique(metric)),
         run    = factor(run, levels = rev(RUNS$run)),          # first run on top
         label  = ifelse(metric == TOTAL_TIME, fmt_hm(value), fmt_num(value)))

p_totals <- ggplot(totals_long, aes(x = value, y = run, fill = run)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), hjust = -0.15,
            size = TXT_SIZE, family = FONT) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  num_axis() +
  ggh4x::facetted_pos_scales(x = list(metric == TOTAL_TIME ~ time_axis(),
                                      metric != TOTAL_TIME ~ num_axis())) +
  scale_fill_manual(values = RUN_COLS, guide = "none") +
  labs(x = NULL, y = NULL,
       title    = "Run totals",
       subtitle = peak_note) +
  theme_pulpo()

# ---- 5. Part (2): main steps, one bar per run ----------------------------------
pool_fun <- function(x, how) if (how == "max") max(x) else sum(x)

step_long <- bind_rows(lapply(seq_len(nrow(STEP_METRICS)), function(i) {
  m <- STEP_METRICS[i, ]
  proc |>
    group_by(run, group) |>
    summarise(value = pool_fun(.data[[m$id]], m$pool),
              n_proc = n(), .groups = "drop") |>
    mutate(metric = m$title)
}))

# Name the pooled row with how many processes it holds
n_other   <- range(step_long$n_proc[step_long$group == OTHER])
other_lab <- if (diff(n_other) == 0) {
  sprintf("%s (%d processes)", OTHER, n_other[1])
} else {
  sprintf("%s (%d–%d processes)", OTHER, n_other[1], n_other[2])
}

step_long <- step_long |>
  mutate(group  = ifelse(group == OTHER, other_lab, group),
         group  = factor(group, levels = c(STEP_ORDER, other_lab)),  # rows, top to bottom
         metric = factor(metric, levels = STEP_METRICS$title),
         run_y  = factor(run, levels = rev(RUNS$run)),              # first run on top
         label  = ifelse(metric == "Realtime", fmt_hm(value), fmt_num(value)))
# Only steps a run executed are present, so rows hold 1-4 bars (no empty slots).

# One cell per step x metric, each with its own x axis (ggh4x). Rows are as
# tall as the number of runs that executed the step.
step_facets <- ggh4x::facet_grid2(group ~ metric, scales = "free",
                                  independent = "x", space = "free_y",
                                  switch = "y")

p_steps <- ggplot(step_long, aes(x = value, y = run_y, fill = run)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = label), hjust = -0.15, size = TXT_SIZE, family = FONT) +
  step_facets +
  num_axis() +
  ggh4x::facetted_pos_scales(x = list(metric == "Realtime" ~ time_axis(),
                                      metric != "Realtime" ~ num_axis())) +
  scale_fill_manual(values = RUN_COLS, breaks = RUNS$run,
                    guide = guide_legend(nrow = 1)) +
  labs(x = NULL, y = NULL,
       title    = "Main pipeline steps",
       subtitle = paste0("Each cell on its own x axis. Realtime: time with at least one task of the step ",
                         "executing, queue waiting excluded\n(parallel tasks counted once; pooled rows summed). ",
                         "Peak RAM: largest single task")) +
  theme_pulpo() +
  theme(legend.position      = "top",
        legend.justification = "left",
        legend.margin        = margin(0, 0, 2, 0),
        axis.text.y          = element_blank(),
        axis.ticks.y         = element_blank(),
        axis.text.x          = element_text(size = rel(0.85)),
        strip.placement      = "outside",
        strip.background.y   = element_blank(),
        strip.text.y.left    = element_text(angle = 0, hjust = 1, margin = margin(r = 4)),
        panel.spacing.y      = unit(5, "pt"),
        panel.spacing.x      = unit(4, "pt"))

# ---- 6. Part (3): per-step resource profile ------------------------------------
strip_lab <- setNames(sprintf("%s (n = %d)", runs$run, runs$samples), as.character(runs$run))

SIZE_RANGE <- c(0.8, 6)   # point diameter range (mm-ish ggplot units)
profile_df <- proc |>
  mutate(label = unname(STEP_LABELS[step]),
         # same mapping as scale_size(range = SIZE_RANGE) (area-proportional),
         # passed to ggrepel so labels keep clear of the big points
         pt_size = SIZE_RANGE[1] + sqrt(rescale(max_cores)) * diff(SIZE_RANGE))

p_profile <- ggplot(profile_df, aes(x = real_h, y = peak_rss_gb)) +
  geom_point(aes(size = max_cores, fill = stage), shape = 21, colour = "black",
             stroke = 0.25, alpha = 0.9) +
  facet_wrap(~ run, nrow = 1, labeller = as_labeller(strip_lab)) +
  # square-root time axis: spreads the hour range, where the large steps sit,
  # while keeping sub-minute processes visible near zero
  scale_x_continuous(transform = "sqrt",
                     breaks = c(1 / 60, 0.5, 1, 4, 8, 12, 24),
                     labels = c("1m", "30m", "1h", "4h", "8h", "12h", "24h"),
                     limits = c(0, 24)) +
  scale_y_log10(breaks = 10^(-2:3),
                labels = c("0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(5e-3, 1500)) +
  scale_size(name = "CPU cores,\nbusiest task", range = SIZE_RANGE,
             breaks = c(1, 5, 10, 20)) +
  scale_fill_manual(values = STAGE_COLS, name = NULL,
                    guide = guide_legend(nrow = 2, byrow = TRUE,
                                         override.aes = list(size = 2.5))) +
  labs(x = "Realtime (square-root scale)",
       y = "Peak RAM, largest task (GB, log)",
       title    = "Per-step resource profile",
       subtitle = paste0("One point per Nextflow process; point area = CPU cores used by its ",
                         "busiest task; colour = pipeline stage (Fig. 2a)")) +
  theme_pulpo() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = rel(1)),
        legend.spacing.x = unit(1, "mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Manual label offsets for crowded spots: x in sqrt(hours) units (axis spans
# 0-4.9), y in log10 units (1 = one decade).
# Dong "Validate outputs" sits beside the NanoComp points, so push it up-left.
LABEL_NUDGE <- tribble(
  ~run,               ~step,                    ~nx,  ~ny,
  "Dong · ONT cDNA",  "VALIDATE_BAMBU_OUTPUTS", -0.7,  0.8
)
label_df <- profile_df |>
  filter(!is.na(label)) |>
  left_join(LABEL_NUDGE |> mutate(run = factor(run, levels = RUNS$run)),
            by = c("run", "step")) |>
  mutate(nx = coalesce(nx, 0), ny = coalesce(ny, 0))

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_profile <- p_profile +
    ggrepel::geom_text_repel(
      data = label_df, aes(label = label, point.size = pt_size),
      nudge_x = label_df$nx, nudge_y = label_df$ny,
      size = TXT_SIZE, family = FONT,
      min.segment.length = 0, segment.size = 0.2, segment.colour = "grey30",
      box.padding = 0.35, point.padding = 0.1, force = 2,
      max.overlaps = Inf, seed = 1)
} else {
  message("Install 'ggrepel' for non-overlapping labels in the profile part.")
  p_profile <- p_profile +
    geom_text(aes(label = label), na.rm = TRUE, vjust = -1.2,
              size = TXT_SIZE, family = FONT)
}

# ---- 7. Part (4): disk I/O per step --------------------------------------------
# Data read vs data written by each main step, one point per run. Equal log
# axes, so the dashed diagonal (written = read) runs at 45 degrees: steps on it
# pass their data through (restranded FASTQ, BAM); steps far below it read a
# lot and write little. "Data read" is Nextflow rchar, so it also counts
# reference and index reads, not only the step's input files.
IO_SHAPES <- c(21, 22, 24, 23, 25)[seq_along(IO_STEPS)]
io_levels <- pretty_step(IO_STEPS)

io_df <- proc |>
  filter(step %in% IO_STEPS) |>
  select(run, step, read_gb, written_gb) |>
  mutate(step_lab = factor(pretty_step(step), levels = io_levels))

IO_LIMITS <- c(0.1, 3000)
io_breaks <- 10^(-1:3)
io_labels <- c("0.1", "1", "10", "100", "1,000")
if (any(unlist(io_df[c("read_gb", "written_gb")]) < IO_LIMITS[1] |
        unlist(io_df[c("read_gb", "written_gb")]) > IO_LIMITS[2])) {
  warning("Some I/O values fall outside IO_LIMITS; widen IO_LIMITS.")
}

p_io <- ggplot(io_df, aes(x = read_gb, y = written_gb)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey55", linewidth = 0.3) +
  annotate("text", x = 6, y = 14, label = "written = read", angle = 45,
           colour = "grey40", size = TXT_SIZE, family = FONT) +
  annotate("text", x = IO_LIMITS[2] / 1.6, y = IO_LIMITS[1] * 1.8,
           label = "reads far more than it writes", hjust = 1,
           colour = "grey40", size = TXT_SIZE, family = FONT) +
  geom_point(aes(fill = run, shape = step_lab), size = 2.2,
             colour = "black", stroke = 0.3) +
  scale_x_log10(breaks = io_breaks, labels = io_labels) +
  scale_y_log10(breaks = io_breaks, labels = io_labels) +
  coord_fixed(xlim = IO_LIMITS, ylim = IO_LIMITS, expand = FALSE) +
  scale_fill_manual(values = RUN_COLS, breaks = RUNS$run, name = "Run",
                    guide = guide_legend(override.aes = list(shape = 21, size = 2.5),
                                         order = 1)) +
  scale_shape_manual(values = setNames(IO_SHAPES, io_levels), name = "Step",
                     guide = guide_legend(override.aes = list(fill = "grey85", size = 2.5),
                                          order = 2)) +
  labs(x = "Data read (GB, log)", y = "Data written (GB, log)",
       title    = "Disk I/O per step",
       subtitle = "Main steps, one point per run; data read includes reference and index reads") +
  theme_pulpo() +
  theme(legend.position = "right", legend.box = "horizontal",
        legend.title = element_text(size = rel(1)),
        legend.justification = "top")

# ---- 8. Assemble and save ------------------------------------------------------
part_plots <- list(totals = p_totals, steps = p_steps,
                   profile = p_profile, io = p_io)
panel_h <- sum(PART_H[PARTS]) + 12      # + main title and subtitle

# The I/O plot has a fixed 1:1 aspect ratio. Aligned with the other rows it
# would shrink every row to its width, so it enters the stack as a free-standing
# element (wrap_elements) and keeps its own width.
stack <- part_plots[PARTS]
# Indent part titles so the A/B/C tag sits to their left, as in Fig. 2/3
TAG_INDENT <- 14   # pt
stack <- lapply(stack, function(p) {
  p + theme(plot.title    = element_text(margin = margin(l = TAG_INDENT, b = 2)),
            plot.subtitle = element_text(margin = margin(l = TAG_INDENT, b = 4)))
})
if ("io" %in% PARTS) stack[["io"]] <- wrap_elements(full = stack[["io"]])

panel <- wrap_plots(stack, ncol = 1) +
  plot_layout(heights = unname(PART_H[PARTS])) +
  plot_annotation(
    title    = "Computational resources",
    subtitle = paste0("Nextflow trace metrics of the final runs (Supplementary Tables 1–5); ",
                      "LongBench: 8 samples per protocol; Dong: 6 samples"),
    tag_levels = "A",                      # parts tagged A, B, C, ... in PARTS order
    theme    = theme_pulpo() + theme(plot.title = element_text(size = rel(1.5)))
  ) &
  theme(plot.tag          = element_text(family = FONT, face = "bold", hjust = 0, vjust = 1,
                                         size = BASE_SIZE * 1.6, colour = "black"),
        plot.tag.position = c(0, 1),         # top-left corner of each part
        plot.tag.location = "plot")

save_fig <- function(plot, name, width_mm, height_mm) {
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(outdir, paste0(name, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm", device = pdf_dev)
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(outdir, paste0(name, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = 600,
         device = png_dev, bg = "white")
  # SVG: svglite keeps text as editable <text> (needs the font installed to
  # display as intended); grDevices::svg turns text into outlined paths.
  svg_dev <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
  if (!requireNamespace("svglite", quietly = TRUE)) {
    message("Install 'svglite' for SVGs with editable text; using grDevices::svg.")
  }
  ggsave(file.path(outdir, paste0(name, ".svg")), plot,
         width = width_mm, height = height_mm, units = "mm", device = svg_dev)
  message("Saved: ", file.path(outdir, name), ".{pdf,png,svg}")
}

save_fig(panel, "resources_panel", 180, panel_h)
for (nm in names(part_plots)) {
  save_fig(part_plots[[nm]], paste0("resources_", nm), 180, PART_H[[nm]])
}

# ---- 9. Source data ------------------------------------------------------------
write.csv(runs,      file.path(outdir, "source_data_resources_runs.csv"),  row.names = FALSE)
write.csv(proc,      file.path(outdir, "source_data_resources_processes.csv"), row.names = FALSE)
write.csv(step_long |> select(-run_y, -label),
          file.path(outdir, "source_data_resources_steps.csv"), row.names = FALSE)
write.csv(io_df, file.path(outdir, "source_data_resources_io.csv"), row.names = FALSE)
