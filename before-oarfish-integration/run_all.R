#!/usr/bin/env Rscript
# =============================================================================
# run_all.R - runs every downstream script in order (Fig. 1 LongBench, Fig. 2 Dong)
#
# Each script runs in its OWN fresh R process (Rscript), so nothing one script
# defines can mask a function in the next one. The console shows one line per
# script (OK / FAILED / SKIPPED + time); everything the script prints goes to
# its log file in downstream/logs/<date_time>/.
#
# Use: Rscript "C:/Users/mbarb/OneDrive/paper/downstream/run_all.R"
#      or open it in RStudio and click Source. Edit the settings in section 0.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Settings -- edit here
# -----------------------------------------------------------------------------
BASE <- "C:/Users/mbarb/OneDrive/paper/downstream"

RUN_FOLDERS   <- c("longbench", "dong")  # which folders, in this order
SKIP          <- character(0)            # e.g. c("dong/01_de_edger.R") to reuse old results
ONLY          <- character(0)            # e.g. c("dong/04_isoformic.R"): run just these
STOP_ON_ERROR <- FALSE  # FALSE: a failed script only skips the scripts that need it
                        # TRUE : stop at the first failure

# scripts in run order, with the scripts they need (same folder) - see README.md
PIPELINE <- list(
  longbench = list(
    "01_sequin_edger.R"             = character(0),
    "02_sirv_ercc.R"                = character(0),
    "03_ujc_report.R"               = character(0),
    "04_ujc_length.R"               = "03_ujc_report.R",
    "05_read_support.R"             = character(0),
    "06_supp_tables_ujc.R"          = "03_ujc_report.R",
    "07_fig1_panels.R"              = c("01_sequin_edger.R", "02_sirv_ercc.R", "03_ujc_report.R",
                                        "04_ujc_length.R", "05_read_support.R")
  ),
  dong = list(
    "01_de_edger.R"                 = character(0),
    "02_sequin_edger.R"             = character(0),
    "03_cis_circos.R"               = "01_de_edger.R",
    "04_isoformic.R"                = c("01_de_edger.R", "03_cis_circos.R"),
    "05_read_support.R"             = character(0),
    "06_fig2_panels.R"              = c("01_de_edger.R", "02_sequin_edger.R", "05_read_support.R"),
    "07_circos_coverage_biotypes.R" = character(0)
  )
)

# -----------------------------------------------------------------------------
# 1. Set-up
# -----------------------------------------------------------------------------
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("Rscript not found at ", rscript)
if (!dir.exists(BASE)) stop("BASE folder not found: ", BASE)

stamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_dir <- file.path(BASE, "logs", stamp)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

steps <- do.call(rbind, lapply(RUN_FOLDERS, function(fd) {
  p <- PIPELINE[[fd]]
  if (is.null(p)) stop("Unknown folder in RUN_FOLDERS: ", fd)
  data.frame(folder = fd, script = names(p), id = file.path(fd, names(p)),
             stringsAsFactors = FALSE)
}))
needs_of <- function(id) {
  fd <- dirname(id)
  file.path(fd, PIPELINE[[fd]][[basename(id)]])
}
bad <- setdiff(c(SKIP, ONLY), steps$id)
if (length(bad)) stop("Not in the pipeline: ", paste(bad, collapse = ", "),
                      " (write them as folder/script.R)")

fmt_time <- function(s) sprintf("%02d:%02d", as.integer(s) %/% 60, as.integer(s) %% 60)
status   <- setNames(rep(NA_character_, nrow(steps)), steps$id)

cat(sprintf("run_all.R  |  R %s  |  logs: %s\n\n", getRversion(), log_dir))

# -----------------------------------------------------------------------------
# 2. Run
# -----------------------------------------------------------------------------
t_all <- Sys.time()
for (i in seq_len(nrow(steps))) {
  id  <- steps$id[i]
  scr <- file.path(BASE, id)
  log <- file.path(log_dir, sub("/", "__", sub("\\.R$", ".log", id)))
  cat(sprintf("[%2d/%d] %-42s ", i, nrow(steps), id))

  ## not selected / skipped on purpose / missing / a needed script did not succeed
  if (length(ONLY) && !id %in% ONLY) { status[id] <- "not selected"; cat("-\n"); next }
  if (id %in% SKIP)                  { status[id] <- "SKIPPED (setting)"; cat("SKIPPED (setting)\n"); next }
  if (!file.exists(scr))             { status[id] <- "FAILED (file missing)"; cat("FAILED: file missing\n") }
  else {
    failed_need <- intersect(needs_of(id), names(status)[grepl("^(FAILED|SKIPPED \\(needs)", status)])
    if (length(failed_need)) {
      status[id] <- "SKIPPED (needs a failed script)"
      cat("SKIPPED: needs", paste(basename(failed_need), collapse = ", "), "\n")
      next
    }
    t0 <- Sys.time()
    owd <- setwd(dirname(scr))
    ## stdout and stderr (messages, warnings, errors) captured together, then
    ## written to the log (the log appears when the script ends)
    out <- tryCatch(
      suppressWarnings(system2(rscript, c("--no-save", "--no-restore", shQuote(scr)),
                               stdout = TRUE, stderr = TRUE)),
      error = function(e) structure(conditionMessage(e), status = 1L),
      finally = setwd(owd))
    writeLines(out, log)
    code <- if (is.null(attr(out, "status"))) 0L else attr(out, "status")
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (identical(as.integer(code), 0L)) {
      status[id] <- "OK"
      cat("OK  ", fmt_time(dt), "\n")
    } else {
      status[id] <- "FAILED"
      cat("FAILED", fmt_time(dt), "\n")
      txt <- if (file.exists(log)) readLines(log, warn = FALSE) else character(0)
      cat("        last lines of", basename(log), ":\n")
      cat(paste0("        | ", tail(txt, 12)), sep = "\n")
    }
  }
  if (STOP_ON_ERROR && grepl("^FAILED", status[id])) {
    cat("\nSTOP_ON_ERROR = TRUE: stopping here.\n")
    break
  }
}

# -----------------------------------------------------------------------------
# 3. Summary
# -----------------------------------------------------------------------------
res <- data.frame(script = names(status), status = ifelse(is.na(status), "not run", status),
                  stringsAsFactors = FALSE)
write.csv(res, file.path(log_dir, "summary.csv"), row.names = FALSE)

## the font each script registered (plot_theme.R prints "Figure font: ...")
fonts <- unlist(lapply(list.files(log_dir, "\\.log$", full.names = TRUE), function(f) {
  l <- grep("^Figure font:", readLines(f, warn = FALSE), value = TRUE)
  if (length(l)) paste0(sub("\\.log$", "", basename(f)), ": ", sub("^Figure font: ", "", l[1]))
}))

cat(sprintf("\nDone in %s.  OK: %d  |  failed: %d  |  skipped: %d\n",
            fmt_time(as.numeric(difftime(Sys.time(), t_all, units = "secs"))),
            sum(res$status == "OK"), sum(grepl("^FAILED", res$status)),
            sum(grepl("^SKIPPED", res$status))))
if (length(fonts)) cat("Figure font per script:\n", paste0("  ", fonts, "\n"), sep = "")
cat("Logs and summary.csv: ", log_dir, "\n", sep = "")
if (any(grepl("^FAILED", res$status)))
  cat("Failed: ", paste(res$script[grepl("^FAILED", res$status)], collapse = ", "),
      "  -> open its .log file\n", sep = "")
