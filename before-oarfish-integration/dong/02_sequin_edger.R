## ===========================================================================
## 02_sequin_edger.R - sequin ground truth for the edgeR workflow (Dong)   (Fig. 2c)
## Dong et al. 2023 (Nat Methods 20:1810): H1975 + sequins Mix A,
## HCC827 + sequins Mix B. Expected log2FC = log2(MixA / MixB), i.e. the
## same direction as the main DE (H1975 / HCC827).
## Same test as 01_de_edger.R: edgeR v4 glmQLFit (robust) -> glmTreat
## (lfc 0.5), FDR 0.05, predFC (prior.count 2) for the display log2FC;
## filter >= 10 reads in >= 3 samples. Library sizes and TMM factors from the
## ratio-1 control sequins (|expected log2FC| < 0.25).
## Truth and scoring identical to longbench/01_sequin_edger.R.
##
## Inputs : dong/bambu_validated/BambuOutput_counts_{gene,transcript}_validated.txt,
##          dong/samplesheet.csv, Anaquin (RnaQuinIsoformMixture)
## Outputs: results/sequins_edger/
##   <level>/sequin_de.csv, sequin_de_results.csv, summary.csv
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(Anaquin)
})

source("C:/Users/mbarb/OneDrive/paper/downstream/dong/00_paths.R")   # data_dir, sheet, seq_dir

## ------------------------------------------------------------- settings ---
tool        <- "edgeR"
outroot     <- seq_dir
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
grp_ref     <- "TKI_sensitive"    # HCC827 + Mix B (denominator)
grp_num     <- "TKI_resistant"    # H1975  + Mix A (numerator)
alpha       <- 0.05               # FDR, as in the main DE
lfc         <- 0.5                # H0 |log2FC| <= lfc, as in the main DE
min_count   <- 10                 # filter, as in the main DE
min_samples <- 3
ctrl_lfc    <- 0.25               # |expected log2FC| below this = ratio-1 control sequin

inputs <- tribble(
  ~level,       ~file,
  "gene",       "BambuOutput_counts_gene_validated.txt",
  "transcript", "BambuOutput_counts_transcript_validated.txt"
)
## full counts only, as in the main DE (unique transcript counts were tested
## once and gave equivalent sequin performance, so they were dropped)

## ---------------------------------------------- expected concentrations ---
data("RnaQuinIsoformMixture", package = "Anaquin")
iso_mix <- as_tibble(RnaQuinIsoformMixture) %>%
  transmute(id = as.character(Name), A = MixA, B = MixB)
gene_mix <- iso_mix %>%
  mutate(id = sub("_[0-9]+$", "", id)) %>%
  group_by(id) %>%
  summarise(A = sum(A), B = sum(B), .groups = "drop")
expected <- list(transcript = iso_mix, gene = gene_mix)

ss <- read.csv(sheet)

## ------------------------------------------------------------ functions ---
read_sequins <- function(file) {
  x <- read.table(file.path(data_dir, "bambu_validated", file), header = TRUE,
                  row.names = 1, sep = "\t", check.names = FALSE)
  x$GENEID <- NULL
  x <- x[grepl("^R[12]_", rownames(x)), , drop = FALSE]
  m <- round(as.matrix(x))
  colnames(m) <- sub("[.]bam$", "", colnames(m))
  m
}

group_of <- function(samples) {
  g <- gsub("-", "_", ss$group[match(samples, ss$sample)])
  stopifnot(!anyNA(g))
  factor(g, levels = c(grp_ref, grp_num))
}

## everything the test needs for input row i of `inputs`
##   m / e      all sequins with an Anaquin concentration
##   ctrl_all   ratio-1 control sequins among m (for normalisation)
##   m2 / ctrl  sequins passing the DE filter + their controls
##   exp_lfc    expected log2FC of m2 (H1975 / HCC827)
prep_input <- function(i) {
  lv  <- inputs$level[i]
  m   <- read_sequins(inputs$file[i])
  e   <- expected[[lv]]
  no_exp <- setdiff(rownames(m), e$id)
  if (length(no_exp))
    message(sprintf("[%s] %d sequins without an Anaquin concentration (dropped): %s",
                    lv, length(no_exp), paste(head(no_exp, 10), collapse = ", ")))
  ids <- intersect(rownames(m), e$id)
  m <- m[ids, , drop = FALSE]
  e <- e[match(ids, e$id), ]
  grp <- group_of(colnames(m))

  exp_all  <- setNames(log2(e$A / e$B), ids)
  ctrl_all <- abs(exp_all) < ctrl_lfc
  keep     <- rowSums(m >= min_count) >= min_samples
  m2       <- m[keep, , drop = FALSE]
  message(sprintf("[%s] %d sequins, %d pass the filter, %d ratio-1 controls, %d truly DE",
                  lv, nrow(m), nrow(m2), sum(ctrl_all[keep]),
                  sum(abs(exp_all[keep]) > lfc)))
  list(level = lv, tag = lv, grp = grp,
       m2 = m2, exp_lfc = exp_all[keep], ctrl = ctrl_all[keep])
}

## r: tibble(id, obs_lfc, disp_lfc, pvalue, padj, mean) from the test
score_calls <- function(r, exp_lfc) {
  r %>%
    mutate(exp_lfc = unname(exp_lfc[id]),
           ratio   = round(exp_lfc),                  # Anaquin ratio group
           truth   = abs(exp_lfc) > lfc,
           called  = !is.na(padj) & padj < alpha,
           outcome = case_when(called & truth & sign(obs_lfc) == sign(exp_lfc) ~ "true positive",
                               called & truth                                  ~ "wrong sign",
                               called & !truth                                 ~ "false positive",
                               !called & truth                                 ~ "false negative",
                               TRUE                                            ~ "true negative"))
}

auc_fn <- function(score, truth) {
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## one row per level
##   sensitivity = TP / truly DE            (wrong-sign calls are not TP)
##   fdr         = (FP + wrong sign) / called
##   precision   = TP / called = 1 - fdr
##   f1          = 2 * precision * sensitivity / (precision + sensitivity)
## `called` is reused below as the NUMBER of calls on purpose (it is a count
## from that line on); the per-sequin logical is only used through `outcome`.
summarise_calls <- function(res) {
  res %>%
    filter(!is.na(pvalue), is.finite(obs_lfc)) %>%
    group_by(tool, level) %>%
    summarise(tested = n(), truly_de = sum(truth), called = sum(called),
              true_pos = sum(outcome == "true positive"),
              false_pos = sum(outcome == "false positive"),
              wrong_sign = sum(outcome == "wrong sign"),
              false_neg = sum(outcome == "false negative"),
              true_neg = sum(outcome == "true negative"),
              sensitivity = true_pos / max(truly_de, 1),
              fdr = if (called > 0) (false_pos + wrong_sign) / called else NA_real_,
              precision = if (called > 0) true_pos / called else NA_real_,
              f1 = if (is.na(precision) || precision + sensitivity == 0) NA_real_ else
                     2 * precision * sensitivity / (precision + sensitivity),
              auc = auc_fn(-log10(pvalue), truth),
              slope = unname(coef(lm(obs_lfc ~ exp_lfc))[2]),
              r2 = summary(lm(obs_lfc ~ exp_lfc))$r.squared,
              rmse_test_lfc = sqrt(mean((obs_lfc - exp_lfc)^2)),
              rmse_display_lfc = sqrt(mean((disp_lfc - exp_lfc)^2, na.rm = TRUE)),
              .groups = "drop")
}

## DGEList with library sizes + TMM factors from control sequins
norm_dge <- function(m, grp, ctrl) {
  y <- DGEList(m, group = grp)
  if (sum(ctrl) >= 3) {
    yc <- normLibSizes(DGEList(m[ctrl, , drop = FALSE]))
    y$samples$lib.size     <- yc$samples$lib.size
    y$samples$norm.factors <- yc$samples$norm.factors
  } else {
    y <- normLibSizes(y)
  }
  y
}

run_edger <- function(m, grp, ctrl) {
  y      <- norm_dge(m, grp, ctrl)
  design <- model.matrix(~ group, data.frame(group = grp))
  cn  <- paste0("group", grp_num)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt  <- topTags(glmTreat(fit, coef = cn, lfc = lfc), n = Inf, sort.by = "none")$table
  shr <- predFC(estimateDisp(y, design), design, prior.count = 2)[, cn]
  eff <- y$samples$lib.size * y$samples$norm.factors
  tibble(id = rownames(tt), obs_lfc = tt$logFC, disp_lfc = unname(shr[rownames(tt)]),
         pvalue = tt$PValue, padj = tt$FDR,
         mean = rowMeans(edgeR::cpm(y)[rownames(tt), , drop = FALSE]) * mean(eff) / 1e6)
}

## ----------------------------------------------------------------- run ----
res_all <- list()
for (i in seq_len(nrow(inputs))) {
  p <- prep_input(i)
  t <- score_calls(run_edger(p$m2, p$grp, p$ctrl), p$exp_lfc) %>%
    mutate(tool = tool, level = p$level, .before = 1)
  outdir <- file.path(outroot, p$tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write_csv(t, file.path(outdir, "sequin_de.csv"))
  res_all[[i]] <- t
}

res_all <- bind_rows(res_all)
summ    <- summarise_calls(res_all)
write_csv(res_all, file.path(outroot, "sequin_de_results.csv"))
write_csv(summ, file.path(outroot, "summary.csv"))
print(dplyr::select(summ, level, tested, truly_de, called, fdr, sensitivity, auc, slope))
message("Done: ", outroot)
