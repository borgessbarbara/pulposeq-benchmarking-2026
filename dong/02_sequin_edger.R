## ===========================================================================
## 02_sequin_edger.R - sequin ground truth for the edgeR workflow (Dong)   (Fig. 2c)
## Dong et al. 2023 (Nat Methods 20:1810): H1975 + sequins Mix A,
## HCC827 + sequins Mix B. Expected log2FC = log2(MixA / MixB), i.e. the
## same direction as the main DE (H1975 / HCC827).
## Same test as 01_de_edger.R: edgeR v4 glmQLFit (robust) -> glmTreat
## (lfc 0.5), FDR 0.05, predFC (prior.count 2) for the display log2FC;
## filter >= 10 reads in >= 3 samples. Library sizes and TMM factors from the
## ratio-1 control sequins (|expected log2FC| < 0.25).
## Truth, scoring, detection rule and diagnostics identical to
## longbench/01_sequin_edger.R (do not edit one without the other).
##
## CHANGES vs previous version (harmonised with longbench/01_sequin_edger.R)
##   * summary.csv keeps every previous column and adds:
##       detected / detected_any (>= 3 reads in >= 2 libraries / >= 1 read),
##       spearman, median_abs_err, n_one_mix_zero, slope_no_zero, r2_no_zero
##       (without sequins with 0 reads in every library of one group: the
##       unshrunk logFC is then ~ +/-10-15 and dominates R2), slope_shrunk /
##       r2_shrunk (predFC, prior.count = 2), truly_de_in_mix,
##       sensitivity_in_mix (TP / every truly-DE sequin in the mix)
##   * sequin_de.csv / sequin_de_results.csv add mean_A, mean_B (reads per
##     library in the Mix A / Mix B samples) and one_mix_zero
##   * sequin_fc_outliers.csv: isoform- or gene-level outliers with every
##     sibling isoform, length, capture efficiency and the gene-level result
##   * identification on chrIS (transcript level): detection, novel models
##     (raw Bambu vs validated, if the GTFs exist), 70% detection limit,
##     reads per library, capture-vs-length rho -> sequin_identification.csv,
##     sequin_per_isoform.csv, sequin_detection_by_conc.csv,
##     novel_models_on_sequins.csv
##
## Inputs : dong/bambu_validated/BambuOutput_counts_{gene,transcript}_validated.txt,
##          dong/bambu{,_validated}/ GTFs (optional), dong/samplesheet.csv,
##          Anaquin (RnaQuinIsoformMixture)
## Outputs: results/sequins_edger/
##   <level>/sequin_de.csv, sequin_de_results.csv, summary.csv,
##   sequin_fc_outliers.csv, sequin_identification.csv, sequin_per_isoform.csv,
##   sequin_detection_by_conc.csv, novel_models_on_sequins.csv
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
mix_of_group <- setNames(c("B", "A"), c(grp_ref, grp_num))
alpha       <- 0.05               # FDR, as in the main DE
lfc         <- 0.5                # H0 |log2FC| <= lfc, as in the main DE
min_count   <- 10                 # filter, as in the main DE
min_samples <- 3
ctrl_lfc    <- 0.25               # |expected log2FC| below this = ratio-1 control sequin
outlier_lfc <- 2                  # |observed - expected log2FC| above this -> sequin_fc_outliers.csv

## identification (as longbench/01_sequin_edger.R)
sequin_seq    <- "chrIS"          # sequin chromosome name in the genome / GTF
det_min_reads <- 3                # detected = >= this many reads ...
det_min_libs  <- 2                # ... in >= this many libraries
det_rule_txt  <- sprintf(">=%d reads in >=%d libraries", det_min_reads, det_min_libs)
lod_p         <- 0.7
conc_breaks   <- c(0, 1, 10, 100, 1000, Inf)
## GTFs / raw counts: edit if your layout differs; a missing file only skips
## that part (novel models, lengths, Bambu-vs-validated detection)
gtf_validated  <- file.path(data_dir, "bambu_validated", "BambuOutput_annotations_validated.gtf")
gtf_raw        <- file.path(data_dir, "bambu", "bambu_novel_transcripts.gtf")
counts_raw_tx  <- file.path(data_dir, "bambu", "BambuOutput_counts_transcript.txt")
sequin_gtf     <- NULL            # optional: full sequin annotation (lengths of zero-read sequins too)

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
stopifnot("Anaquin mixture: non-positive concentration" = all(iso_mix$A > 0 & iso_mix$B > 0))
## sequin length from the Anaquin mixture table when it has one (all 164,
## no GTF needed); otherwise from the GTFs below
anaquin_len <- if ("Length" %in% names(RnaQuinIsoformMixture))
  setNames(as.numeric(RnaQuinIsoformMixture$Length), as.character(RnaQuinIsoformMixture$Name)) else NULL

ss <- read.csv(sheet)

## ------------------------------------------------------------ functions ---
## sequin rows of one Bambu count table, rounded per sample
read_counts <- function(path) {
  x <- read.table(path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  x$GENEID <- NULL
  x <- x[grepl("^R[12]_", rownames(x)), , drop = FALSE]
  m <- round(as.matrix(x))
  colnames(m) <- sub("[.]bam$", "", colnames(m))
  m
}
read_sequins <- function(file) read_counts(file.path(data_dir, "bambu_validated", file))

## rows for `ids`, 0 for anything absent (the validated table only lists sequins with reads)
take <- function(m, ids) {
  out <- matrix(0, length(ids), ncol(m), dimnames = list(ids, colnames(m)))
  hit <- intersect(ids, rownames(m))
  out[hit, ] <- m[hit, ]
  out
}

detect <- function(x, min_reads = det_min_reads, min_libs = det_min_libs)
  rowSums(x >= min_reads) >= min_libs

group_of <- function(samples) {
  g <- gsub("-", "_", ss$group[match(samples, ss$sample)])
  stopifnot(!anyNA(g))
  factor(g, levels = c(grp_ref, grp_num))
}
mix_of <- function(samples) unname(mix_of_group[as.character(group_of(samples))])   # "A" / "B"

## everything the test needs for input row i of `inputs`
##   m_mix      every sequin in the mix at this level (0 if absent)
##   ctrl_all   ratio-1 control sequins (for normalisation)
##   m2 / ctrl  sequins passing the DE filter + their controls
##   exp_lfc    expected log2FC of m2 (H1975 / HCC827 = Mix A / Mix B)
prep_input <- function(i) {
  lv  <- inputs$level[i]
  m   <- read_sequins(inputs$file[i])
  e   <- expected[[lv]]
  no_exp <- setdiff(rownames(m), e$id)
  if (length(no_exp))
    message(sprintf("[%s] %d sequins without an Anaquin concentration (dropped): %s",
                    lv, length(no_exp), paste(head(no_exp, 10), collapse = ", ")))
  m_mix <- take(m, e$id)
  grp   <- group_of(colnames(m_mix))
  mix   <- mix_of(colnames(m_mix))

  exp_all  <- setNames(log2(e$A / e$B), e$id)
  ctrl_all <- abs(exp_all) < ctrl_lfc
  det_all  <- detect(m_mix)
  keep     <- rowSums(m_mix >= min_count) >= min_samples
  m2       <- m_mix[keep, , drop = FALSE]
  message(sprintf("[%s] %d sequins in the mix, %d detected (%s), %d pass the DE filter, %d ratio-1 controls, %d truly DE",
                  lv, nrow(m_mix), sum(det_all), det_rule_txt, nrow(m2), sum(ctrl_all[keep]),
                  sum(abs(exp_all[keep]) > lfc)))
  list(level = lv, tag = lv, grp = grp, mix = mix,
       m2 = m2, exp_lfc = exp_all[keep], ctrl = ctrl_all[keep],
       mix_reads = tibble(level = lv, id = e$id, exp_lfc = unname(exp_all),
                          mean_A = rowMeans(m_mix[, mix == "A", drop = FALSE]),
                          mean_B = rowMeans(m_mix[, mix == "B", drop = FALSE]),
                          detected = unname(det_all)),
       detected = sum(det_all), detected_any = sum(rowSums(m_mix) >= 1),
       truly_de_in_mix = sum(abs(exp_all) > lfc))
}

## r: tibble(id, obs_lfc, disp_lfc, pvalue, padj, mean, ...) from the test
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

r2_slope <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 3) return(c(slope = NA_real_, r2 = NA_real_))
  f <- lm(y[ok] ~ x[ok])
  c(slope = unname(coef(f)[2]), r2 = summary(f)$r.squared)
}

## one row per level
##   sensitivity = TP / truly DE            (wrong-sign calls are not TP)
##   fdr         = (FP + wrong sign) / called
##   precision   = TP / called = 1 - fdr
##   f1          = 2 * precision * sensitivity / (precision + sensitivity)
## `called` is reused below as the NUMBER of calls on purpose (it is a count
## from that line on); the per-sequin logical is only used through `outcome`.
## `per_level`: tibble(level, detected, detected_any, truly_de_in_mix)
summarise_calls <- function(res, per_level) {
  res %>%
    filter(!is.na(pvalue), is.finite(obs_lfc)) %>%
    group_by(tool, level) %>%
    summarise(fit_all = list(r2_slope(obs_lfc, exp_lfc)),
              fit_nz  = list(r2_slope(obs_lfc[!one_mix_zero], exp_lfc[!one_mix_zero])),
              fit_shr = list(r2_slope(disp_lfc, exp_lfc)),
              spearman       = cor(obs_lfc, exp_lfc, method = "spearman"),
              median_abs_err = median(abs(obs_lfc - exp_lfc)),
              n_one_mix_zero = sum(one_mix_zero),
              tested = n(), truly_de = sum(truth), called = sum(called),
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
              rmse_test_lfc = sqrt(mean((obs_lfc - exp_lfc)^2)),
              rmse_display_lfc = sqrt(mean((disp_lfc - exp_lfc)^2, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(slope         = map_dbl(fit_all, "slope"), r2        = map_dbl(fit_all, "r2"),
           slope_no_zero = map_dbl(fit_nz,  "slope"), r2_no_zero = map_dbl(fit_nz,  "r2"),
           slope_shrunk  = map_dbl(fit_shr, "slope"), r2_shrunk = map_dbl(fit_shr, "r2")) %>%
    dplyr::select(-fit_all, -fit_nz, -fit_shr) %>%
    left_join(per_level, by = "level") %>%
    mutate(sensitivity_in_mix = true_pos / pmax(truly_de_in_mix, 1)) %>%
    dplyr::relocate(detection_rule, detected, detected_any, .after = level)
}

## DGEList with library sizes + TMM factors from control sequins
norm_dge <- function(m, grp, ctrl) {
  y <- DGEList(m, group = grp)
  if (sum(ctrl) >= 3) {
    yc <- normLibSizes(DGEList(m[ctrl, , drop = FALSE]))
    y$samples$lib.size     <- yc$samples$lib.size
    y$samples$norm.factors <- yc$samples$norm.factors
  } else {
    message("<3 ratio-1 sequins; TMM on all sequins")
    y <- normLibSizes(y)
  }
  y
}

run_edger <- function(m, grp, mix, ctrl) {
  y      <- norm_dge(m, grp, ctrl)
  design <- model.matrix(~ group, data.frame(group = grp))
  cn  <- paste0("group", grp_num)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt  <- topTags(glmTreat(fit, coef = cn, lfc = lfc), n = Inf, sort.by = "none")$table
  shr <- predFC(estimateDisp(y, design), design, prior.count = 2)[, cn]
  eff <- y$samples$lib.size * y$samples$norm.factors
  ids <- rownames(tt)
  tibble(id = ids, obs_lfc = tt$logFC, disp_lfc = unname(shr[ids]),
         pvalue = tt$PValue, padj = tt$FDR,
         mean = rowMeans(edgeR::cpm(y)[ids, , drop = FALSE]) * mean(eff) / 1e6,
         ## raw reads per library in each mix; 0 in every library of one mix
         ## makes the unshrunk logFC extreme (triangles in Fig. 2c)
         mean_A = rowMeans(m[ids, mix == "A", drop = FALSE]),
         mean_B = rowMeans(m[ids, mix == "B", drop = FALSE]),
         one_mix_zero = mean_A == 0 | mean_B == 0)
}

## ----------------------------------------------------------------- run ----
res_all <- list(); per_level <- list(); mix_reads <- list()
for (i in seq_len(nrow(inputs))) {
  p <- prep_input(i)
  t <- score_calls(run_edger(p$m2, p$grp, p$mix, p$ctrl), p$exp_lfc) %>%
    mutate(tool = tool, level = p$level, .before = 1)
  outdir <- file.path(outroot, p$tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write_csv(t, file.path(outdir, "sequin_de.csv"))
  res_all[[i]]   <- t
  mix_reads[[i]] <- p$mix_reads
  per_level[[i]] <- tibble(level = p$level, detection_rule = det_rule_txt,
                           detected = p$detected, detected_any = p$detected_any,
                           truly_de_in_mix = p$truly_de_in_mix)
}

res_all   <- bind_rows(res_all)
mix_reads <- bind_rows(mix_reads)
summ      <- summarise_calls(res_all, bind_rows(per_level))
write_csv(res_all, file.path(outroot, "sequin_de_results.csv"))
write_csv(summ, file.path(outroot, "summary.csv"))
print(dplyr::select(summ, level, detected, tested, truly_de, called, wrong_sign, fdr,
                    sensitivity, sensitivity_in_mix, auc, slope, r2,
                    n_one_mix_zero, r2_no_zero, spearman), width = Inf)

## ---------------------------------------------- identification on chrIS ---
read_gtf_seq <- function(f, seq) {
  g <- read.table(f, sep = "\t", comment.char = "#", quote = "",
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attr"),
                  colClasses = c("character", "character", "character", "integer", "integer",
                                 "character", "character", "character", "character"))
  g %>%
    filter(feature == "exon", seqname == seq) %>%
    mutate(transcript_id = str_match(attr, 'transcript_id "([^"]+)"')[, 2],
           gene_id       = str_match(attr, 'gene_id "([^"]+)"')[, 2])
}
tx_spans <- function(ex) {
  ex %>%
    group_by(transcript_id, gene_id, strand) %>%
    summarise(exonic_len = sum(end - start + 1L), exons = n(),
              start = min(start), end = max(end), .groups = "drop")
}
prf <- function(tp, fn, fp) {
  p  <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  r  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  f1 <- if (is.na(p) || is.na(r) || p + r == 0) NA_real_ else 2 * p * r / (p + r)
  c(precision = p, recall = r, f1 = f1)
}

## detection (positives = all 164 sequin isoforms in the mix)
m_val <- read_sequins(inputs$file[inputs$level == "transcript"])
sm    <- take(m_val, iso_mix$id)
det   <- detect(sm)
tp    <- sum(det); fn <- sum(!det)
if (fn > 0)
  message(sprintf("sequin isoforms not detected (%s; %d of %d): %s", det_rule_txt, fn,
                  length(det), paste(head(names(det)[!det], 30), collapse = ", ")))
det_raw <- if (file.exists(counts_raw_tx)) sum(detect(take(read_counts(counts_raw_tx), iso_mix$id))) else NA_integer_

## dose-response: per library, >= 1 read vs that library's mix concentration
mix  <- mix_of(colnames(sm))
cmat <- vapply(mix, function(mx) if (mx == "A") iso_mix$A else iso_mix$B, numeric(nrow(iso_mix)))
hit  <- as.vector(sm >= 1)
lc   <- log2(as.vector(cmat))
lfit <- suppressWarnings(glm(hit ~ lc, family = binomial))
l_b0 <- unname(coef(lfit)[1]); l_b1 <- unname(coef(lfit)[2])
lod  <- if (is.finite(l_b1) && l_b1 > 0) 2^((qlogis(lod_p) - l_b0) / l_b1) else NA_real_
by_conc <- tibble(conc = as.vector(cmat), hit = hit) %>%
  mutate(band = cut(conc, conc_breaks, right = FALSE, dig.lab = 5)) %>%
  group_by(band) %>%
  summarise(pairs = n(), frac_libraries_with_reads = mean(hit), .groups = "drop")

## lengths: optional full sequin GTF, else the validated GTF (sequins with reads only)
len_src <- c(sequin_gtf, gtf_validated)
len_src <- len_src[file.exists(len_src)]
val_ex  <- if (file.exists(gtf_validated)) read_gtf_seq(gtf_validated, sequin_seq) else NULL
seq_len_tbl <- if (length(len_src))
  map_dfr(len_src, ~ tx_spans(filter(read_gtf_seq(.x, sequin_seq), !grepl("^BambuTx", transcript_id)))) %>%
    distinct(transcript_id, .keep_all = TRUE) else
  tibble(transcript_id = character(), exonic_len = integer())
if (!length(len_src))
  message("no GTF found: novel models skipped", if (is.null(anaquin_len)) " and no sequin lengths" else "",
          " (edit `gtf_validated`)")

## capture efficiency (as longbench/01): reads per amol/ul, centred on the median detected isoform
rA  <- rowMeans(sm[, mix == "A", drop = FALSE])
rB  <- rowMeans(sm[, mix == "B", drop = FALSE])
cap <- log2((rA / iso_mix$A + rB / iso_mix$B) / 2)
cap[!det | !is.finite(cap)] <- NA
cap <- cap - median(cap, na.rm = TRUE)
len <- if (!is.null(anaquin_len)) unname(anaquin_len[iso_mix$id]) else
  seq_len_tbl$exonic_len[match(iso_mix$id, seq_len_tbl$transcript_id)]
cap_rho <- suppressWarnings(cor(cap, log2(len), method = "spearman", use = "complete.obs"))

per_iso <- tibble(id = iso_mix$id, gene = sub("_[0-9]+$", "", iso_mix$id), len = len,
                  conc_A = iso_mix$A, conc_B = iso_mix$B, exp_lfc = log2(iso_mix$A / iso_mix$B),
                  reads_A = rA, reads_B = rB, libs_with_reads = rowSums(sm >= 1),
                  detected = unname(det), detected_any = rowSums(sm) >= 1, capture_log2 = cap)

## novel models on chrIS, raw Bambu vs validated
novel_template <- tibble(stage = character(), transcript_id = character(), gene_id = character(),
                         strand = character(), exonic_len = integer(), exons = integer(),
                         start = integer(), end = integer())
novel <- list()
for (st in c("raw", "validated")) {
  f <- if (st == "raw") gtf_raw else gtf_validated
  novel[[st]] <- if (file.exists(f))
    mutate(tx_spans(filter(read_gtf_seq(f, sequin_seq), grepl("^BambuTx", transcript_id))), stage = st, .before = 1) else NULL
}
n_nov <- function(st) if (is.null(novel[[st]])) NA_integer_ else nrow(novel[[st]])

ident <- tibble(annotated = nrow(iso_mix), detection_rule = det_rule_txt,
                detected_bambu = det_raw, detected = tp, not_detected = fn,
                detected_any = sum(rowSums(sm) >= 1),
                sequin_reads_per_library = mean(colSums(sm)),
                lod70_amol_ul = lod, lod70_logistic_slope = l_b1,
                capture_length_rho = cap_rho)
for (st in c("raw", "validated")) {
  nv <- n_nov(st)
  s  <- prf(tp, fn, if (is.na(nv)) 0 else nv)
  ident[[paste0("novel_", st)]]     <- nv
  ident[[paste0("precision_", st)]] <- if (is.na(nv)) NA_real_ else unname(s["precision"])
  ident[[paste0("f1_", st)]]        <- if (is.na(nv)) NA_real_ else unname(s["f1"])
  ident[["recall"]]                 <- unname(s["recall"])
}

write_csv(ident,   file.path(outroot, "sequin_identification.csv"))
write_csv(per_iso, file.path(outroot, "sequin_per_isoform.csv"))
write_csv(by_conc, file.path(outroot, "sequin_detection_by_conc.csv"))
## unname(): a NAMED list would be read by bind_rows() as columns, not tables
write_csv(bind_rows(c(list(novel_template), unname(novel))),
          file.path(outroot, "novel_models_on_sequins.csv"))
print(ident, width = Inf)
print(by_conc)

## ---- fold-change outliers and their sibling isoforms (transcript + gene level) ----
## isoform level: reads assigned to the wrong sibling (error depends on sibling
##   abundances, which differ between mixes); gene level: the gene total is a
##   capture-weighted sum of its isoforms (see longbench/01_sequin_edger.R)
tx_res <- res_all %>% filter(level == "transcript") %>% mutate(gene = sub("_[0-9]+$", "", id))
out_genes <- bind_rows(
  tx_res %>% filter(abs(obs_lfc - exp_lfc) > outlier_lfc | outcome %in% c("wrong sign", "false positive")),
  res_all %>% filter(level == "gene",
                     abs(obs_lfc - exp_lfc) > outlier_lfc | outcome %in% c("wrong sign", "false positive")) %>%
    mutate(gene = id)) %>%
  distinct(gene)
gene_res <- res_all %>% filter(level == "gene") %>%
  transmute(gene = id, gene_exp_lfc = exp_lfc, gene_obs_lfc = obs_lfc, gene_outcome = outcome)
outliers <- mix_reads %>%
  filter(level == "transcript") %>%
  mutate(gene = sub("_[0-9]+$", "", id)) %>%
  semi_join(out_genes, by = "gene") %>%
  left_join(dplyr::select(tx_res, id, obs_lfc, padj, outcome), by = "id") %>%
  left_join(gene_res, by = "gene") %>%
  left_join(dplyr::select(per_iso, id, len, capture_log2), by = "id") %>%
  mutate(tested = !is.na(outcome)) %>%
  dplyr::select(gene, id, len, capture_log2, exp_lfc, obs_lfc, padj, outcome, tested,
                mean_A, mean_B, gene_exp_lfc, gene_obs_lfc, gene_outcome) %>%
  arrange(gene, id)
write_csv(outliers, file.path(outroot, "sequin_fc_outliers.csv"))
if (nrow(outliers)) print(outliers, n = Inf, width = Inf)

message("Done: ", outroot)
