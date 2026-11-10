## ===========================================================================
## 01_sequin_edger.R - LongBench sequin ground truth with edgeR      (Fig. 1d)
##   oarfish version of before-oarfish-integration/longbench/01_sequin_edger.R:
##   same test, truth and scoring; the counts are oarfish's, with the bootstrap
##   replicates used as edgeR intends (downstream/helpers/oarfish_counts.R):
##     transcript level  edgeR::catchOarfish()     -> counts / RTA overdispersion
##     gene level        edgeR::catchOarfishGene() -> counts / gene overdispersion
##                       (primary), and undivided summed counts (sensitivity check)
##   The expression filter is applied to the undivided read counts; edgeR then
##   tests the divided counts (not rounded). Needs edgeR >= 4.99.12.
##
## Test (as dong/02_sequin_edger.R and dong/01_de_edger.R):
##   edgeR v4 glmQLFit(robust = TRUE) -> glmTreat(lfc = 0.5), FDR 0.05;
##   display log2FC = predFC(prior.count = 2); filter >= min_count reads in
##   >= min_samples samples; library sizes AND TMM factors from the ratio-1
##   control sequins (|expected log2FC| < 0.25), so the normalisation never
##   sees the sequins that are meant to change.
##   Contrast: Mix B (SCLC-A) vs Mix A (expected log2FC = log2(B / A)).
##
## Truth and scoring (identical to dong/02_sequin_edger.R)
##   truly DE  = |expected log2FC| > lfc
##   called    = padj < alpha
##   true positive = called, truly DE, correct sign; wrong-sign calls are
##   false discoveries; precision = TP / called = 1 - FDR; F1 from precision
##   and sensitivity. AUC ranks -log10(p) against truly DE.
##
## Identification (no DE tool involved): every sequin isoform in the mix
## (Anaquin, 164) is a positive; detected = >= 1 read (rounded per sample,
## summed) in the curated oarfish table; also reported before curation, for
## every model oarfish quantified and for Bambu's own table. FP = novel BambuTx
## models on chrIS, counted in Bambu's GTF (raw) and after pulposeq curation.
##
## Inputs : <run>/transcriptome/oarfish/ (catchOarfish), <run>/transcriptome/curated/
##          counts + annotation, <run>/transcriptome/bambu/, data/samplesheet_*.csv, Anaquin
## Outputs: spike-ins/results/sequins/
##   edger/<platform>/<level>/<count_mode>/sequin_de.csv
##   edger/sequin_de_results.csv           every sequin, every protocol / level / count mode
##   edger/fold_change_summary_edger.csv   one row per protocol x level x count mode
##                                         (count_mode "divided" = primary analysis)
##   sequin_identification.csv             identification on chrIS
##   novel_models_on_sequins.csv           novel models on chrIS, raw vs validated
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(Anaquin)
})

source("C:/Users/mbarb/OneDrive/paper/downstream/longbench/00_paths.R")
source(file.path(shared_helper_dir, "oarfish_counts.R"))   # read_oarfish()
plat_labels <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio")
seqroot <- seq_dir
outroot <- file.path(seqroot, "edger")
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------- settings ---
## platforms: 00_paths.R
sheets    <- c(ont_cdna = "samplesheet_ont.csv",
               ont_drs  = "samplesheet_ont_drs.csv",
               pb       = "samplesheet_pb.csv")
levels_to_run <- c("gene", "transcript")
alpha       <- 0.05   # LongBench DE settings
lfc         <- 0.5    # LongBench DE settings
min_count   <- 10     # LongBench DE settings
min_samples <- 2      # LongBench DE settings
ctrl_lfc    <- 0.25   # |expected log2FC| below this = ratio-1 control sequin
## counts edgeR receives: "divided" = counts / RTA overdispersion (primary);
## "undivided" = summed oarfish counts, gene level only, as a sensitivity check
count_modes <- list(transcript = "divided", gene = c("divided", "undivided"))

## identification (no DE tool involved: reads per sequin + novel models on chrIS)
sequin_seq       <- "chrIS"   # sequin chromosome name in the genome / GTF
sequin_min_reads <- 1         # an annotated sequin isoform counts as detected with >= this many reads
gtf_stage <- list(
  raw       = setNames(vapply(platforms, run_file, "", key = "bambu_novel_gtf"), platforms),
  validated = setNames(vapply(platforms, run_file, "", key = "validated_gtf"), platforms))
tool        <- "edgeR"

## ---------------------------------------------- expected concentrations ---
data("RnaQuinIsoformMixture", package = "Anaquin")
iso_mix <- as_tibble(RnaQuinIsoformMixture) %>%
  transmute(id = as.character(Name), A = MixA, B = MixB)
gene_mix <- iso_mix %>%
  mutate(gene = sub("_[0-9]+$", "", id)) %>%
  group_by(gene) %>%
  summarise(A = sum(A), B = sum(B), .groups = "drop") %>%
  transmute(id = gene, A, B)
expected <- list(transcript = iso_mix, gene = gene_mix)

## -------------------------------------------------------------- inputs ----
## sequin rows of the curated oarfish counts: undivided reads, divided counts
## and the RTA overdispersion (catchOarfish / catchOarfishGene)
read_sequins <- function(plat, lv) {
  x <- read_oarfish(run_dir[[plat]], lv, cache_dir = oarfish_cache)
  s <- grepl("^R[12]_", rownames(x$counts))
  list(reads   = x$counts[s, , drop = FALSE],
       divided = x$divided[s, , drop = FALSE],
       od      = x$overdispersion[s])
}

sample_mix <- function(samples, sheet) {
  ss  <- read.csv(file.path(data_dir, sheet))
  grp <- gsub("-", "_", ss$group[match(samples, ss$sample)])
  stopifnot(!anyNA(grp))
  setNames(factor(ifelse(grp == "SCLC_A", "B", "A"), levels = c("A", "B")), samples)
}

auc_fn <- function(score, truth) {
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---------------------------------------------------------- scoring -------
## same definitions as dong/02_sequin_edger.R (do not edit one without the other)
score_calls <- function(t) {
  t %>%
    mutate(truth   = abs(exp_lfc) > lfc,
           called  = !is.na(padj) & padj < alpha,
           outcome = case_when(called & truth & sign(obs_lfc) == sign(exp_lfc) ~ "true positive",
                               called & truth                                  ~ "wrong sign",
                               called & !truth                                 ~ "false positive",
                               !called & truth                                 ~ "false negative",
                               TRUE                                            ~ "true negative"))
}

summarise_calls <- function(t) {
  fit <- lm(obs_lfc ~ exp_lfc, data = t)
  t %>%
    summarise(tested      = n(),
              truly_de    = sum(truth),
              n_called    = sum(called),
              true_pos    = sum(outcome == "true positive"),
              false_pos   = sum(outcome == "false positive"),
              wrong_sign  = sum(outcome == "wrong sign"),
              false_neg   = sum(outcome == "false negative"),
              true_neg    = sum(outcome == "true negative"),
              sensitivity = true_pos / max(truly_de, 1),
              fdr         = if (n_called > 0) (false_pos + wrong_sign) / n_called else NA_real_,
              precision   = if (n_called > 0) true_pos / n_called else NA_real_,
              f1          = if (is.na(precision) || precision + sensitivity == 0) NA_real_ else
                              2 * precision * sensitivity / (precision + sensitivity),
              auc         = auc_fn(-log10(pvalue), truth),
              slope       = unname(coef(fit)[2]),
              intercept   = unname(coef(fit)[1]),
              r2          = summary(fit)$r.squared,
              rmse        = sqrt(mean((obs_lfc - exp_lfc)^2)))
}

## ------------------------------------------------------------- edgeR ------
## cpm() is always called as edgeR::cpm(): a plain cpm() defined by another
## script in the same R session would otherwise mask it
## DGEList with library sizes + TMM factors from the control sequins
norm_dge <- function(m, mix, ctrl, lab) {
  y <- DGEList(m, group = mix)
  if (sum(ctrl) >= 3) {
    yc <- normLibSizes(DGEList(m[ctrl, , drop = FALSE]))
    y$samples$lib.size     <- yc$samples$lib.size
    y$samples$norm.factors <- yc$samples$norm.factors
  } else {
    message(sprintf("[%s] <3 ratio-1 sequins; TMM on all sequins", lab))
    y <- normLibSizes(y)
  }
  y
}

run_edger <- function(m, mix, ctrl, lab) {
  y      <- norm_dge(m, mix, ctrl, lab)
  design <- model.matrix(~ mix)
  cn     <- "mixB"                                         # Mix B vs Mix A
  fit    <- glmQLFit(y, design, robust = TRUE)
  tt     <- topTags(glmTreat(fit, coef = cn, lfc = lfc), n = Inf, sort.by = "none")$table
  shr    <- predFC(estimateDisp(y, design), design, prior.count = 2)[, cn]
  eff    <- y$samples$lib.size * y$samples$norm.factors
  tibble(id = rownames(tt), obs_lfc = tt$logFC, disp_lfc = unname(shr[rownames(tt)]),
         pvalue = tt$PValue, padj = tt$FDR,
         mean = rowMeans(edgeR::cpm(y)[rownames(tt), , drop = FALSE]) * mean(eff) / 1e6)
}

## ----------------------------------------------------------------- run ----
res_all <- list(); fc_all <- list()

for (plat in platforms) {
  for (lv in levels_to_run) {
    x   <- read_sequins(plat, lv)
    r   <- x$reads
    mix <- sample_mix(colnames(r), sheets[[plat]])
    e   <- expected[[lv]]
    ids <- intersect(rownames(r), e$id)
    e   <- e[match(ids, e$id), ]

    exp_all  <- setNames(log2(e$B / e$A), ids)
    ctrl_all <- abs(exp_all) < ctrl_lfc
    detected <- sum(rowSums(r[ids, , drop = FALSE]) > 0)
    ## filter on reads (undivided), so every count mode tests the same sequins
    keep     <- rowSums(r[ids, , drop = FALSE] >= min_count) >= min_samples

    for (mode in count_modes[[lv]]) {
      lab    <- paste(plat_labels[[plat]], lv, mode, "(edgeR)")
      outdir <- file.path(outroot, plat, lv, mode)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      m  <- if (mode == "divided") x$divided else x$reads
      m2 <- m[ids, , drop = FALSE][keep, , drop = FALSE]
      message(sprintf("[%s] %d sequins, %d pass the filter, %d ratio-1 controls, %d truly DE",
                      lab, length(ids), nrow(m2), sum(ctrl_all[keep]), sum(abs(exp_all[keep]) > lfc)))

      ## ---- DE: Mix B vs Mix A ----
      t <- run_edger(m2, mix, ctrl_all[keep], lab) %>%
        mutate(exp_lfc = unname(exp_all[id]), ratio = round(exp_lfc),
               overdispersion = unname(x$od[id])) %>%
        filter(!is.na(pvalue), is.finite(obs_lfc)) %>%
        score_calls() %>%
        mutate(tool = tool, platform = plat, level = lv, count_mode = mode, .before = 1)
      write_csv(t, file.path(outdir, "sequin_de.csv"))
      res_all[[length(res_all) + 1]] <- t

      fc_all[[length(fc_all) + 1]] <- summarise_calls(t) %>%
        mutate(tool = tool, platform = plat, level = lv, count_mode = mode,
               detected = detected,
               median_overdispersion = median(x$od[ids][keep]), .before = 1)
    }
  }
}

res_all <- bind_rows(res_all)
fc_all  <- bind_rows(fc_all)
write_csv(res_all, file.path(outroot, "sequin_de_results.csv"))
write_csv(fc_all,  file.path(outroot, "fold_change_summary_edger.csv"))
print(dplyr::select(fc_all, platform, level, count_mode, tested, truly_de, n_called, wrong_sign,
                    sensitivity, precision, f1, auc, slope, r2), n = Inf)

## ---------------------------------------------- identification on chrIS ---
## exon records of one GTF on the sequin chromosome
read_gtf_seq <- function(f, seq) {
  if (!file.exists(f))
    stop("GTF not found: ", f, "\n  edit `gtf_stage` at the top of the script")
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

## one row per transcript (exon_len first: summarised start / end would mask it)
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

ident_all <- list(); novel_seq_all <- list()
for (plat in platforms) {
  prot <- plat_labels[[plat]]

  ## positives = ALL sequin isoforms in the mix (Anaquin, 164); one missing
  ## from a count table had no reads there. Reads are rounded per sample,
  ## then summed - the same rule as `detected` in fold_change_summary.csv.
  sequin_reads <- function(f) {
    x <- read.table(f, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    x$GENEID <- NULL
    obs <- rowSums(round(as.matrix(x[grepl("^R[12]_", rownames(x)), , drop = FALSE])))
    r <- setNames(numeric(nrow(iso_mix)), iso_mix$id)
    hit <- intersect(names(obs), iso_mix$id)
    r[hit] <- obs[hit]
    list(reads = r, in_table = length(hit))
  }
  ## curated oarfish table (the counts the DE analysis uses, undivided)
  val <- sequin_reads(run_file(plat, "oarfish_tx"))
  seq_reads <- val$reads
  det <- seq_reads >= sequin_min_reads
  tp  <- sum(det); fn <- sum(!det)
  if (fn > 0)
    message(sprintf("[%s] sequin isoforms not detected (%d of %d): %s", prot, fn, length(det),
                    paste(head(names(seq_reads)[!det], 20), collapse = ", ")))

  ## before curation: every model oarfish quantified, and Bambu's own table
  all_oar <- read_oarfish(run_dir[[plat]], "transcript", curated = FALSE, cache_dir = oarfish_cache)$counts
  oar_obs <- rowSums(round(all_oar[grepl("^R[12]_", rownames(all_oar)), , drop = FALSE]))
  det_oar_all <- sum(oar_obs[intersect(names(oar_obs), iso_mix$id)] >= sequin_min_reads)
  f_raw <- run_file(plat, "bambu_counts_tx")
  raw   <- if (file.exists(f_raw)) sequin_reads(f_raw) else NULL
  det_raw <- if (is.null(raw)) NA_integer_ else sum(raw$reads >= sequin_min_reads)
  if (det_oar_all > tp)
    message(sprintf("[%s] %d sequin isoforms have oarfish reads but are not in the curated set: %s",
                    prot, det_oar_all - tp,
                    paste(head(setdiff(names(oar_obs)[oar_obs >= sequin_min_reads],
                                       names(seq_reads)[det]), 20), collapse = ", ")))

  ## annotated sequins on chrIS (from the validated GTF) - also checks the name
  val_ex <- read_gtf_seq(gtf_stage$validated[[plat]], sequin_seq)
  ann    <- tx_spans(filter(val_ex, !grepl("^BambuTx", transcript_id)))
  if (!nrow(ann))
    stop(sprintf("[%s] no annotated transcripts on '%s' in the validated GTF - set `sequin_seq`",
                 prot, sequin_seq))

  row <- tibble(platform = plat, protocol = prot,
                annotated = length(seq_reads),                 # 164 (Anaquin)
                in_validated_table = val$in_table,             # curated oarfish table
                in_bambu_table = if (is.null(raw)) NA_integer_ else raw$in_table,
                detected_oarfish_all_models = det_oar_all,     # >= 1 oarfish read before curation
                detected_bambu = det_raw,                      # >= 1 Bambu read (Bambu's raw table)
                detected = tp, not_detected = fn)              # curated oarfish table
  for (st in c("raw", "validated")) {
    ex  <- if (st == "validated") val_ex else read_gtf_seq(gtf_stage$raw[[plat]], sequin_seq)
    nov <- tx_spans(filter(ex, grepl("^BambuTx", transcript_id)))
    if (nrow(nov)) {
      ## overlapping annotated sequins, same vs opposite strand
      ov <- function(i, same) {
        hit <- ann$start <= nov$end[i] & ann$end >= nov$start[i] &
               (if (same) ann$strand == nov$strand[i] else ann$strand != nov$strand[i])
        paste(ann$transcript_id[hit], collapse = ";")
      }
      nov <- nov %>%
        mutate(same_strand_overlap     = map_chr(seq_len(n()), ov, same = TRUE),
               opposite_strand_overlap = map_chr(seq_len(n()), ov, same = FALSE),
               candidate_mirror        = opposite_strand_overlap != "" & same_strand_overlap == "")
      message(sprintf("[%s] %s: %d novel model(s) on %s (%d candidate mirror(s))", prot, st,
                      nrow(nov), sequin_seq, sum(nov$candidate_mirror)))
      novel_seq_all[[length(novel_seq_all) + 1]] <-
        mutate(nov, platform = plat, protocol = prot, stage = st, .before = 1)
    }
    s <- prf(tp, fn, nrow(nov))
    row[[paste0("novel_", st)]]     <- nrow(nov)
    row[[paste0("precision_", st)]] <- unname(s["precision"])
    row[[paste0("f1_", st)]]        <- unname(s["f1"])
    row[["recall"]]                 <- unname(s["recall"])
  }
  ident_all[[plat]] <- row
}
ident_all     <- bind_rows(ident_all)
novel_seq_all <- bind_rows(novel_seq_all)
write_csv(ident_all,     file.path(seqroot, "sequin_identification.csv"))
write_csv(novel_seq_all, file.path(seqroot, "novel_models_on_sequins.csv"))   # empty if none
print(ident_all, width = Inf)

message("Done: ", outroot)
