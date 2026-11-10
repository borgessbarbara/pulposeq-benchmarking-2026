## ===========================================================================
## 02_sirv_ercc.R - SIRV-Set 4 benchmark per protocol (transcript level)  (Fig. 1b, c)
##
## Expected amounts: Lexogen SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx
##   SIRV E0     69 isoforms, all 60 amol/ul (equimolar)
##   long SIRVs  15 transcripts, 4-12 kb, also 60 amol/ul
##   ERCC        92 transcripts, one mix (SIRV ERCC Mix 1), ~1e6-fold
##               concentration range; every sample received the same mix
##
## SIRV annotation given to pulposeq: Lexogen "correct" (c) annotation =
## exactly the 69 E0 isoforms. The design sheet flags each isoform as part of
## the correct (c, 69), insufficient (i, 44) and/or over-annotated (o, 100)
## annotation; the 31 isoforms that are only in "o" (SIRV104, SIRV110, ...)
## are not in the mix AND not in the GTF used here, so they cannot receive
## reads. The script detects this from the count table and then reports the
## decoy test as not applicable instead of as a passed specificity check.
##
## Per protocol (ONT cDNA, ONT DRS, PacBio):
##   1. SIRV E0 isoform accuracy: log2(observed / expected) per isoform and
##      sample; within 2-fold, lost (< 1/8 expected), unstable across samples
##      (SD of log2 ratio > 1, e.g. EM splits like SIRV701 / SIRV705)
##   2. reads on over-annotation-only SIRV isoforms - only if they are in the
##      annotation (otherwise NA, see above)
##   3. novel models (BambuTx*) on spike-in sequences, at two stages, as in
##      the dong sequin check (sequin_mirror_isoforms.md):
##        raw        Bambu discovery (bambu/bambu_novel_transcripts.gtf)
##        validated  pulposeq output (bambu_validated/BambuOutput_annotations_validated.gtf)
##      With the SIRV loci correctly annotated, every novel SIRV model is a
##      false positive; ERCCs are single-exon, so any novel ERCC model is an
##      artefact too
##   4. SIRV identification: precision / recall / F1, for raw and validated
##        positives = 69 E0 + 15 long SIRVs (all in the mix)
##        TP = positives with >= sirv_min_reads reads (summed over samples)
##        FN = positives below that
##        FP = novel SIRV models at that stage (+ reads on over-annotation-only
##             isoforms, when those are in the annotation)
##   5. length bias: E0 + long SIRVs are all equimolar, so reads per molecule
##      vs transcript length (191 nt - 12 kb) shows the length bias directly
##   6. ERCC: dose response (log2 CPM vs log2 input, linear fit), detection limit
##      (70% detection, logistic fit as in Anaquin) and detection above fixed
##      concentrations. No F1: every annotated ERCC is in the mix (no negative
##      set) and all samples share one mix (no differential truth), so recall
##      would only reflect sequencing depth.
##
## Normalisation: CPM within each family (counts: a full-length read = one
## molecule). SIRV and ERCC input did not change between DRS libraries.
## Inputs : <run>/bambu_validated/BambuOutput_counts_transcript_validated.txt,
##          <run>/bambu{,_validated}/ GTFs, data/SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx
## Outputs: results/spikeins/sirv_ercc/
##   sirv_ercc_summary.csv, length_bias.csv, ercc_per_transcript.csv,
##   novel_models_on_spikeins.csv, <platform>/sirv_e0_per_isoform.csv
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("C:/Users/mbarb/OneDrive/paper/downstream/longbench/00_paths.R")   # platforms, base
plat_labels <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio")
outroot <- sirv_dir
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)

lex_file  <- file.path(data_dir, "SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx")

## Novel models are counted at two stages (same layout as the dong check,
## sequin_mirror_isoforms.md):
##   raw        <protocol>/bambu/bambu_novel_transcripts.gtf          (Bambu discovery)
##   validated  <protocol>/bambu_validated/BambuOutput_annotations_validated.gtf
##              (what pulposeq delivers; novel = transcript_id BambuTx*)
gtf_stage <- list(
  raw       = setNames(file.path(base, platforms, "bambu", "bambu_novel_transcripts.gtf"), platforms),
  validated = setNames(file.path(base, platforms, "bambu_validated",
                                 "BambuOutput_annotations_validated.gtf"), platforms))

sirv_annotation <- "correct (c): 69 E0 isoforms"   # recorded in the summary table
sirv_min_reads  <- 1          # a SIRV isoform counts as detected with >= this many reads
ercc_conc_cuts  <- c(1, 11.71875)                  # amol/ul; detection reported above these

## ----------------------------------------------------- expected amounts ---
## readxl drops empty leading columns, so columns are located relative to
## the ID column (found by its content: SIRV101 / SIRV4001 / ERCC-00002),
## not by fixed position
id_re <- "^(SIRV\\d{3,}|ERCC-\\d+)$"
read_sheet <- function(sheet, conc_off, len_off, flag_off = NULL) {
  d  <- readxl::read_excel(lex_file, sheet = sheet, col_names = FALSE,
                           .name_repair = "minimal")
  ic <- which.max(vapply(d, function(x) sum(grepl(id_re, as.character(x))), numeric(1)))
  out <- tibble(id   = as.character(d[[ic]]),
                conc = suppressWarnings(as.numeric(d[[ic + conc_off]])),
                len  = suppressWarnings(as.numeric(d[[ic + len_off]])))
  if (!is.null(flag_off))                               # annotation flags c / i / o
    out <- mutate(out,
                  in_correct      = !is.na(suppressWarnings(as.numeric(d[[ic + flag_off[1]]]))),
                  in_insufficient = !is.na(suppressWarnings(as.numeric(d[[ic + flag_off[2]]]))),
                  in_over         = !is.na(suppressWarnings(as.numeric(d[[ic + flag_off[3]]]))))
  filter(out, grepl(id_re, id))
}
sirv_all <- read_sheet("SIRV isoform", conc_off = 18, len_off = 21, flag_off = c(2, 3, 4))
e0       <- filter(sirv_all, !is.na(conc))                         # 69 in the mix
over_only <- filter(sirv_all, is.na(conc), in_over, !in_correct)$id # 31, "o" annotation only
long     <- read_sheet("longSIRV", conc_off = 3, len_off = 6)
ercc     <- read_sheet("ERCC",     conc_off = 4, len_off = 7) %>% filter(!is.na(conc))
message(sprintf("Lexogen: %d E0 isoforms, %d over-annotation-only, %d long SIRVs, %d ERCC",
                nrow(e0), length(over_only), nrow(long), nrow(ercc)))
## sanity: counts, E0 + long SIRVs at 60 amol/ul, long SIRVs 4-12 kb, the
## correct annotation is exactly the E0 mix
stopifnot("Lexogen sheet not read as expected" =
            nrow(e0) == 69 && length(over_only) == 31 && nrow(long) == 15 && nrow(ercc) == 92 &&
            all(e0$conc == 60) && all(long$conc == 60) &&
            setequal(e0$id, filter(sirv_all, in_correct)$id) &&
            all(long$len > 3000 & long$len < 13000) &&
            all(e0$len > 100 & e0$len < 5000) && !anyNA(ercc$len))
positives <- c(e0$id, long$id)                                      # 84 SIRV molecules in the mix

## ------------------------------------------------------------- helpers ----
read_tx <- function(plat) {
  f <- file.path(base, plat, "bambu_validated", "BambuOutput_counts_transcript_validated.txt")
  x <- read.table(f, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  x$GENEID <- NULL
  m <- as.matrix(x[grepl("^(SIRV|ERCC)", rownames(x)), , drop = FALSE])
  colnames(m) <- sub("[.]bam$", "", colnames(m))
  m
}

## novel Bambu models (BambuTx*) on SIRV / ERCC sequences in one GTF.
## Transcripts are rebuilt from exon records, so it works whether or not the
## GTF has "transcript" lines.
read_novel_spikein <- function(f, plat, stage) {
  if (!file.exists(f))
    stop(sprintf("%s GTF not found for %s: %s\n  edit `gtf_stage` at the top of the script",
                 stage, plat, f))
  g <- read.table(f, sep = "\t", comment.char = "#", quote = "",
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attr"),
                  colClasses = c("character", "character", "character", "integer", "integer",
                                 "character", "character", "character", "character"))
  g %>%
    filter(feature == "exon", grepl("^(SIRV|ERCC)", seqname)) %>%
    mutate(transcript_id = str_match(attr, 'transcript_id "([^"]+)"')[, 2],
           gene_id       = str_match(attr, 'gene_id "([^"]+)"')[, 2]) %>%
    filter(grepl("^BambuTx", transcript_id)) %>%
    group_by(transcript_id, gene_id, seqname, strand) %>%
    ## exonic_len first: once `start` / `end` are summarised they would mask
    ## the per-exon columns (same pitfall as the old sequin scoring bug)
    summarise(exonic_len = sum(end - start + 1L), exons = n(),
              start = min(start), end = max(end), .groups = "drop") %>%
    transmute(protocol = plat_labels[[plat]], stage,
              family = if_else(grepl("^SIRV", seqname), "SIRV", "ERCC"),
              seqname, start, end, strand, exons, exonic_len, gene_id, transcript_id)
}

prf <- function(tp, fn, fp) {
  p  <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  r  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  f1 <- if (is.na(p) || is.na(r) || p + r == 0) NA_real_ else 2 * p * r / (p + r)
  c(precision = p, recall = r, f1 = f1)
}

## rows for `ids`, 0 for anything absent from the matrix (= not detected)
take <- function(m, ids) {
  out <- matrix(0, length(ids), ncol(m), dimnames = list(ids, colnames(m)))
  hit <- intersect(ids, rownames(m))
  out[hit, ] <- m[hit, ]
  out
}
cpm_mat <- function(m) sweep(m, 2, pmax(colSums(m), 1), "/") * 1e6   # not `cpm`: would mask edgeR::cpm

## ----------------------------------------------------------------- run ----
len_all <- list(); ercc_all <- list(); summ <- list(); novel_all <- list()

for (plat in platforms) {
  prot   <- plat_labels[[plat]]
  outdir <- file.path(outroot, plat)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  m <- read_tx(plat)

  ## ---- 1. SIRV E0 isoform accuracy ----
  e0_m <- take(m, e0$id)
  e0_l <- as_tibble(log2((cpm_mat(e0_m) + 1) / (1e6 / nrow(e0))), rownames = "id") %>%
    pivot_longer(-id, names_to = "sample", values_to = "log2_ratio") %>%
    mutate(protocol = prot, gene = sub("^(SIRV\\d)\\d+$", "\\1", id))
  e0_iso <- e0_l %>%
    group_by(protocol, gene, id) %>%
    summarise(median = median(log2_ratio), sd = sd(log2_ratio), .groups = "drop")
  write_csv(e0_iso, file.path(outdir, "sirv_e0_per_isoform.csv"))

  ## ---- 2. reads on over-annotation-only SIRV isoforms (if annotated) ----
  ## Bambu reports every annotated transcript, even with 0 reads, so an
  ## over-annotation-only ID missing from the table means it was not in the GTF.
  ## (Assumes the validated table keeps zero-count annotated transcripts; with
  ## the correct annotation used here this is moot - SIRV104 etc. are not in
  ## the Lexogen GTF at all.)
  over_in_annot <- intersect(over_only, rownames(m))
  decoy_testable <- length(over_in_annot) > 0
  if (!decoy_testable)
    message(sprintf("[%s] none of the %d over-annotation-only SIRV isoforms is in the annotation: decoy test not applicable",
                    prot, length(over_only)))
  over_reads    <- if (decoy_testable) sum(m[over_in_annot, ]) else NA_real_
  over_with_rds <- if (decoy_testable) sum(rowSums(m[over_in_annot, , drop = FALSE]) >= sirv_min_reads) else NA_integer_
  sirv_tot      <- sum(take(m, c(positives, over_in_annot)))

  ## ---- 3. novel models on spike-in sequences: raw Bambu vs validated ----
  nov <- bind_rows(read_novel_spikein(gtf_stage$raw[[plat]],       plat, "raw"),
                   read_novel_spikein(gtf_stage$validated[[plat]], plat, "validated"))
  novel_all[[plat]] <- nov
  n_nov <- function(fam, st) sum(nov$family == fam & nov$stage == st)
  for (st in c("raw", "validated")) {
    x <- filter(nov, stage == st)
    if (nrow(x))
      message(sprintf("[%s] %s: novel models on spike-ins: %s", prot, st,
                      paste(sprintf("%s %s:%d-%d %s, %d exon(s), %d nt", x$transcript_id, x$seqname,
                                    x$start, x$end, x$strand, x$exons, x$exonic_len), collapse = "; ")))
  }

  ## ---- 4. SIRV identification: precision / recall / F1, per stage ----
  ## TP / FN from the validated counts (annotated SIRVs are not curated away);
  ## FP = novel SIRV models at that stage (+ reads on over-annotation-only
  ## isoforms when those are annotated)
  pos_det <- rowSums(take(m, positives)) >= sirv_min_reads
  tp <- sum(pos_det)
  fn <- sum(!pos_det)
  decoy_fp <- if (decoy_testable) over_with_rds else 0L
  fp_raw   <- n_nov("SIRV", "raw")       + decoy_fp
  fp_val   <- n_nov("SIRV", "validated") + decoy_fp
  s_raw    <- prf(tp, fn, fp_raw)
  s_val    <- prf(tp, fn, fp_val)
  if (fn > 0)
    message(sprintf("[%s] SIRV isoforms in the mix not detected: %s", prot,
                    paste(positives[!pos_det], collapse = ", ")))

  ## ---- 5. length bias: E0 + long SIRVs, all equimolar ----
  eq      <- bind_rows(mutate(e0, set = "SIRV E0"), mutate(long, set = "long SIRV"))
  eq_m    <- take(m, eq$id)
  eq_rel  <- log2((rowMeans(cpm_mat(eq_m)) + 1) / (1e6 / nrow(eq)))
  len_t   <- tibble(protocol = prot, id = eq$id, set = eq$set, len = eq$len,
                    log2_ratio = unname(eq_rel[eq$id]))
  len_all[[plat]] <- len_t
  len_fit <- lm(log2_ratio ~ log2(len), data = len_t)

  ## ---- 6. ERCC ----
  er_m   <- take(m, ercc$id)
  er_cpm <- cpm_mat(er_m)
  det    <- rowMeans(er_m > 0)                              # fraction of samples detecting
  er_t   <- tibble(protocol = prot, id = ercc$id, conc = ercc$conc,
                   mean_count = rowMeans(er_m), detected = det,
                   log2_cpm = log2(rowMeans(er_cpm)))
  ercc_all[[plat]] <- er_t

  ok  <- is.finite(er_t$log2_cpm)
  efit <- if (sum(ok) >= 3) lm(log2_cpm ~ log2(conc), data = er_t[ok, ]) else NULL
  lfit <- suppressWarnings(glm(cbind(round(det * ncol(m)), ncol(m) - round(det * ncol(m))) ~
                                 log2(conc), family = binomial, data = er_t))
  lod  <- 2^((qlogis(0.7) - coef(lfit)[1]) / coef(lfit)[2])  # conc at 70% detection

  ## ERCC detection above fixed concentrations (any sample)
  ercc_above <- map_dfc(ercc_conc_cuts, function(cc) {
    s <- er_t$conc >= cc
    tibble(!!sprintf("ercc_detected_ge_%s", format(cc, trim = TRUE)) := sum(det[s] > 0),
           !!sprintf("ercc_total_ge_%s",    format(cc, trim = TRUE)) := sum(s))
  })

  summ[[plat]] <- bind_cols(tibble(
    protocol = prot,
    sirv_annotation      = sirv_annotation,
    ## SIRV E0 quantification
    e0_detected          = sum(rowSums(e0_m) >= sirv_min_reads),            # same rule as sirv_tp
    e0_median_abs_log2   = median(abs(e0_iso$median)),
    e0_within_2fold      = mean(abs(e0_iso$median) <= 1),
    e0_lost              = sum(e0_iso$median < -3),
    e0_unstable          = sum(e0_iso$sd > 1),
    long_detected        = sum(rowSums(take(m, long$id)) >= sirv_min_reads), # e0 + long = sirv_tp
    ## decoy test (NA when the over-annotation isoforms are not in the GTF)
    decoy_test           = if (decoy_testable) "applicable" else "not applicable (not in annotation)",
    over_only_in_annotation   = length(over_in_annot),
    over_only_with_reads      = over_with_rds,
    over_only_read_pct        = if (decoy_testable) 100 * over_reads / max(sirv_tot, 1) else NA_real_,
    ## assembly: novel models on spike-in sequences, Bambu vs after validation
    novel_sirv_raw       = n_nov("SIRV", "raw"),
    novel_sirv_validated = n_nov("SIRV", "validated"),
    novel_ercc_raw       = n_nov("ERCC", "raw"),
    novel_ercc_validated = n_nov("ERCC", "validated"),
    ## SIRV identification (positives = 69 E0 + 15 long SIRVs)
    sirv_tp              = tp,
    sirv_fn              = fn,
    sirv_fp_raw          = fp_raw,
    sirv_fp_validated    = fp_val,
    sirv_recall          = unname(s_val["recall"]),
    sirv_precision_raw   = unname(s_raw["precision"]),
    sirv_f1_raw          = unname(s_raw["f1"]),
    sirv_precision       = unname(s_val["precision"]),     # validated = what pulposeq delivers
    sirv_f1              = unname(s_val["f1"]),
    ## length bias
    length_slope         = unname(coef(len_fit)[2]),   # log2 ratio per doubling of length
    long_vs_e0_log2      = mean(len_t$log2_ratio[len_t$set == "long SIRV"]) -
                           mean(len_t$log2_ratio[len_t$set == "SIRV E0"]),
    ## ERCC
    ercc_detected        = sum(det > 0),
    ercc_total           = nrow(ercc),
    ercc_slope           = if (is.null(efit)) NA_real_ else unname(coef(efit)[2]),
    ercc_r2              = if (is.null(efit)) NA_real_ else summary(efit)$r.squared,
    ercc_lod70_amol_ul   = unname(lod)), ercc_above)
}

## ------------------------------------------------------------- outputs ----
summ      <- bind_rows(summ)
len_all   <- bind_rows(len_all)
ercc_all  <- bind_rows(ercc_all)
novel_all <- bind_rows(novel_all)
write_csv(summ,      file.path(outroot, "sirv_ercc_summary.csv"))
write_csv(len_all,   file.path(outroot, "length_bias.csv"))
write_csv(ercc_all,  file.path(outroot, "ercc_per_transcript.csv"))
write_csv(novel_all, file.path(outroot, "novel_models_on_spikeins.csv"))
print(summ, width = Inf)

message("Done: ", outroot)
