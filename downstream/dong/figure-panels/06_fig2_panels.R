## ===========================================================================
## 06_fig2_panels.R - Fig. 2 a-g: the Dong panels drawn with ggplot2
## H1975 (TKI-resistant) vs HCC827 (TKI-sensitive), ONT PCR-cDNA, Bambu full
## counts, 3 independent cultures per line (Dong et al. 2023 Nat Methods).
##
## Panel  object  content
##   a     pA      read-level support across biotypes
##   b     pB      novel transcripts by biotype and GffCompare class code
##   c     pID     sequin identification, vertical bars with values on top: recall
##                 (detected / expected), precision, F1 (false positive = novel model
##                 on a sequin; Bambu value before validation in brackets); as Fig. 1b
##   d     pC      sequin ground truth, isoform level: observed vs expected log2FC,
##                 DE precision / sensitivity / F1 top-left, sequins tested +
##                 R2 / slope bottom-right
##   e     pD      sample-sample Spearman correlation, isoform level (edgeR logCPM)
##   f     pE      volcano, gene level, Up / Down totals; marker genes labelled
##                 (only those tested AND DE)
##   g     pF      DE isoforms of the marker genes (Bambu gene assignment), rows
##                 grouped by gene and ranked by log2FC; z-score of logCPM per
##                 isoform; last column = isoform DE call
##   h, j-l        DE isoform profiles: 04_isoformic.R, section 7e
##                 (isoformic/figures/fig2_panels/)
##   i             EGFR locus: 04_isoformic.R (isoformic/figures/locus_*/EGFR.*)
##   m             cis-acting pairs circos: 03_cis_circos.R (paper_cis_circos.*)
##   n             cis-pair heatmap: 03_cis_circos.R (paper_cis_heatmap.*)
## (v2: new panel c; the old c-m are now d-n. Sequin identification numbers
## live only in panel c.)
##   ED Fig. 1     every event + cis-pair profile: 04_isoformic.R
##
## Layout = the final Fig. 2 (manuscript/Figure2.pdf, an A3 page, 297 x 420 mm):
## every panel is saved at the size it has on that page, with the font sizes,
## legend positions, key sizes and panel-border widths measured there
## (titles / subtitles 10 pt, legend and tick text 8 pt; see each panel), and the
## orientation as printed (vertical heatmap in f). No bold text. Geometry
## (points, lines, tiles) as drawn by make_fig_dong2.R (theme base size 6).
##
## Needs (run first): 01_de_edger.R, 02_sequin_edger.R, 05_read_support.R
## Outputs
##   figures/fig2_panels/fig2{a..g}_*.{pdf,svg,png}   (svg = editable text)
##   results/figure_tables/supp_novel_biotype_classcode_dong.csv   (Fig. 2b numbers)
##   results/figure_tables/sequin_identification_dong.csv           (Fig. 2c numbers)
## Needs: tidyverse, ggrepel, RColorBrewer, svglite, showtext, systemfonts
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

source("C:/Users/mbarb/OneDrive/pulposeq-oarfish/pulposeq-benchmarking-2026/downstream/dong/00_paths.R")
source(file.path(helper_dir, "plot_theme.R"))       # font_family, base_theme, group_cols, npg_cols
source(file.path(helper_dir, "figure_helpers.R"))   # save_panel, keep_level, read_support_plot
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------- settings ---
alpha    <- 0.05                  # as in 01_de_edger.R / 02_sequin_edger.R
lfc      <- 0.5
line_lab <- c(TKI_sensitive = "HCC827", TKI_resistant = "H1975")
lv_lab   <- c(gene = "gene", transcript = "isoform")
markers  <- c(
  # EGFR driver locus
  "EGFR",
  # LUAD lineage identity
  "NAPSA", "CEACAM5", "KRT7", "NKX2-1", "SFTPB", "MUC1", "SLC2A1", "TOP2A", "TPX2",
  "MALAT1", "CYTOR", "MIR4435-2HG", "DARS1-AS1", "HAGLR",
  # OXPHOS isoform switching
  "NDUFA12", "COX6B1", "NDUFA6",
  # hypoxia (VHL-HIF axis)
  "ELOC", "MIR210HG",
  # lncRNA regulation
  "UBE2I", "UBE2I-AS1", "SNHG17",
  # translation / survival
  "TPT1",
  # cis-pair partner (with BambuTx946_k)
  "LSP1P5"
)                                 # panels f (labels) and g (heatmap)
hm_only_de <- TRUE                # panel g: TRUE = DE isoforms only; FALSE = every tested isoform
max_iso_per_gene <- Inf           # panel g: top isoforms per gene by padj (Inf = all)
z_lim      <- 2                   # panel g: colour scale limits (row z-score)
y_cap      <- 9.5                 # panel d: |observed log2FC| drawn at most this far
cor_level  <- "transcript"        # panel e: "gene" or "transcript" (isoform)
## panel g: sample columns, left to right, as in the final figure (H1975 first;
## the DE call is the last column); other names = HCC827 then H1975 as in the data
hm_samples <- c("H1975_2", "H1975_1", "H1975_3", "HCC827_2", "HCC827_1", "HCC827_3")

dir_cols <- c(Down = "#1170ff", Up = "darkred")
## DE heatmap colours (as in the DE scripts): Moreland cool-warm
hm_cols  <- colorRampPalette(c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2",
                               "#F2CBB7", "#EE8468", "#B40426"), space = "Lab")(100)
cor_cols <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(50)
cc_cols  <- c("junction match (j)"      = npg_cols[["red"]],          # as Fig. 1 / the old Fig. 1e
              "contains reference (k)"  = npg_cols[["salmon"]],
              "retained intron (m, n)"  = npg_cols[["brown"]],
              "sense overlap (i, o, y)" = npg_cols[["beige"]],
              "antisense (i, x)"        = npg_cols[["greyblue"]],
              "intergenic (u)"          = npg_cols[["mint"]])

## ------------------------------------------------ figure-wide styling -----
## base theme: line widths / spacing of the original panels (base size 6);
## every text size is set explicitly from the final figure
gg35 <- utils::packageVersion("ggplot2") >= "3.5.0"
theme_set(theme_bw(base_family = font_family, base_size = 6) + base_theme +
  theme(text            = element_text(face = "plain", lineheight = 1.2),   # no bold; line spacing as in the figure
        plot.title      = element_text(size = 10, face = "plain"),
        plot.subtitle   = element_text(size = 10, colour = "grey30", lineheight = 1),
        plot.title.position = "plot",                            # titles flush with the plot's left edge
        strip.text      = element_text(size = 8, face = "plain", margin = margin(1.5, 1.5, 1.5, 1.5)),
        axis.text       = element_text(size = 8),
        axis.title      = element_text(size = 10, face = "plain"),
        legend.text     = element_text(size = 8),
        legend.title    = element_text(size = 8, face = "plain"),
        legend.key.size = unit(4, "mm"),
        legend.position = "right",
        legend.margin   = margin(0, 0, 0, 0),
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        plot.margin     = margin(0.5, 0.5, 0.5, 0.5, "mm")))
update_geom_defaults("text", list(family = font_family, size = 8 / .pt))
pt2lw  <- function(pt) pt / (.pt * 0.75)                          # line width in pt -> ggplot linewidth
border <- function(pt) theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = pt2lw(pt)))
## legend beside the panel: vertical position (0 = bottom, 1 = top), as placed in the figure
leg_right <- function(y) theme(legend.position = "right", legend.justification = c(0, y))

## small boxed labels with a leader line (7.5 pt, as in the figure)
lab_repel <- function(data, mapping, ...) {
  geom_label_repel(data = data, mapping = mapping, size = 7.5 / .pt, family = font_family,
                   colour = "black", fill = "white", label.size = 0.15,
                   label.r = unit(0.05, "lines"), label.padding = unit(0.08, "lines"),
                   box.padding = 0.3, point.padding = 0.15, min.segment.length = 0,
                   max.overlaps = Inf, segment.colour = "black", segment.size = 0.25,
                   show.legend = FALSE, seed = 1, ...)
}

## ================================================================ a =====
## read-level support (05_read_support.R; same design as Fig. 1f)
rs <- read_csv(file.path(rs_dir, "read_support_all_runs_long.csv"), show_col_types = FALSE)
pA <- read_support_plot(rs, label_fs = 6, n_sep = "\n") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 8, colour = "grey20"),
        axis.ticks.y = element_line()) +
  border(0.5) + leg_right(0.69)

## ================================================================ b =====
## novel transcripts by biotype x class code (validated pulposeq metadata;
## the class-code grouping of Fig. 1 / the old Fig. 1e)
cc_group <- function(cls) case_when(
  grepl("\\(j\\)", cls)                   ~ "junction match (j)",
  grepl("\\(k\\)", cls)                   ~ "contains reference (k)",
  grepl("retained intron", cls)           ~ "retained intron (m, n)",
  grepl("^antisense", cls)                ~ "antisense (i, x)",
  grepl("\\(u\\)", cls)                   ~ "intergenic (u)",
  TRUE                                    ~ "sense overlap (i, o, y)")
bt_lab <- c(novel_lncRNA = "lncRNA", novel_protein_coding = "protein-coding",
            novel_non_coding = "non-coding")

dn <- read_csv(run_path("novel_meta"),
               show_col_types = FALSE)
dnov <- dn %>%
  transmute(biotype = factor(pulposeq_transcript_biotype, levels = names(bt_lab)),
            cls     = factor(cc_group(classification), levels = names(cc_cols)))
if (anyNA(dnov$biotype))
  message("  ! ", sum(is.na(dnov$biotype)), " novel transcripts with an unexpected biotype: ",
          paste(unique(dn$pulposeq_transcript_biotype[is.na(dnov$biotype)]), collapse = ", "))

dnov %>%
  dplyr::count(biotype, cls, .drop = FALSE) %>%
  mutate(`% of novel transcripts` = 100 * n / sum(n)) %>%
  group_by(biotype) %>%
  mutate(`% within biotype` = 100 * n / na_if(sum(n), 0)) %>%   # NA when a biotype has no novel transcripts
  ungroup() %>%
  transmute(Dataset = "Dong", Protocol = "Dong ONT cDNA", Biotype = bt_lab[as.character(biotype)],
            `Class code (GffCompare)` = cls, `Novel transcripts` = n,
            `% of novel transcripts`, `% within biotype`) %>%
  write_csv(file.path(tab_dir, "supp_novel_biotype_classcode_dong.csv"))

b_tot <- dplyr::count(dnov, biotype)
b_lab <- setNames(sprintf("%s (%s)", bt_lab[as.character(b_tot$biotype)], scales::comma(b_tot$n)),
                  as.character(b_tot$biotype))
pB <- ggplot(dplyr::count(dnov, biotype, cls), aes(biotype, n, fill = cls)) +
  geom_col(width = 0.75) +
  geom_text(data = b_tot, aes(x = biotype, y = n, label = scales::comma(n)),
            inherit.aes = FALSE, vjust = -0.4, size = 7 / .pt) +
  scale_fill_manual(values = cc_cols, drop = FALSE) +
  scale_x_discrete(labels = b_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = NULL, fill = "GffCompare class code", title = "Novel transcripts") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  border(0.5) + leg_right(1)

## ---- spike-in identification panel (Fig. 1b; the same code draws Fig. 2c) ----
## One row per spike-in set x protocol: expected (isoforms in the mix), detected
## (TP), fp_raw / fp_val = novel models on that set's sequences, Bambu output /
## after pulposeq validation (every spike-in is annotated, so any novel model on
## one is a false positive). Recall = TP / expected; precision = TP / (TP + FP).
## (v3: vertical bars with the values on top, instead of lollipops)
id_metrics <- function(id) {
  id %>% mutate(recall        = detected / expected,
                precision     = detected / (detected + fp_val),
                f1            = 2 * precision * recall / (precision + recall),
                precision_raw = detected / (detected + fp_raw),
                f1_raw        = 2 * precision_raw * recall / (precision_raw + recall))
}
## Vertical bars, one facet per metric (recall, precision, F1); x = spike-in
## set, one bar per protocol, the value printed on top of each bar (recall as
## detected / expected). Bars and values are after pulposeq validation; where
## the Bambu output (before validation) differs, its value follows in brackets.
id_bars <- function(id, set_levels, prot_levels, cols, label_pt = 6, x_angle = 35) {
  m_lev <- c("Recall", "Precision", "F1")
  m_lab <- c(Recall = "Recall (detected / expected)", Precision = "Precision", F1 = "F1")
  d <- id_metrics(id) %>%
    mutate(set = factor(set, levels = set_levels), protocol = factor(protocol, levels = prot_levels))
  with_raw <- function(v, raw) ifelse(abs(raw - v) >= 0.005, sprintf("%.2f (%.2f)", v, raw), sprintf("%.2f", v))
  long <- bind_rows(
    d %>% transmute(set, protocol, metric = "Recall", value = recall,
                    label = sprintf("%d/%d", as.integer(detected), as.integer(expected))),
    d %>% transmute(set, protocol, metric = "Precision", value = precision,
                    label = with_raw(precision, precision_raw)),
    d %>% transmute(set, protocol, metric = "F1", value = f1, label = with_raw(f1, f1_raw))) %>%
    mutate(metric = factor(metric, levels = m_lev))
  dodge <- position_dodge(width = 0.85)
  ggplot(long, aes(set, value, fill = protocol)) +
    geom_col(position = dodge, width = 0.8) +
    geom_text(aes(label = label, group = protocol), position = dodge, angle = 90, hjust = -0.12,
              vjust = 0.5, size = label_pt / .pt, family = font_family, colour = "grey15") +
    facet_grid(. ~ metric, labeller = labeller(metric = m_lab)) +
    ## headroom above 1 for the labels; the axis stops at 1
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.5", "0.75", "1"), expand = c(0, 0)) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(),
          axis.text.x = if (x_angle == 0) element_text() else element_text(angle = x_angle, hjust = 1, vjust = 1))
}

## ================================================================ c =====
## sequin identification (02_sequin_edger.R, sequin_identification.csv): of the
## 164 sequin isoforms in the mix, detected = >= 3 reads in >= 2 libraries;
## false positive = a novel model on a sequin (Bambu output vs after pulposeq
## validation). Same panel design as Fig. 1b.
sq_id <- read_csv(file.path(seq_dir, "sequin_identification.csv"), show_col_types = FALSE)
if ("count_mode" %in% names(sq_id)) sq_id <- filter(sq_id, count_mode == "divided")
idf <- tibble(set = "Sequins", protocol = "ONT cDNA",
              detected = sq_id$detected[1], expected = sq_id$annotated[1],
              fp_raw = sq_id$novel_raw[1], fp_val = sq_id$novel_validated[1])
write_csv(id_metrics(idf), file.path(tab_dir, "sequin_identification_dong.csv"))
pID <- id_bars(idf, "Sequins", "ONT cDNA", c("ONT cDNA" = "grey30"), label_pt = 7, x_angle = 0) +
  guides(fill = "none") +                          # one protocol: no colour key
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Sequin identification",
       subtitle = paste0("Detected: >= 3 reads in >= 2 libraries\n",
                         "False positive: a novel model on a sequin\n",
                         "Precision, F1: after validation (before, if different)")) +
  border(0.5)

## ================================================================ d =====
## sequin DE ground truth, isoform level (02_sequin_edger.R: same test as the DE)
out_lev  <- c("true positive", "true negative", "false positive", "false negative", "wrong sign")
out_cols <- c("true positive" = "grey15", "true negative" = "grey75", "false positive" = "#D55E00",
              "false negative" = "#E69F00", "wrong sign" = "#CC79A7")
out_labs <- c(out_lev[-5], "wrong sign\n(opposite DE direction)")

## primary analysis: counts divided by their RTA overdispersion (count_mode "divided")
sq <- read_csv(file.path(seq_dir, "sequin_de_results.csv"), show_col_types = FALSE) %>%
  filter(level == "transcript", count_mode == "divided",
         !is.na(pvalue), is.finite(obs_lfc)) %>%   # as summarise_calls()
  mutate(outcome = factor(outcome, levels = out_lev),
         capped  = abs(obs_lfc) > y_cap,
         y       = pmax(pmin(obs_lfc, y_cap), -y_cap))
## top-left: DE calls among the TESTED sequins (the points); bottom-right: how
## many sequins could be tested (>= 10 reads in >= 3 libraries), then R2 / slope;
## that corner is empty by design (expected down, observed strongly up).
## Sequins detected / in the mix: panel c.
sq_s <- read_csv(file.path(seq_dir, "summary.csv"), show_col_types = FALSE) %>%
  filter(level == "transcript", count_mode == "divided") %>%
  mutate(lab_top = sprintf("precision = %.2f\nsensitivity = %.2f\nF1 = %.2f",
                           precision, sensitivity, f1),
         lab_bot = sprintf("%d tested\nR\u00b2 = %.2f\nslope = %.2f",
                           as.integer(tested), r2, slope))

pC <- ggplot(sq, aes(exp_lfc, y)) +
  annotate("rect", xmin = -lfc, xmax = lfc, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_point(aes(colour = outcome, shape = capped), size = 1.0, alpha = 0.9, show.legend = TRUE,
             position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  geom_text(data = sq_s, aes(x = -Inf, y = Inf, label = lab_top),
            hjust = -0.05, vjust = 1.1, lineheight = 1.46, size = 7 / .pt) +
  geom_text(data = sq_s, aes(x = Inf, y = -Inf, label = lab_bot),
            hjust = 1.05, vjust = -0.3, lineheight = 1.46, size = 7 / .pt) +
  scale_colour_manual(values = out_cols, labels = out_labs, drop = FALSE, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), breaks = "TRUE",
                     labels = sprintf("|observed log2FC| > %.1f", y_cap), name = NULL) +
  ## bottom room: the tested + R2 / slope block sits below the data
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.25))) +
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = 2)),
         shape  = guide_legend(order = 2, override.aes = list(size = 2, shape = 2))) +
  labs(x = "Expected log2FC\n(Mix A / Mix B)", y = "Observed log2FC (edgeR)",
       title = "Sequin ground truth",
       ## glmTreat TESTS against |log2FC| <= 0.5; it is not a fold-change filter
       ## three lines: the second line no longer runs past the panel edge
       subtitle = paste0("Mix A (H1975) vs Mix B (HCC827)\nedgeR glmTreat (H0 |log2FC| <= ", lfc,
                         "), FDR < 0.05\nshaded: not truly DE")) +
  theme(legend.spacing.y = unit(0, "mm")) + leg_right(1)

## ================================================================ e =====
## sample-sample Spearman correlation (01_de_edger.R: cor(logCPM, method =
## "spearman") over the isoforms it tested), samples ordered by line
rho_m <- as.matrix(read.csv(file.path(de_dir, cor_level, "sample_spearman.csv"),
                            row.names = 1, check.names = FALSE))
ord <- c(sort(grep(paste0("^", line_lab[["TKI_sensitive"]]), colnames(rho_m), value = TRUE)),
         sort(grep(paste0("^", line_lab[["TKI_resistant"]]), colnames(rho_m), value = TRUE)))
rho_l <- as_tibble(rho_m, rownames = "s1") %>%
  pivot_longer(-s1, names_to = "s2", values_to = "rho") %>%
  mutate(s1 = factor(s1, levels = rev(ord)), s2 = factor(s2, levels = ord))

pD <- ggplot(rho_l, aes(s2, s1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradientn(colours = cor_cols, name = "Spearman \u03c1") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = "Sample correlation",
       subtitle = paste0(str_to_sentence(lv_lab[[cor_level]]), "-level expression (logCPM)")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.key.height = unit(6, "mm")) +                    # colour bar 5 x 6 = 30 mm: as tall as
                                                                # the 6 tiles (8 mm ran past the panel)
  leg_right(1)

## ================================================================ f =====
## volcano, gene level; only markers that were tested and are DE are labelled
res <- read_csv(file.path(de_dir, "all_results.csv"), show_col_types = FALSE)

vol <- res %>%
  filter(level == "gene") %>%
  mutate(symbol = sub("_ENSG.*$", "", feature),               # "LSP1P5_ENSG..." -> "LSP1P5"
         call   = factor(call, levels = c("Up", "Down", "NS")),
         y      = -log10(pmax(padj, 1e-300))) %>%
  arrange(dplyr::desc(call))                                 # NS drawn first
mk_l <- filter(vol, symbol %in% markers, call != "NS")
not_tested <- setdiff(markers, vol$symbol)
not_de     <- setdiff(filter(vol, symbol %in% markers, call == "NS")$symbol, mk_l$symbol)
if (length(not_tested))
  message("marker genes not tested at gene level (not labelled): ", paste(not_tested, collapse = ", "))
if (length(not_de))
  message("marker genes not DE at gene level (not labelled): ", paste(not_de, collapse = ", "))
n_up <- sum(vol$call == "Up"); n_dn <- sum(vol$call == "Down")

pE <- ggplot(vol, aes(log2FC, y, colour = call)) +
  geom_point(data = filter(vol, call == "NS"), size = 0.7, alpha = 0.4, stroke = 0) +
  geom_point(data = filter(vol, call != "NS"), size = 1.0, alpha = 0.6, stroke = 0) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(alpha), linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_point(data = mk_l, shape = 21, fill = NA, colour = "black", size = 1.6, stroke = 0.5) +
  lab_repel(mk_l, aes(label = symbol)) +
  annotate("text", x = -Inf, y = Inf, label = paste("Down:", scales::comma(n_dn)),
           hjust = -0.15, vjust = 1.6, colour = dir_cols[["Down"]], fontface = "plain",
           size = 10 / .pt) +
  annotate("text", x = Inf, y = Inf, label = paste("Up:", scales::comma(n_up)),
           hjust = 1.15, vjust = 1.6, colour = dir_cols[["Up"]], fontface = "plain",
           size = 10 / .pt) +
  scale_colour_manual(values = c(dir_cols[c("Up", "Down")], NS = "grey80"),
                      labels = c(Up = "upregulated", Down = "downregulated", NS = "not significant"),
                      name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "log2FC, H1975 / HCC827 (predFC)", y = "-log10(FDR)",
       title = "Differential Gene Expression",
       subtitle = paste0("FDR ", alpha, ", |log2FC| > ", lfc)) +
  leg_right(0.66)

## ================================================================ g =====
## DE isoforms of the marker genes. Gene = Bambu's gene assignment (GENEID of
## the transcript count table, i.e. what the DE used), named via the
## reference metadata. Rows = isoforms grouped by gene (marker order) and
## ranked by log2FC; columns = samples + the isoform DE call.
iso_val <- read_tsv(run_path("counts_tx"),
                    show_col_types = FALSE, col_select = c(TXNAME, GENEID))
xt   <- readRDS(file.path(de_dir, "transcript", "logcpm.rds"))
grp  <- setNames(as.character(xt$group), colnames(xt$mat))
smp  <- c(names(grp)[grp == "TKI_sensitive"], names(grp)[grp == "TKI_resistant"])
refm <- read_csv(run_path("ann_meta"),
                 col_select = c(gene_id_version, gene_name), show_col_types = FALSE) %>%
  distinct(gene_id_version, .keep_all = TRUE)

iso <- res %>%
  filter(level == "transcript", feature %in% rownames(xt$mat)) %>%
  left_join(dplyr::select(iso_val, id = TXNAME, gene_id = GENEID), by = "id") %>%
  left_join(refm, by = c(gene_id = "gene_id_version")) %>%
  filter(gene_name %in% markers, !hm_only_de | call != "NS") %>%
  mutate(gene = factor(gene_name, levels = intersect(markers, gene_name))) %>%
  group_by(gene) %>%
  { if (is.finite(max_iso_per_gene))
      slice_min(., padj, n = max_iso_per_gene, with_ties = FALSE) else . } %>%
  ungroup() %>%
  arrange(gene, dplyr::desc(log2FC_test))
message(sprintf("panel g: %d isoforms from %d marker genes%s", nrow(iso), n_distinct(iso$gene),
                if (hm_only_de) " (DE isoforms only)" else ""))
no_iso <- setdiff(markers, as.character(iso$gene))
if (length(no_iso))
  message("marker genes with no ", if (hm_only_de) "DE " else "", "isoform in panel g: ",
          paste(no_iso, collapse = ", "))

z <- t(scale(t(xt$mat[iso$feature, smp, drop = FALSE])))
## columns, right to left: DE call (x = 0), HCC827, a gap of 0.3 column, H1975
smp_lr <- if (setequal(hm_samples, smp)) hm_samples else rev(smp)       # left to right
smp_rl <- rev(smp_lr)
x_pos  <- setNames(-(seq_along(smp_rl) + 0.3 * (grp[smp_rl] == "TKI_resistant")), smp_rl)
hm <- as_tibble(z, rownames = "feature") %>%
  pivot_longer(-feature, names_to = "sample", values_to = "z") %>%
  left_join(dplyr::select(iso, feature, gene), by = "feature") %>%
  mutate(feature = factor(feature, levels = rev(unique(iso$feature))),   # first isoform on top
         x = unname(x_pos[sample]),
         z = pmax(pmin(z, z_lim), -z_lim))
de_col <- iso %>%
  transmute(feature = factor(feature, levels = levels(hm$feature)), gene,
            call = factor(call, levels = c("Up", "Down", "NS")))
## "DE" under the DE column (horizontal; the sample names are vertical)
de_lab <- tibble(gene = factor(tail(levels(iso$gene)[levels(iso$gene) %in% iso$gene], 1),
                               levels = levels(iso$gene)), x = 0, label = "DE")

pF <- ggplot(hm, aes(x, feature)) +
  geom_tile(aes(fill = z)) +
  geom_point(data = de_col, aes(0, feature, colour = call), shape = 15, size = 1.3) +
  geom_text(data = de_lab, aes(x, -Inf, label = label), inherit.aes = FALSE,
            vjust = 1.4, size = 8.5 / .pt, colour = "grey30") +
  coord_cartesian(clip = "off") +
  scale_fill_gradientn(colours = hm_cols, limits = c(-z_lim, z_lim), name = "z-score (logCPM)") +
  scale_colour_manual(values = c(dir_cols[c("Up", "Down")], NS = "grey80"),
                      labels = c(Up = "Up in H1975", Down = "Down in H1975", NS = "not DE"),
                      name = NULL, drop = TRUE) +
  scale_x_continuous(breaks = unname(x_pos), labels = names(x_pos),
                     expand = expansion(add = 0.5)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(gene ~ ., scales = "free_y", space = "free_y") +
  guides(fill = guide_colourbar(order = 1),
         colour = guide_legend(order = 2, override.aes = list(size = 2.5))) +
  labs(x = NULL, y = NULL, title = "Differential Transcript Expression",
       subtitle = paste(c(if (is.finite(max_iso_per_gene)) paste("Top", max_iso_per_gene),
                          if (hm_only_de) "DE" else "Tested",
                          "isoforms per loci, ranked by log2FC"), collapse = " ")) +
  theme(strip.text.y.right = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 8.5, colour = "grey10"),
        strip.background = element_blank(),
        strip.clip = "off",
        panel.spacing.y = unit(0.84, "mm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  leg_right(1)

## ===================================================== single panels =====
## one file per panel; size (width, height, mm) = the panel's extent on the
## final Fig. 2 page (A3) + the 0.5 mm plot margin on each side
## v2: new panel c (size chosen for it: place it and adjust); d = the old c,
## 4.1 mm taller for the third subtitle line; e-g = the old d-f, same sizes.
## Old files fig2c_sequin_isoform.*, fig2d_sample_correlation.*, ... are not all
## overwritten (new names): delete the old ones from figures/fig2_panels/.
panels <- list(
  fig2a_read_support          = list(p = pA,  size = c(129.3,  69.8)),
  fig2b_novel_biotype         = list(p = pB,  size = c( 60.9,  74.2)),
  fig2c_sequin_identification = list(p = pID, size = c( 75.1,  60.0)),
  fig2d_sequin_isoform        = list(p = pC,  size = c( 75.1,  71.0)),
  fig2e_sample_correlation    = list(p = pD,  size = c( 68.6,  55.2)),
  fig2f_volcano               = list(p = pE,  size = c(111.3,  85.2)),
  fig2g_isoform_heatmap       = list(p = pF,  size = c( 80.0, 215.8))
)
for (nm in names(panels)) save_panel(panels[[nm]]$p, nm, panels[[nm]]$size, panel_dir)
message("Done: ", length(panels), " panels in ", panel_dir)
