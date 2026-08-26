#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(forcats)
  library(scales)
  library(stringr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1]]
}
input_dir <- value_after("--input-dir")
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260742)

branch_levels <- c(
  "MetaBAT2-MOCK1-only", "MetaBAT2-multisample",
  "SemiBin2-self-supervised", "VAMB-taxonomy-free", "TaxVAMB-Kraken2"
)
branch_labels <- c(
  "MetaBAT2\nMOCK1", "MetaBAT2\n2-sample",
  "SemiBin2\nself-sup.", "VAMB\nno taxonomy", "TaxVAMB\nKraken2"
)
branch_colors <- c(
  "MetaBAT2-MOCK1-only" = "#56B4E9",
  "MetaBAT2-multisample" = "#0072B2",
  "SemiBin2-self-supervised" = "#009E73",
  "VAMB-taxonomy-free" = "#E69F00",
  "TaxVAMB-Kraken2" = "#CC79A7"
)
tier_colors <- c(
  "HQ proxy" = "#009E73", "MQ proxy" = "#56B4E9",
  "Below MQ proxy" = "#BDBDBD", "QC pass" = "#0072B2", "QC fail" = "#D55E00"
)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
      axis.text = element_text(color = "black"), legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot"
    )
}
save_pub <- function(plot, file_base, width = 190, height = 125, dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}

bins <- read_tsv(file.path(input_dir, "bin-quality-truth-audit.tsv"), show_col_types = FALSE) |>
  mutate(
    Branch = factor(Branch, levels = branch_levels),
    TruthProxyTier = factor(TruthProxyTier, levels = c("Below MQ proxy", "MQ proxy", "HQ proxy")),
    QCMinimumPass = as.logical(QCMinimumPass)
  )
summary <- read_tsv(file.path(input_dir, "binner-summary.tsv"), show_col_types = FALSE) |>
  mutate(Branch = factor(Branch, levels = branch_levels))

truth_yield <- bins |>
  count(Branch, TruthProxyTier, name = "Bins")
p_truth <- ggplot(truth_yield, aes(x = Branch, y = Bins, fill = TruthProxyTier)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_text(aes(label = if_else(Bins > 0, as.character(Bins), "")),
            position = position_stack(vjust = 0.5), size = 3.1) +
  scale_x_discrete(labels = branch_labels) +
  scale_fill_manual(values = tier_colors, drop = FALSE) +
  labs(title = "Known-truth audit", x = NULL, y = "Candidate bins", fill = NULL) +
  theme_pub() + theme(axis.text.x = element_text(size = 9), legend.position = "bottom")

qc_yield <- bins |>
  mutate(QC = if_else(QCMinimumPass, "QC pass", "QC fail")) |>
  count(Branch, QC, name = "Bins")
p_qc <- ggplot(qc_yield, aes(x = Branch, y = Bins, fill = QC)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_text(aes(label = if_else(Bins > 0, as.character(Bins), "")),
            position = position_stack(vjust = 0.5), size = 3.1) +
  scale_x_discrete(labels = branch_labels) +
  scale_fill_manual(values = tier_colors, drop = FALSE) +
  labs(title = "Reference-free minimum QC", subtitle = "CheckM2 >=50%, <10%; GUNC pass",
       x = NULL, y = "Candidate bins", fill = NULL) +
  theme_pub() + theme(axis.text.x = element_text(size = 9), legend.position = "bottom")
p_yield <- (p_truth | p_qc) +
  plot_annotation(
    title = "Binner yield depends on the quality definition",
    subtitle = "Known mock truth is used only after binning"
  )
save_pub(p_yield, file.path(output_dir, "42-binner-quality-yield"), width = 235, height = 135)

p_scatter <- ggplot(
  bins,
  aes(x = DominantGenomeRecoveryPct, y = AlignedContaminationProxyPct,
      color = Branch, shape = QCMinimumPass, size = BinBp)
) +
  geom_vline(xintercept = c(50, 90), linetype = c(3, 2), color = "#666666", linewidth = 0.45) +
  geom_hline(yintercept = c(10, 5), linetype = c(3, 2), color = "#666666", linewidth = 0.45) +
  geom_point(alpha = 0.76, stroke = 0.45) +
  scale_color_manual(values = branch_colors, labels = branch_labels) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), na.translate = FALSE,
                     labels = c(`FALSE` = "QC fail", `TRUE` = "QC pass")) +
  scale_size_continuous(range = c(1.5, 6), labels = label_number(scale = 1e-6, suffix = " MB")) +
  coord_cartesian(xlim = c(0, 105), ylim = c(0, 100)) +
  labs(
    title = "Recovery and purity expose different binning failures",
    subtitle = "Coordinates are from the blinded post-hoc mock audit",
    x = "Dominant-genome recovery (%)", y = "Aligned contamination proxy (%)",
    color = NULL, shape = "CheckM2 + GUNC", size = "Bin size"
  ) +
  theme_pub() +
  theme(panel.grid.major = element_line(color = "#EEEEEE"), legend.position = "right")
save_pub(p_scatter, file.path(output_dir, "42-recovery-purity"), width = 205, height = 145)

single_multi <- read_tsv(file.path(input_dir, "single-vs-multisample.tsv"), show_col_types = FALSE)
single_multi_long <- tibble(
  Design = factor(c("MOCK1 only", "Multisample"), levels = c("MOCK1 only", "Multisample")),
  Shared = rep(single_multi$SharedBinnedContigs, 2),
  Unique = c(single_multi$SingleOnlyContigs, single_multi$MultiOnlyContigs)
) |>
  pivot_longer(c(Shared, Unique), names_to = "Status", values_to = "Contigs")
p_design <- ggplot(single_multi_long, aes(x = Design, y = Contigs, fill = Status)) +
  geom_col(width = 0.64, color = "white", linewidth = 0.25) +
  geom_text(
    data = filter(single_multi_long, Status == "Shared"),
    aes(label = Contigs), position = position_stack(vjust = 0.5), size = 3.5
  ) +
  geom_text(
    data = filter(single_multi_long, Status == "Unique"),
    aes(y = pmax(Contigs + 260, 360), label = paste0("Unique: ", Contigs)),
    size = 3.25
  ) +
  scale_fill_manual(values = c(Shared = "#0072B2", Unique = "#E69F00")) +
  labs(
    title = "Adding a second depth column changes bin membership",
    subtitle = sprintf("Adjusted Rand index on %s jointly binned contigs: %.3f",
                       comma(single_multi$AdjustedRandContigs),
                       single_multi$AdjustedRandOnSharedBinnedContigs),
    x = NULL, y = "Binned contigs", fill = NULL
  ) +
  theme_pub() + theme(legend.position = "bottom")
save_pub(p_design, file.path(output_dir, "42-single-vs-multisample"), width = 175, height = 125)

tax <- read_tsv(file.path(input_dir, "taxonomy-summary.tsv"), show_col_types = FALSE) |>
  mutate(
    Annotation = if_else(KrakenStatus == "U" | DeepestRank == "unclassified",
                         "Unclassified", str_to_title(DeepestRank)),
    Annotation = factor(
      Annotation,
      levels = c("Unclassified", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    )
  ) |>
  group_by(Annotation) |>
  summarise(Contigs = sum(Contigs), .groups = "drop")
tax_colors <- c(
  Unclassified = "#BDBDBD", Domain = "#999999", Phylum = "#56B4E9",
  Class = "#0072B2", Order = "#009E73", Family = "#E69F00",
  Genus = "#D55E00", Species = "#CC79A7"
)
p_tax <- ggplot(tax, aes(x = "Kraken2 Standard-8", y = Contigs, fill = Annotation)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = if_else(Contigs >= 200, comma(Contigs), "")),
    position = position_stack(vjust = 0.5), size = 3.2
  ) +
  scale_fill_manual(values = tax_colors, drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "TaxVAMB receives partial, not perfect, taxonomy",
    subtitle = "Deepest contiguous canonical rank for 10,203 contigs",
    x = NULL, y = "Contigs", fill = "Deepest rank"
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_tax, file.path(output_dir, "42-taxonomy-coverage"), width = 180, height = 130)
