#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(forcats)
  library(scales)
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
set.seed(20260745)

pal <- c("Article44-selected" = "#0072B2", "Article42-QC-pass" = "#E69F00")
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#EEEEEE"),
      axis.text = element_text(color = "black"), legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA), plot.title.position = "plot"
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

threshold <- read_tsv(file.path(input_dir, "threshold-summary.tsv"), show_col_types = FALSE) |>
  mutate(Branch = factor(Branch, levels = c("Species 95% ANI", "Near-clone 99.9% ANI"))) |>
  select(Branch, Representatives, GenomesRemoved) |>
  pivot_longer(-Branch, names_to = "Fate", values_to = "Genomes") |>
  mutate(Fate = recode(Fate, Representatives = "Representatives", GenomesRemoved = "Clustered duplicates"))
p_yield <- ggplot(threshold, aes(x = Branch, y = Genomes, fill = Fate)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = Genomes), position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Representatives" = "#0072B2", "Clustered duplicates" = "#999999")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "ANI threshold changes the biological unit represented",
    subtitle = "All 124 genomes pass the same CheckM2 and GUNC input gate",
    x = NULL, y = "Genomes", fill = NULL
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_yield, file.path(output_dir, "45-dereplication-yield"), width = 175, height = 120)

membership <- read_tsv(file.path(input_dir, "cluster-membership.tsv.gz"), show_col_types = FALSE) |>
  filter(Branch == "Species 95% ANI") |>
  mutate(IsRepresentative = as.logical(IsRepresentative))
p_quality <- ggplot(membership, aes(x = Completeness, y = Contamination)) +
  annotate("rect", xmin = 90, xmax = Inf, ymin = -Inf, ymax = 5, fill = "#009E73", alpha = 0.07) +
  geom_vline(xintercept = c(50, 90), linetype = c(3, 2), color = "#555555", linewidth = 0.45) +
  geom_hline(yintercept = c(10, 5), linetype = c(3, 2), color = "#555555", linewidth = 0.45) +
  geom_point(aes(color = SourceStage, shape = IsRepresentative, size = ClusterSize), alpha = 0.78, stroke = 0.55) +
  scale_color_manual(
    values = pal,
    labels = c(
      "Article44-selected" = "Article 44 selected",
      "Article42-QC-pass" = "Article 42 QC pass"
    )
  ) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "Cluster member", `TRUE` = "Representative")) +
  scale_size_continuous(range = c(2, 7), breaks = pretty_breaks(n = 4)) +
  coord_cartesian(xlim = c(48, 102), ylim = c(0, max(10, max(membership$Contamination) * 1.08))) +
  labs(
    title = "Representative choice combines quality and cluster centrality",
    subtitle = "A representative is not automatically a complete MIMAG high-quality genome",
    x = "CheckM2 completeness (%)", y = "CheckM2 contamination (%)",
    color = "Input stage", shape = NULL, size = "Cluster size"
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_quality, file.path(output_dir, "45-representative-quality"), width = 205, height = 140)

pairs <- read_tsv(file.path(input_dir, "pairwise-ani.tsv.gz"), show_col_types = FALSE) |>
  filter(Branch == "Near-clone 99.9% ANI") |>
  mutate(SameCluster = as.logical(SameCluster))
p_ani <- ggplot(pairs, aes(x = ANIPct, y = AlignmentFractionPct, color = SameCluster)) +
  geom_vline(xintercept = 99.9, linetype = 2, color = "#555555") +
  geom_hline(yintercept = 30, linetype = 2, color = "#555555") +
  geom_point(alpha = 0.62, size = 2) +
  scale_color_manual(values = c(`FALSE` = "#999999", `TRUE` = "#009E73"),
                     labels = c(`FALSE` = "Split at 99.9% ANI", `TRUE` = "Same near-clone cluster")) +
  coord_cartesian(xlim = c(99.45, 100.01), ylim = c(0, 101)) +
  labs(
    title = "ANI and aligned fraction must be interpreted together",
    subtitle = "The 99.9% sensitivity branch separates four comparisons below the ANI cutoff",
    x = "fastANI (%)", y = "Alignment fraction (%)", color = NULL
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_ani, file.path(output_dir, "45-ani-alignment-audit"), width = 190, height = 130)

retention <- read_tsv(file.path(input_dir, "source-retention.tsv"), show_col_types = FALSE) |>
  filter(Branch == "Species 95% ANI") |>
  select(SourceBranch, InputGenomes, Representatives) |>
  pivot_longer(-SourceBranch, names_to = "Set", values_to = "Genomes") |>
  mutate(
    Set = recode(Set, InputGenomes = "Input genomes", Representatives = "Representatives"),
    SourceBranch = fct_reorder(SourceBranch, Genomes, .fun = max)
  )
p_source <- ggplot(retention, aes(y = SourceBranch, x = Genomes, fill = Set)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  geom_text(aes(label = Genomes), position = position_dodge(width = 0.72), hjust = -0.15, size = 3.2) +
  scale_fill_manual(values = c("Input genomes" = "#999999", "Representatives" = "#0072B2")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Dereplication removes reconstruction redundancy",
    subtitle = "Source counts are an audit, not a ranking of binning algorithms",
    x = "Genomes", y = NULL, fill = NULL
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_source, file.path(output_dir, "45-source-retention"), width = 205, height = 130)
