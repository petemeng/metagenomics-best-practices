#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
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
set.seed(20260741)

pal_pub <- c(
  MOCK1 = "#0072B2", MOCK2 = "#D55E00",
  `No concordant alignment` = "#BDBDBD",
  `Unique concordant alignment` = "#009E73",
  `Multiple concordant alignments` = "#CC79A7"
)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot"
    )
}
save_pub <- function(plot, file_base, width = 180, height = 125, dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}

fate <- read_tsv(file.path(input_dir, "mapping-fate-long.tsv"), show_col_types = FALSE) |>
  mutate(
    Sample = factor(Sample, levels = c("MOCK1", "MOCK2")),
    ReadFate = factor(
      ReadFate,
      levels = c(
        "No concordant alignment",
        "Unique concordant alignment",
        "Multiple concordant alignments"
      )
    )
  )
p_fate <- ggplot(fate, aes(x = Sample, y = Percent, fill = ReadFate)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = if_else(Percent >= 3, sprintf("%.1f%%", Percent), "")),
    position = position_stack(vjust = 0.5), size = 3.5, color = "black"
  ) +
  scale_fill_manual(
    values = pal_pub, drop = FALSE,
    labels = c("No concordant", "Unique concordant", "Multiple concordant")
  ) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Paired-read alignment fate",
    subtitle = "Bowtie2 concordant categories conserve the input-pair ledger",
    x = NULL, y = "Input read pairs", fill = NULL
  ) +
  theme_pub() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 10))
save_pub(p_fate, file.path(output_dir, "41-mapping-fate"), width = 175, height = 120)

depth <- read_tsv(file.path(input_dir, "contig-depth-wide.tsv.gz"), show_col_types = FALSE)
correlation <- read_tsv(file.path(input_dir, "depth-correlation.tsv"), show_col_types = FALSE)
p_depth <- ggplot(
  depth,
  aes(x = MOCK1MeanDepth + 0.01, y = MOCK2MeanDepth + 0.01, color = log10(LengthBp))
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#4D4D4D", linewidth = 0.45) +
  geom_point(alpha = 0.38, size = 0.85) +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  scale_y_log10(labels = label_number(accuracy = 0.01)) +
  scale_color_viridis_c(option = "C", end = 0.92) +
  coord_equal() +
  labs(
    title = "Co-abundance is informative but not identical across samples",
    subtitle = sprintf(
      "18,354 shared contigs; Pearson r = %.3f and Spearman rho = %.3f\non log1p depth",
      correlation$PearsonLog1p, correlation$SpearmanLog1p
    ),
    x = "MOCK1 JGI mean depth (+0.01)",
    y = "MOCK2 JGI mean depth (+0.01)",
    color = expression(log[10] * " contig length")
  ) +
  theme_pub() +
  theme(legend.position = "right", panel.grid.major = element_line(color = "#EEEEEE"))
save_pub(p_depth, file.path(output_dir, "41-depth-concordance"), width = 175, height = 145)

long <- read_tsv(file.path(input_dir, "contig-depth-long.tsv.gz"), show_col_types = FALSE) |>
  mutate(Sample = factor(Sample, levels = c("MOCK1", "MOCK2")))
p_breadth <- ggplot(long, aes(x = JgiMeanDepth + 0.01, y = BreadthPct, color = Sample)) +
  geom_point(alpha = 0.20, size = 0.55) +
  geom_hline(yintercept = 50, linetype = 3, color = "#666666", linewidth = 0.45) +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 100), labels = label_percent(scale = 1)) +
  scale_color_manual(values = pal_pub) +
  facet_wrap(vars(Sample), nrow = 1) +
  labs(
    title = "Mean depth and breadth answer different questions",
    subtitle = "Breadth uses primary mapped records; JGI depth additionally applies\n97% end-to-end identity",
    x = "JGI mean depth (+0.01)", y = "Contig breadth", color = "Sample"
  ) +
  theme_pub() +
  theme(legend.position = "none", panel.grid.major.x = element_line(color = "#EEEEEE"))
save_pub(p_breadth, file.path(output_dir, "41-depth-breadth"), width = 190, height = 120)
