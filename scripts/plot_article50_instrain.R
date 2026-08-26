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
set.seed(20260750)

pal <- c("MOCK1" = "#0072B2", "MOCK2" = "#E69F00",
         "Same-strain call" = "#009E73", "Below one or both gates" = "#CC79A7")
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EEEEEE"),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot"
    )
}
save_pub <- function(plot, file_base, width = 190, height = 130, dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}

profile <- read_tsv(file.path(input_dir, "profile-genome-long.tsv.gz"), show_col_types = FALSE) |>
  mutate(Sample = factor(Sample, levels = c("MOCK1", "MOCK2")))

diversity <- profile |>
  filter(PresentBreadth50AtMinCov, is.finite(NucleotideDiversity),
         NucleotideDiversity > 0)
label_diversity <- diversity |>
  group_by(SGB) |>
  summarise(MaxDiversity = max(NucleotideDiversity), .groups = "drop") |>
  slice_max(MaxDiversity, n = 6, with_ties = FALSE) |>
  inner_join(filter(diversity, Sample == "MOCK2"), by = "SGB")
p_diversity <- ggplot(
  diversity,
  aes(x = Sample, y = NucleotideDiversity, group = SGB)
) +
  geom_line(color = "#B8B8B8", linewidth = 0.55) +
  geom_point(aes(color = Sample), size = 2.8, alpha = 0.9) +
  geom_text(data = label_diversity, aes(label = SGB), nudge_x = 0.08,
            size = 2.8, hjust = 0, check_overlap = TRUE, color = "#333333") +
  scale_color_manual(values = pal[c("MOCK1", "MOCK2")]) +
  scale_y_log10(labels = label_scientific()) +
  coord_cartesian(xlim = c(0.85, 2.28), clip = "off") +
  labs(
    title = "Within-population diversity is genome specific",
    subtitle = "Observed-depth estimate; lines join the same SGB across samples",
    x = NULL, y = "Nucleotide diversity (log scale)", color = NULL
  ) +
  theme_pub() +
  theme(legend.position = "top", plot.margin = margin(5.5, 34, 5.5, 5.5))
save_pub(p_diversity, file.path(output_dir, "50-nucleotide-diversity"), 190, 135)

top_snv <- profile |>
  filter(PresentBreadth50AtMinCov, is.finite(SNVsPerMbpConsidered)) |>
  slice_max(SNVsPerMbpConsidered, n = 6, with_ties = FALSE)
p_snv <- ggplot(
  filter(profile, PresentBreadth50AtMinCov, is.finite(SNVsPerMbpConsidered)),
  aes(x = Coverage, y = SNVsPerMbpConsidered, color = Sample)
) +
  geom_point(size = 3, alpha = 0.82) +
  geom_text(data = top_snv, aes(label = SGB), nudge_y = 0.035,
            size = 2.8, show.legend = FALSE, check_overlap = TRUE) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_color_manual(values = pal[c("MOCK1", "MOCK2")]) +
  labs(
    title = "SNV burden must be read together with depth",
    subtitle = "Denominator: genome bases reaching the prespecified 5x minimum coverage",
    x = "Mean coverage (x, log scale)",
    y = "Within-sample SNVs per Mbp assessed (log scale)", color = NULL
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_snv, file.path(output_dir, "50-snv-density"), 190, 135)

detection_label <- profile |>
  filter(!PresentBreadth50AtMinCov | BreadthAtMinCovPct < 70) |>
  group_by(SGB) |>
  slice_min(BreadthAtMinCovPct, n = 1, with_ties = FALSE) |>
  ungroup()
p_detection <- ggplot(
  profile,
  aes(x = Coverage, y = BreadthAtMinCovPct, color = Sample)
) +
  geom_hline(yintercept = 50, linetype = 2, color = "#555555") +
  geom_point(aes(size = pmax(SNVsPerMbpConsidered, 0)), alpha = 0.82) +
  geom_text(data = detection_label, aes(label = SGB), nudge_y = -2.1,
            size = 2.8, show.legend = FALSE, check_overlap = TRUE) +
  scale_x_log10(labels = label_number()) +
  scale_size_continuous(range = c(2.2, 6), labels = label_number()) +
  scale_color_manual(values = pal[c("MOCK1", "MOCK2")]) +
  coord_cartesian(ylim = c(0, 102)) +
  labs(
    title = "Database-mode comparison applies a breadth gate",
    subtitle = "Dashed line: at least 50% of the genome reaches 5x coverage",
    x = "Mean coverage (x, log scale)", y = "Genome at >=5x coverage (%)",
    color = NULL, size = "SNVs per\nMbp assessed"
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_detection, file.path(output_dir, "50-profile-detection-audit"), 195, 140)

comparison <- read_tsv(
  file.path(input_dir, "pairwise-genome-comparison.tsv"), show_col_types = FALSE
) |>
  filter(Compared, is.finite(PopANIPct), is.finite(PercentGenomeCompared)) |>
  mutate(
    PopulationDifferencesPerMbp = pmax((100 - PopANIPct) * 10000, 0.05),
    Call = if_else(SameStrainRule, "Same-strain call", "Below one or both gates")
  )
label_compare <- comparison |>
  arrange(desc(PopulationDifferencesPerMbp), PercentGenomeCompared) |>
  slice_head(n = 8)
p_compare <- ggplot(
  comparison,
  aes(x = PercentGenomeCompared, y = PopulationDifferencesPerMbp, color = Call)
) +
  geom_vline(xintercept = 50, linetype = 2, color = "#555555") +
  geom_hline(yintercept = 10, linetype = 2, color = "#555555") +
  geom_point(aes(size = ComparedBases), alpha = 0.88) +
  geom_text(data = label_compare, aes(label = SGB), nudge_y = 0.08,
            size = 2.8, show.legend = FALSE, check_overlap = TRUE) +
  scale_y_log10(labels = label_number()) +
  scale_size_continuous(range = c(2.5, 6), labels = label_number(scale = 1e-6,
                                                                suffix = " M")) +
  scale_color_manual(values = pal[c("Same-strain call", "Below one or both gates")]) +
  coord_cartesian(xlim = c(0, 102)) +
  labs(
    title = "A strain call requires identity and comparable genome support",
    subtitle = "Operational gates: popANI >=99.999% and genome compared >=50%",
    x = "Genome compared (%)",
    y = "Population differences per Mbp (log scale)",
    color = NULL, size = "Compared bases"
  ) +
  guides(
    size = guide_legend(order = 1, nrow = 1),
    color = guide_legend(order = 2, nrow = 1)
  ) +
  theme_pub() +
  theme(
    legend.position = "top", legend.box = "vertical",
    legend.justification = "left"
  )
save_pub(p_compare, file.path(output_dir, "50-popani-overlap"), 205, 155)
