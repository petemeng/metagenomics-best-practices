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
set.seed(20260748)

pal <- c("MOCK1" = "#0072B2", "MOCK2" = "#E69F00",
         "Catalog assigned" = "#0072B2", "Unmapped" = "#999999")
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

coverage <- read_tsv(file.path(input_dir, "coverm-long.tsv.gz"), show_col_types = FALSE)
primary <- coverage |>
  filter(Branch == "Primary 95% identity") |>
  group_by(SGB) |>
  mutate(MeanAbundance = mean(RelativeAbundancePct)) |>
  ungroup() |>
  mutate(SGB = fct_reorder(SGB, MeanAbundance))

p_heatmap <- ggplot(primary, aes(x = Sample, y = SGB, fill = RelativeAbundancePct)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", RelativeAbundancePct),
                color = RelativeAbundancePct >= 4), size = 2.8) +
  scale_color_manual(values = c(`FALSE` = "#111111", `TRUE` = "white"), guide = "none") +
  scale_fill_gradientn(colors = c("#F7FBFF", "#6BAED6", "#08306B"),
                       trans = log1p_trans(), breaks = c(1, 2, 5, 10)) +
  labs(
    title = "One non-redundant catalog, one abundance coordinate system",
    subtitle = "CoverM relative abundance after the 95% identity filter",
    x = NULL, y = "Species-level genome bin", fill = "Relative\nabundance (%)"
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_heatmap, file.path(output_dir, "48-mag-abundance-heatmap"), width = 190, height = 185)

low_breadth <- primary |> filter(CoveredFractionPct < 90)
p_depth <- ggplot(primary, aes(x = MeanDepth, y = CoveredFractionPct, color = Sample)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#555555") +
  geom_hline(yintercept = c(50, 90), linetype = c(2, 3), color = "#555555") +
  geom_point(size = 3, alpha = 0.82) +
  geom_text(data = low_breadth, aes(label = SGB), nudge_y = -2.2, size = 3,
            show.legend = FALSE, check_overlap = TRUE) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20), labels = label_number()) +
  scale_color_manual(values = pal[c("MOCK1", "MOCK2")]) +
  coord_cartesian(ylim = c(45, 101)) +
  labs(
    title = "Breadth distinguishes genome-wide support from local pileups",
    subtitle = "Detection requires mean depth >=1x and breadth >=50%; 90% marks broad support",
    x = "Mean depth (x, log scale)", y = "Covered fraction (%)", color = NULL
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_depth, file.path(output_dir, "48-breadth-depth-audit"), width = 190, height = 130)

sensitivity <- read_tsv(file.path(input_dir, "stringency-sensitivity.tsv"), show_col_types = FALSE)
top_sensitivity <- sensitivity |>
  slice_max(abs(DeltaRelativeAbundancePctPoints), n = 1, with_ties = FALSE)
p_sensitivity <- ggplot(
  sensitivity,
  aes(x = PrimaryRelativeAbundancePct, y = StrictRelativeAbundancePct, color = Sample)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555") +
  geom_point(size = 3, alpha = 0.82) +
  geom_text(data = top_sensitivity, aes(label = SGB), nudge_y = 0.12,
            size = 3, show.legend = FALSE) +
  scale_x_log10(breaks = c(1, 2, 5, 10), labels = label_number()) +
  scale_y_log10(breaks = c(1, 2, 5, 10), labels = label_number()) +
  scale_color_manual(values = pal[c("MOCK1", "MOCK2")]) +
  coord_equal() +
  labs(
    title = "Catalog abundance is stable at a stricter identity filter",
    subtitle = "Each point is one SGB; dashed line: no change",
    x = "Relative abundance at 95% identity (%)",
    y = "Relative abundance at 97% identity (%)", color = NULL
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_sensitivity, file.path(output_dir, "48-identity-sensitivity"), width = 185, height = 145)

capture <- read_tsv(file.path(input_dir, "sample-capture-summary.tsv"), show_col_types = FALSE) |>
  transmute(
    Sample,
    Branch = factor(Branch, levels = c("Primary 95% identity", "Strict 97% identity")),
    `Catalog assigned` = CatalogRelativeAbundancePct,
    Unmapped = UnmappedPct
  ) |>
  pivot_longer(c(`Catalog assigned`, Unmapped), names_to = "Read fate", values_to = "Percent") |>
  mutate(`Read fate` = factor(`Read fate`, levels = c("Unmapped", "Catalog assigned")))
p_capture <- ggplot(capture, aes(x = Branch, y = Percent, fill = `Read fate`)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.1f", Percent)), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.4) +
  facet_wrap(~ Sample, nrow = 1) +
  scale_fill_manual(values = pal[c("Catalog assigned", "Unmapped")]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 100.5)) +
  labs(
    title = "Catalog capture is not the same as total metagenome coverage",
    subtitle = "Unmapped includes organisms or sequence not represented by the 24-SGB catalog",
    x = NULL, y = "Reads (%)", fill = NULL
  ) +
  theme_pub() + theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))
save_pub(p_capture, file.path(output_dir, "48-catalog-capture"), width = 190, height = 125)
