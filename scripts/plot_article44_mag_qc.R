#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(forcats)
  library(scales)
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
set.seed(20260744)

tier_levels <- c("High quality", "Medium quality", "Low/failed")
tier_colors <- c("High quality" = "#009E73", "Medium quality" = "#0072B2", "Low/failed" = "#D55E00")
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

qc <- read_tsv(file.path(input_dir, "mag-quality-summary.tsv"), show_col_types = FALSE) |>
  mutate(
    MIMAGQuality = factor(MIMAGQuality, levels = tier_levels),
    GUNCPass = as.logical(GUNCPass),
    Label = sprintf("MAG %02d", row_number())
  )
p_quality <- ggplot(
  qc,
  aes(x = CheckM2Completeness, y = CheckM2Contamination, color = MIMAGQuality,
      shape = GUNCPass, size = BinBp)
) +
  annotate("rect", xmin = 90, xmax = Inf, ymin = -Inf, ymax = 5, fill = "#009E73", alpha = 0.07) +
  geom_vline(xintercept = c(50, 90), linetype = c(3, 2), color = "#555555", linewidth = 0.45) +
  geom_hline(yintercept = c(10, 5), linetype = c(3, 2), color = "#555555", linewidth = 0.45) +
  geom_point(alpha = 0.86, stroke = 0.55) +
  scale_color_manual(values = tier_colors) +
  scale_shape_manual(values = c(`FALSE` = 4, `TRUE` = 16), labels = c(`FALSE` = "GUNC fail", `TRUE` = "GUNC pass")) +
  scale_size_continuous(range = c(2.2, 6.5), labels = label_number(scale = 1e-6, suffix = " MB")) +
  coord_cartesian(xlim = c(0, 102), ylim = c(0, max(12, max(qc$CheckM2Contamination, na.rm = TRUE) * 1.12))) +
  labs(
    title = "MIMAG quality requires more than the CheckM2 quadrant",
    subtitle = "High quality additionally requires complete 5S/16S/23S rRNA and >=18 tRNA isotypes",
    x = "CheckM2 completeness (%)", y = "CheckM2 contamination (%)",
    color = "MIMAG quality", shape = NULL, size = "Genome size"
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_quality, file.path(output_dir, "44-quality-landscape"), width = 205, height = 145)

requirements <- qc |>
  transmute(
    MAG, Label,
    `Completeness >90%` = CheckM2Completeness > 90,
    `Contamination <5%` = CheckM2Contamination < 5,
    `GUNC pass` = GUNCPass,
    `Complete 5S` = as.logical(Complete5S),
    `Complete 16S` = as.logical(Complete16S),
    `Complete 23S` = as.logical(Complete23S),
    `>=18 tRNA isotypes` = TRNAIsotypes >= 18
  ) |>
  pivot_longer(-c(MAG, Label), names_to = "Requirement", values_to = "Pass") |>
  mutate(
    Requirement = factor(
      Requirement,
      levels = c("Completeness >90%", "Contamination <5%", "GUNC pass", "Complete 5S", "Complete 16S", "Complete 23S", ">=18 tRNA isotypes")
    ),
    Label = factor(Label, levels = rev(unique(Label)))
  )
p_requirements <- ggplot(requirements, aes(x = Requirement, y = Label, fill = Pass)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = if_else(Pass, "PASS", "FAIL")), size = 3.0) +
  scale_fill_manual(values = c(`FALSE` = "#E69F00", `TRUE` = "#009E73"), guide = "none") +
  labs(
    title = "Every high-quality requirement is checked explicitly",
    subtitle = "Rows are selected refinement candidates; labels replace taxonomy until GTDB-Tk",
    x = NULL, y = NULL
  ) +
  theme_pub() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_pub(p_requirements, file.path(output_dir, "44-mimag-requirements"), width = 210, height = max(110, 8 * nrow(qc) + 55))

graph <- read_tsv(file.path(input_dir, "assembly-graph-audit.tsv"), show_col_types = FALSE) |>
  left_join(qc |> select(MAG, MIMAGQuality, Label), by = "MAG")
p_components <- ggplot(graph, aes(y = fct_reorder(Label, K141Components), x = K141Components, fill = MIMAGQuality)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tier_colors) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(title = "Within-bin graph fragmentation", x = "k141 graph components", y = NULL, fill = NULL) +
  theme_pub() + theme(legend.position = "none")
p_boundary <- ggplot(graph, aes(x = K141BoundaryEdges, y = PairedBoundaryPct, color = MIMAGQuality, size = MAGContigs)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = tier_colors) +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    title = "Boundary evidence for manual review",
    subtitle = "Not used to assign MIMAG quality",
    x = "k141 boundary edges", y = "Paired-read boundary weight (%)", color = NULL, size = "Contigs"
  ) +
  theme_pub() + theme(legend.position = "right")
p_graph <- (p_components | p_boundary) + plot_layout(widths = c(1, 1.08)) +
  plot_annotation(title = "Assembly graph and paired-read links are a separate audit layer")
save_pub(p_graph, file.path(output_dir, "44-assembly-graph-audit"), width = 220, height = 150)

comparison <- qc |>
  select(MAG, MIMAGQuality, CheckM1StrainHeterogeneityPct,
         CheckM1Completeness, CheckM2Completeness, CheckM1Contamination, CheckM2Contamination) |>
  pivot_longer(
    cols = matches("CheckM[12](Completeness|Contamination)"),
    names_to = c("Tool", "Metric"), names_pattern = "(CheckM[12])(Completeness|Contamination)",
    values_to = "Value"
  ) |>
  pivot_wider(names_from = Tool, values_from = Value)
p_checkm <- ggplot(comparison, aes(x = CheckM1, y = CheckM2, color = MIMAGQuality, size = CheckM1StrainHeterogeneityPct)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "#555555") +
  geom_point(alpha = 0.82) +
  facet_wrap(~Metric, scales = "free") +
  scale_color_manual(values = tier_colors) +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    title = "CheckM1 is retained as an orthogonal strain audit",
    subtitle = "CheckM2 defines the primary completeness and contamination values",
    x = "CheckM1 estimate (%)", y = "CheckM2 estimate (%)", color = NULL,
    size = "CheckM1 strain\nheterogeneity (%)"
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_checkm, file.path(output_dir, "44-checkm-audit"), width = 205, height = 135)
