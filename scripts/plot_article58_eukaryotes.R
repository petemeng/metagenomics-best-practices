#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: plot_article58_eukaryotes.R FROZEN_DIR FIGURE_DIR")
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- args[[2]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260758)

pal_pub <- c(
  "Blue" = "#0072B2", "Orange" = "#D55E00",
  "Green" = "#009E73", "Sky" = "#56B4E9",
  "Yellow" = "#E69F00", "Purple" = "#CC79A7",
  "Gray" = "#7A7A7A", "Light" = "#E9EEF2",
  "Dark" = "#253238", "White" = "#FFFFFF"
)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) + theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#ECEFF1", linewidth = 0.3),
    axis.text = element_text(color = "black"),
    legend.key = element_blank(),
    strip.background = element_rect(fill = "#EEF3F5", color = NA),
    plot.title.position = "plot",
    plot.caption = element_text(color = "#455A64", hjust = 0)
  )
}
save_pub <- function(plot, stem, width = 205, height = 145, dpi = 350) {
  base <- file.path(figure_dir, stem)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}

evidence <- read_tsv(file.path(input_dir, "eukdetect-species-evidence.tsv"), show_col_types = FALSE) %>%
  mutate(Name = factor(Name, levels = rev(Name[order(TotalMarkerReads)])))
p_markers <- ggplot(evidence, aes(x = ObservedMarkers, y = Name)) +
  geom_segment(aes(x = 0, xend = ObservedMarkers, yend = Name), color = pal_pub[["Light"]], linewidth = 2.2) +
  geom_point(aes(size = TotalMarkerReads, fill = PercentIdentity), shape = 21, color = pal_pub[["Dark"]]) +
  scale_fill_gradient(low = pal_pub[["Sky"]], high = pal_pub[["Blue"]], name = "Identity (%)") +
  scale_size_continuous(range = c(3, 9), labels = label_number(), name = "Marker reads") +
  labs(x = "Observed marker genes", y = NULL, title = "A  Marker evidence") +
  theme_pub()
p_rpksb <- ggplot(evidence, aes(x = RPKSB, y = Name)) +
  geom_col(fill = pal_pub[["Orange"]], width = 0.65) +
  geom_text(data = filter(evidence, RPKSB >= 0.8),
            aes(label = number(RPKSB, accuracy = 0.001)), hjust = 1.08,
            size = 3.2, color = "white", fontface = "bold") +
  geom_text(data = filter(evidence, RPKSB < 0.8),
            aes(label = number(RPKSB, accuracy = 0.001)), hjust = -0.12,
            size = 3.2, color = pal_pub[["Dark"]], fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "RPKSB", y = NULL, title = "B  Normalized signal") +
  theme_pub() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p_evidence <- p_markers + p_rpksb + plot_layout(widths = c(1.35, 1)) + plot_annotation(
  title = "EukDetect2 evidence in the Zymo D6300 mock",
  caption = "Species-level calls pass the built-in >=2-marker and >=4-read gate; RPKSB is not cell abundance."
)
save_pub(p_evidence, "58-eukdetect-marker-evidence", 225, 130)

composition <- read_tsv(file.path(input_dir, "eukdetect-composition.tsv"), show_col_types = FALSE) %>%
  group_by(Denominator, Component) %>% summarise(Percent = sum(Percent), .groups = "drop") %>%
  mutate(
    Denominator = factor(Denominator, levels = c("Expected total DNA", "Expected eukaryote-only", "EukDetect RelEuk")),
    Component = factor(Component, levels = c("Other organisms", "Other eukaryotes", "Cryptococcus neoformans complex", "Saccharomyces cerevisiae"))
  )
component_colors <- c(
  "Saccharomyces cerevisiae" = pal_pub[["Blue"]],
  "Cryptococcus neoformans complex" = pal_pub[["Orange"]],
  "Other eukaryotes" = pal_pub[["Purple"]],
  "Other organisms" = pal_pub[["Light"]]
)
p_composition <- ggplot(composition, aes(x = Denominator, y = Percent, fill = Component)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = component_colors, drop = FALSE) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = c(0, 0)) +
  labs(
    x = NULL, y = "Percent", fill = NULL,
    title = "The denominator changes the biological statement",
    caption = "RelEuk is normalized only among detected eukaryotes; it is not a fraction of the whole microbial community."
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pub() + theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "bottom", legend.text = element_text(size = 9.5)
  )
save_pub(p_composition, "58-eukdetect-abundance-denominator", 215, 155)

reference_metrics <- read_tsv(file.path(input_dir, "eukrep-reference-metrics.tsv"), show_col_types = FALSE) %>%
  filter(Metric %in% c("Sensitivity", "Specificity")) %>%
  mutate(Mode = factor(Mode, levels = c("Strict", "Balanced", "Lenient")))
p_reference <- ggplot(reference_metrics, aes(x = FragmentLength / 1000, y = Percent, color = Mode, shape = Metric)) +
  geom_hline(yintercept = 90, linetype = 3, color = pal_pub[["Gray"]]) +
  geom_line(aes(group = interaction(Mode, Metric)), linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Strict" = pal_pub[["Blue"]], "Balanced" = pal_pub[["Green"]], "Lenient" = pal_pub[["Orange"]])) +
  scale_shape_manual(values = c("Sensitivity" = 16, "Specificity" = 17)) +
  scale_x_continuous(breaks = c(3, 5, 10, 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x = "Reference fragment length (kb)", y = "Performance (%)", color = "Mode", shape = "Metric",
    title = "EukRep performance depends on fragment length and stringency",
    caption = "Deterministic benchmark: 80 evenly spaced fragments per species from 10 official Zymo reference genomes."
  ) +
  theme_pub() + theme(legend.position = "bottom")
save_pub(p_reference, "58-eukrep-reference-benchmark", 205, 145)

assembly_calls <- read_tsv(file.path(input_dir, "eukrep-assembly-calls.tsv"), show_col_types = FALSE) %>%
  count(Mode, Truth, Prediction, name = "Contigs") %>%
  mutate(
    Mode = factor(Mode, levels = c("Strict", "Balanced", "Lenient")),
    Truth = factor(Truth, levels = c("Eukaryote", "Prokaryote", "Unresolved"))
  )
p_assembly <- ggplot(assembly_calls, aes(x = Truth, y = Contigs, fill = Prediction)) +
  geom_col(position = "stack", width = 0.68) +
  facet_wrap(~Mode, nrow = 1) +
  scale_fill_manual(values = c("Eukaryote" = pal_pub[["Orange"]], "Prokaryote" = pal_pub[["Blue"]], "Unclassified" = pal_pub[["Gray"]])) +
  scale_y_continuous(labels = label_number()) +
  labs(
    x = "Reference-supported contig truth", y = "Contigs >=3 kb", fill = "EukRep call",
    title = "Assembly-level calls require reference and mode audits",
    caption = "Truth: >=95% identity and >=80% query coverage to an official Zymo v2 reference.\nOther contigs remain unresolved."
  ) +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
save_pub(p_assembly, "58-eukrep-assembly-audit", 220, 145)

ladder <- tibble::tibble(
  Level = factor(c("Marker reads", "Eukaryotic contigs", "Genome bin", "Activity", "Colonization"),
                 levels = rev(c("Marker reads", "Eukaryotic contigs", "Genome bin", "Activity", "Colonization"))),
  Evidence = c(
    "Taxon detected",
    "Domain-classified sequence",
    "Coherent genome reconstruction",
    "RNA / protein / growth evidence",
    "Persistence or experimental support"
  ),
  Ceiling = 1:5
)
p_ladder <- ggplot(ladder, aes(x = Ceiling, y = Level)) +
  geom_segment(aes(x = 0, xend = Ceiling, yend = Level), color = pal_pub[["Light"]], linewidth = 6, lineend = "round") +
  geom_point(aes(fill = Ceiling), shape = 21, size = 5, color = pal_pub[["Dark"]]) +
  geom_text(aes(label = Evidence), hjust = -0.06, size = 3.5) +
  scale_fill_gradient(low = pal_pub[["Sky"]], high = pal_pub[["Orange"]], guide = "none") +
  scale_x_continuous(limits = c(0, 8.4), breaks = 1:5, labels = c("Detection", "Sequence", "Genome", "Activity", "Persistence")) +
  labs(
    x = "Maximum claim supported by the evidence", y = NULL,
    title = "DNA evidence does not establish viability or colonization",
    caption = "Stop the biological claim at the weakest independently supported layer."
  ) +
  theme_pub() + theme(panel.grid = element_blank(), axis.ticks.y = element_blank())
save_pub(p_ladder, "58-eukaryote-evidence-ladder", 220, 135)
