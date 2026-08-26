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
if (length(args) != 2) {
  stop("Usage: plot_article56_host_evidence.R SUMMARY_DIR FIGURE_DIR")
}
summary_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- args[[2]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260756)

pal_pub <- c(
  "Blue" = "#0072B2",
  "Orange" = "#D55E00",
  "Green" = "#009E73",
  "Sky" = "#56B4E9",
  "Yellow" = "#E69F00",
  "Purple" = "#CC79A7",
  "Gray" = "#8A9499",
  "Light" = "#E9EEF2",
  "Dark" = "#253238",
  "White" = "#FFFFFF"
)
scale_color_pub <- function(...) scale_color_manual(values = pal_pub, ...)
scale_fill_pub <- function(...) scale_fill_manual(values = pal_pub, ...)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#ECEFF1", linewidth = 0.3),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot",
      plot.caption = element_text(color = "#455A64", hjust = 0),
      plot.subtitle = element_text(color = "#455A64")
    )
}
save_pub <- function(plot, file_base, width = 210, height = 145,
                     units = "mm", dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = units, device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = units, dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = units, dpi = dpi, compression = "lzw", bg = "white")
  invisible(plot)
}
wrap_text <- function(x, width = 30) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), "")
}

evidence <- read_tsv(
  file.path(summary_dir, "evidence-hierarchy.tsv"),
  show_col_types = FALSE
) %>%
  arrange(Rank) %>%
  mutate(
    EvidenceLabel = wrap_text(Evidence, 25),
    EvidenceLabel = factor(EvidenceLabel, levels = rev(EvidenceLabel)),
    SourceClass = if_else(
      grepl("Not integrated", iPHoPComponent),
      "Orthogonal to iPHoP",
      "Available to iPHoP"
    )
  )

p_ladder <- ggplot(
  evidence,
  aes(y = EvidenceLabel, x = Directness, color = SourceClass)
) +
  geom_segment(
    aes(x = 0.55, xend = Directness, yend = EvidenceLabel),
    linewidth = 1.05, color = "#D9E0E4"
  ) +
  geom_point(size = 4.2) +
  geom_text(
    aes(label = paste0("Tier ", Rank)),
    x = 0.35, hjust = 0, color = pal_pub[["Dark"]], size = 3.2
  ) +
  scale_color_manual(
    values = c(
      "Available to iPHoP" = pal_pub[["Blue"]],
      "Orthogonal to iPHoP" = pal_pub[["Orange"]]
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 5.25), breaks = 1:5,
    labels = c("Predictive", "", "", "", "Direct linkage")
  ) +
  labs(
    title = "Virus-host evidence has a claim ceiling",
    subtitle = "Rank reflects evidentiary directness, not a calibrated posterior probability",
    x = "Evidence directness", y = NULL,
    caption = paste0(
      "An integrated prophage requires host flanks and validated junctions.\n",
      "Co-binning in the same MAG is not sufficient."
    )
  ) +
  theme_pub() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.margin = margin(10, 18, 10, 10)
  )
save_pub(
  p_ladder,
  file.path(figure_dir, "56-evidence-ladder"),
  width = 225, height = 150
)

scope <- read_tsv(
  file.path(summary_dir, "benchmark-scope.tsv"),
  show_col_types = FALSE
)
lifestyle <- scope %>%
  filter(Category %in% c(
    "Temperate stratum",
    "Virulent stratum",
    "Not in either displayed lifestyle stratum"
  )) %>%
  mutate(
    CategoryLabel = recode(
      Category,
      "Temperate stratum" = "Temperate",
      "Virulent stratum" = "Virulent",
      "Not in either displayed lifestyle stratum" = "Neither"
    ),
    CategoryLabel = factor(
      CategoryLabel,
      levels = c("Neither", "Virulent", "Temperate")
    ),
    Share = Count / sum(Count),
    FillKey = recode(
      CategoryLabel,
      "Temperate" = "Blue",
      "Virulent" = "Orange",
      "Neither" = "Gray"
    )
  )
p_lifestyle <- ggplot(lifestyle, aes(x = 1, y = Count, fill = FillKey)) +
  geom_col(width = 0.58, color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = paste0(CategoryLabel, "\n", comma(Count), " (", percent(Share, accuracy = 0.1), ")")),
    position = position_stack(vjust = 0.5), color = "white", fontface = "bold",
    size = 3.3, lineheight = 0.95
  ) +
  coord_flip() +
  scale_fill_pub(guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  labs(
    title = "A. Held-out test composition",
    subtitle = "1,870 viruses spanning 170 host genera",
    x = NULL, y = "Viral genomes"
  ) +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

scale_data <- scope %>%
  filter(Category %in% c(
    "All test viruses",
    "High-quality prokaryotic-virus genomes",
    "Eukaryotic-virus genomes"
  )) %>%
  mutate(
    DatasetLabel = recode(
      Category,
      "All test viruses" = "Held-out test",
      "High-quality prokaryotic-virus genomes" = "IMG/VR application",
      "Eukaryotic-virus genomes" = "Euk-virus negative control"
    ),
    DatasetLabel = factor(
      DatasetLabel,
      levels = c("Held-out test", "Euk-virus negative control", "IMG/VR application")
    ),
    FillKey = recode(
      DatasetLabel,
      "Held-out test" = "Blue",
      "Euk-virus negative control" = "Orange",
      "IMG/VR application" = "Green"
    )
  )
p_scale <- ggplot(scale_data, aes(DatasetLabel, Count, fill = FillKey)) +
  geom_col(width = 0.67) +
  geom_text(aes(label = comma(Count)), hjust = -0.12, fontface = "bold", size = 3.7) +
  coord_flip() +
  scale_fill_pub(guide = "none") +
  scale_y_log10(
    breaks = c(1e3, 1e4, 1e5),
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "B. Benchmark and stress-test scale",
    subtitle = "Log scale; bars are different evaluation sets",
    x = NULL, y = "Viral genomes (log10 scale)"
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

p_scope <- p_lifestyle / p_scale +
  plot_annotation(
    title = "The iPHoP paper separates calibration, deployment, and negative controls",
    caption = paste0(
      "The 258-virus category is derived as 1,870 - 949 - 663 and is labelled\n",
      "only as outside the two displayed lifestyle strata."
    )
  )
save_pub(
  p_scope,
  file.path(figure_dir, "56-iphop-benchmark-scope"),
  width = 230, height = 180
)

confidence <- read_tsv(
  file.path(summary_dir, "confidence-contract.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    Use = factor(
      Use,
      levels = c(
        "Exploratory candidate list",
        "Default high-confidence screen",
        "Stringent sensitivity analysis"
      )
    ),
    FillKey = c("Orange", "Blue", "Green")
  )
p_confidence <- ggplot(
  confidence,
  aes(MinimumScore, NominalMaximumFDRPct, color = FillKey)
) +
  geom_line(color = "#B0BEC5", linewidth = 1) +
  geom_point(size = 5) +
  geom_text(
    aes(label = paste0("Score >=", MinimumScore, "\nFDR <=", NominalMaximumFDRPct, "%")),
    nudge_y = c(2.2, 2.2, 2.2), fontface = "bold", size = 3.5,
    color = pal_pub[["Dark"]]
  ) +
  geom_hline(
    yintercept = c(5, 10), linetype = "dashed",
    color = c(pal_pub[["Green"]], pal_pub[["Blue"]])
  ) +
  scale_color_pub(guide = "none") +
  scale_x_continuous(breaks = c(75, 90, 95), limits = c(72, 98)) +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 25), limits = c(0, 30)) +
  labs(
    title = "A confidence score is an empirical calibration contract",
    subtitle = "The nominal FDR is approximately 100 - score",
    x = "Minimum iPHoP confidence score", y = "Nominal maximum FDR (%)",
    caption = paste0(
      "The paper warns that calibration depends on benchmark composition.\n",
      "A score is not experimental confirmation of an individual link."
    )
  ) +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())
save_pub(
  p_confidence,
  file.path(figure_dir, "56-confidence-contract"),
  width = 205, height = 145
)

negative <- read_tsv(
  file.path(summary_dir, "negative-control.tsv"),
  show_col_types = FALSE
)
counts <- negative %>%
  filter(Metric %in% c(
    "Eukaryotic viruses tested",
    "Erroneous prokaryotic-host predictions"
  )) %>%
  mutate(
    Label = recode(
      Metric,
      "Eukaryotic viruses tested" = "Euk-virus controls",
      "Erroneous prokaryotic-host predictions" = "False prokaryotic-host calls"
    ),
    Label = factor(Label, levels = rev(c("Euk-virus controls", "False prokaryotic-host calls"))),
    FillKey = if_else(grepl("False", Label), "Orange", "Gray")
  )
p_negative_count <- ggplot(counts, aes(Label, Count, fill = FillKey)) +
  geom_col(width = 0.66) +
  geom_text(
    aes(label = if_else(
      grepl("False", Label),
      paste0(comma(Count), " (12.5%)"),
      comma(Count)
    )),
    hjust = -0.12, fontface = "bold", size = 3.8
  ) +
  coord_flip() +
  scale_fill_pub(guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "A. Domain-mismatched negative control",
    x = NULL, y = "Viral genomes"
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

diagnostics <- negative %>%
  filter(Metric %in% c(
    "Errors originating from k-mer comparison",
    "Errors with iPHoP score below 90",
    "Riboviria among errors",
    "Monodnaviria among errors"
  )) %>%
  mutate(
    Label = recode(
      Metric,
      "Errors originating from k-mer comparison" = "k-mer-derived",
      "Errors with iPHoP score below 90" = "Score below 90",
      "Riboviria among errors" = "Riboviria",
      "Monodnaviria among errors" = "Monodnaviria"
    ),
    Label = factor(Label, levels = rev(c(
      "k-mer-derived", "Score below 90", "Riboviria", "Monodnaviria"
    ))),
    FillKey = if_else(Label %in% c("k-mer-derived", "Score below 90"), "Orange", "Purple")
  )
p_negative_diag <- ggplot(diagnostics, aes(Label, Percent, fill = FillKey)) +
  geom_col(width = 0.66) +
  geom_text(
    aes(label = paste0(round(Percent, 1), "%")),
    hjust = -0.12, fontface = "bold", size = 3.7
  ) +
  coord_flip() +
  scale_fill_pub(guide = "none") +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100)) +
  labs(
    title = "B. Diagnostics among the 1,018 errors",
    subtitle = "Bars overlap and must not be summed",
    x = NULL, y = "Share of erroneous predictions (%)"
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

p_negative <- p_negative_count / p_negative_diag +
  plot_annotation(
    title = "A negative control reveals the main failure mode",
    caption = paste0(
      "RefSeq r214 eukaryotic viruses should not receive bacterial or archaeal hosts.\n",
      "Domain screening and score >=90 remove many, but not all, risks."
    )
  )
save_pub(
  p_negative,
  file.path(figure_dir, "56-negative-control"),
  width = 225, height = 180
)

ceiling <- read_tsv(
  file.path(summary_dir, "claim-ceiling.tsv"),
  show_col_types = FALSE
)
evidence_levels <- c(
  "Integrated prophage + host flanks",
  "CRISPR spacer match",
  "Long-read / Hi-C linkage",
  "iPHoP score >=90",
  "k-mer / composition alone",
  "Co-abundance alone"
)
claim_levels <- c(
  "Host domain",
  "Host family",
  "Host genus",
  "Host species/strain",
  "Active infection now",
  "Causal ecological effect"
)
ceiling <- ceiling %>%
  mutate(
    EvidenceLabel = factor(wrap_text(Evidence, 24), levels = rev(wrap_text(evidence_levels, 24))),
    ClaimLabel = factor(wrap_text(Claim, 18), levels = wrap_text(claim_levels, 18)),
    TileLabel = recode(Status, "Allowed" = "A", "Conditional" = "C", "Avoid" = "X")
  )
p_ceiling <- ggplot(ceiling, aes(ClaimLabel, EvidenceLabel, fill = Status)) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(aes(label = TileLabel), fontface = "bold", size = 4) +
  scale_fill_manual(
    values = c(
      "Allowed" = "#B8E0D2",
      "Conditional" = "#F6D186",
      "Avoid" = "#E7A6A1"
    ),
    breaks = c("Allowed", "Conditional", "Avoid"),
    name = "Reporting status"
  ) +
  labs(
    title = "Evidence type limits what can be claimed",
    subtitle = "A = allowed, C = conditional, X = avoid",
    x = NULL, y = NULL,
    caption = paste0(
      "Even direct physical linkage does not by itself prove active infection or causality.\n",
      "Taxonomic resolution also depends on the linked host sequence."
    )
  ) +
  theme_pub() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9.5),
    legend.position = "bottom"
  )
save_pub(
  p_ceiling,
  file.path(figure_dir, "56-claim-ceiling"),
  width = 245, height = 175
)

write_tsv(
  tibble(
    Package = c("R", "ggplot2", "dplyr", "readr", "tidyr", "scales", "patchwork"),
    Version = c(
      as.character(getRversion()),
      as.character(packageVersion("ggplot2")),
      as.character(packageVersion("dplyr")),
      as.character(packageVersion("readr")),
      as.character(packageVersion("tidyr")),
      as.character(packageVersion("scales")),
      as.character(packageVersion("patchwork"))
    )
  ),
  file.path(summary_dir, "plot-software-versions.tsv")
)
