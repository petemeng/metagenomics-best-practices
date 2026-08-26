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
  stop("Usage: plot_article55_virus_abundance.R SUMMARY_DIR FIGURE_DIR")
}
summary_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- args[[2]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260755)

pal_pub <- c(
  "Blue" = "#0072B2",
  "Orange" = "#D55E00",
  "Green" = "#009E73",
  "Sky" = "#56B4E9",
  "Yellow" = "#E69F00",
  "Purple" = "#CC79A7",
  "Gray" = "#999999",
  "Light" = "#E9EEF2",
  "Dark" = "#253238"
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
      plot.caption = element_text(color = "#455A64", hjust = 0)
    )
}
save_pub <- function(plot, file_base, width = 205, height = 145,
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
compact_label <- function(x) {
  vapply(x, function(value) {
    if (is.na(value)) return(NA_character_)
    if (abs(value) >= 1e9) return(paste0(format(round(value / 1e9, 1), trim = TRUE), "B"))
    if (abs(value) >= 1e6) return(paste0(format(round(value / 1e6, 1), trim = TRUE), "M"))
    if (abs(value) >= 1e3) return(paste0(format(round(value / 1e3, 1), trim = TRUE), "k"))
    format(signif(value, 3), scientific = FALSE, trim = TRUE)
  }, character(1))
}
pretty_phage <- function(x) {
  recode(
    x,
    "vB_Eco_SLUR29" = "SLUR29",
    "vB_Eco_mar002J2" = "J2",
    "vB_Vpa_sm033" = "SM033",
    "vB_Eco_mar003J3" = "J3",
    "vB_EcoS_swan01" = "SWAN",
    "vB_Eco_mar001J1" = "J1",
    "vB_Eco_mar005P1" = "P1",
    "vB_VpaS_sm032" = "SM032",
    "phix174" = "phiX174",
    .default = x
  )
}

abundance <- read_tsv(
  file.path(summary_dir, "illumina-abundance-evidence.tsv"),
  show_col_types = FALSE
) %>%
  filter(Library != "pool") %>%
  mutate(
    PhageLabel = pretty_phage(PhageID),
    LibraryLabel = factor(
      paste("Library", Library),
      levels = paste("Library", 1:3)
    ),
    Present = Present == "True"
  )
phage_order <- abundance %>%
  group_by(PhageID, PhageLabel) %>%
  summarise(InputGenomeCopies = first(InputGenomeCopies), .groups = "drop") %>%
  arrange(InputGenomeCopies) %>%
  pull(PhageLabel)
abundance <- abundance %>%
  mutate(PhageLabel = factor(PhageLabel, levels = phage_order))

p_breadth <- ggplot(abundance, aes(LibraryLabel, PhageLabel, fill = BreadthPct)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_point(
    data = filter(abundance, Present),
    shape = 21, fill = "white", color = pal_pub[["Dark"]], size = 2.3, stroke = 0.45
  ) +
  scale_fill_gradientn(
    colors = c("#F4F7F9", pal_pub[["Sky"]], pal_pub[["Blue"]]),
    values = rescale(c(0, 75, 100)), limits = c(0, 100),
    breaks = c(0, 50, 75, 100), name = "Breadth (%)"
  ) +
  labs(
    title = "Genome breadth, not depth alone, defines detection",
    subtitle = "Open circles pass mean depth >=1x and breadth >=75%",
    x = NULL, y = NULL,
    caption = "Three unamplified Illumina library replicates from the same 15-phage mock pool."
  ) +
  theme_pub() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(size = 9.2, face = "italic")
  )
save_pub(p_breadth, file.path(figure_dir, "55-breadth-heatmap"), width = 180, height = 175)

thresholds <- read_tsv(
  file.path(summary_dir, "breadth-threshold-sensitivity.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    DatasetLabel = recode(
      Dataset,
      "lib1_illumina" = "Library 1",
      "lib2_illumina" = "Library 2",
      "lib3_illumina" = "Library 3",
      "pooled_illumina" = "Pooled"
    ),
    DatasetLabel = factor(
      DatasetLabel,
      levels = c("Library 1", "Library 2", "Library 3", "Pooled")
    )
  )
p_threshold <- ggplot(
  thresholds,
  aes(BreadthThresholdPct, PhagesDetected, color = DatasetLabel)
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.8) +
  geom_vline(xintercept = 75, linetype = "dashed", color = pal_pub[["Dark"]]) +
  annotate("text", x = 75, y = 8.45, label = "Primary gate", hjust = -0.08,
           color = pal_pub[["Dark"]], size = 3.2) +
  scale_color_manual(
    values = c(
      "Library 1" = pal_pub[["Blue"]],
      "Library 2" = pal_pub[["Orange"]],
      "Library 3" = pal_pub[["Green"]],
      "Pooled" = pal_pub[["Purple"]]
    ),
    name = NULL
  ) +
  scale_x_continuous(breaks = c(50, 70, 75, 90)) +
  scale_y_continuous(breaks = 8:12, limits = c(8, 12)) +
  labs(
    title = "A. Detection counts depend on the breadth gate",
    x = "Minimum genome breadth (%)", y = "Phages detected"
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

prevalence <- read_tsv(
  file.path(summary_dir, "replicate-prevalence.tsv"),
  show_col_types = FALSE
) %>%
  count(ReplicatesPresent, name = "Phages") %>%
  complete(ReplicatesPresent = 0:3, fill = list(Phages = 0)) %>%
  mutate(FillKey = if_else(ReplicatesPresent == 3, "Green", "Orange"))
p_prevalence <- ggplot(prevalence, aes(factor(ReplicatesPresent), Phages, fill = FillKey)) +
  geom_col(width = 0.68, color = "white") +
  geom_text(aes(label = Phages), vjust = -0.35, fontface = "bold", size = 4) +
  scale_fill_pub(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "B. Technical-replicate occupancy",
    subtitle = "Not cohort prevalence; 75% breadth / 1x depth",
    x = "Libraries with detection (of 3)", y = "Phages"
  ) +
  theme_pub()

p_sensitivity <- p_threshold + p_prevalence +
  plot_layout(widths = c(1.25, 0.85)) +
  plot_annotation(
    title = "Threshold choice and replicate handling change the observed set",
    caption = "Pooling rescues coverage but removes replicate-level uncertainty."
  )
save_pub(p_sensitivity, file.path(figure_dir, "55-threshold-prevalence"), width = 225, height = 145)

depth <- read_tsv(
  file.path(summary_dir, "input-versus-observed-depth.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    PhageLabel = pretty_phage(PhageID),
    Risk = case_when(
      AbundanceInterpretation == "Duplicate-reference multi-mapping risk" ~ "Identical references",
      AbundanceInterpretation == "Same-vOTU cross-mapping risk" ~ "Same vOTU",
      ReplicatesPresent == 0 ~ "Not detected",
      TRUE ~ "Other references"
    ),
    Label = if_else(
      Risk != "Other references" | PhageID %in% c("S-RSM4", "vB_VpaS_sm032"),
      PhageLabel, ""
    )
  )
p_input_depth <- ggplot(
  depth,
  aes(InputGenomeCopies, MeanDepthAcrossLibrariesX, color = Risk)
) +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_text(
    aes(label = Label), hjust = -0.1, vjust = -0.4,
    size = 3.0, check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_x_log10(labels = compact_label) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.02),
                     labels = compact_label) +
  scale_color_manual(
    values = c(
      "Identical references" = pal_pub[["Orange"]],
      "Same vOTU" = pal_pub[["Purple"]],
      "Not detected" = pal_pub[["Gray"]],
      "Other references" = pal_pub[["Blue"]]
    ),
    name = NULL
  ) +
  labs(
    title = "A. Input and sequence depth diverge",
    x = "Input genome copies", y = "Mean depth across libraries (x)"
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

j12 <- abundance %>%
  filter(PhageID %in% c("vB_Eco_mar001J1", "vB_Eco_mar002J2")) %>%
  mutate(PhageLabel = pretty_phage(PhageID))
p_j12 <- ggplot(j12, aes(LibraryLabel, MeanDepthX, color = PhageLabel, group = PhageLabel)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3.0) +
  scale_color_manual(
    values = c("J1" = pal_pub[["Orange"]], "J2" = pal_pub[["Blue"]]),
    name = NULL
  ) +
  scale_y_continuous(labels = compact_label) +
  labs(
    title = "B. J1 and J2 receive the same reads",
    subtitle = "100% ANI and 100% aligned span",
    x = NULL, y = "Published mean depth (x)"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

p_multimap <- p_input_depth + p_j12 +
  plot_layout(widths = c(1.25, 0.85)) +
  plot_annotation(
    title = "Redundant references and ambiguous mapping can duplicate abundance",
    caption = "The study used BBMap minid=0.90 with ambiguous=all. Collapse the catalog to vOTUs before final mapping; never sum depths across members."
  )
save_pub(p_multimap, file.path(figure_dir, "55-abundance-multimapping"), width = 235, height = 155)

pairs <- read_tsv(
  file.path(summary_dir, "votu-pairwise-threshold-audit.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    SameVOTU = SameVOTU == "True",
    PairLabel = case_when(
      PhageA == "vB_Eco_mar001J1" & PhageB == "vB_Eco_mar002J2" ~ "J1 / J2",
      SameVOTU ~ paste(pretty_phage(PhageA), pretty_phage(PhageB), sep = " / "),
      TRUE ~ ""
    ),
    FillKey = if_else(SameVOTU, "Green", "Gray")
  )
p_votu <- ggplot(pairs, aes(ANIPct, ShorterAlignmentFractionPct, fill = FillKey)) +
  annotate("rect", xmin = 95, xmax = Inf, ymin = 85, ymax = Inf,
           fill = pal_pub[["Green"]], alpha = 0.08) +
  geom_vline(xintercept = 95, linetype = "dashed", color = pal_pub[["Dark"]]) +
  geom_hline(yintercept = 85, linetype = "dashed", color = pal_pub[["Dark"]]) +
  geom_point(shape = 21, size = 3.4, color = "white", stroke = 0.5) +
  geom_text(aes(label = PairLabel), hjust = -0.08, vjust = -0.45,
            size = 2.8, check_overlap = TRUE) +
  scale_fill_pub(guide = "none") +
  coord_cartesian(xlim = c(74, 101), ylim = c(40, 102), clip = "off") +
  labs(
    title = "A. vOTUs require two gates",
    subtitle = "Seven alignable non-self pairs in the 15-reference catalog",
    x = "Average nucleotide identity (%)",
    y = "Alignment fraction of shorter genome (%)"
  ) +
  theme_pub()

taxonomy <- read_tsv(
  file.path(summary_dir, "taxonomy-votu-ledger.tsv"),
  show_col_types = FALSE
) %>%
  count(geNomadFamily, name = "References") %>%
  mutate(
    geNomadFamily = reorder(geNomadFamily, References),
    FillKey = if_else(geNomadFamily == "Unclassified at family", "Gray", "Blue")
  )
p_taxonomy <- ggplot(taxonomy, aes(References, geNomadFamily, fill = FillKey)) +
  geom_col(width = 0.68, color = "white") +
  geom_text(aes(label = References), hjust = -0.25, fontface = "bold", size = 3.7) +
  scale_fill_pub(guide = "none") +
  scale_x_continuous(breaks = 0:4, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "B. Unresolved ranks stay explicit",
    subtitle = "geNomad DB v1.9 / ICTV MSL39",
    x = "Reference genomes", y = NULL
  ) +
  theme_pub()

p_catalog <- p_votu + p_taxonomy +
  plot_layout(widths = c(1.35, 0.85)) +
  plot_annotation(
    title = "A 15-reference catalog becomes 12 vOTUs",
    caption = "J1, J2, SWAN, and SLUR29 form one four-member vOTU. Taxonomy release and vOTU identity are separate provenance fields."
  )
save_pub(p_catalog, file.path(figure_dir, "55-votu-taxonomy"), width = 235, height = 155)

lifecycle <- read_tsv(
  file.path(summary_dir, "lifecycle-evidence-ledger.tsv"),
  show_col_types = FALSE
) %>%
  filter(ConfirmedLifecycle %in% c("Temperate", "Virulent")) %>%
  mutate(
    PhageLabel = pretty_phage(PhageID),
    Confirmed = recode(ConfirmedLifecycle, "Temperate" = "Lysogenic", "Virulent" = "Lytic")
  ) %>%
  select(PhageLabel, Confirmed, DeepPL = DeepPLPrediction, PhaTYP = PhaTYPPrediction) %>%
  pivot_longer(-PhageLabel, names_to = "Evidence", values_to = "Call") %>%
  mutate(
    Evidence = factor(Evidence, levels = c("Confirmed", "DeepPL", "PhaTYP")),
    Call = factor(Call, levels = c("Lysogenic", "Lytic"))
  )
lifecycle_order <- lifecycle %>%
  filter(Evidence == "Confirmed") %>%
  arrange(Call, PhageLabel) %>%
  pull(PhageLabel)
lifecycle <- lifecycle %>%
  mutate(PhageLabel = factor(PhageLabel, levels = rev(lifecycle_order)))
p_lifecycle <- ggplot(lifecycle, aes(Evidence, PhageLabel, fill = Call)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("Lysogenic" = pal_pub[["Green"]], "Lytic" = pal_pub[["Orange"]]),
    name = "Lifecycle class"
  ) +
  labs(
    title = "A. Predictions remain predictions",
    subtitle = "Twelve phages with experimentally supported labels",
    x = NULL, y = NULL
  ) +
  theme_pub() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

evidence_map <- read_tsv(
  file.path(summary_dir, "state-lifecycle-evidence-map.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    Row = rev(seq_len(n())),
    Observation = wrap_text(Observation, 21),
    Supports = wrap_text(Supports, 22),
    DoesNotProve = wrap_text(DoesNotProve, 22)
  )
p_state <- ggplot(evidence_map) +
  geom_tile(aes(0.12, Row), width = 0.18, height = 0.72,
            fill = pal_pub[["Sky"]], color = "white") +
  geom_text(aes(0.25, Row, label = Observation), hjust = 0,
            fontface = "bold", size = 3.0, lineheight = 0.92) +
  geom_text(aes(1.55, Row, label = Supports), hjust = 0,
            size = 2.65, lineheight = 0.92) +
  geom_text(aes(2.85, Row, label = DoesNotProve), hjust = 0,
            size = 2.65, lineheight = 0.92) +
  annotate(
    "text", x = c(0.25, 1.55, 2.85), y = max(evidence_map$Row) + 0.72,
    label = c("Observation", "Supports", "Does not prove"),
    hjust = 0, fontface = "bold", color = pal_pub[["Dark"]], size = 3.4
  ) +
  coord_cartesian(xlim = c(0, 4.15), ylim = c(0.45, 4.95), clip = "off") +
  labs(
    title = "B. Physical state and lifecycle answer different questions"
  ) +
  theme_void(base_size = 12) +
  theme(plot.title.position = "plot", plot.margin = margin(8, 12, 8, 8))

p_life <- p_lifecycle + p_state +
  plot_layout(widths = c(0.8, 1.45)) +
  plot_annotation(
    title = "Lifestyle claims require an evidence ledger",
    caption = "A free virion can come from a temperate phage; an integrated sequence does not by itself prove inducibility."
  )
save_pub(p_life, file.path(figure_dir, "55-lifecycle-evidence"), width = 250, height = 170)

write_tsv(
  tibble(
    Package = c("R", "ggplot2", "dplyr", "readr", "tidyr", "scales", "patchwork"),
    Version = c(
      paste(R.version$major, R.version$minor, sep = "."),
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

message("Article 55 figures written to ", normalizePath(figure_dir))
