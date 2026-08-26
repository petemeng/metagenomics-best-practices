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
  stop("Usage: plot_article57_plasmids.R <summary-dir> <figure-dir>")
}
summary_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- args[[2]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260757)

pal_pub <- c(
  Blue = "#0072B2", Orange = "#D55E00", Green = "#009E73",
  Sky = "#56B4E9", Yellow = "#E69F00", Purple = "#CC79A7",
  Gray = "#7A7A7A", Light = "#E9EEF2", Dark = "#253238"
)

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#ECEFF1", linewidth = 0.3),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      strip.text = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "#4D4D4D"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_pub <- function(plot, stem, width = 190, height = 130) {
  base <- file.path(figure_dir, stem)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = 350, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = 350, compression = "lzw", bg = "white")
}

# Figure 1: reference-labelled replicon audit
confusion <- read_tsv(file.path(summary_dir, "reference-confusion.tsv"), show_col_types = FALSE) %>%
  mutate(
    ReferenceLabel = factor(ReferenceLabel, levels = c("Other replicon", "Plasmid")),
    geNomadCall = factor(geNomadCall, levels = c("Not plasmid", "Plasmid")),
    Cell = case_when(
      ReferenceLabel == "Plasmid" & geNomadCall == "Plasmid" ~ "Concordant",
      ReferenceLabel == "Other replicon" & geNomadCall == "Not plasmid" ~ "Concordant",
      TRUE ~ "Discordant"
    )
  )
metrics <- read_tsv(file.path(summary_dir, "reference-metrics.tsv"), show_col_types = FALSE) %>%
  mutate(Metric = factor(Metric, levels = rev(c("Sensitivity", "Precision", "Specificity"))))

p_conf <- ggplot(confusion, aes(geNomadCall, ReferenceLabel, fill = Cell)) +
  geom_tile(color = "white", linewidth = 1.3) +
  geom_text(aes(label = comma(Count)), size = 5, fontface = "bold") +
  scale_fill_manual(values = c(Concordant = pal_pub[["Green"]], Discordant = pal_pub[["Orange"]])) +
  labs(x = "geNomad call", y = "GenBank sequence label", title = "Confusion matrix", fill = NULL) +
  coord_fixed() + theme_pub() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

p_metric <- ggplot(metrics, aes(Estimate, Metric)) +
  geom_errorbarh(aes(xmin = Wilson95Lower, xmax = Wilson95Upper), height = 0.14,
                 linewidth = 0.7, color = pal_pub[["Dark"]]) +
  geom_point(size = 3.3, color = pal_pub[["Blue"]]) +
  geom_text(aes(label = percent(Estimate, accuracy = 0.1)), nudge_y = 0.22,
            size = 3.6, fontface = "bold") +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1.03), breaks = seq(0, 1, 0.25)) +
  labs(x = "Estimate with Wilson 95% CI", y = NULL,
       title = "Wilson intervals") +
  theme_pub() + theme(panel.grid.major.y = element_blank())

reference_audit <- p_conf + p_metric +
  plot_layout(widths = c(0.9, 1.25)) +
  plot_annotation(
    title = "Length-matched reference audit",
    subtitle = "Performance is conditional on 43 labelled plasmids and 43 matched other sequences"
  )
save_pub(reference_audit,
         "57-reference-plasmid-audit", 205, 112)

# Figure 2: metagenome candidate evidence
candidates <- read_tsv(file.path(summary_dir, "coassembly-plasmid-candidates.tsv"), show_col_types = FALSE) %>%
  mutate(
    EvidenceContext = factor(
      EvidenceContext,
      levels = c(
        "ARG + plasmid call + conjugation marker",
        "ARG + plasmid call",
        "Plasmid call + conjugation marker",
        "Plasmid call only"
      )
    ),
    ReferenceSupport = factor(
      ReferenceSupport,
      levels = c(
        "Reference-plasmid supported",
        "Complete-cellular-reference conflict",
        "Other-reference-label conflict",
        "No high-coverage reference support"
      )
    )
  )
context_colors <- c(
  "ARG + plasmid call + conjugation marker" = pal_pub[["Purple"]],
  "ARG + plasmid call" = pal_pub[["Orange"]],
  "Plasmid call + conjugation marker" = pal_pub[["Blue"]],
  "Plasmid call only" = pal_pub[["Gray"]]
)
support_shapes <- c(
  "Reference-plasmid supported" = 16,
  "Complete-cellular-reference conflict" = 17,
  "Other-reference-label conflict" = 15,
  "No high-coverage reference support" = 1
)

p_candidate <- ggplot(candidates, aes(LengthBp, PlasmidScore, color = EvidenceContext,
                                       shape = ReferenceSupport)) +
  geom_hline(yintercept = 0.7, linetype = 2, color = "#666666") +
  geom_point(size = 2.8, alpha = 0.88, stroke = 0.8) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(limits = c(0.68, 1.005), breaks = c(0.7, 0.8, 0.9, 1.0)) +
  scale_color_manual(
    values = context_colors,
    labels = c(
      "ARG + plasmid call + conjugation marker" = "ARG + transfer marker",
      "ARG + plasmid call" = "ARG detected",
      "Plasmid call + conjugation marker" = "Transfer marker present",
      "Plasmid call only" = "Plasmid call only"
    ),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = support_shapes,
    labels = c(
      "Reference-plasmid supported" = "Reference plasmid",
      "Complete-cellular-reference conflict" = "Complete cellular conflict",
      "Other-reference-label conflict" = "Other reference conflict",
      "No high-coverage reference support" = "No high-coverage match"
    ),
    drop = TRUE
  ) +
  labs(
    x = "Candidate length (bp, log scale)", y = "geNomad plasmid score",
    color = "Evidence on the contig", shape = "Exact-reference audit",
    title = "A plasmid score is only the first evidence layer",
    subtitle = "Independent CARD calls and transfer-related markers determine the claim ceiling"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2))
save_pub(p_candidate, "57-coassembly-plasmid-evidence", 205, 145)

# Figure 3: claimable context for all primary CARD calls
arg_context <- read_tsv(file.path(summary_dir, "arg-context-summary.tsv"), show_col_types = FALSE) %>%
  mutate(
    EvidenceContext = factor(
      EvidenceContext,
      levels = rev(c(
        "ARG + plasmid call + conjugation marker",
        "ARG + plasmid call",
        "ARG on non-plasmid contig"
      ))
    ),
    Label = paste0(ARGCalls, "  (", percent(Fraction, accuracy = 0.1), ")")
  )
p_arg <- ggplot(arg_context, aes(ARGCalls, EvidenceContext, fill = EvidenceContext)) +
  geom_col(width = 0.66, show.legend = FALSE) +
  geom_text(aes(label = Label), hjust = -0.08, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(
    "ARG + plasmid call + conjugation marker" = pal_pub[["Purple"]],
    "ARG + plasmid call" = pal_pub[["Orange"]],
    "ARG on non-plasmid contig" = pal_pub[["Gray"]]
  ), drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.40)), breaks = pretty_breaks()) +
  labs(x = "CARD Perfect/Strict calls", y = NULL,
       title = "None of 34 primary ARG calls linked to a plasmid candidate",
       subtitle = "Counts refer to gene calls; co-location is not observed transfer") +
  theme_pub() + theme(panel.grid.major.y = element_blank())
save_pub(p_arg, "57-arg-mobility-context", 195, 108)

# Figure 4: real USA300 positive-control plasmids
usa <- read_tsv(file.path(summary_dir, "usa300-replicon-audit.tsv"), show_col_types = FALSE) %>%
  filter(ReferenceLabel == "Plasmid") %>%
  arrange(LengthBp) %>%
  mutate(DisplayName = factor(DisplayName, levels = DisplayName))
features <- read_tsv(file.path(summary_dir, "usa300-plasmid-features.tsv"), show_col_types = FALSE) %>%
  filter(FeatureType != "None") %>%
  mutate(
    PositionPercent = 100 * (as.numeric(Start) + as.numeric(End)) / 2 / as.numeric(RepliconLengthBp),
    DisplayName = factor(DisplayName, levels = levels(usa$DisplayName)),
    FeatureDisplay = case_when(
      grepl("mupA", Feature, fixed = TRUE) ~ "mupA",
      TRUE ~ Feature
    )
  )

p_usa_score <- ggplot(usa, aes(PlasmidScore, DisplayName)) +
  geom_segment(aes(x = 0.7, xend = PlasmidScore, yend = DisplayName), linewidth = 2.8,
               color = pal_pub[["Light"]]) +
  geom_point(aes(size = CARDPrimaryARGCount), color = pal_pub[["Blue"]], alpha = 0.95) +
  scale_size_continuous(range = c(3, 8), breaks = 0:3) +
  scale_x_continuous(limits = c(0.68, 1.01), breaks = c(0.7, 0.8, 0.9, 1.0)) +
  labs(x = "geNomad plasmid score", y = NULL, size = "CARD ARG calls",
       title = "Classification and ARG calls") +
  theme_pub() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())

p_usa_map <- ggplot() +
  geom_segment(data = usa, aes(x = 0, xend = 100, y = DisplayName, yend = DisplayName),
               linewidth = 4, color = pal_pub[["Light"]], lineend = "round") +
  geom_point(data = features, aes(PositionPercent, DisplayName, color = FeatureType), size = 3.3) +
  geom_text(
    data = features %>% filter(FeatureType == "CARD Perfect/Strict ARG"),
    aes(PositionPercent, DisplayName, label = FeatureDisplay), angle = 35, hjust = -0.05,
    vjust = -0.45, size = 3.1, check_overlap = TRUE
  ) +
  scale_color_manual(values = c(
    "CARD Perfect/Strict ARG" = pal_pub[["Orange"]],
    "Conjugation-related marker" = pal_pub[["Purple"]]
  )) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105), breaks = seq(0, 100, 25)) +
  labs(x = "Position along reference plasmid", y = NULL, color = NULL,
       title = "Feature coordinates") +
  theme_pub() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())

usa_panel <- p_usa_score + p_usa_map +
  plot_layout(widths = c(0.9, 1.25), guides = "collect") +
  plot_annotation(title = "USA300 reference plasmid positive control") &
  theme(legend.position = "bottom")
save_pub(usa_panel,
         "57-usa300-positive-control", 215, 130)

# Figure 5: evidence ladder and claim ceiling
ladder <- read_tsv(file.path(summary_dir, "mobility-evidence-ladder.tsv"), show_col_types = FALSE) %>%
  mutate(
    Rank = factor(Rank, levels = rev(sort(unique(Rank)))),
    EvidenceShort = recode(
      Evidence,
      "Predicted plasmid + ARG + transfer marker" = "Plasmid + ARG + transfer marker",
      "Predicted plasmid + ARG" = "Plasmid + ARG"
    ),
    EvidenceLabel = paste0("Level ", Rank, " · ", EvidenceShort),
    ClaimWrapped = vapply(MaximumClaim, function(x) paste(strwrap(x, width = 34), collapse = "\n"), character(1)),
    Strength = 7 - as.numeric(as.character(Rank))
  )
p_ladder <- ggplot(ladder, aes(Strength, Rank)) +
  geom_segment(aes(x = 1, xend = Strength, yend = Rank), linewidth = 4.5,
               color = pal_pub[["Light"]], lineend = "round") +
  geom_point(aes(color = Strength), size = 5) +
  geom_text(aes(x = 0.72, label = EvidenceLabel), hjust = 1, size = 3.6, fontface = "bold") +
  geom_text(aes(x = 6.35, label = ClaimWrapped), hjust = 0, size = 3.3, lineheight = 0.95) +
  scale_color_gradient(low = pal_pub[["Gray"]], high = pal_pub[["Green"]], guide = "none") +
  scale_x_continuous(limits = c(-11.2, 13.2), breaks = NULL) +
  labs(x = NULL, y = NULL, title = "Evidence strength sets the maximum defensible mobility claim",
       subtitle = "Sequence co-location prioritizes candidates; only transfer assays demonstrate transfer") +
  theme_pub() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank())
save_pub(p_ladder, "57-mobility-evidence-ladder", 220, 148)

write_tsv(
  tibble(
    Package = c("R", "ggplot2", "dplyr", "readr", "tidyr", "scales", "patchwork"),
    Version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      as.character(packageVersion("ggplot2")), as.character(packageVersion("dplyr")),
      as.character(packageVersion("readr")), as.character(packageVersion("tidyr")),
      as.character(packageVersion("scales")), as.character(packageVersion("patchwork"))
    )
  ),
  file.path(summary_dir, "plot-software-versions.tsv")
)
