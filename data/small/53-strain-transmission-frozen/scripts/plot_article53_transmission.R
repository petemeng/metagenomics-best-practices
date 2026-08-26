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
value_after <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) stop(paste("Missing", flag))
  args[[position + 1]]
}
work_dir <- normalizePath(value_after("--work-dir"), mustWork = TRUE)
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260753)

pair_pal <- c("Matched mother" = "#D55E00", "Time-matched other mother" = "#0072B2")
test_pal <- c("Control hold-out" = "#0072B2", "rCDI / FMT" = "#D55E00")
outcome_pal <- c("Resolved" = "#009E73", "Failed" = "#D55E00")
source_pal <- c("Donor" = "#0072B2", "Self" = "#E69F00", "Both" = "#CC79A7", "Unique" = "#999999")

pal_pub <- c("Primary" = "#0072B2", "Contrast" = "#D55E00", "Neutral" = "#999999")
scale_color_pub <- function(...) scale_color_manual(values = pal_pub, ...)
scale_fill_pub <- function(...) scale_fill_manual(values = pal_pub, ...)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.3),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}
save_pub <- function(plot, file_base, width = 200, height = 145,
                     units = "mm", dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = units, device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = units, dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = units, dpi = dpi, compression = "lzw", bg = "white")
  invisible(plot)
}

summary_dir <- file.path(work_dir, "summary")
write_tsv(
  tibble(
    Software = c("R", "ggplot2", "dplyr", "readr", "tidyr", "scales", "patchwork"),
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

# Species-level negative control: each infant is compared with every mother
# available at the same nominal time point.
pairs <- read_tsv(
  file.path(summary_dir, "mother-infant-pairwise-sharing.tsv"),
  show_col_types = FALSE
) %>%
  filter(ThresholdPctExclusive == 0.1) %>%
  mutate(
    InfantLabel = paste0("Pair ", InfantPair, " · ", TimePoint),
    Comparison = if_else(MatchedPair, "Matched mother", "Time-matched other mother")
  )
ranks <- read_tsv(
  file.path(summary_dir, "mother-infant-rank-sensitivity.tsv"),
  show_col_types = FALSE
) %>%
  filter(ThresholdPctExclusive == 0.1) %>%
  transmute(
    Infant,
    RankLabel = if_else(
      CandidateMothers > 1,
      paste0("rank ", MatchedRank, "/", CandidateMothers),
      "no negative control"
    )
  )
label_order <- pairs %>%
  distinct(Infant, InfantPair, TimePoint, InfantLabel) %>%
  arrange(InfantPair, TimePoint) %>%
  pull(InfantLabel)
pairs <- pairs %>%
  mutate(InfantLabel = factor(InfantLabel, levels = label_order)) %>%
  left_join(ranks, by = "Infant")
p_species <- ggplot(pairs, aes(InfantLabel, JaccardSimilarity)) +
  geom_point(
    data = pairs %>% filter(!MatchedPair),
    aes(color = Comparison), size = 2.2, alpha = 0.72,
    position = position_jitter(width = 0.13, height = 0, seed = 20260753)
  ) +
  geom_point(
    data = pairs %>% filter(MatchedPair),
    aes(color = Comparison), size = 3.5, shape = 18
  ) +
  geom_text(
    data = pairs %>% filter(MatchedPair),
    aes(label = RankLabel), nudge_y = 0.016, size = 3.05,
    color = "#4D4D4D", check_overlap = TRUE
  ) +
  scale_color_manual(values = pair_pal, name = NULL) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.03, 0.14))
  ) +
  labs(
    title = "Species sharing does not identify the transmitting mother",
    subtitle = "Named MetaPhlAn2 species present at >0.1%; comparisons are time-point matched",
    x = "Infant sample", y = "Species-presence Jaccard similarity",
    caption = "Only 3 of 7 infants with an unrelated-mother control rank their matched mother first."
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )
save_pub(p_species, file.path(output_dir, "53-species-sharing-negative-control"), width = 220, height = 155)

# Published marker-SNV divergence, checked against the pinned article XML.
strain <- read_tsv(
  file.path(summary_dir, "published-mother-infant-strain-evidence.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    SpeciesLabel = recode(
      Species,
      "Bifidobacterium bifidum" = "B. bifidum",
      "Coprococcus comes" = "C. comes",
      "Ruminococcus bromii" = "R. bromii"
    ),
    SpeciesLabel = factor(SpeciesLabel, levels = c("B. bifidum", "C. comes", "R. bromii"))
  )
strain_long <- strain %>%
  transmute(
    SpeciesLabel,
    Matched = IntraPairDivergencePct,
    `Closest unrelated` = ClosestOtherDivergencePct,
    ClosestOtherQualifier
  ) %>%
  pivot_longer(c(Matched, `Closest unrelated`), names_to = "Comparison", values_to = "DivergencePct") %>%
  mutate(
    Label = case_when(
      Comparison == "Closest unrelated" & ClosestOtherQualifier == "at least" ~ paste0(">=", number(DivergencePct, accuracy = 0.01), "%"),
      TRUE ~ paste0(number(DivergencePct, accuracy = 0.01), "%")
    )
  )
p_strain <- ggplot(strain_long, aes(SpeciesLabel, DivergencePct, group = SpeciesLabel)) +
  geom_line(color = "#777777", linewidth = 0.75) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 3.7) +
  geom_text(aes(label = Label, color = Comparison), nudge_x = 0.12, hjust = 0, size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = c("Matched" = "#D55E00", "Closest unrelated" = "#0072B2"), name = NULL) +
  scale_shape_manual(values = c("Matched" = 18, "Closest unrelated" = 16), name = NULL) +
  scale_y_log10(
    breaks = c(0.04, 0.07, 0.13, 0.3, 0.6, 1, 1.6),
    labels = function(x) paste0(number(x, accuracy = 0.01), "%"),
    limits = c(0.03, 2.25),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Marker-SNV evidence separates matched pairs from unrelated strains",
    subtitle = "Published StrainPhlAn divergence; PanPhlAn independently corroborated gene content",
    x = NULL, y = "Marker divergence (log scale)"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", plot.margin = margin(6, 26, 6, 6))
save_pub(p_strain, file.path(output_dir, "53-mother-infant-strain-evidence"), width = 195, height = 145)

# Negative-control performance at progressively finer taxonomic resolution.
classifier <- read_tsv(
  file.path(summary_dir, "relatedness-classifier-performance.tsv"),
  show_col_types = FALSE
) %>%
  mutate(TaxonomicLevel = factor(TaxonomicLevel, levels = c("Family", "Genus", "Species", "Strain"))) %>%
  pivot_longer(c(AUROC, AUPR), names_to = "Metric", values_to = "Performance") %>%
  mutate(Metric = factor(Metric, levels = c("AUROC", "AUPR")))
p_classifier <- ggplot(
  classifier,
  aes(TaxonomicLevel, Performance, color = TestSet, group = TestSet)
) +
  geom_hline(yintercept = 0.5, color = "#BBBBBB", linetype = "dotted", linewidth = 0.45) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.12)) +
  geom_point(size = 3.0, position = position_dodge(width = 0.12)) +
  geom_text(
    aes(label = number(Performance, accuracy = 0.01)),
    size = 2.8, nudge_y = 0.035, show.legend = FALSE
  ) +
  facet_wrap(~Metric, nrow = 1) +
  scale_color_manual(values = test_pal, name = "Test set") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Strain sharing best separates related from unrelated sample pairs",
    subtitle = "Control-trained logistic models evaluated on held-out controls and rCDI/FMT pairs",
    x = "Resolution of shared-taxon evidence", y = "Performance"
  ) +
  theme_pub() +
  theme(legend.position = "bottom")
save_pub(p_classifier, file.path(output_dir, "53-relatedness-classifier"), width = 210, height = 145)

# Case-wise recipient persistence versus donor engraftment.
fmt_long <- read_tsv(
  file.path(summary_dir, "fmt-casewise-sharing-long.tsv"),
  show_col_types = FALSE
) %>%
  group_by(Case, Resolution) %>%
  filter(n_distinct(Comparison) == 2) %>%
  ungroup() %>%
  mutate(
    Comparison = factor(Comparison, levels = c("Pre-FMT / post-FMT", "Donor / post-FMT")),
    Resolution = factor(Resolution, levels = c("Species", "Strain"))
  )
fmt_summary <- read_tsv(
  file.path(summary_dir, "fmt-casewise-sharing-summary.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    Resolution,
    Label = paste0("medians: ", number(PrePostMedian), " vs ", number(DonorPostMedian)),
    x = 1.5,
    y = if_else(Resolution == "Strain", 33, 103)
  )
p_fmt <- ggplot(fmt_long, aes(Comparison, SharedTaxa, group = Case, color = FMTOutcome)) +
  geom_line(alpha = 0.38, linewidth = 0.55) +
  geom_point(alpha = 0.78, size = 2.0) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.2) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", shape = 18, size = 3.6) +
  geom_text(
    data = fmt_summary,
    aes(x = x, y = y, label = Label), inherit.aes = FALSE,
    size = 3.15, fontface = "bold"
  ) +
  facet_wrap(~Resolution, scales = "free_y", nrow = 1) +
  scale_color_manual(values = outcome_pal, name = "FMT outcome") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
  labs(
    title = "Donor/post-FMT pairs share more taxa in most cases",
    subtitle = "Twenty-five cases with both comparisons; black diamonds and lines show medians",
    x = "Pair compared", y = "Shared taxa"
  ) +
  theme_pub() +
  theme(legend.position = "bottom")
save_pub(p_fmt, file.path(output_dir, "53-fmt-casewise-sharing"), width = 210, height = 150)

# Authors' source labels for competing-strain events. These are event counts,
# not independent-patient replicates.
sources <- read_tsv(
  file.path(summary_dir, "fmt-source-event-summary.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    Source = factor(Source, levels = c("Unique", "Both", "Self", "Donor")),
    FMTOutcome = factor(FMTOutcome, levels = c("Resolved", "Failed")),
    Label = if_else(EventFraction >= 0.055, paste0(Events, "\n", percent(EventFraction, accuracy = 0.1)), "")
  )
p_source <- ggplot(sources, aes(FMTOutcome, EventFraction, fill = Source)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.35) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3.1) +
  scale_fill_manual(values = source_pal, name = "Post-FMT source") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "Post-FMT competing-strain events retain distinct source histories",
    subtitle = "408 author-classified events: 314 in resolved cases and 94 in failed cases",
    x = "Clinical outcome", y = "Fraction of strain events",
    caption = "Event counts are descriptive and are not independent patient-level replicates."
  ) +
  theme_pub() +
  theme(legend.position = "bottom")
save_pub(p_source, file.path(output_dir, "53-fmt-source-events"), width = 185, height = 145)
