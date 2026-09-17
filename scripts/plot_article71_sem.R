#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1L]]
}

prepared_dir <- normalizePath(
  value_after("--prepared-dir", value_after("--input-dir")),
  mustWork = TRUE
)
analysis_dir <- normalizePath(
  value_after("--analysis-dir", value_after("--input-dir")),
  mustWork = TRUE
)
figure_dir <- value_after("--figure-dir")
if (is.null(figure_dir)) stop("Missing argument: --figure-dir")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- normalizePath(figure_dir, mustWork = TRUE)

set.seed(20260771)

pal_pub <- c(
  "Exposure" = "#D55E00",
  "Microbiome" = "#2A9D8F",
  "Phenotype" = "#4C78A8",
  "CD" = "#E45756",
  "UC" = "#59A14F",
  "Control" = "#4C78A8",
  "Neutral" = "#9AA5B1",
  "Warning" = "#C43C39"
)

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E9ECEF", linewidth = 0.3),
      axis.text = element_text(color = "#263238"),
      axis.title = element_text(color = "#263238"),
      legend.key = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "#EDF2F4", color = NA),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(color = "#546E7A"),
      plot.caption = element_text(color = "#546E7A", hjust = 0),
      plot.title.position = "plot"
    )
}

panel_label <- function(plot, label) {
  plot + labs(tag = label) + theme(plot.tag = element_text(face = "bold", size = 12))
}

save_pub <- function(plot, stem, width = 10.5, height = 5.8) {
  base <- file.path(figure_dir, stem)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".png"), plot, width = width, height = height,
         dpi = 360, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height,
         dpi = 360, compression = "lzw", bg = "white")
}

read_prepared <- function(name) {
  read.delim(file.path(prepared_dir, name), check.names = FALSE)
}

read_analysis <- function(name) {
  read.delim(file.path(analysis_dir, name), check.names = FALSE)
}

empty_canvas <- function(xlim = c(0, 10), ylim = c(0, 10)) {
  ggplot() +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    theme_void()
}

box_layer <- function(plot, xmin, xmax, ymin, ymax, label, fill, color = "#607D8B") {
  plot +
    annotate(
      "rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = fill, color = color, linewidth = 0.7
    ) +
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, size = 3.5, lineheight = 0.95)
}

arrow_layer <- function(plot, x, xend, y, yend, color = "#607D8B") {
  plot + annotate(
    "segment", x = x, xend = xend, y = y, yend = yend,
    color = color, linewidth = 0.7,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  )
}

# Figure 71.1: original study-data schematic; no publisher artwork is copied.
study_anchor <- function() {
  p <- empty_canvas()
  p <- box_layer(p, 0.3, 2.4, 6.2, 8.0, "PRISM\n155 participants", "#EAF2F8")
  p <- box_layer(p, 0.3, 2.4, 2.2, 4.0, "Validation\n65 participants", "#EAF2F8")
  p <- box_layer(p, 3.5, 6.2, 6.2, 8.0, "Pinned metadata.tsv\nclinical variables", "#F4F6F7")
  p <- box_layer(p, 3.5, 6.2, 2.2, 4.0, "Pinned genera.tsv\n11,720 genus columns", "#F4F6F7")
  p <- box_layer(p, 7.2, 9.7, 6.9, 8.4, "Antibiotic\nYes / No", "#FDEBD0", pal_pub[["Exposure"]])
  p <- box_layer(p, 7.2, 9.7, 4.3, 5.8, "Shannon entropy\ncomplete profile", "#D5F5E3", pal_pub[["Microbiome"]])
  p <- box_layer(p, 7.2, 9.7, 1.7, 3.2, "Fecal calprotectin\nlog1p transformed", "#D6EAF8", pal_pub[["Phenotype"]])
  p <- arrow_layer(p, 2.4, 3.5, 7.1, 7.1)
  p <- arrow_layer(p, 2.4, 3.5, 3.1, 3.1)
  p <- arrow_layer(p, 6.2, 7.2, 7.1, 7.6)
  p <- arrow_layer(p, 6.2, 7.2, 3.1, 5.0)
  p <- arrow_layer(p, 6.2, 7.2, 7.1, 2.5)
  p <- p +
    annotate("text", x = 5.0, y = 9.2,
             label = "Public processed tables → tutorial analysis nodes",
             fontface = "bold", size = 4.4) +
    annotate("text", x = 5.0, y = 0.55,
             label = "Original schematic based on Franzosa et al. study design; publisher artwork not reproduced",
             color = "#546E7A", size = 3.1)
  save_pub(p, "71-study-data-anchor", width = 10.5, height = 5.4)
}

data_positivity <- function() {
  attrition <- read_prepared("sample-attrition.tsv")
  overlap <- read_analysis("antibiotic-overlap-by-diagnosis.tsv")
  comparison <- read_prepared("prism-complete-case-comparison.tsv")

  attrition$Stage <- factor(attrition$Stage, levels = rev(attrition$Stage))
  p1 <- ggplot(attrition, aes(Subjects, Stage)) +
    geom_col(fill = "#8CB6C6", width = 0.62) +
    geom_text(aes(label = Subjects), hjust = -0.25, size = 3.2) +
    coord_cartesian(xlim = c(0, 235), clip = "off") +
    labs(title = "Complete-case losses are visible", x = "Participants", y = NULL) +
    theme_pub()

  overlap_long <- rbind(
    data.frame(Diagnosis = overlap$Diagnosis, Exposure = "Unexposed", N = overlap$Unexposed),
    data.frame(Diagnosis = overlap$Diagnosis, Exposure = "Exposed", N = overlap$Exposed)
  )
  overlap_long$Diagnosis <- factor(
    overlap_long$Diagnosis,
    levels = c("Control", "CD", "UC")
  )
  p2 <- ggplot(overlap_long, aes(Diagnosis, N, fill = Exposure)) +
    geom_col(width = 0.66) +
    geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 3.1) +
    scale_fill_manual(values = c("Unexposed" = "#B0BEC5", "Exposed" = pal_pub[["Exposure"]])) +
    labs(title = "No exposed control participant", x = NULL, y = "Participants") +
    theme_pub()

  comparison$Variable <- factor(comparison$Variable, levels = rev(comparison$Variable))
  comparison$AbsSMD <- abs(comparison$StandardizedDifference)
  comparison$MissingLabel <- sprintf("Excluded missing %.0f%%", comparison$ExcludedMissingPct)
  p3 <- ggplot(comparison, aes(AbsSMD, Variable)) +
    geom_vline(xintercept = 0.2, linetype = 2, color = pal_pub[["Warning"]]) +
    geom_segment(aes(x = 0, xend = AbsSMD, yend = Variable), color = "#B0BEC5") +
    geom_point(aes(color = AbsSMD >= 0.2), size = 2.6) +
    geom_text(aes(label = MissingLabel), x = 0.58, hjust = 1, size = 2.55, color = "#607D8B") +
    scale_color_manual(values = c("TRUE" = pal_pub[["Warning"]], "FALSE" = pal_pub[["Phenotype"]]), guide = "none") +
    coord_cartesian(xlim = c(0, 0.62), clip = "off") +
    labs(title = "Included and excluded PRISM participants differ", x = "Absolute standardized difference", y = NULL) +
    theme_pub(base_size = 9.5)

  combined <- panel_label(p1, "A") | panel_label(p2, "B") | panel_label(p3, "C")
  save_pub(combined + plot_layout(widths = c(1.05, 0.85, 1.35)),
           "71-data-positivity", width = 13.2, height = 5.6)
}

prespecified_dag <- function() {
  p1 <- empty_canvas()
  p1 <- box_layer(p1, 0.4, 2.5, 5.4, 7.1, "Antibiotic\nmetadata", "#FDEBD0", pal_pub[["Exposure"]])
  p1 <- box_layer(p1, 3.9, 6.1, 5.4, 7.1, "Shannon\ndiversity", "#D5F5E3", pal_pub[["Microbiome"]])
  p1 <- box_layer(p1, 7.5, 9.7, 5.4, 7.1, "log1p fecal\ncalprotectin", "#D6EAF8", pal_pub[["Phenotype"]])
  p1 <- box_layer(p1, 3.7, 6.3, 1.2, 3.1, "Diagnosis + age +\nthree medications", "#ECEFF1")
  p1 <- arrow_layer(p1, 2.5, 3.9, 6.25, 6.25, pal_pub[["Exposure"]])
  p1 <- arrow_layer(p1, 6.1, 7.5, 6.25, 6.25, pal_pub[["Microbiome"]])
  p1 <- arrow_layer(p1, 2.4, 7.5, 5.65, 5.65, pal_pub[["Exposure"]])
  p1 <- arrow_layer(p1, 4.5, 1.5, 3.1, 5.4)
  p1 <- arrow_layer(p1, 5.0, 5.0, 3.1, 5.4)
  p1 <- arrow_layer(p1, 5.5, 8.6, 3.1, 5.4)
  p1 <- p1 + annotate("text", 5, 8.5, label = "Prespecified conditional-association graph", fontface = "bold", size = 4.2)

  timeline <- data.frame(
    Node = c(
      "Treatment\ninitiation", "Antibiotic\nmetadata",
      "Stool genus\nprofile", "Fecal\ncalprotectin"
    ),
    X = c(1.0, 3.8, 6.4, 9.0),
    Status = c("unknown", "recorded", "recorded", "recorded")
  )
  p2 <- ggplot(timeline, aes(X, 0)) +
    geom_hline(yintercept = 0, color = "#B0BEC5", linewidth = 0.7) +
    geom_point(aes(color = Status), size = 3.4) +
    geom_text(aes(label = Node), y = 0.34, size = 2.9, lineheight = 0.9) +
    geom_text(aes(label = Status), y = -0.28, size = 2.8, color = "#607D8B") +
    scale_color_manual(values = c("unknown" = "#E0A800", "recorded" = "#90A4AE"), guide = "none") +
    coord_cartesian(xlim = c(0, 10), ylim = c(-0.8, 0.9), clip = "off") +
    labs(title = "A path arrow is not a time stamp") +
    theme_void() + theme(plot.title = element_text(face = "bold"))

  save_pub(panel_label(p1, "A") | panel_label(p2, "B"),
           "71-prespecified-dag", width = 11.2, height = 4.8)
}

local_paths <- function() {
  effects <- read_analysis("path-effect-summary.tsv")
  effects <- subset(effects, Model == "Primary Shannon path")
  get_effect <- function(name) effects[effects$Effect == name, ]
  a <- get_effect("A"); b <- get_effect("B"); direct <- get_effect("Direct")

  p1 <- empty_canvas()
  p1 <- box_layer(p1, 0.4, 2.5, 4.7, 6.4, "Antibiotic", "#FDEBD0", pal_pub[["Exposure"]])
  p1 <- box_layer(p1, 3.9, 6.1, 4.7, 6.4, "Shannon", "#D5F5E3", pal_pub[["Microbiome"]])
  p1 <- box_layer(p1, 7.5, 9.7, 4.7, 6.4, "Calprotectin", "#D6EAF8", pal_pub[["Phenotype"]])
  p1 <- arrow_layer(p1, 2.5, 3.9, 5.55, 5.55, pal_pub[["Exposure"]])
  p1 <- arrow_layer(p1, 6.1, 7.5, 5.55, 5.55, pal_pub[["Microbiome"]])
  p1 <- arrow_layer(p1, 2.4, 7.5, 4.9, 4.9, pal_pub[["Exposure"]])
  p1 <- p1 +
    annotate("text", 3.2, 6.2, label = sprintf("a = %.2f\n95%% boot [%.2f, %.2f]", a$Estimate, a$CILower, a$CIUpper), size = 3.0, color = pal_pub[["Exposure"]]) +
    annotate("text", 6.8, 6.2, label = sprintf("b = %.2f\n95%% boot [%.2f, %.2f]", b$Estimate, b$CILower, b$CIUpper), size = 3.0, color = pal_pub[["Microbiome"]]) +
    annotate("text", 5.0, 4.15, label = sprintf("c' = %.2f [%.2f, %.2f]", direct$Estimate, direct$CILower, direct$CIUpper), size = 3.0, color = pal_pub[["Exposure"]]) +
    annotate("text", 5.0, 8.2, label = "Only the exposure-to-diversity path is precise", fontface = "bold", size = 4.0)

  coefficients <- read_analysis("local-path-coefficients-hc3.tsv")
  coefficients <- subset(
    coefficients,
    (Model == "Microbiome node" & Term %in% c("Antibiotic", "CD", "UC")) |
      (Model == "Phenotype node" & Term %in% c("Antibiotic", "ShannonZ", "CD", "UC"))
  )
  coefficients$Label <- paste(
    ifelse(coefficients$Model == "Microbiome node", "Microbiome", "Phenotype"),
    coefficients$Term,
    sep = ": "
  )
  coefficients$Label <- factor(coefficients$Label, levels = rev(coefficients$Label))
  coefficients$Color <- ifelse(coefficients$Term == "ShannonZ", "Microbiome",
                               ifelse(coefficients$Term == "Antibiotic", "Exposure", "Neutral"))
  p2 <- ggplot(coefficients, aes(Estimate, Label, color = Color)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#90A4AE") +
    geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.18) +
    geom_point(size = 2.4) +
    scale_color_manual(values = pal_pub, guide = "none") +
    labs(title = "Local equations, one coefficient at a time", x = "HC3 coefficient", y = NULL) +
    theme_pub()

  save_pub(panel_label(p1, "A") | panel_label(p2, "B"),
           "71-local-paths", width = 11.6, height = 5.2)
}

path_decomposition <- function() {
  effects <- subset(read_analysis("path-effect-summary.tsv"), Model == "Primary Shannon path")
  effects <- subset(effects, Effect %in% c("Direct", "Indirect", "Total"))
  effects$Effect <- factor(effects$Effect, levels = c("Total", "Indirect", "Direct"))
  p1 <- ggplot(effects, aes(Estimate, Effect, color = Effect)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#90A4AE") +
    geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.18) +
    geom_point(size = 2.8) +
    scale_color_manual(values = c("Direct" = pal_pub[["Exposure"]], "Indirect" = pal_pub[["Microbiome"]], "Total" = "#B565A7"), guide = "none") +
    labs(title = "Direct and indirect estimates point in opposite directions", x = "Path effect (PRISM SD units)", y = NULL) +
    theme_pub()

  draws <- read_analysis("sem-path-bootstrap.tsv.gz")
  density_data <- stack(draws[c("Direct", "Indirect", "Total")])
  names(density_data) <- c("Value", "Effect")
  p2 <- ggplot(density_data, aes(Value, color = Effect, fill = Effect)) +
    geom_density(alpha = 0.10, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = 2, color = "#90A4AE") +
    scale_color_manual(values = c("Direct" = pal_pub[["Exposure"]], "Indirect" = pal_pub[["Microbiome"]], "Total" = "#B565A7")) +
    scale_fill_manual(values = c("Direct" = pal_pub[["Exposure"]], "Indirect" = pal_pub[["Microbiome"]], "Total" = "#B565A7")) +
    labs(title = "5,000 full-path bootstrap refits", x = "Bootstrap path effect", y = "Density", color = NULL, fill = NULL) +
    theme_pub()
  save_pub(panel_label(p1, "A") | panel_label(p2, "B"),
           "71-path-decomposition", width = 11.2, height = 4.8)
}

model_fit_dsep <- function() {
  fits <- read_analysis("sem-fit-comparison.tsv")
  fits$Short <- c("Partial", "Microbiome-only", "Reverse")
  p1 <- ggplot(fits, aes(Short, DeltaAIC, fill = Short)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = sprintf("AIC %.1f", AIC)), y = 0.06, vjust = -0.5, size = 3.0) +
    scale_fill_manual(values = c("Partial" = "#F4B183", "Microbiome-only" = pal_pub[["Microbiome"]], "Reverse" = "#B0BEC5"), guide = "none") +
    labs(title = "Direction is not selected by AIC", x = NULL, y = expression(Delta * "AIC")) +
    theme_pub()

  constrained <- fits[fits$Short == "Microbiome-only", ]
  p2 <- empty_canvas()
  p2 <- box_layer(p2, 0.8, 4.4, 4.4, 7.0, "Partial graph\n0 d-sep claims\nSaturated", "#FDEBD0")
  p2 <- box_layer(p2, 5.6, 9.2, 4.4, 7.0,
                  sprintf("Microbiome-only\n1 d-sep claim\nFisher C = %.2f\nP = %.3f", constrained$FisherC, constrained$FisherP),
                  "#D5F5E3")
  p2 <- p2 + annotate("text", 5, 8.3, label = "Global fit exists only when the graph omits a path", fontface = "bold", size = 3.9)
  save_pub(panel_label(p1, "A") | panel_label(p2, "B"),
           "71-model-fit-dsep", width = 11.2, height = 4.8)
}

overlap_influence <- function() {
  overlap <- read_analysis("antibiotic-overlap-by-diagnosis.tsv")
  p1 <- ggplot(overlap, aes(Diagnosis, Exposed, fill = Diagnosis)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = paste0(Exposed, "/", Exposed + Unexposed)), vjust = -0.35, size = 3.2) +
    scale_fill_manual(values = pal_pub, guide = "none") +
    coord_cartesian(ylim = c(0, max(overlap$Exposed) + 2), clip = "off") +
    labs(title = "Structural positivity fails in controls", x = NULL, y = "Antibiotic-exposed participants") +
    theme_pub()

  leave <- read_analysis("leave-one-out-paths.tsv")
  leave <- leave[order(leave$Indirect), ]
  leave$Index <- seq_len(nrow(leave))
  p2 <- ggplot(leave, aes(Index, Indirect, color = factor(OmittedAntibiotic))) +
    geom_hline(yintercept = median(leave$Indirect), linetype = 2, color = "#90A4AE") +
    geom_line(color = pal_pub[["Microbiome"]], linewidth = 0.7) +
    geom_point(size = 1.8) +
    scale_color_manual(values = c("0" = pal_pub[["Microbiome"]], "1" = pal_pub[["Exposure"]]), labels = c("0" = "Unexposed omitted", "1" = "Exposed omitted")) +
    labs(title = "Leave-one-out checks single-participant influence", x = "Omitted participant, ordered by estimate", y = "Indirect path", color = NULL) +
    theme_pub()
  save_pub(panel_label(p1, "A") | panel_label(p2, "B"),
           "71-overlap-influence", width = 11.4, height = 4.9)
}

transport_sensitivity <- function() {
  transport <- read_analysis("outcome-path-transport.tsv")
  transport$Label <- ifelse(grepl("PRISM", transport$CohortModel), "PRISM (n=90)", "Validation (n=38)")
  transport$Label <- factor(transport$Label, levels = rev(transport$Label))
  p1 <- ggplot(transport, aes(Estimate, Label, color = Label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#90A4AE") +
    geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.16) +
    geom_point(size = 2.6) +
    scale_color_manual(values = c("PRISM (n=90)" = pal_pub[["Exposure"]], "Validation (n=38)" = pal_pub[["Phenotype"]]), guide = "none") +
    labs(title = "The outcome path changes sign in validation", x = "Shannon-to-calprotectin coefficient", y = NULL) +
    theme_pub()

  effects <- read_analysis("path-effect-summary.tsv")
  effects <- subset(effects, Effect == "Indirect")
  effects$Label <- c("Shannon", "Faecalibacterium")
  effects$Label <- factor(effects$Label, levels = rev(effects$Label))
  p2 <- ggplot(effects, aes(Estimate, Label, color = Label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#90A4AE") +
    geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.16) +
    geom_point(size = 2.6) +
    scale_color_manual(values = c("Shannon" = pal_pub[["Microbiome"]], "Faecalibacterium" = "#B565A7"), guide = "none") +
    labs(title = "Changing the mediator changes the estimate", x = "Antibiotic-to-microbiome-to-calprotectin path", y = NULL) +
    theme_pub()
  save_pub(panel_label(p1, "A") | panel_label(p2, "B"),
           "71-transport-sensitivity", width = 11.3, height = 4.6)
}

power_positive_control <- function() {
  power <- read_analysis("path-power-simulation.tsv")
  p1 <- ggplot(power, aes(SampleSize, Power, color = Target, fill = Target)) +
    geom_hline(yintercept = 0.8, linetype = 2, color = "#607D8B") +
    geom_ribbon(aes(ymin = CILower, ymax = CIUpper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("B path (HC3)" = pal_pub[["Phenotype"]], "Joint a and b (HC3)" = pal_pub[["Microbiome"]])) +
    scale_fill_manual(values = c("B path (HC3)" = pal_pub[["Phenotype"]], "Joint a and b (HC3)" = pal_pub[["Microbiome"]])) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Observed b magnitude needs about 500 participants", subtitle = "1,000 simulations per point; exposure fraction fixed at 13/90 with overlap in every diagnosis stratum", x = "Sample size", y = "Detection probability", color = NULL, fill = NULL) +
    theme_pub()

  positive <- read_analysis("positive-control-paths.tsv")
  positive <- subset(positive, Effect %in% c("A", "B", "Indirect"))
  positive$Effect <- factor(positive$Effect, levels = c("Indirect", "B", "A"))
  p2 <- ggplot(positive, aes(Estimate, Effect)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#90A4AE") +
    geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.16, color = pal_pub[["Microbiome"]]) +
    geom_point(color = pal_pub[["Microbiome"]], size = 2.7) +
    geom_point(aes(x = TrueValue), shape = 4, stroke = 1.2, size = 3.0, color = "black") +
    labs(title = "A positive control shows what detection looks like", subtitle = "Crosses: generating truth; circles: estimate with 95% bootstrap interval", x = "Path effect", y = NULL) +
    theme_pub()
  save_pub(panel_label(p1, "A") | panel_label(p2, "B"),
           "71-power-positive-control", width = 12.0, height = 5.0)
}

mediation_rho <- function() {
  sensitivity <- read_analysis("mediation-rho-sensitivity.tsv")
  p <- ggplot(sensitivity, aes(Rho, Indirect)) +
    geom_hline(yintercept = 0, linetype = 2, color = "#607D8B") +
    geom_vline(xintercept = 0, linetype = 3, color = "#90A4AE") +
    geom_ribbon(aes(ymin = CILower, ymax = CIUpper), fill = pal_pub[["Microbiome"]], alpha = 0.16) +
    geom_line(color = pal_pub[["Microbiome"]], linewidth = 0.9) +
    annotate("segment", x = -0.15, xend = -0.15, y = -0.17, yend = 0,
             color = pal_pub[["Warning"]], linetype = 2) +
    annotate("text", x = -0.15, y = -0.22, label = "Point estimate crosses 0 at rho = -0.15", hjust = 0, size = 3.2, color = pal_pub[["Warning"]]) +
    coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.8, 1.2)) +
    labs(title = "Residual-correlation sensitivity under a hypothetical causal model", subtitle = "The interval already crosses zero at rho = 0", x = expression(rho~"between mediator and outcome errors"), y = "Indirect effect") +
    theme_pub()
  save_pub(p, "71-mediation-rho-sensitivity", width = 8.4, height = 5.0)
}

study_anchor()
data_positivity()
prespecified_dag()
local_paths()
path_decomposition()
model_fit_dsep()
overlap_influence()
transport_sensitivity()
power_positive_control()
mediation_rho()

cat("figures\t", figure_dir, "\n", sep = "")
