#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    paste(
      "Usage: render_article07_figures.R",
      "<result-dir> <index-evidence.tsv> <figure-dir>"
    ),
    call. = FALSE
  )
}

result_dir <- normalizePath(args[[1]], mustWork = TRUE)
index_evidence_path <- normalizePath(args[[2]], mustWork = TRUE)
figure_dir <- args[[3]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260719)

classification <- utils::read.delim(
  file.path(result_dir, "contaminant-classification.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_burden <- utils::read.delim(
  file.path(result_dir, "sample-contamination-burden.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
index_evidence <- utils::read.delim(
  index_evidence_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

stopifnot(
  nrow(classification) == 1951L,
  nrow(sample_burden) == 569L,
  nrow(index_evidence) == 3L,
  sum(classification$GlobalCombinedContaminant) == 56L,
  all(c(
    "DisplayLabel",
    "SourceType",
    "LowerPercent",
    "PointPercent",
    "UpperPercent"
  ) %in% colnames(index_evidence))
)

pal_pub <- c(
  "#0072B2",
  "#D55E00",
  "#009E73",
  "#CC79A7",
  "#E69F00",
  "#56B4E9",
  "#F0E442",
  "#999999"
)

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(
    base_size = base_size,
    base_family = "sans"
  ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.3
      ),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.3
      ),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8),
        color = NA
      )
    )
}

save_pub <- function(
  plot,
  file_base,
  width = 174,
  height = 118,
  units = "mm",
  dpi = 350
) {
  ggplot2::ggsave(
    paste0(file_base, ".pdf"),
    plot,
    width = width,
    height = height,
    units = units,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    paste0(file_base, ".png"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".tiff"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  invisible(plot)
}

sample_burden$SampleClass <- factor(
  sample_burden$SampleClass,
  levels = c("Biological sample", "Negative control")
)
sample_colors <- c(
  "Biological sample" = pal_pub[[1]],
  "Negative control" = pal_pub[[2]]
)

p_library <- ggplot2::ggplot(
  sample_burden,
  ggplot2::aes(
    x = SampleClass,
    y = LibrarySize,
    fill = SampleClass
  )
) +
  ggplot2::geom_violin(
    width = 0.72,
    trim = FALSE,
    alpha = 0.24,
    color = NA
  ) +
  ggplot2::geom_boxplot(
    width = 0.23,
    outlier.shape = NA,
    alpha = 0.76,
    linewidth = 0.42
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.12,
      height = 0,
      seed = 20260719
    ),
    shape = 21,
    size = 1.25,
    stroke = 0.18,
    alpha = 0.45
  ) +
  ggplot2::scale_fill_manual(
    values = sample_colors,
    guide = "none"
  ) +
  ggplot2::scale_y_log10(
    breaks = c(
      100,
      300,
      1000,
      3000,
      10000,
      30000
    ),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::labs(
    title = "Negative controls remain analytically informative",
    subtitle = paste0(
      sum(sample_burden$SampleClass == "Biological sample"),
      " biological samples · ",
      sum(sample_burden$SampleClass == "Negative control"),
      " controls · 6 plates"
    ),
    x = NULL,
    y = "Library size (reads; log scale)"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = ggplot2::element_text(
      color = "grey35",
      size = 9.5
    ),
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(face = "bold")
  )

save_pub(
  p_library,
  file.path(figure_dir, "07-control-library-size"),
  width = 164,
  height = 112
)

classification$Classification <- factor(
  ifelse(
    classification$GlobalCombinedContaminant,
    "Likely contaminant",
    "Not classified"
  ),
  levels = c("Not classified", "Likely contaminant")
)
classification_colors <- c(
  "Not classified" = "#B8B8B8",
  "Likely contaminant" = pal_pub[[2]]
)
label_features <- intersect(
  c("Seq30", "Seq175", "Seq3"),
  classification$FeatureID
)

p_prevalence <- ggplot2::ggplot(
  classification,
  ggplot2::aes(
    x = NegativePrevalence,
    y = BiologicalPrevalence,
    color = Classification,
    size = TotalReads
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(alpha = 0.62) +
  ggplot2::geom_text(
    data = classification[
      classification$FeatureID %in% label_features,
      ,
      drop = FALSE
    ],
    ggplot2::aes(label = FeatureID),
    size = 3.0,
    color = "grey10",
    nudge_y = 0.035,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = classification_colors,
    name = "Combined classifier"
  ) +
  ggplot2::scale_size_continuous(
    trans = "log10",
    range = c(0.8, 4.4),
    breaks = c(10, 1000, 100000, 1000000),
    labels = scales::label_number(big.mark = ","),
    name = "Total reads"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::labs(
    title = "Controls and concentration identify a candidate branch",
    subtitle = paste0(
      sum(classification$GlobalCombinedContaminant),
      " of ",
      nrow(classification),
      " ASVs classified at threshold = 0.10"
    ),
    x = "Prevalence in negative controls",
    y = "Prevalence in biological samples"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = ggplot2::element_text(
      color = "grey35",
      size = 9.3
    ),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "right",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

save_pub(
  p_prevalence,
  file.path(figure_dir, "07-contaminant-prevalence"),
  width = 180,
  height = 128
)

biological_rows <- (
  sample_burden$SampleClass == "Biological sample"
)
burden_correlation <- stats::cor.test(
  sample_burden$ContaminantFraction[biological_rows],
  sample_burden$DNAConcentration[biological_rows],
  method = "spearman",
  exact = FALSE
)
rho_label <- sprintf(
  "Biological samples: Spearman rho = %.2f",
  unname(burden_correlation$estimate)
)

p_burden <- ggplot2::ggplot(
  sample_burden,
  ggplot2::aes(
    x = DNAConcentration,
    y = ContaminantFraction,
    color = SampleClass
  )
) +
  ggplot2::geom_point(alpha = 0.62, size = 1.8) +
  ggplot2::scale_color_manual(
    values = sample_colors,
    name = "Sample class"
  ) +
  ggplot2::scale_x_log10(
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.0001
    ),
    breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
    labels = scales::percent_format(accuracy = 0.01),
    limits = c(0, 1)
  ) +
  ggplot2::annotate(
    "label",
    x = Inf,
    y = Inf,
    label = rho_label,
    hjust = 1.04,
    vjust = 1.25,
    size = 3.1,
    label.size = 0.25,
    fill = scales::alpha("white", 0.86),
    color = "grey15"
  ) +
  ggplot2::labs(
    title = "Candidate burden rises as DNA concentration falls",
    subtitle = paste0(
      "Burden is the read fraction assigned to ",
      sum(classification$GlobalCombinedContaminant),
      " globally classified ASVs"
    ),
    x = "DNA concentration (PicoGreen intensity; log scale)",
    y = "Candidate contaminant read fraction"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = ggplot2::element_text(
      color = "grey35",
      size = 9.3
    ),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

save_pub(
  p_burden,
  file.path(figure_dir, "07-contaminant-burden"),
  width = 176,
  height = 118
)

index_evidence$DisplayLabel <- factor(
  index_evidence$DisplayLabel,
  levels = rev(index_evidence$DisplayLabel)
)
source_colors <- c(
  "Peer-reviewed study" = pal_pub[[1]],
  "Vendor guidance" = pal_pub[[2]],
  "Vendor worked example" = pal_pub[[3]]
)

p_index <- ggplot2::ggplot(
  index_evidence,
  ggplot2::aes(
    y = DisplayLabel,
    color = SourceType
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = LowerPercent,
      xend = UpperPercent,
      yend = DisplayLabel
    ),
    linewidth = 4.4,
    alpha = 0.55,
    lineend = "round",
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = index_evidence[
      !is.na(index_evidence$PointPercent),
      ,
      drop = FALSE
    ],
    ggplot2::aes(x = PointPercent),
    shape = 21,
    fill = "white",
    size = 3.4,
    stroke = 1.1
  ) +
  ggplot2::scale_color_manual(
    values = source_colors,
    labels = c(
      "Peer-reviewed study" = "Study",
      "Vendor guidance" = "Vendor guidance",
      "Vendor worked example" = "Vendor example"
    ),
    name = "Source"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 6.4),
    breaks = 0:6,
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::labs(
    title = "Reported index hopping varies by context",
    subtitle = "Ranges are not directly comparable or pass/fail thresholds",
    x = "Reported index-hopped reads",
    y = NULL
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = ggplot2::element_text(
      color = "grey35",
      size = 9.3
    ),
    axis.title.x = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 9.3
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    panel.grid.major.y = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

save_pub(
  p_index,
  file.path(figure_dir, "07-index-hopping-evidence"),
  width = 178,
  height = 112
)

message("Rendered four Article 07 figures in PDF, PNG and LZW TIFF.")
