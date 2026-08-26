#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    paste(
      "Usage: render_article08_figures.R",
      "<result-dir> <figure-dir>"
    ),
    call. = FALSE
  )
}

result_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- args[[2]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260719)

read_tsv <- function(filename) {
  utils::read.delim(
    file.path(result_dir, filename),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

benchmark <- read_tsv("platform-benchmark-audit.tsv")
span_survival <- read_tsv("span-survival.tsv")

stopifnot(
  nrow(benchmark) == 3L,
  nrow(span_survival) == 51L,
  all(c(
    "MeanReadLengthBp",
    "ReadLengthSDBp",
    "MeanMappedIdentityPercent",
    "AssemblyN50Bp",
    "RecoveredFullGenomes"
  ) %in% colnames(benchmark))
)

platform_levels <- c("Illumina", "ONT", "PacBio")
platform_labels <- c(
  Illumina = "Illumina HiSeq 3000",
  ONT = "ONT MinION R9",
  PacBio = "PacBio Sequel II CCS"
)
platform_colors <- c(
  Illumina = "#0072B2",
  ONT = "#D55E00",
  PacBio = "#009E73"
)

benchmark$PlatformKey <- factor(
  benchmark$PlatformKey,
  levels = platform_levels
)
span_survival$PlatformKey <- factor(
  span_survival$PlatformKey,
  levels = platform_levels
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

format_compact <- function(x) {
  labels <- ifelse(
    x >= 1e6,
    sprintf("%.1fM", x / 1e6),
    ifelse(
      x >= 1e3,
      sprintf("%.1fk", x / 1e3),
      format(
        round(x),
        scientific = FALSE,
        trim = TRUE
      )
    )
  )
  sub("\\.0([kM])$", "\\1", labels)
}

geometry_labels <- c(
  Illumina = "Illumina PE",
  ONT = "ONT R9 (historical)",
  PacBio = "PacBio CCS"
)
benchmark$GeometryLabel <- unname(
  geometry_labels[as.character(benchmark$PlatformKey)]
)
benchmark$LabelX <- c(
  Illumina = 225,
  ONT = 4100,
  PacBio = 8200
)[as.character(benchmark$PlatformKey)]
benchmark$LabelY <- c(
  Illumina = 99.15,
  ONT = 90.25,
  PacBio = 98.85
)[as.character(benchmark$PlatformKey)]

p_geometry <- ggplot2::ggplot(
  benchmark,
  ggplot2::aes(
    x = MeanReadLengthBp,
    y = MeanMappedIdentityPercent,
    color = PlatformKey
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = pmax(MeanReadLengthBp - ReadLengthSDBp, 1),
      xend = MeanReadLengthBp + ReadLengthSDBp,
      yend = MeanMappedIdentityPercent
    ),
    linewidth = 1.2,
    alpha = 0.45
  ) +
  ggplot2::geom_point(
    size = 4.8,
    alpha = 0.9
  ) +
  ggplot2::geom_label(
    ggplot2::aes(
      x = LabelX,
      y = LabelY,
      label = GeometryLabel
    ),
    size = 3.05,
    label.size = 0.22,
    fill = scales::alpha("white", 0.9),
    show.legend = FALSE
  ) +
  ggplot2::scale_x_log10(
    breaks = c(100, 500, 1000, 5000, 10000),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(86.5, 100),
    breaks = c(88, 90, 92, 94, 96, 98, 100),
    labels = scales::percent_format(scale = 1, accuracy = 1)
  ) +
  ggplot2::scale_color_manual(
    values = platform_colors,
    labels = platform_labels,
    name = "Platform"
  ) +
  ggplot2::labs(
    title = "Read length and mapped identity are separate axes",
    subtitle = "Same 71-strain MOCK1 DNA; full-run post-QC metrics",
    x = "Mean read length (bp; log scale)",
    y = "Mean end-to-end mapped identity"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = ggplot2::element_text(
      color = "grey35",
      size = 9.4
    ),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

save_pub(
  p_geometry,
  file.path(figure_dir, "08-read-geometry"),
  width = 180,
  height = 122
)

p_span <- ggplot2::ggplot(
  span_survival,
  ggplot2::aes(
    x = ThresholdBp,
    y = FractionSpanning,
    color = PlatformKey,
    group = PlatformKey
  )
) +
  ggplot2::geom_step(
    linewidth = 1.15,
    direction = "hv"
  ) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::scale_x_log10(
    breaks = c(50, 100, 150, 500, 1000, 5000, 10000, 20000),
    labels = c("50", "100", "150", "500", "1k", "5k", "10k", "20k")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::scale_color_manual(
    values = platform_colors,
    labels = platform_labels,
    name = "Platform"
  ) +
  ggplot2::labs(
    title = "Kilobase-scale span is a property of the read distribution",
    subtitle = "Deterministic 5,000-read prefixes; not full-run yield estimates",
    x = "Minimum interval spanned by one read (bp; log scale)",
    y = "Reads at or above threshold"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = ggplot2::element_text(
      color = "grey35",
      size = 9.4
    ),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

save_pub(
  p_span,
  file.path(figure_dir, "08-span-survival"),
  width = 180,
  height = 122
)

assembly_long <- rbind(
  data.frame(
    PlatformKey = benchmark$PlatformKey,
    Metric = "Assembly N50",
    Value = benchmark$AssemblyN50Bp,
    Display = format_compact(benchmark$AssemblyN50Bp),
    stringsAsFactors = FALSE
  ),
  data.frame(
    PlatformKey = benchmark$PlatformKey,
    Metric = "Recovered full genomes",
    Value = benchmark$RecoveredFullGenomes,
    Display = as.character(benchmark$RecoveredFullGenomes),
    stringsAsFactors = FALSE
  )
)
assembly_long$Metric <- factor(
  assembly_long$Metric,
  levels = c("Assembly N50", "Recovered full genomes")
)

p_assembly <- ggplot2::ggplot(
  assembly_long,
  ggplot2::aes(
    x = PlatformKey,
    y = Value,
    fill = PlatformKey
  )
) +
  ggplot2::geom_col(
    width = 0.66,
    alpha = 0.88
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = Display),
    vjust = -0.35,
    size = 3.5,
    fontface = "bold"
  ) +
  ggplot2::facet_wrap(
    ~Metric,
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      Illumina = "Illumina",
      ONT = "ONT R9",
      PacBio = "PacBio CCS"
    )
  ) +
  ggplot2::scale_y_continuous(
    labels = format_compact,
    expand = ggplot2::expansion(mult = c(0, 0.16))
  ) +
  ggplot2::scale_fill_manual(
    values = platform_colors,
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Read geometry propagates into assembly outcomes",
    subtitle = "Same MOCK1; platform-specific assemblers and unequal read budgets",
    x = NULL,
    y = "Reported value"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = ggplot2::element_text(
      color = "grey35",
      size = 9.4
    ),
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      angle = 0
    ),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.spacing = grid::unit(10, "pt")
  )

save_pub(
  p_assembly,
  file.path(figure_dir, "08-assembly-impact"),
  width = 182,
  height = 116
)

message("Rendered three Article 08 figures in PDF, PNG and LZW TIFF.")
