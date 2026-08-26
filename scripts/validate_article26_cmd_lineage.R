#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260725L
set.seed(primary_seed)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(args[[i]], 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("project-root", "input-dir", "output-dir", "figure-dir", "chapter")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c("ggplot2", "patchwork", "scales", "jsonlite", "digest")
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
input_dir <- normalizePath(args[["input-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
chapter_path <- normalizePath(args[["chapter"]], mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

font_cache_dir <- file.path(tempdir(), "article26-fontconfig-cache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

validation_log <- file.path(output_dir, "validation.log")
writeLines(
  c(
    "Article 26 curatedMetagenomicData lineage validation",
    paste0("StartedUTC\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Seed\t", primary_seed)
  ),
  validation_log
)
log_msg <- function(...) {
  line <- paste0(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    "\t", paste0(..., collapse = "")
  )
  cat(line, "\n")
  cat(line, "\n", file = validation_log, append = TRUE)
  invisible(line)
}

checks <- data.frame(
  Category = character(), CheckID = character(), Status = character(),
  Detail = character(), stringsAsFactors = FALSE
)
add_check <- function(category, check_id, passed, detail) {
  checks <<- rbind(
    checks,
    data.frame(
      Category = category,
      CheckID = check_id,
      Status = if (isTRUE(passed)) "PASS" else "FAIL",
      Detail = paste(detail, collapse = "; "),
      stringsAsFactors = FALSE
    )
  )
  invisible(passed)
}
sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
read_tsv <- function(name) {
  utils::read.delim(
    file.path(input_dir, name), check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
}
read_tsv_gz <- function(name) {
  utils::read.delim(
    gzfile(file.path(input_dir, name)), check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
}

verify_checksum_manifest <- function(directory) {
  manifest_path <- file.path(directory, "file-checksums.sha256")
  lines <- readLines(manifest_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  expected_files <- character()
  for (line in lines) {
    pieces <- strsplit(line, "[[:space:]]+", perl = TRUE)[[1L]]
    expected <- pieces[[1L]]
    relative <- paste(pieces[-1L], collapse = " ")
    path <- file.path(directory, relative)
    observed <- if (file.exists(path)) sha256_file(path) else "missing"
    add_check("Frozen input", paste0("sha256-", relative), identical(observed, expected), observed)
    expected_files <- c(expected_files, relative)
  }
  payloads <- sort(setdiff(basename(list.files(directory, full.names = TRUE)), "file-checksums.sha256"))
  add_check(
    "Frozen input", "checksum-manifest-complete",
    identical(payloads, sort(expected_files)),
    paste0("payloads=", length(payloads), "; entries=", length(expected_files))
  )
  length(lines)
}

checksum_entries <- verify_checksum_manifest(input_dir)
notice_path <- file.path(project_root, "data", "small", "26-data-NOTICE.txt")
notice <- readLines(notice_path, warn = FALSE)
add_check("Frozen input", "notice-present", length(notice) >= 28L, notice_path)
for (token in c(
  "22,588 samples, 93 studies and 141 columns",
  "22,710 samples from 94 cohorts",
  "177 IDs occur in two studies",
  "MetaPhlAn3 and CHOCOPhlAn 201901",
  "counts=TRUE multiplies percentage"
)) {
  add_check(
    "Frozen input", paste0("notice-", gsub("[^a-z0-9]+", "-", tolower(token))),
    any(grepl(token, notice, fixed = TRUE)), token
  )
}

contract <- read_tsv("analysis-contract.tsv")
resource_catalog <- read_tsv("resource-catalog.tsv")
release_summary <- read_tsv("resource-release-summary.tsv")
query_resolution <- read_tsv("asnicar-query-resolution.tsv")
sample_catalog <- read_tsv_gz("sample-catalog.tsv.gz")
sample_id_collisions <- read_tsv("sample-id-collisions.tsv")
metadata_completeness <- read_tsv("metadata-completeness.tsv")
attrition <- read_tsv("query-attrition.tsv")
selection_audit <- read_tsv_gz("sample-selection-audit.tsv.gz")
selected_metadata <- read_tsv_gz("selected-sample-metadata.tsv.gz")
merged_profile <- read_tsv_gz("merged-species-relative-abundance.tsv.gz")
merged_contract <- read_tsv("merged-object-contract.tsv")
counts_boundary <- read_tsv("counts-boundary.tsv")
lineage_cards <- read_tsv("lineage-cards.tsv")
compatibility <- read_tsv("merge-compatibility.tsv")
resource_manifest <- read_tsv("resource-manifest.tsv")

contract_value <- function(item) contract$Value[match(item, contract$Item)]
expected_contract <- c(
  package_version = "3.12.0",
  bioconductor_release = "3.19",
  package_git_commit = "c5711e9",
  package_snapshot_samples = "22588",
  package_snapshot_studies = "93",
  metadata_columns = "141",
  globally_colliding_sample_ids = "177",
  resource_titles = "732",
  resource_types = "6",
  study_type_latest_resources = "558",
  study_type_pairs_with_two_releases = "174",
  query_stool_profiles = "319",
  query_independent_units = "261",
  merge_function = "curatedMetagenomicData::mergeData",
  merge_assay = "relative_abundance",
  merge_unit = "percent",
  seed = as.character(primary_seed),
  paper_collection_samples = "22710",
  paper_collection_cohorts = "94"
)
add_check("Contract", "contract-row-count", nrow(contract) == 24L, nrow(contract))
for (item in names(expected_contract)) {
  observed <- contract_value(item)
  add_check("Contract", paste0("contract-", gsub("_", "-", item)), identical(observed, expected_contract[[item]]), observed)
}
add_check(
  "Contract", "paper-package-scope-separated",
  grepl("HeQ_2017: 122 samples", contract_value("paper_only_addition"), fixed = TRUE),
  contract_value("paper_only_addition")
)

expected_types <- c(
  "gene_families", "marker_abundance", "marker_presence",
  "pathway_abundance", "pathway_coverage", "relative_abundance"
)
add_check("Resource catalog", "resource-rows", nrow(resource_catalog) == 732L, nrow(resource_catalog))
add_check("Resource catalog", "resource-studies", length(unique(resource_catalog$Study)) == 93L, length(unique(resource_catalog$Study)))
add_check("Resource catalog", "resource-types", setequal(unique(resource_catalog$DataType), expected_types), paste(sort(unique(resource_catalog$DataType)), collapse = "/"))
add_check("Resource catalog", "resources-per-type", all(table(resource_catalog$DataType) == 122L), paste(table(resource_catalog$DataType), collapse = "/"))
add_check("Resource catalog", "latest-resources", sum(resource_catalog$LatestForStudyType) == 558L, sum(resource_catalog$LatestForStudyType))

pair_catalog <- unique(resource_catalog[, c("Study", "DataType", "ReleaseCountForStudyType")])
release_pair_counts <- table(pair_catalog$ReleaseCountForStudyType)
add_check("Resource catalog", "study-type-pairs", nrow(pair_catalog) == 558L, nrow(pair_catalog))
add_check("Resource catalog", "single-release-pairs", unname(release_pair_counts[["1"]]) == 384L, release_pair_counts[["1"]])
add_check("Resource catalog", "double-release-pairs", unname(release_pair_counts[["2"]]) == 174L, release_pair_counts[["2"]])
date_counts <- table(resource_catalog$ReleaseDate)
expected_date_counts <- c(
  `2021-03-31` = 498L, `2021-04-02` = 18L, `2021-10-14` = 174L,
  `2022-04-13` = 24L, `2022-10-19` = 18L
)
add_check("Resource catalog", "release-date-counts", identical(as.integer(date_counts[names(expected_date_counts)]), unname(expected_date_counts)), paste(date_counts, collapse = "/"))
add_check("Resource catalog", "release-summary-rows", nrow(release_summary) == 30L, nrow(release_summary))
add_check("Resource catalog", "asnicar-query-rows", nrow(query_resolution) == 12L, nrow(query_resolution))
add_check("Resource catalog", "asnicar-latest-six", sum(query_resolution$LatestForStudyType) == 6L, sum(query_resolution$LatestForStudyType))
add_check("Resource catalog", "asnicar-old-superseded", all(query_resolution$UnqualifiedQuerySelection[query_resolution$ReleaseDate == "2021-03-31"] == "Superseded"), "2021-03-31")
add_check("Resource catalog", "asnicar-new-selected", all(query_resolution$UnqualifiedQuerySelection[query_resolution$ReleaseDate == "2021-10-14"] == "Selected"), "2021-10-14")

add_check("Sample catalog", "sample-rows", nrow(sample_catalog) == 22588L, nrow(sample_catalog))
add_check("Sample catalog", "sample-studies", length(unique(sample_catalog$study_name)) == 93L, length(unique(sample_catalog$study_name)))
add_check("Sample catalog", "study-sample-key-unique", !anyDuplicated(sample_catalog$StudySampleKey), sum(duplicated(sample_catalog$StudySampleKey)))
add_check("Sample catalog", "raw-sample-id-collisions", length(unique(sample_catalog$sample_id[duplicated(sample_catalog$sample_id)])) == 177L, length(unique(sample_catalog$sample_id[duplicated(sample_catalog$sample_id)])))
add_check("Sample catalog", "collision-table-rows", nrow(sample_id_collisions) == 177L, nrow(sample_id_collisions))
add_check("Sample catalog", "collision-occurrences-two", all(sample_id_collisions$Occurrences == 2L), paste(range(sample_id_collisions$Occurrences), collapse = "/"))
add_check("Sample catalog", "collision-qualified", all(grepl("\\.", sample_id_collisions$QualifiedIDs)), "study-qualified IDs")
add_check("Sample catalog", "metadata-completeness-rows", nrow(metadata_completeness) == 36L, nrow(metadata_completeness))
lookup_completeness <- function(context, field) {
  hit <- metadata_completeness$Context == context & metadata_completeness$Field == field
  metadata_completeness$Completeness[hit]
}
add_check("Sample catalog", "age-completeness", abs(lookup_completeness("Package snapshot", "age") - 13553 / 22588) < 1e-15, lookup_completeness("Package snapshot", "age"))
add_check("Sample catalog", "accession-completeness", abs(lookup_completeness("Package snapshot", "NCBI_accession") - 17968 / 22588) < 1e-15, lookup_completeness("Package snapshot", "NCBI_accession"))
add_check("Sample catalog", "disease-complete", identical(lookup_completeness("Package snapshot", "disease"), 1), lookup_completeness("Package snapshot", "disease"))

expected_attrition <- c(22588L, 928L, 319L, 319L, 261L)
add_check("Query", "attrition-rows", nrow(attrition) == 5L, nrow(attrition))
add_check("Query", "attrition-counts", identical(as.integer(attrition$Samples), expected_attrition), paste(attrition$Samples, collapse = "/"))
add_check("Query", "selection-pool-rows", nrow(selection_audit) == 319L, nrow(selection_audit))
add_check("Query", "selection-representatives", sum(selection_audit$Representative) == 261L, sum(selection_audit$Representative))
add_check("Query", "selected-metadata-rows", nrow(selected_metadata) == 261L, nrow(selected_metadata))
add_check("Query", "stool-only", all(selection_audit$body_site == "stool"), paste(unique(selection_audit$body_site), collapse = "/"))
add_check("Query", "selected-subject-unique", !anyDuplicated(paste(selected_metadata$study_name, selected_metadata$subject_id, sep = "::")), "study x subject")
study_selected_counts <- table(factor(selected_metadata$study_name, levels = c("AsnicarF_2017", "HMP_2012", "ZellerG_2014")))
add_check("Query", "selected-study-counts", identical(as.integer(study_selected_counts), c(10L, 95L, 156L)), paste(study_selected_counts, collapse = "/"))
group_max <- tapply(selection_audit$NumberReadsNumeric, paste(selection_audit$study_name, selection_audit$subject_id, sep = "::"), max, na.rm = TRUE)
selected_max <- selection_audit$NumberReadsNumeric[selection_audit$Representative]
names(selected_max) <- paste(selection_audit$study_name[selection_audit$Representative], selection_audit$subject_id[selection_audit$Representative], sep = "::")
add_check("Query", "representative-highest-depth", all(selected_max[names(group_max)] == group_max), "maximum number_reads")

profile_sample_ids <- names(merged_profile)[-(1:2)]
profile_values <- data.matrix(merged_profile[, profile_sample_ids, drop = FALSE])
add_check("Merged profile", "merged-feature-rows", nrow(merged_profile) == 923L, nrow(merged_profile))
add_check("Merged profile", "merged-sample-columns", length(profile_sample_ids) == 261L, length(profile_sample_ids))
add_check("Merged profile", "merged-sample-order", identical(profile_sample_ids, selected_metadata$sample_id), paste(head(profile_sample_ids), collapse = "/"))
add_check("Merged profile", "merged-feature-unique", !anyDuplicated(merged_profile$Lineage), sum(duplicated(merged_profile$Lineage)))
add_check("Merged profile", "merged-finite", all(is.finite(profile_values)), sum(!is.finite(profile_values)))
add_check("Merged profile", "merged-nonnegative", min(profile_values) >= 0, min(profile_values))
add_check("Merged profile", "merged-percent-bound", max(profile_values) <= 100 + 1e-8, max(profile_values))
add_check("Merged profile", "object-contract-rows", nrow(merged_contract) == 4L, nrow(merged_contract))
add_check("Merged profile", "object-feature-counts", identical(as.integer(merged_contract$Features), c(298L, 740L, 652L, 923L)), paste(merged_contract$Features, collapse = "/"))
add_check("Merged profile", "object-sample-counts", identical(as.integer(merged_contract$Samples), c(10L, 95L, 156L, 261L)), paste(merged_contract$Samples, collapse = "/"))
add_check("Merged profile", "object-assay", all(merged_contract$Assay == "relative_abundance"), paste(unique(merged_contract$Assay), collapse = "/"))
add_check("Merged profile", "object-unit-percent", all(merged_contract$Unit == "percent"), paste(unique(merged_contract$Unit), collapse = "/"))
add_check("Merged profile", "column-sums-near-percent", min(colSums(profile_values)) > 98 & max(colSums(profile_values)) < 101, paste(range(colSums(profile_values)), collapse = "/"))

add_check("Counts boundary", "counts-rows", nrow(counts_boundary) == 24L, nrow(counts_boundary))
add_check("Counts boundary", "counts-not-raw", all(!counts_boundary$ExactRawFeatureCountsRecovered), paste(unique(counts_boundary$ExactRawFeatureCountsRecovered), collapse = "/"))
add_check("Counts boundary", "rounding-difference-present", any(counts_boundary$PseudoCountMinusNumberReads != 0), paste(range(counts_boundary$PseudoCountMinusNumberReads), collapse = "/"))
add_check("Counts boundary", "relative-sums-audited", all(counts_boundary$RelativeAbundanceSumPct > 94 & counts_boundary$RelativeAbundanceSumPct <= 100.001), paste(range(counts_boundary$RelativeAbundanceSumPct), collapse = "/"))
add_check("Counts boundary", "positive-features-conserved", all(counts_boundary$PositiveRelativeFeatures == counts_boundary$PositivePseudoCountFeatures), "positive features")

required_lineage_fields <- c(
  "Study", "Accession", "RawReadsAvailable", "ProfileSource", "Profiler",
  "ProfilerVersion", "DatabaseRelease", "FeatureType", "Unit", "Normalization",
  "OriginalPublication", "ReprocessedOrAuthorProvided", "ProfileFileSHA256"
)
add_check("Lineage", "lineage-card-fields", identical(names(lineage_cards), required_lineage_fields), paste(names(lineage_cards), collapse = "/"))
add_check("Lineage", "lineage-card-rows", nrow(lineage_cards) == 5L, nrow(lineage_cards))
add_check("Lineage", "lineage-hashes", all(nchar(lineage_cards$ProfileFileSHA256) == 64L), paste(nchar(lineage_cards$ProfileFileSHA256), collapse = "/"))
add_check("Lineage", "cmd-taxonomy-metaphlan3", all(grepl("3.x", lineage_cards$ProfilerVersion[lineage_cards$Profiler == "MetaPhlAn" & lineage_cards$Study != "MeslierE_2022_MOCK1"], fixed = TRUE)), "MetaPhlAn3")
add_check("Lineage", "cmd-taxonomy-database", all(lineage_cards$DatabaseRelease[lineage_cards$Profiler == "MetaPhlAn" & lineage_cards$Study != "MeslierE_2022_MOCK1"] == "CHOCOPhlAn 201901"), "CHOCOPhlAn 201901")
add_check("Lineage", "humann-exact-release-not-invented", grepl("exact releases not encoded", lineage_cards$DatabaseRelease[lineage_cards$Profiler == "HUMAnN"], fixed = TRUE), lineage_cards$DatabaseRelease[lineage_cards$Profiler == "HUMAnN"])
add_check("Lineage", "metaphlan4-separated", lineage_cards$ProfilerVersion[lineage_cards$Study == "MeslierE_2022_MOCK1"] == "4.2.5", lineage_cards$ProfilerVersion[lineage_cards$Study == "MeslierE_2022_MOCK1"])
add_check("Lineage", "compatibility-rows", nrow(compatibility) == 7L, nrow(compatibility))
add_check("Lineage", "conditional-merge-count", sum(compatibility$Verdict == "Conditional merge") == 2L, sum(compatibility$Verdict == "Conditional merge"))
add_check("Lineage", "duplicate-release-blocked", compatibility$Verdict[compatibility$Comparison == "Asnicar old + latest tax"] == "Choose one release", compatibility$Verdict[compatibility$Comparison == "Asnicar old + latest tax"])
add_check("Lineage", "tax-function-separated", compatibility$Verdict[compatibility$Comparison == "Asnicar tax + pathway"] == "Separate assays", compatibility$Verdict[compatibility$Comparison == "Asnicar tax + pathway"])
add_check("Lineage", "cmd3-cmd4-separated", compatibility$Verdict[compatibility$Comparison == "cMD3 tax + MetaPhlAn4 tax"] == "Reprocess or harmonize", compatibility$Verdict[compatibility$Comparison == "cMD3 tax + MetaPhlAn4 tax"])
add_check("Lineage", "resource-manifest-rows", nrow(resource_manifest) == 7L, nrow(resource_manifest))
add_check("Lineage", "figure-license", resource_manifest$License[resource_manifest$Resource == "Manghi et al. Figure 1"] == "CC BY 4.0", resource_manifest$License[resource_manifest$Resource == "Manghi et al. Figure 1"])

anchor_path <- file.path(figure_dir, "26-cmd3-fig1-original.png")
add_check("Anchor figure", "anchor-present", file.exists(anchor_path), anchor_path)
add_check("Anchor figure", "anchor-sha256", file.exists(anchor_path) && identical(sha256_file(anchor_path), "559c95bd7e5c35c99e853bc16bc6c8d9739a9d7555992bf58efef3fd5b77b7c1"), if (file.exists(anchor_path)) sha256_file(anchor_path) else "missing")
add_check("Anchor figure", "anchor-bytes", file.exists(anchor_path) && file.info(anchor_path)$size == 478436, if (file.exists(anchor_path)) file.info(anchor_path)$size else -1)

chapter_lines <- readLines(chapter_path, warn = FALSE)
chapter_text <- paste(chapter_lines, collapse = "\n")
required_chapter_tokens <- c(
  "draft: false", "eval: true", "freeze: auto", "expected_images: 5",
  "## 这一步对应论文里的哪张图", "## 理论：", "## 准备工作",
  "## 可复制代码", "## 审计与升级", "## 出版级美化",
  "## 常见坑", "## 这段 Methods 怎么写", "## 换成你自己的数据怎么做",
  "## 参考", "returnSamples", "mergeData", "counts = TRUE",
  "22,588", "22,710", "HeQ_2017", "MetaPhlAn 3", "HUMAnN 3",
  "CHOCOPhlAn 201901", "data lineage", "study_name", "sample_id",
  "26-resource-release-audit.png", "26-metadata-completeness.png",
  "26-query-attrition.png", "26-lineage-compatibility.png",
  "[@manghi2026cmd3]", "[@pasolli2017cmd]", "[@beghini2021biobakery]"
)
for (token in required_chapter_tokens) {
  add_check(
    "Chapter", paste0("chapter-token-", gsub("[^a-z0-9]+", "-", tolower(token))),
    grepl(token, chapter_text, fixed = TRUE), token
  )
}
for (banned in c(
  "Planned chapter", "Do not publish", "本篇可独立跑通",
  "这体现全系列", "作者代码通常长这样", "（即本文）"
)) {
  add_check(
    "Chapter", paste0("chapter-banned-", gsub("[^a-z0-9]+", "-", tolower(banned))),
    !grepl(banned, chapter_text, fixed = TRUE), banned
  )
}
add_check("Chapter", "chapter-single-source-free", !grepl('source("R/theme_pub.R")', chapter_text, fixed = TRUE), "inline plotting functions")
add_check("Chapter", "chapter-citation-anchor", grepl("10.1038/s41467-025-66888-1", chapter_text, fixed = TRUE), "cMD3 DOI")
add_check("Chapter", "chapter-original-figure-hash", grepl("559c95bd7e5c35c99e853bc16bc6c8d9739a9d7555992bf58efef3fd5b77b7c1", chapter_text, fixed = TRUE), "anchor SHA-256")

pal_pub <- c(
  blue = "#0072B2", orange = "#E69F00", green = "#009E73",
  red = "#D55E00", purple = "#7A5195", sky = "#56B4E9",
  yellow = "#F0E442", grey = "#7A7A7A", light = "#E8EEF3"
)
theme_pub <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = base_size),
      plot.caption = ggplot2::element_text(color = "grey35", size = base_size - 1, hjust = 0),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}
save_pub <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")), plot,
    device = grDevices::cairo_pdf, width = width, height = height, units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "in", dpi = 350, bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".tiff")), plot,
    device = "tiff", compression = "lzw", width = width, height = height,
    units = "in", dpi = 350, bg = "white"
  )
}

release_totals <- aggregate(Resources ~ ReleaseDate, release_summary, sum)
release_totals$StudySnapshots <- release_totals$Resources / 6
release_totals$ReleaseDate <- factor(release_totals$ReleaseDate, levels = sort(unique(release_totals$ReleaseDate)))
p_release <- ggplot2::ggplot(release_totals, ggplot2::aes(ReleaseDate, StudySnapshots)) +
  ggplot2::geom_col(fill = pal_pub[["blue"]], width = 0.72) +
  ggplot2::geom_text(ggplot2::aes(label = StudySnapshots), vjust = -0.35, fontface = "bold", size = 3.2) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(
    title = "A. Date-stamped study snapshots",
    subtitle = "Each snapshot contributes six assay resources",
    x = "Resource date", y = "Study snapshots"
  ) +
  theme_pub(10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")

multiplicity <- as.data.frame(table(pair_catalog$ReleaseCountForStudyType), stringsAsFactors = FALSE)
names(multiplicity) <- c("Releases", "StudyAssayPairs")
multiplicity$Label <- ifelse(multiplicity$Releases == "1", "One release", "Two releases")
p_multiplicity <- ggplot2::ggplot(multiplicity, ggplot2::aes(Label, StudyAssayPairs, fill = Label)) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::geom_text(ggplot2::aes(label = StudyAssayPairs), vjust = -0.4, fontface = "bold", size = 3.5) +
  ggplot2::scale_fill_manual(values = c("One release" = pal_pub[["green"]], "Two releases" = pal_pub[["orange"]])) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(
    title = "B. Automatic latest-release choice",
    subtitle = "174 study-assay pairs expose two snapshots",
    x = NULL, y = "Study-assay pairs"
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "none")

resource_plot <- (p_release + p_multiplicity) +
  patchwork::plot_annotation(
    title = "Resource discovery is not resource identity",
    subtitle = "curatedMetagenomicData 3.12.0: 732 titles = 122 study snapshots x 6 assay types",
    caption = "An unqualified query selects the latest date within each study-assay pair; record the selected title before analysis."
  )
save_pub(resource_plot, "26-resource-release-audit", 10.5, 5.2)

field_labels <- c(
  subject_id = "Subject ID", body_site = "Body site", study_condition = "Study condition",
  disease = "Disease", age = "Age", gender = "Gender", country = "Country",
  number_reads = "Read count", NCBI_accession = "NCBI accession",
  sequencing_platform = "Sequencing platform", antibiotics_current_use = "Current antibiotics",
  days_from_first_collection = "Days from first collection"
)
metadata_completeness$FieldLabel <- unname(field_labels[metadata_completeness$Field])
metadata_completeness$FieldLabel <- factor(metadata_completeness$FieldLabel, levels = rev(unname(field_labels)))
metadata_completeness$Context <- factor(
  metadata_completeness$Context,
  levels = c("Package snapshot", "Three named studies", "Selected stool profiles")
)
completeness_plot <- ggplot2::ggplot(
  metadata_completeness,
  ggplot2::aes(Context, FieldLabel, fill = Completeness)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(Completeness, accuracy = 1)),
    color = ifelse(metadata_completeness$Completeness >= 0.62, "white", "black"),
    size = 3.1, fontface = "bold"
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7FBFF", pal_pub[["sky"]], pal_pub[["blue"]]),
    limits = c(0, 1), labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Metadata completeness changes with the analytic cohort",
    subtitle = "A field can be available in the package yet sparse after a biological filter",
    x = NULL, y = NULL, fill = "Complete",
    caption = "Completeness means non-missing and non-empty; it does not guarantee correct coding or causal interpretability."
  ) +
  theme_pub(10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
save_pub(completeness_plot, "26-metadata-completeness", 9.4, 6.3)

attrition$Step <- factor(attrition$Step, levels = rev(attrition$Step))
p_attrition <- ggplot2::ggplot(attrition, ggplot2::aes(Samples, Step)) +
  ggplot2::geom_segment(ggplot2::aes(x = 1, xend = Samples, yend = Step), color = "grey75", linewidth = 1.1) +
  ggplot2::geom_point(size = 4.2, color = pal_pub[["blue"]]) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(Samples)), hjust = -0.18,
    fontface = "bold", size = 3.2
  ) +
  ggplot2::scale_x_log10(
    breaks = c(10, 100, 1000, 10000),
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0.02, 0.2))
  ) +
  ggplot2::labs(
    title = "A. Query attrition",
    subtitle = "Filters are part of the data lineage",
    x = "Samples (log scale)", y = NULL
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "none")

candidate_counts <- as.data.frame(table(selection_audit$study_name), stringsAsFactors = FALSE)
names(candidate_counts) <- c("Study", "Samples")
candidate_counts$Stage <- "Stool profiles"
representative_counts <- as.data.frame(table(selected_metadata$study_name), stringsAsFactors = FALSE)
names(representative_counts) <- c("Study", "Samples")
representative_counts$Stage <- "One per study-subject"
study_counts <- rbind(candidate_counts, representative_counts)
study_counts$Study <- factor(study_counts$Study, levels = c("AsnicarF_2017", "HMP_2012", "ZellerG_2014"))
p_study <- ggplot2::ggplot(study_counts, ggplot2::aes(Study, Samples, fill = Stage)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = Samples),
    position = ggplot2::position_dodge(width = 0.72), vjust = -0.35,
    size = 3.1, fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = c(
    "Stool profiles" = pal_pub[["sky"]],
    "One per study-subject" = pal_pub[["green"]]
  )) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
  ggplot2::labs(
    title = "B. Repeated samples are explicit",
    subtitle = "Highest read depth, then sample ID, fixes the representative",
    x = NULL, y = "Samples", fill = NULL
  ) +
  theme_pub(10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

query_plot <- (p_attrition + p_study) +
  patchwork::plot_annotation(
    title = "From 22,588 catalog rows to 261 independent stool units",
    caption = "The three-study example retains AsnicarF_2017, HMP_2012 and ZellerG_2014; the final profile has 923 aligned species features."
  )
save_pub(query_plot, "26-query-attrition", 11.2, 5.8)

criterion_columns <- c(
  "SameFeatureType", "SameProfilerLineage", "SameDatabaseRelease",
  "SameUnitAndDenominator", "DisjointBiologicalUnits"
)
criterion_labels <- c(
  SameFeatureType = "Feature type", SameProfilerLineage = "Profiler lineage",
  SameDatabaseRelease = "Database release", SameUnitAndDenominator = "Unit / denominator",
  DisjointBiologicalUnits = "Biological units"
)
compatibility_long <- do.call(
  rbind,
  lapply(seq_len(nrow(compatibility)), function(i) {
    data.frame(
      Comparison = compatibility$Comparison[[i]],
      Criterion = criterion_columns,
      Match = as.logical(unlist(compatibility[i, criterion_columns, drop = FALSE], use.names = FALSE)),
      stringsAsFactors = FALSE
    )
  })
)
compatibility_long$Criterion <- factor(
  unname(criterion_labels[compatibility_long$Criterion]),
  levels = unname(criterion_labels)
)
compatibility_long$Comparison <- factor(
  compatibility_long$Comparison,
  levels = rev(compatibility$Comparison)
)
compatibility_long$Status <- ifelse(compatibility_long$Match, "Matched", "Mismatch")
compatibility_long$Mark <- ifelse(compatibility_long$Match, "PASS", "FAIL")
compatibility_plot <- ggplot2::ggplot(
  compatibility_long,
  ggplot2::aes(Criterion, Comparison, fill = Status)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.9) +
  ggplot2::geom_text(ggplot2::aes(label = Mark), color = "white", size = 2.7, fontface = "bold") +
  ggplot2::scale_fill_manual(values = c("Matched" = pal_pub[["green"]], "Mismatch" = pal_pub[["red"]])) +
  ggplot2::labs(
    title = "A shared assay name is not a complete merge contract",
    subtitle = "Every row must be reviewed against feature, software, database, unit and biological identity",
    x = NULL, y = NULL, fill = NULL,
    caption = "Matched contracts permit alignment; they do not remove study effects. A mismatch requires separation, harmonization or reprocessing."
  ) +
  theme_pub(10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )
save_pub(compatibility_plot, "26-lineage-compatibility", 10.8, 6.2)

figure_stems <- c(
  "26-resource-release-audit", "26-metadata-completeness",
  "26-query-attrition", "26-lineage-compatibility"
)
figure_audit <- do.call(
  rbind,
  lapply(figure_stems, function(stem) {
    do.call(
      rbind,
      lapply(c("pdf", "png", "tiff"), function(extension) {
        path <- file.path(figure_dir, paste0(stem, ".", extension))
        data.frame(
          Figure = stem,
          Format = extension,
          Exists = file.exists(path),
          Bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
          SHA256 = if (file.exists(path)) sha256_file(path) else NA_character_,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
add_check("Figures", "figure-files", nrow(figure_audit) == 12L && all(figure_audit$Exists), paste0(sum(figure_audit$Exists), "/12"))
add_check("Figures", "figure-nonempty", all(figure_audit$Bytes > 10000), paste(range(figure_audit$Bytes), collapse = "/"))

resource_audit <- checks[checks$Category == "Resource catalog", , drop = FALSE]
metadata_audit <- checks[checks$Category %in% c("Sample catalog", "Query"), , drop = FALSE]
lineage_audit <- checks[checks$Category %in% c("Merged profile", "Counts boundary", "Lineage"), , drop = FALSE]
write_tsv(resource_audit, file.path(output_dir, "resource-audit.tsv"))
write_tsv(metadata_audit, file.path(output_dir, "metadata-query-audit.tsv"))
write_tsv(lineage_audit, file.path(output_dir, "lineage-merge-audit.tsv"))
write_tsv(figure_audit, file.path(output_dir, "figure-audit.tsv"))
write_tsv(checks, file.path(output_dir, "validation-checks.tsv"))

failures <- checks[checks$Status != "PASS", , drop = FALSE]
summary <- list(
  status = if (nrow(failures) == 0L) "passed" else "failed",
  article = 26L,
  seed = primary_seed,
  package_snapshot = list(samples = 22588L, studies = 93L, metadata_columns = 141L),
  paper_collection = list(samples = 22710L, cohorts = 94L, paper_only_study = "HeQ_2017", paper_only_samples = 122L),
  resource_catalog = list(titles = 732L, latest_resources = 558L, double_release_pairs = 174L),
  query = list(stool_profiles = 319L, independent_units = 261L, merged_species = 923L),
  sample_id_collisions = 177L,
  checks = nrow(checks),
  passed = sum(checks$Status == "PASS"),
  failed = nrow(failures),
  checksum_entries = checksum_entries,
  generated_figure_files = nrow(figure_audit),
  generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
jsonlite::write_json(
  summary, file.path(output_dir, "validation-summary.json"),
  auto_unbox = TRUE, pretty = TRUE
)

log_msg(
  "checks=", nrow(checks), "; passed=", sum(checks$Status == "PASS"),
  "; failed=", nrow(failures)
)
if (nrow(failures) > 0L) {
  print(failures, row.names = FALSE)
  stop("Article 26 validation failed.", call. = FALSE)
}
log_msg("Article 26 validation passed.")
cat("Article 26 validation passed:", nrow(checks), "checks.\n")
