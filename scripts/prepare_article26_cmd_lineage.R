#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, digits = 17)
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
required <- c(
  "hmp-rda", "zeller-rda", "asnicar-rds", "asnicar-pathway-rda",
  "metaphlan4-profile", "output-dir", "notice"
)
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c(
  "curatedMetagenomicData", "TreeSummarizedExperiment",
  "SummarizedExperiment", "S4Vectors", "digest"
)
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

paths <- lapply(args[required[seq_len(5L)]], normalizePath, mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
notice_path <- normalizePath(args[["notice"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(notice_path), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
load_one <- function(path, expected_name) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!identical(loaded, expected_name)) {
    stop("Unexpected object in ", basename(path), ": ", paste(loaded, collapse = ", "), call. = FALSE)
  }
  env[[expected_name]]
}
present <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}
accession_tokens <- function(x) {
  tokens <- trimws(unlist(strsplit(as.character(x[present(x)]), ";", fixed = TRUE)))
  sort(unique(tokens[nzchar(tokens)]))
}

expected_hashes <- c(
  `hmp-rda` = "f42713b0d59d15f65172ead2375c2ea2c1b6309b44cb2b29231edb190e56ab23",
  `zeller-rda` = "d8e0f3fd00b2339b1aa929197ca0869c43990ff885a04fc675e70d4aff5604b2",
  `asnicar-rds` = "2952774730bff2af9e13c9c40058320aed524dc2dc7408b44dc6e697c06564b2",
  `asnicar-pathway-rda` = "ead7c78c075fec92a7d641b731594e068b2ba2a47479151d081c338f615af121",
  `metaphlan4-profile` = "f081ac9e89439523b625dcf5712bb3305d47c0da8c0122552854895b26973c0d"
)
observed_hashes <- vapply(paths, sha256_file, character(1L))
if (!identical(unname(observed_hashes), unname(expected_hashes))) {
  mismatch <- names(expected_hashes)[observed_hashes != expected_hashes]
  stop("Source checksum mismatch: ", paste(mismatch, collapse = ", "), call. = FALSE)
}
if (!identical(as.character(utils::packageVersion("curatedMetagenomicData")), "3.12.0")) {
  stop("curatedMetagenomicData 3.12.0 is required.", call. = FALSE)
}

cmd_namespace <- asNamespace("curatedMetagenomicData")
resource_titles <- get("resourceTitles", envir = cmd_namespace)
title_parts <- strsplit(resource_titles, "\\.")
if (!all(lengths(title_parts) == 3L)) {
  stop("Every curatedMetagenomicData resource title must contain date.study.type.", call. = FALSE)
}
resource_catalog <- data.frame(
  ResourceTitle = resource_titles,
  ReleaseDate = vapply(title_parts, `[[`, character(1L), 1L),
  Study = vapply(title_parts, `[[`, character(1L), 2L),
  DataType = vapply(title_parts, `[[`, character(1L), 3L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
resource_catalog <- resource_catalog[order(
  resource_catalog$Study, resource_catalog$DataType, resource_catalog$ReleaseDate
), , drop = FALSE]
group_key <- paste(resource_catalog$Study, resource_catalog$DataType, sep = "::")
release_count <- ave(resource_catalog$ReleaseDate, group_key, FUN = function(x) length(unique(x)))
latest_release <- ave(resource_catalog$ReleaseDate, group_key, FUN = max)
resource_catalog$ReleaseCountForStudyType <- as.integer(release_count)
resource_catalog$LatestForStudyType <- resource_catalog$ReleaseDate == latest_release
resource_catalog$UnqualifiedQuerySelection <- ifelse(
  resource_catalog$LatestForStudyType, "Selected", "Superseded"
)

stopifnot(
  nrow(resource_catalog) == 732L,
  length(unique(resource_catalog$Study)) == 93L,
  length(unique(resource_catalog$DataType)) == 6L,
  all(table(resource_catalog$DataType) == 122L),
  sum(resource_catalog$LatestForStudyType) == 558L,
  sum(resource_catalog$ReleaseCountForStudyType == 2L) == 348L
)

release_summary <- aggregate(
  ResourceTitle ~ ReleaseDate + DataType,
  resource_catalog,
  length
)
names(release_summary)[names(release_summary) == "ResourceTitle"] <- "Resources"
release_summary$Studies <- release_summary$Resources
release_summary <- release_summary[order(release_summary$ReleaseDate, release_summary$DataType), ]

query_example <- resource_catalog[
  resource_catalog$Study == "AsnicarF_2017",
  c(
    "ResourceTitle", "ReleaseDate", "Study", "DataType",
    "ReleaseCountForStudyType", "LatestForStudyType", "UnqualifiedQuerySelection"
  ),
  drop = FALSE
]
stopifnot(nrow(query_example) == 12L, sum(query_example$LatestForStudyType) == 6L)

metadata_all <- curatedMetagenomicData::sampleMetadata
stopifnot(
  is.data.frame(metadata_all),
  identical(dim(metadata_all), c(22588L, 141L)),
  length(unique(metadata_all$study_name)) == 93L
)

sample_id_frequency <- table(metadata_all$sample_id)
colliding_sample_ids <- names(sample_id_frequency)[sample_id_frequency > 1L]
stopifnot(length(colliding_sample_ids) == 177L, all(sample_id_frequency[colliding_sample_ids] == 2L))
sample_id_collisions <- do.call(
  rbind,
  lapply(colliding_sample_ids, function(sample_id) {
    hit <- metadata_all[metadata_all$sample_id == sample_id, c("study_name", "sample_id"), drop = FALSE]
    data.frame(
      SampleID = sample_id,
      Occurrences = nrow(hit),
      Studies = paste(sort(unique(hit$study_name)), collapse = " | "),
      QualifiedIDs = paste(
        paste(hit$sample_id, hit$study_name, sep = "."),
        collapse = " | "
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
sample_id_collisions <- sample_id_collisions[order(sample_id_collisions$SampleID), , drop = FALSE]

catalog_fields <- c(
  "study_name", "sample_id", "subject_id", "body_site", "body_subsite",
  "study_condition", "disease", "age", "infant_age", "age_category",
  "gender", "country", "non_westernized", "sequencing_platform",
  "DNA_extraction_kit", "PMID", "number_reads", "number_bases",
  "NCBI_accession", "antibiotics_current_use", "days_from_first_collection"
)
sample_catalog <- metadata_all[, catalog_fields, drop = FALSE]
sample_catalog$StudySampleKey <- paste(sample_catalog$study_name, sample_catalog$sample_id, sep = "::")
stopifnot(!anyDuplicated(sample_catalog$StudySampleKey))

selected_studies <- c("AsnicarF_2017", "HMP_2012", "ZellerG_2014")
in_studies <- metadata_all$study_name %in% selected_studies
is_stool <- in_studies & metadata_all$body_site == "stool"
has_ids <- is_stool & present(metadata_all$NCBI_accession) & present(metadata_all$subject_id)

selection_pool <- metadata_all[has_ids, , drop = FALSE]
selection_pool$NumberReadsNumeric <- suppressWarnings(as.numeric(selection_pool$number_reads))
read_sort <- selection_pool$NumberReadsNumeric
read_sort[!is.finite(read_sort)] <- -1
selection_order <- order(
  match(selection_pool$study_name, selected_studies),
  selection_pool$subject_id,
  -read_sort,
  selection_pool$sample_id
)
selection_pool <- selection_pool[selection_order, , drop = FALSE]
selection_key <- paste(selection_pool$study_name, selection_pool$subject_id, sep = "::")
selection_pool$Representative <- !duplicated(selection_key)
selection_pool$SelectionRule <- ifelse(
  selection_pool$Representative,
  "highest finite number_reads within study x subject; lexical sample_id tie-break",
  "lower number_reads or lexical tie-break within study x subject"
)
selected_metadata <- selection_pool[selection_pool$Representative, , drop = FALSE]

stopifnot(
  sum(in_studies) == 928L,
  sum(is_stool) == 319L,
  sum(has_ids) == 319L,
  nrow(selected_metadata) == 261L,
  identical(
    as.integer(table(factor(selected_metadata$study_name, levels = selected_studies))),
    c(10L, 95L, 156L)
  ),
  !anyDuplicated(paste(selected_metadata$study_name, selected_metadata$subject_id, sep = "::"))
)

attrition <- data.frame(
  Step = factor(
    c(
      "Package snapshot", "Three named studies", "Stool profiles",
      "Traceable sample and subject IDs", "One sample per study and subject"
    ),
    levels = c(
      "Package snapshot", "Three named studies", "Stool profiles",
      "Traceable sample and subject IDs", "One sample per study and subject"
    )
  ),
  Samples = c(nrow(metadata_all), sum(in_studies), sum(is_stool), sum(has_ids), nrow(selected_metadata)),
  Studies = c(
    length(unique(metadata_all$study_name)),
    length(unique(metadata_all$study_name[in_studies])),
    length(unique(metadata_all$study_name[is_stool])),
    length(unique(metadata_all$study_name[has_ids])),
    length(unique(selected_metadata$study_name))
  ),
  IndependentStudySubjects = c(
    length(unique(paste(metadata_all$study_name, metadata_all$subject_id, sep = "::"))),
    length(unique(paste(metadata_all$study_name[in_studies], metadata_all$subject_id[in_studies], sep = "::"))),
    length(unique(paste(metadata_all$study_name[is_stool], metadata_all$subject_id[is_stool], sep = "::"))),
    length(unique(paste(metadata_all$study_name[has_ids], metadata_all$subject_id[has_ids], sep = "::"))),
    nrow(selected_metadata)
  ),
  Rule = c(
    "curatedMetagenomicData 3.12.0 sampleMetadata",
    paste(selected_studies, collapse = " | "),
    "body_site == stool",
    "non-missing NCBI_accession and subject_id",
    "maximum number_reads; lexical sample_id tie-break"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
attrition$Step <- as.character(attrition$Step)

completeness_fields <- c(
  "subject_id", "body_site", "study_condition", "disease", "age", "gender",
  "country", "number_reads", "NCBI_accession", "sequencing_platform",
  "antibiotics_current_use", "days_from_first_collection"
)
contexts <- list(
  `Package snapshot` = metadata_all,
  `Three named studies` = metadata_all[in_studies, , drop = FALSE],
  `Selected stool profiles` = selection_pool
)
metadata_completeness <- do.call(
  rbind,
  lapply(names(contexts), function(context_name) {
    x <- contexts[[context_name]]
    data.frame(
      Context = context_name,
      Field = completeness_fields,
      Complete = vapply(x[completeness_fields], function(z) sum(present(z)), integer(1L)),
      Total = nrow(x),
      Completeness = vapply(x[completeness_fields], function(z) mean(present(z)), numeric(1L)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
rownames(metadata_completeness) <- NULL

hmp_matrix <- load_one(paths[["hmp-rda"]], "2021-03-31.HMP_2012.relative_abundance")
zeller_matrix <- load_one(paths[["zeller-rda"]], "2021-03-31.ZellerG_2014.relative_abundance")
asnicar_tse <- readRDS(paths[["asnicar-rds"]])
asnicar_pathway <- load_one(
  paths[["asnicar-pathway-rda"]],
  "2021-10-14.AsnicarF_2017.pathway_abundance"
)

stopifnot(
  is.matrix(hmp_matrix), identical(dim(hmp_matrix), c(1188L, 748L)),
  is.matrix(zeller_matrix), identical(dim(zeller_matrix), c(1019L, 156L)),
  inherits(asnicar_tse, "TreeSummarizedExperiment"), identical(dim(asnicar_tse), c(298L, 24L)),
  is.matrix(asnicar_pathway), identical(dim(asnicar_pathway), c(11173L, 24L))
)

phylogenetic_tree <- get("phylogeneticTree", envir = cmd_namespace)
row_data_long <- get("rowDataLong", envir = cmd_namespace)
make_relative_tse <- function(source_matrix, study_name) {
  keep_tips <- intersect(rownames(source_matrix), phylogenetic_tree[["tip.label"]])
  study_metadata <- metadata_all[metadata_all$study_name == study_name, , drop = FALSE]
  study_metadata <- study_metadata[match(colnames(source_matrix), study_metadata$sample_id), , drop = FALSE]
  stopifnot(!anyNA(study_metadata$sample_id), identical(study_metadata$sample_id, colnames(source_matrix)))
  rownames(study_metadata) <- study_metadata$sample_id
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = S4Vectors::SimpleList(
      relative_abundance = source_matrix[keep_tips, , drop = FALSE]
    ),
    colData = S4Vectors::DataFrame(study_metadata),
    rowTree = phylogenetic_tree
  )
  keep_ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  SummarizedExperiment::rowData(tse) <- S4Vectors::DataFrame(
    row_data_long[keep_tips, keep_ranks, drop = FALSE]
  )
  tse
}

hmp_tse <- make_relative_tse(hmp_matrix, "HMP_2012")
zeller_tse <- make_relative_tse(zeller_matrix, "ZellerG_2014")
study_objects <- list(
  AsnicarF_2017 = asnicar_tse,
  HMP_2012 = hmp_tse,
  ZellerG_2014 = zeller_tse
)
for (study_name in names(study_objects)) {
  expected_ids <- selected_metadata$sample_id[selected_metadata$study_name == study_name]
  missing_ids <- setdiff(expected_ids, colnames(study_objects[[study_name]]))
  if (length(missing_ids) > 0L) {
    stop(study_name, " selected IDs missing from profile: ", paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  study_objects[[study_name]] <- study_objects[[study_name]][, expected_ids, drop = FALSE]
}

merged_tse <- curatedMetagenomicData::mergeData(study_objects)
merged_matrix <- SummarizedExperiment::assay(merged_tse, "relative_abundance")
merged_coldata <- as.data.frame(SummarizedExperiment::colData(merged_tse))
stopifnot(
  ncol(merged_tse) == 261L,
  nrow(merged_tse) > max(vapply(study_objects, nrow, integer(1L))),
  identical(colnames(merged_tse), unlist(lapply(study_objects, colnames), use.names = FALSE)),
  identical(merged_coldata$study_name, rep(names(study_objects), vapply(study_objects, ncol, integer(1L)))),
  all(is.finite(merged_matrix)),
  min(merged_matrix) >= 0,
  max(merged_matrix) <= 100 + 1e-8
)

selected_export_fields <- c(
  "study_name", "sample_id", "subject_id", "body_site", "body_subsite",
  "study_condition", "disease", "age", "infant_age", "age_category",
  "gender", "country", "sequencing_platform", "PMID", "number_reads",
  "number_bases", "NCBI_accession", "antibiotics_current_use",
  "days_from_first_collection", "NumberReadsNumeric", "Representative", "SelectionRule"
)
selection_export_fields <- c(
  "study_name", "sample_id", "subject_id", "body_site", "body_subsite",
  "study_condition", "disease", "age", "gender", "number_reads",
  "NCBI_accession", "NumberReadsNumeric", "Representative", "SelectionRule"
)

merged_profile <- data.frame(
  Lineage = rownames(merged_matrix),
  Species = as.character(SummarizedExperiment::rowData(merged_tse)$species),
  as.data.frame(merged_matrix, check.names = FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

profile_contract <- do.call(
  rbind,
  c(
    lapply(names(study_objects), function(study_name) {
      object <- study_objects[[study_name]]
      matrix <- SummarizedExperiment::assay(object, "relative_abundance")
      data.frame(
        Object = paste0(study_name, " selected relative_abundance"),
        Class = class(object)[[1L]],
        Features = nrow(object),
        Samples = ncol(object),
        Assay = "relative_abundance",
        Unit = "percent",
        MinimumColumnSum = min(colSums(matrix)),
        MedianColumnSum = stats::median(colSums(matrix)),
        MaximumColumnSum = max(colSums(matrix)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }),
    list(data.frame(
      Object = "Merged selected relative_abundance",
      Class = class(merged_tse)[[1L]],
      Features = nrow(merged_tse),
      Samples = ncol(merged_tse),
      Assay = "relative_abundance",
      Unit = "percent",
      MinimumColumnSum = min(colSums(merged_matrix)),
      MedianColumnSum = stats::median(colSums(merged_matrix)),
      MaximumColumnSum = max(colSums(merged_matrix)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  )
)

asnicar_meta <- as.data.frame(SummarizedExperiment::colData(asnicar_tse))
asnicar_ra <- SummarizedExperiment::assay(asnicar_tse, "relative_abundance")
counts_boundary <- do.call(
  rbind,
  lapply(seq_len(ncol(asnicar_ra)), function(i) {
    reads <- suppressWarnings(as.numeric(asnicar_meta$number_reads[[i]]))
    pseudo <- round(asnicar_ra[, i] * reads / 100)
    data.frame(
      Study = "AsnicarF_2017",
      SampleID = colnames(asnicar_ra)[[i]],
      NumberReads = reads,
      RelativeAbundanceSumPct = sum(asnicar_ra[, i]),
      PseudoCountSum = sum(pseudo),
      PseudoCountMinusNumberReads = sum(pseudo) - reads,
      PositiveRelativeFeatures = sum(asnicar_ra[, i] > 0),
      PositivePseudoCountFeatures = sum(pseudo > 0),
      ExactRawFeatureCountsRecovered = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)

metadata_by_study <- lapply(selected_studies, function(study_name) {
  metadata_all[metadata_all$study_name == study_name, , drop = FALSE]
})
names(metadata_by_study) <- selected_studies
accession_summary <- lapply(metadata_by_study, function(x) {
  list(
    sample_records = nrow(x),
    populated = sum(present(x$NCBI_accession)),
    run_tokens = length(accession_tokens(x$NCBI_accession))
  )
})

lineage_cards <- data.frame(
  Study = c("HMP_2012", "ZellerG_2014", "AsnicarF_2017", "AsnicarF_2017", "MeslierE_2022_MOCK1"),
  Accession = c(
    sprintf("%d SRA run tokens across %d sample records", accession_summary$HMP_2012$run_tokens, accession_summary$HMP_2012$sample_records),
    sprintf("%d ENA/SRA run tokens across %d sample records", accession_summary$ZellerG_2014$run_tokens, accession_summary$ZellerG_2014$sample_records),
    sprintf("SRR4052021-SRR4052044; %d run tokens", accession_summary$AsnicarF_2017$run_tokens),
    sprintf("SRR4052021-SRR4052044; %d run tokens", accession_summary$AsnicarF_2017$run_tokens),
    "PRJEB52977 / SAMEA14435832 / ERR9765746"
  ),
  RawReadsAvailable = c(
    sprintf("Yes; NCBI_accession populated for %d/%d sample records", accession_summary$HMP_2012$populated, accession_summary$HMP_2012$sample_records),
    sprintf("Yes; NCBI_accession populated for %d/%d sample records", accession_summary$ZellerG_2014$populated, accession_summary$ZellerG_2014$sample_records),
    sprintf("Yes; NCBI_accession populated for %d/%d sample records", accession_summary$AsnicarF_2017$populated, accession_summary$AsnicarF_2017$sample_records),
    sprintf("Yes; NCBI_accession populated for %d/%d sample records", accession_summary$AsnicarF_2017$populated, accession_summary$AsnicarF_2017$sample_records),
    "Yes; checksum-locked paired FASTQ subset"
  ),
  ProfileSource = c(
    "2021-03-31.HMP_2012.relative_abundance",
    "2021-03-31.ZellerG_2014.relative_abundance",
    "2021-10-14.AsnicarF_2017.relative_abundance (EH7091)",
    "2021-10-14.AsnicarF_2017.pathway_abundance (EH7089)",
    "Article 15 frozen profile-all.tsv / species-profile.tsv"
  ),
  Profiler = c("MetaPhlAn", "MetaPhlAn", "MetaPhlAn", "HUMAnN", "MetaPhlAn"),
  ProfilerVersion = c(
    "3.x; exact patch not encoded in resource title",
    "3.x; exact patch not encoded in resource title",
    "3.x; exact patch not encoded in resource title",
    "3.x; exact patch not encoded in resource title",
    "4.2.5"
  ),
  DatabaseRelease = c(
    "CHOCOPhlAn 201901", "CHOCOPhlAn 201901", "CHOCOPhlAn 201901",
    "UniRef90 and MetaCyc; exact releases not encoded in resource title",
    "mpa_vJan26_CHOCOPhlAnSGB_202605"
  ),
  FeatureType = c(
    "Taxonomic lineage; species rows used",
    "Taxonomic lineage; species rows used",
    "Taxonomic lineage; species rows used",
    "MetaCyc pathway; community and stratified rows coexist",
    "Taxonomic lineage / SGB; species rows shown"
  ),
  Unit = c(
    "Relative abundance percent", "Relative abundance percent", "Relative abundance percent",
    "Relative abundance proportion with multiple row scopes", "Relative abundance percent"
  ),
  Normalization = c(
    rep("MetaPhlAn compositional profile; align one rank before analysis", 3L),
    "HUMAnN relative abundance; split community and stratified rows before closure",
    "MetaPhlAn rel_ab_w_read_stats; UNCLASSIFIED retained in the primary profile"
  ),
  OriginalPublication = c(
    "PMID 22699609; DOI 10.1038/nature11234",
    "PMID 25432777; DOI 10.15252/msb.20145645",
    "PMID 28144631; DOI 10.1128/mSystems.00164-16",
    "PMID 28144631; DOI 10.1128/mSystems.00164-16",
    "DOI 10.1038/s41597-022-01762-z"
  ),
  ReprocessedOrAuthorProvided = c(
    rep("Uniform cMD3 / bioBakery3 reprocessing of public raw reads", 4L),
    "Local checksum-locked MetaPhlAn4 reprocessing; not a cMD resource"
  ),
  ProfileFileSHA256 = unname(observed_hashes[c(
    "hmp-rda", "zeller-rda", "asnicar-rds", "asnicar-pathway-rda", "metaphlan4-profile"
  )]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

compatibility <- data.frame(
  Comparison = c(
    "HMP tax + Zeller tax",
    "HMP tax + Asnicar latest tax",
    "Asnicar old + latest tax",
    "Asnicar tax + pathway",
    "cMD3 tax + MetaPhlAn4 tax",
    "cMD percent + counts=TRUE",
    "Repeated HMP stool + subject model"
  ),
  SameFeatureType = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  SameProfilerLineage = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  SameDatabaseRelease = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  SameUnitAndDenominator = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
  DisjointBiologicalUnits = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
  Verdict = c(
    "Conditional merge", "Conditional merge", "Choose one release",
    "Separate assays", "Reprocess or harmonize", "Derived, not independent",
    "Deduplicate or model subject"
  ),
  RequiredAction = c(
    "Align species IDs, preserve study labels, and model study",
    "Align species IDs, preserve study labels, and model study",
    "Never concatenate duplicate snapshots of the same samples",
    "Keep taxonomy and function in separate assay spaces",
    "Do not equate MetaPhlAn3 species with MetaPhlAn4 SGB features",
    "Record the reconstruction; do not call values observed taxon reads",
    "Use one sample per study x subject or a repeated-measures model"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

analysis_contract <- data.frame(
  Item = c(
    "package_version", "bioconductor_release", "package_git_commit",
    "package_snapshot_samples", "package_snapshot_studies", "metadata_columns",
    "globally_colliding_sample_ids",
    "resource_titles", "resource_types", "study_type_snapshots",
    "study_type_latest_resources", "study_type_pairs_with_two_releases",
    "query_studies", "query_body_site", "subject_selection_rule",
    "query_stool_profiles", "query_independent_units", "merge_function",
    "merge_assay", "merge_unit", "seed", "paper_collection_samples",
    "paper_collection_cohorts", "paper_only_addition"
  ),
  Value = c(
    "3.12.0", "3.19", "c5711e9", "22588", "93", "141", "177", "732", "6", "732",
    "558", "174", paste(selected_studies, collapse = " | "), "stool",
    "highest finite number_reads within study x subject; lexical sample_id tie-break",
    "319", "261", "curatedMetagenomicData::mergeData", "relative_abundance",
    "percent", as.character(primary_seed), "22710", "94",
    "HeQ_2017: 122 samples used in the paper but absent from the package snapshot"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

resource_manifest <- data.frame(
  Resource = c(
    "curatedMetagenomicData package snapshot", "HMP_2012 relative abundance",
    "ZellerG_2014 relative abundance", "AsnicarF_2017 relative abundance",
    "AsnicarF_2017 pathway abundance", "MetaPhlAn4 comparison profile",
    "Manghi et al. Figure 1"
  ),
  Version = c(
    "3.12.0 / Bioconductor 3.19 / git c5711e9", "2021-03-31", "2021-03-31",
    "2021-10-14 / EH7091", "2021-10-14 / EH7089", "MetaPhlAn 4.2.5 / vJan26_202605",
    "Nature Communications 17:196 (2026)"
  ),
  Source = c(
    "Bioconductor curatedMetagenomicData",
    "ExperimentHub curatedMetagenomicData/HMP_2012",
    "ExperimentHub curatedMetagenomicData/ZellerG_2014",
    "ExperimentHub curatedMetagenomicData/AsnicarF_2017",
    "ExperimentHub curatedMetagenomicData/AsnicarF_2017",
    "Article 15 checksum-locked local run",
    "https://doi.org/10.1038/s41467-025-66888-1"
  ),
  SHA256 = c(
    "package-data-object-recorded-by-derived-file-checksums",
    observed_hashes[["hmp-rda"]], observed_hashes[["zeller-rda"]],
    observed_hashes[["asnicar-rds"]], observed_hashes[["asnicar-pathway-rda"]],
    observed_hashes[["metaphlan4-profile"]],
    "559c95bd7e5c35c99e853bc16bc6c8d9739a9d7555992bf58efef3fd5b77b7c1"
  ),
  License = c(
    rep("Artistic-2.0 package plus original-study terms", 5L),
    "Derived tutorial evidence; source study terms apply", "CC BY 4.0"
  ),
  InterpretationBoundary = c(
    "Package snapshot has 22588 samples from 93 studies; it is not the 22710-sample paper analysis set",
    "MetaPhlAn3 lineage and 201901 database; not interchangeable with current MetaPhlAn4 SGBs",
    "MetaPhlAn3 lineage and 201901 database; not interchangeable with current MetaPhlAn4 SGBs",
    "Latest resource is selected automatically from two date-stamped snapshots",
    "Community and stratified pathway rows must be separated",
    "Comparison lineage only; do not merge directly with cMD3 features",
    "Original resource overview; tutorial figures are independent audits"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write_tsv(resource_catalog, file.path(output_dir, "resource-catalog.tsv"))
write_tsv(release_summary, file.path(output_dir, "resource-release-summary.tsv"))
write_tsv(query_example, file.path(output_dir, "asnicar-query-resolution.tsv"))
write_tsv_gz(sample_catalog, file.path(output_dir, "sample-catalog.tsv.gz"))
write_tsv(sample_id_collisions, file.path(output_dir, "sample-id-collisions.tsv"))
write_tsv(metadata_completeness, file.path(output_dir, "metadata-completeness.tsv"))
write_tsv(attrition, file.path(output_dir, "query-attrition.tsv"))
write_tsv_gz(
  selection_pool[, selection_export_fields, drop = FALSE],
  file.path(output_dir, "sample-selection-audit.tsv.gz")
)
write_tsv_gz(
  selected_metadata[, selected_export_fields, drop = FALSE],
  file.path(output_dir, "selected-sample-metadata.tsv.gz")
)
write_tsv_gz(merged_profile, file.path(output_dir, "merged-species-relative-abundance.tsv.gz"))
write_tsv(profile_contract, file.path(output_dir, "merged-object-contract.tsv"))
write_tsv(counts_boundary, file.path(output_dir, "counts-boundary.tsv"))
write_tsv(lineage_cards, file.path(output_dir, "lineage-cards.tsv"))
write_tsv(compatibility, file.path(output_dir, "merge-compatibility.tsv"))
write_tsv(analysis_contract, file.path(output_dir, "analysis-contract.tsv"))
write_tsv(resource_manifest, file.path(output_dir, "resource-manifest.tsv"))

notice_lines <- c(
  "Article 26 data notice — curatedMetagenomicData deep use and lineage cards",
  "",
  "Source and license",
  "- The package snapshot is curatedMetagenomicData 3.12.0 from Bioconductor 3.19 (Artistic-2.0).",
  "- Human profiles are uniformly reprocessed cMD3 resources; original-study terms remain applicable.",
  "- The original Figure 1 is from Manghi et al. 2026 under CC BY 4.0.",
  "",
  "Snapshot boundary",
  "- Package sampleMetadata contains 22,588 samples, 93 studies and 141 columns.",
  "- The paper analysis contains 22,710 samples from 94 cohorts because HeQ_2017 (122 samples) was added for the paper but is absent from this package snapshot.",
  "- Resource discovery contains 732 date-stamped titles: 122 per each of six data types.",
  "- sample_id is not globally unique: 177 IDs occur in two studies and must be qualified by study_name when merged.",
  "- An unqualified data query automatically keeps the latest date per study x data type; dryrun output can still show both releases.",
  "",
  "Profile lineage boundary",
  "- cMD3 taxonomy was generated with MetaPhlAn3 and CHOCOPhlAn 201901; cMD3 function was generated with HUMAnN3.",
  "- Exact MetaPhlAn3/HUMAnN3 patch versions and every functional database release are not encoded in ExperimentHub resource titles; the lineage card records this as unavailable rather than inventing a value.",
  "- cMD3 features must not be directly equated with the MetaPhlAn 4.2.5 / mpa_vJan26_CHOCOPhlAnSGB_202605 feature universe used in Article 15.",
  "",
  "Sampling and unit boundary",
  "- The worked query uses AsnicarF_2017, HMP_2012 and ZellerG_2014 stool samples with traceable sample and subject IDs.",
  "- One sample per study x subject is selected deterministically by maximum finite number_reads, then lexical sample_id.",
  "- relative_abundance is a percentage composition, not absolute microbial load.",
  "- counts=TRUE multiplies percentage by total reads and rounds; it does not recover observed feature-level read counts.",
  "- Taxonomy and function remain separate assay spaces even when they describe the same biological samples.",
  "",
  "Frozen outputs",
  "- resource-catalog.tsv and sample-catalog.tsv.gz preserve resource and sample discovery surfaces.",
  "- lineage-cards.tsv records study, accession, raw-read availability, profiler, version, database, feature, unit, normalization, publication and processing origin.",
  "- merged-species-relative-abundance.tsv.gz is the audited output of curatedMetagenomicData::mergeData on 261 independent stool units.",
  "- file-checksums.sha256 makes all files in data/small/26-cmd-lineage immutable for rendering."
)
writeLines(notice_lines, notice_path, useBytes = TRUE)

payloads <- sort(setdiff(list.files(output_dir, full.names = TRUE), file.path(output_dir, "file-checksums.sha256")))
checksum_lines <- vapply(
  payloads,
  function(path) paste(sha256_file(path), basename(path), sep = "  "),
  character(1L)
)
writeLines(checksum_lines, file.path(output_dir, "file-checksums.sha256"), useBytes = TRUE)

cat("Article 26 frozen data prepared.\n")
cat("Resources:", nrow(resource_catalog), "titles;", sum(resource_catalog$LatestForStudyType), "latest study-type resources.\n")
cat("Metadata:", nrow(metadata_all), "samples x", ncol(metadata_all), "fields from", length(unique(metadata_all$study_name)), "studies.\n")
cat("Query:", nrow(selection_pool), "stool profiles ->", nrow(selected_metadata), "study-subject units.\n")
cat("Merged object:", nrow(merged_tse), "species features x", ncol(merged_tse), "samples.\n")
cat("Checksummed payloads:", length(payloads), "\n")
