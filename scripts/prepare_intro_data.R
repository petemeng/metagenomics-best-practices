args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: prepare_intro_data.R <meta_all.tsv> ",
    "<cohort_summary.tsv> <summary.json>"
  )
}

input_path <- args[[1]]
cohort_path <- args[[2]]
summary_path <- args[[3]]

metadata <- utils::read.delim(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

required_columns <- c("Sample_ID", "Study", "Group", "Library_Size")
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(metadata$Sample_ID)) {
  stop("Sample_ID values must be unique")
}

discovery_studies <- c("FR-CRC", "AT-CRC", "CN-CRC", "US-CRC", "DE-CRC")
validation_studies <- c("IT-CRC", "IT-CRC-2", "JP-CRC")
study_order <- c(discovery_studies, validation_studies)

crc <- metadata[
  metadata$Study %in% study_order & metadata$Group %in% c("CRC", "CTR"),
  ,
  drop = FALSE
]
crc$Role <- ifelse(
  crc$Study %in% discovery_studies,
  "Discovery / meta-analysis",
  "Independent validation"
)

counts <- as.data.frame(
  stats::xtabs(~ Study + Group, data = crc),
  stringsAsFactors = FALSE
)
counts <- counts[counts$Freq > 0L, , drop = FALSE]
wide <- reshape(
  counts,
  idvar = "Study",
  timevar = "Group",
  direction = "wide"
)
names(wide) <- sub("^Freq\\.", "", names(wide))
for (column in c("CRC", "CTR")) {
  if (!column %in% names(wide)) {
    wide[[column]] <- 0L
  }
}
wide$Role <- ifelse(
  wide$Study %in% discovery_studies,
  "Discovery / meta-analysis",
  "Independent validation"
)
wide$Total <- wide$CRC + wide$CTR
wide <- wide[match(study_order, wide$Study), c("Study", "Role", "CRC", "CTR", "Total")]
names(wide)[names(wide) == "CTR"] <- "Control"

if (anyNA(wide$Study)) {
  stop("One or more declared CRC studies are absent from meta_all.tsv")
}

expected_study_totals <- c(
  "FR-CRC" = 114L,
  "AT-CRC" = 109L,
  "CN-CRC" = 128L,
  "US-CRC" = 104L,
  "DE-CRC" = 120L,
  "IT-CRC" = 53L,
  "IT-CRC-2" = 60L,
  "JP-CRC" = 80L
)
stopifnot(
  identical(as.integer(wide$Total), unname(expected_study_totals[wide$Study]))
)

dir.create(dirname(cohort_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  wide,
  cohort_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

source_samples <- nrow(metadata)
source_studies <- length(unique(metadata$Study))
crc_cohorts <- length(unique(crc$Study))
crc_case_control_samples <- nrow(crc)
discovery_samples <- sum(crc$Role == "Discovery / meta-analysis")
validation_samples <- sum(crc$Role == "Independent validation")
crc_cases <- sum(crc$Group == "CRC")
controls <- sum(crc$Group == "CTR")

json <- sprintf(
  paste0(
    "{\n",
    "  \"source_samples\": %d,\n",
    "  \"source_studies\": %d,\n",
    "  \"crc_cohorts\": %d,\n",
    "  \"crc_case_control_samples\": %d,\n",
    "  \"discovery_samples\": %d,\n",
    "  \"validation_samples\": %d,\n",
    "  \"crc_cases\": %d,\n",
    "  \"controls\": %d\n",
    "}\n"
  ),
  source_samples,
  source_studies,
  crc_cohorts,
  crc_case_control_samples,
  discovery_samples,
  validation_samples,
  crc_cases,
  controls
)
writeLines(json, summary_path, useBytes = TRUE)

stopifnot(
  source_samples == 1892L,
  source_studies == 14L,
  crc_cohorts == 8L,
  crc_case_control_samples == 768L,
  discovery_samples == 575L,
  validation_samples == 193L
)

cat(json)
