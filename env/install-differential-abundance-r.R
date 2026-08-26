#!/usr/bin/env Rscript

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 900,
  stringsAsFactors = FALSE
)

required_base <- c(BiocManager = "1.30.27", remotes = "2.5.0", digest = "0.6.39")
for (package in names(required_base)) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

install_cran_version <- function(package, version) {
  installed <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
  if (is.na(installed) ||
      utils::compareVersion(installed, version) != 0L) {
    remotes::install_version(
      package, version = version, dependencies = FALSE,
      upgrade = "never", quiet = FALSE
    )
  }
  observed <- as.character(utils::packageVersion(package))
  if (utils::compareVersion(observed, version) != 0L) {
    stop(package, " version mismatch: expected ", version, ", observed ", observed)
  }
}

# Source-installed CRAN packages are pinned because they are not represented
# in the conda explicit lock generated from the validated environment.
cran_versions <- c(
  zigg = "0.0.2",
  quadprog = "1.5-8",
  Rfast = "2.1.5.2",
  directlabels = "2026.4.23",
  collapse = "2.1.7",
  nanonext = "1.10.1",
  mirai = "2.7.2",
  ggnewscale = "0.5.2",
  patchwork = "1.3.2"
)
for (package in names(cran_versions)) {
  install_cran_version(package, cran_versions[[package]])
}

if (!requireNamespace("ALDEx2", quietly = TRUE) ||
    as.character(utils::packageVersion("ALDEx2")) != "1.42.0") {
  BiocManager::install(
    "ALDEx2", version = "3.22", update = FALSE,
    ask = FALSE, Ncpus = 1
  )
}

maaslin_commit <- "3a194ece449ef249354df394b58bfe3e6f951ca3"
maaslin_url <- paste0(
  "https://codeload.github.com/biobakery/maaslin3/tar.gz/",
  maaslin_commit
)
maaslin_sha256 <- "358f35e094e03026c8d694d731c0d1e1e6c9a15dec3c466649b9db9329ca1f07"
archive <- tempfile(fileext = ".tar.gz")
utils::download.file(maaslin_url, archive, mode = "wb", quiet = FALSE)
observed_sha256 <- digest::digest(
  file = archive, algo = "sha256", serialize = FALSE
)
if (!identical(observed_sha256, maaslin_sha256)) {
  stop("MaAsLin3 source SHA-256 mismatch: ", observed_sha256)
}
if (!requireNamespace("maaslin3", quietly = TRUE) ||
    as.character(utils::packageVersion("maaslin3")) != "1.5.3") {
  remotes::install_local(
    archive, dependencies = FALSE, upgrade = "never",
    force = TRUE, quiet = FALSE
  )
}

expected <- c(
  maaslin3 = "1.5.3",
  ANCOMBC = "2.12.0",
  ALDEx2 = "1.42.0",
  collapse = "2.1.7",
  mirai = "2.7.2",
  nanonext = "1.10.1",
  ggnewscale = "0.5.2",
  patchwork = "1.3.2"
)
observed <- vapply(
  names(expected),
  function(package) as.character(utils::packageVersion(package)),
  character(1L)
)
if (!identical(unname(observed), unname(expected))) {
  stop(
    "Final R package mismatch: ",
    paste(names(expected), expected, observed, sep = "=", collapse = "; ")
  )
}
writeLines(
  paste(names(observed), observed, sep = "\t"),
  con = stdout()
)
