#!/usr/bin/env Rscript

# 00_inspect_reference_datasets.R
# Purpose:
#   Download and inspect GEO supplementary files BEFORE any differential analysis.
#   This prevents assumptions about file structure, identifiers, or count type.
#
# Datasets:
#   NEC bulk references:
#   GSE108621, GSE172027, GSE154230, GSE134234, GSE268650
#
#   FPT bulk references:
#   GSE182638, GSE247883, GSE255459, GSE319384, GSE287284
#
#   Generic-stress comparator:
#   GSE104664
#
#   scRNA/human-neuron validation datasets (GSE287439, GSE152988)
#   are intentionally handled in separate workflows.
#
# Output:
#   results/reference_inspection/inspection_summary.csv
#   results/reference_inspection/inspection_log.txt
#
# Please send both files back after running this script.

# -------------------------------------------------------------------------
# 0. Install/load minimal dependencies
# -------------------------------------------------------------------------

cran_pkgs <- c("data.table", "readr", "dplyr", "tibble")
bioc_pkgs <- c("GEOquery")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(readr)
  library(dplyr)
  library(tibble)
})

# -------------------------------------------------------------------------
# 1. Configuration
# -------------------------------------------------------------------------

accessions <- c(
  "GSE108621",
  "GSE172027",
  "GSE154230",
  "GSE134234",
  "GSE268650",
  "GSE182638",
  "GSE247883",
  "GSE255459",
  "GSE319384",
  "GSE287284",
  "GSE104664"
)

base_raw <- "data/inspection"
result_dir <- "results/reference_inspection"

dir.create(base_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(result_dir, "inspection_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# -------------------------------------------------------------------------
# 2. File helpers
# -------------------------------------------------------------------------

is_table_candidate <- function(path) {
  x <- tolower(basename(path))
  grepl("\\.(csv|txt|tsv)(\\.gz)?$", x)
}

safe_fread <- function(path) {
  tryCatch(
    {
      data.table::fread(path, data.table = FALSE, check.names = FALSE)
    },
    error = function(e) {
      message("fread failed for ", basename(path), ": ", conditionMessage(e))
      NULL
    }
  )
}

integer_fraction <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(abs(x - round(x)) < 1e-8)
}

numeric_fraction <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  mean(is.finite(y))
}

summarize_table <- function(acc, path, tbl) {

  if (is.null(tbl) || nrow(tbl) == 0 || ncol(tbl) == 0) {
    return(tibble(
      accession = acc,
      file = basename(path),
      rows = NA_integer_,
      cols = NA_integer_,
      numeric_cols = NA_integer_,
      integer_like_numeric_cols = NA_integer_,
      median_integer_fraction = NA_real_,
      all_values_nonnegative = NA,
      likely_matrix_type = "UNREADABLE_OR_EMPTY"
    ))
  }

  num_frac <- vapply(tbl, numeric_fraction, numeric(1))
  numeric_cols <- which(num_frac > 0.95)

  int_frac <- if (length(numeric_cols) > 0) {
    vapply(tbl[numeric_cols], integer_fraction, numeric(1))
  } else {
    numeric(0)
  }

  nonnegative <- NA
  if (length(numeric_cols) > 0) {
    vals <- unlist(lapply(tbl[numeric_cols], function(z) {
      suppressWarnings(as.numeric(z))
    }), use.names = FALSE)
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0) nonnegative <- all(vals >= 0)
  }

  med_int <- if (length(int_frac) > 0) median(int_frac, na.rm = TRUE) else NA_real_

  matrix_type <- dplyr::case_when(
    length(numeric_cols) < 2 ~ "NOT_EXPRESSION_MATRIX_OR_NEEDS_MANUAL_REVIEW",
    is.finite(med_int) && med_int >= 0.999 ~ "INTEGER_LIKE_COUNTS_CANDIDATE",
    is.finite(med_int) && med_int < 0.999 ~ "NONINTEGER_NORMALIZED_CANDIDATE",
    TRUE ~ "UNKNOWN"
  )

  tibble(
    accession = acc,
    file = basename(path),
    rows = nrow(tbl),
    cols = ncol(tbl),
    numeric_cols = length(numeric_cols),
    integer_like_numeric_cols = sum(int_frac >= 0.999, na.rm = TRUE),
    median_integer_fraction = med_int,
    all_values_nonnegative = nonnegative,
    likely_matrix_type = matrix_type
  )
}

# -------------------------------------------------------------------------
# 3. Download and inspect
# -------------------------------------------------------------------------

all_summary <- list()

for (acc in accessions) {

  log_msg("")
  log_msg("============================================================")
  log_msg("ACCESSION: ", acc)
  log_msg("============================================================")

  acc_dir <- file.path(base_raw, acc)
  dir.create(acc_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg("[1] Downloading GEO supplementary files...")

  tryCatch(
    {
      GEOquery::getGEOSuppFiles(
        GEO = acc,
        makeDirectory = FALSE,
        baseDir = acc_dir,
        fetch_files = TRUE
      )
    },
    error = function(e) {
      log_msg("DOWNLOAD ERROR: ", conditionMessage(e))
    }
  )

  files <- list.files(
    acc_dir,
    recursive = TRUE,
    full.names = TRUE
  )

  log_msg("[2] Files found: ", length(files))

  if (length(files) == 0) {
    all_summary[[acc]] <- tibble(
      accession = acc,
      file = NA_character_,
      rows = NA_integer_,
      cols = NA_integer_,
      numeric_cols = NA_integer_,
      integer_like_numeric_cols = NA_integer_,
      median_integer_fraction = NA_real_,
      all_values_nonnegative = NA,
      likely_matrix_type = "NO_FILES"
    )
    next
  }

  for (f in files) {
    info <- file.info(f)
    log_msg("  - ", basename(f), " | ",
            round(info$size / 1024^2, 3), " MB")
  }

  table_files <- files[vapply(files, is_table_candidate, logical(1))]

  log_msg("[3] Candidate tabular files: ", length(table_files))

  if (length(table_files) == 0) {
    all_summary[[acc]] <- tibble(
      accession = acc,
      file = NA_character_,
      rows = NA_integer_,
      cols = NA_integer_,
      numeric_cols = NA_integer_,
      integer_like_numeric_cols = NA_integer_,
      median_integer_fraction = NA_real_,
      all_values_nonnegative = NA,
      likely_matrix_type = "NO_TABULAR_SUPPLEMENT"
    )
    next
  }

  for (f in table_files) {

    log_msg("")
    log_msg("FILE: ", basename(f))

    tbl <- safe_fread(f)

    s <- summarize_table(acc, f, tbl)
    all_summary[[paste(acc, basename(f), sep = "__")]] <- s

    if (is.null(tbl)) {
      log_msg("  Could not read table.")
      next
    }

    log_msg("  Dimensions: ", nrow(tbl), " rows x ", ncol(tbl), " columns")
    log_msg("  Column names:")
    log_msg("    ", paste(names(tbl), collapse = " | "))

    log_msg("  First 3 rows:")
    preview <- utils::capture.output(print(utils::head(tbl, 3)))
    for (line in preview) log_msg("    ", line)

    log_msg("  Numeric columns: ", s$numeric_cols)
    log_msg("  Integer-like numeric columns: ", s$integer_like_numeric_cols)
    log_msg("  Median integer fraction: ",
            ifelse(is.na(s$median_integer_fraction), "NA",
                   round(s$median_integer_fraction, 5)))
    log_msg("  Nonnegative numeric values: ", s$all_values_nonnegative)
    log_msg("  Preliminary type: ", s$likely_matrix_type)

    # Approximate library sizes only for integer-like candidate matrices.
    if (identical(s$likely_matrix_type, "INTEGER_LIKE_COUNTS_CANDIDATE")) {

      num_frac <- vapply(tbl, numeric_fraction, numeric(1))
      num_idx <- which(num_frac > 0.95)

      if (length(num_idx) >= 2) {
        num_mat <- suppressWarnings(
          as.matrix(data.frame(lapply(tbl[num_idx], as.numeric),
                               check.names = FALSE))
        )

        lib <- colSums(num_mat, na.rm = TRUE)

        log_msg("  Approximate numeric-column sums:")
        for (i in seq_along(lib)) {
          log_msg("    ", names(tbl)[num_idx[i]], ": ", format(lib[i], scientific = FALSE))
        }
      }
    }
  }
}

# -------------------------------------------------------------------------
# 4. Save summary
# -------------------------------------------------------------------------

summary_tbl <- bind_rows(all_summary)

readr::write_csv(
  summary_tbl,
  file.path(result_dir, "inspection_summary.csv")
)

log_msg("")
log_msg("============================================================")
log_msg("DONE")
log_msg("============================================================")
log_msg("Please send these two files:")
log_msg("1. ", file.path(result_dir, "inspection_summary.csv"))
log_msg("2. ", file.path(result_dir, "inspection_log.txt"))

cat("\nInspection complete.\n")
cat("Outputs:\n")
cat("  ", file.path(result_dir, "inspection_summary.csv"), "\n")
cat("  ", file.path(result_dir, "inspection_log.txt"), "\n")
