#!/usr/bin/env Rscript

# GSE255459 — FPT cross-inducer calibration
#
# Verified GEO design:
#   BT549, HS578, SUM159
#   DMSO, Erastin, RSL3
#   biological duplicates per condition
#
# Analysis principle:
#   - analyze each cell line separately
#   - estimate Erastin-vs-DMSO and RSL3-vs-DMSO
#   - shrink effects with ashr and retain lfsr
#   - compute direction-aware e for the locked FPT dictionary
#   - aggregate inducer evidence within cell line
#   - aggregate cell-line evidence within study
#
# This script intentionally does NOT merge this dataset with any other study.

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ashr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

source("R/evidence_functions.R")

dir.create("data/raw/GSE255459", recursive = TRUE, showWarnings = FALSE)
dir.create("results/GSE255459", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Download supplementary count matrices
# -------------------------------------------------------------------------

GEOquery::getGEOSuppFiles(
  "GSE255459",
  makeDirectory = FALSE,
  baseDir = "data/raw/GSE255459"
)

files <- list.files(
  "data/raw/GSE255459",
  pattern = "Count_matrix.*\\.csv\\.gz$",
  full.names = TRUE
)

if (length(files) != 2) {
  stop("Expected two GSE255459 count-matrix files; found: ", length(files))
}

# -------------------------------------------------------------------------
# 2. Robustly extract verified sample columns
# -------------------------------------------------------------------------

verified_samples <- tibble::tribble(
  ~sample,            ~cell_line, ~condition,
  "BT549_DMSO_a",     "BT549",    "DMSO",
  "BT549_DMSO_b",     "BT549",    "DMSO",
  "BT549_Erastin_a",  "BT549",    "Erastin",
  "BT549_Erastin_b",  "BT549",    "Erastin",
  "BT549_RSL3_a",     "BT549",    "RSL3",
  "BT549_RSL3_b",     "BT549",    "RSL3",
  "HS578_DMSO_a",     "HS578",    "DMSO",
  "HS578_DMSO_b",     "HS578",    "DMSO",
  "HS578_Erastin_a",  "HS578",    "Erastin",
  "HS578_Erastin_b",  "HS578",    "Erastin",
  "HS578_RSL3_a",     "HS578",    "RSL3",
  "HS578_RSL3_b",     "HS578",    "RSL3",
  "SUM159_DMSO_a",    "SUM159",   "DMSO",
  "SUM159_DMSO_b",    "SUM159",   "DMSO",
  "SUM159_Erastin_a", "SUM159",   "Erastin",
  "SUM159_Erastin_b", "SUM159",   "Erastin",
  "SUM159_RSL3_a",    "SUM159",   "RSL3",
  "SUM159_RSL3_b",    "SUM159",   "RSL3"
)

clean_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

read_geo_count_file <- function(path, verified_sample_names) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  original <- names(x)
  cleaned <- clean_name(original)

  target_clean <- clean_name(verified_sample_names)
  sample_idx <- match(target_clean, cleaned)
  sample_idx <- sample_idx[!is.na(sample_idx)]

  if (length(sample_idx) == 0) {
    stop("No verified GEO sample columns found in ", basename(path))
  }

  # Choose the first non-sample column as the gene identifier.
  non_sample_idx <- setdiff(seq_along(original), sample_idx)
  if (length(non_sample_idx) == 0) {
    stop("No gene identifier column detected in ", basename(path))
  }

  gene_col <- non_sample_idx[1]

  out <- x[, c(gene_col, sample_idx), drop = FALSE]
  names(out)[1] <- "gene"

  names(out)[-1] <- verified_sample_names[
    match(clean_name(names(out)[-1]), target_clean)
  ]

  out
}

parts <- lapply(files, read_geo_count_file,
                verified_sample_names = verified_samples$sample)

counts_tbl <- Reduce(
  function(a, b) full_join(a, b, by = "gene"),
  parts
)

present_samples <- intersect(verified_samples$sample, names(counts_tbl))
missing_samples <- setdiff(verified_samples$sample, present_samples)

if (length(missing_samples) > 0) {
  stop("Missing verified samples: ", paste(missing_samples, collapse = ", "))
}

counts_tbl <- counts_tbl %>%
  select(gene, all_of(verified_samples$sample))

# Handle duplicate identifiers explicitly and reproducibly.
counts_tbl <- counts_tbl %>%
  group_by(gene) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

count_mat <- as.matrix(counts_tbl[, -1])
rownames(count_mat) <- counts_tbl$gene
storage.mode(count_mat) <- "numeric"

if (any(count_mat < 0, na.rm = TRUE)) {
  stop("Negative values detected; matrix does not behave like raw counts.")
}

integer_fraction <- mean(abs(count_mat - round(count_mat)) < 1e-8, na.rm = TRUE)
if (integer_fraction < 0.999) {
  stop("Matrix is not integer-like enough for DESeq2. Integer fraction = ",
       round(integer_fraction, 4))
}

count_mat <- round(count_mat)

# Save locked sample metadata.
write_csv(verified_samples, "results/GSE255459/sample_metadata.csv")

# -------------------------------------------------------------------------
# 3. Basic QC
# -------------------------------------------------------------------------

qc <- tibble(
  sample = colnames(count_mat),
  library_size = colSums(count_mat),
  detected_genes = colSums(count_mat > 0)
) %>%
  left_join(verified_samples, by = "sample")

write_csv(qc, "results/GSE255459/library_qc.csv")

# -------------------------------------------------------------------------
# 4. Differential analysis per cell line
# -------------------------------------------------------------------------

fpt_dict <- read_csv("config/fpt_gene_dictionary.csv", show_col_types = FALSE) %>%
  filter(direction %in% c(-1, 1))

all_effects <- list()
all_dict_e <- list()

for (cl in unique(verified_samples$cell_line)) {

  meta <- verified_samples %>%
    filter(cell_line == cl) %>%
    mutate(
      condition = factor(condition, levels = c("DMSO", "Erastin", "RSL3"))
    )

  sub_counts <- count_mat[, meta$sample, drop = FALSE]

  filtered <- prefilter_counts(
    sub_counts,
    group = meta$condition,
    min_count = 10L
  )

  coldata <- as.data.frame(meta)
  rownames(coldata) <- coldata$sample

  dds <- fit_deseq2(filtered, coldata, ~ condition)

  contrasts <- list(
    Erastin_vs_DMSO = c("condition", "Erastin", "DMSO"),
    RSL3_vs_DMSO    = c("condition", "RSL3", "DMSO")
  )

  for (contrast_id in names(contrasts)) {
    eff <- ashr_effect_from_contrast(dds, contrasts[[contrast_id]]) %>%
      mutate(
        accession = "GSE255459",
        cell_line = cl,
        contrast_id = contrast_id
      )

    all_effects[[paste(cl, contrast_id, sep = "__")]] <- eff

    scored <- score_dictionary_genes(eff, fpt_dict) %>%
      mutate(
        accession = "GSE255459",
        cell_line = cl,
        contrast_id = contrast_id
      )

    all_dict_e[[paste(cl, contrast_id, sep = "__")]] <- scored
  }
}

effect_tbl <- bind_rows(all_effects)
dict_e_tbl <- bind_rows(all_dict_e)

write_csv(
  effect_tbl,
  "results/GSE255459/per_contrast_effects_all_genes.csv"
)

write_csv(
  dict_e_tbl,
  "results/GSE255459/per_contrast_FPT_dictionary_evidence.csv"
)

# -------------------------------------------------------------------------
# 5. Hierarchical within-study aggregation
# -------------------------------------------------------------------------

cell_line_e <- dict_e_tbl %>%
  group_by(gene, cell_line) %>%
  summarise(
    e_cell_line = median(e, na.rm = TRUE),
    .groups = "drop"
  )

study_e <- cell_line_e %>%
  group_by(gene) %>%
  summarise(
    e_study = median(e_cell_line, na.rm = TRUE),
    supporting_cell_lines = sum(is.finite(e_cell_line)),
    .groups = "drop"
  ) %>%
  left_join(
    fpt_dict %>%
      select(gene, stage, role, direction,
             evidence_tier, mechanistic_prior_M),
    by = "gene"
  )

write_csv(
  cell_line_e,
  "results/GSE255459/per_cell_line_FPT_evidence.csv"
)

write_csv(
  study_e,
  "results/GSE255459/FPT_study_evidence.csv"
)

# -------------------------------------------------------------------------
# 6. Session info
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = "results/GSE255459/sessionInfo.txt"
)

message("GSE255459 analysis completed.")
