#!/usr/bin/env Rscript

# GSE182638 — primary rescue-validated ferroptosis calibration
#
# Verified design:
#   Human MM1R and MM1S multiple-myeloma cell lines
#   Untreated / RSL3 / Fer-1 pre-treatment + RSL3
#   n = 3 biological replicates per condition per cell line
#
# Treatment:
#   RSL3: 5 uM for 3 h
#   Rescue: 2 uM Ferrostatin-1 pre-treatment for 2 h, then 5 uM RSL3 for 3 h
#
# Analysis:
#   - download raw supplementary count matrix from GEO
#   - recover sample metadata from GEO
#   - analyze MM1R and MM1S separately
#   - block by biological replicate/passage: ~ replicate + condition
#   - induction: RSL3 - Untreated
#   - reversal: RSL3 - (RSL3 + Fer-1)
#   - DESeq2 effect + SE -> ashr posterior beta + lfsr
#   - compute rescue-validated evidence for the mechanistic gene universe
#   - derive an ALL-GENE rescue-concordant candidate table for later FPT-ESR
#
# Run from repository root:
#   source("ferroptosis-only/scripts/01_GSE182638_primary_rescue.R")

required_pkgs <- c(
  "DESeq2", "ashr", "dplyr", "tidyr", "readr", "tibble",
  "data.table", "ggplot2", "pheatmap", "AnnotationDbi", "org.Hs.eg.db"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ashr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

source("R/evidence_functions.R")

accession <- "GSE182638"
raw_dir <- file.path("ferroptosis-only", "data", "raw", accession)
result_dir <- file.path("ferroptosis-only", "results", accession)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Locked GEO metadata
# -------------------------------------------------------------------------

message("[1/9] Loading locked GEO metadata...")

meta <- tibble::tribble(
  ~sample,      ~title,                                           ~cell_line, ~condition,  ~replicate, ~passage, ~verified_short_name,
  "GSM5534024", "MM1R cells treated with Fer-1 and RSL3, rep 1", "MM1R",     "RSL3_Fer1", "1",       "P4",    "MM1R-T3-P4-FR",
  "GSM5534025", "MM1R cells treated with RSL3, rep 1",           "MM1R",     "RSL3",      "1",       "P4",    "MM1R-T3-P4-RSL3",
  "GSM5534026", "Untreated MM1R cells, rep 1",                   "MM1R",     "Untreated", "1",       "P4",    "MM1R-T3-P4-UT",
  "GSM5534027", "MM1R cells treated with Fer-1 and RSL3, rep 2", "MM1R",     "RSL3_Fer1", "2",       "P6",    "MM1R-T3-P6-FR",
  "GSM5534028", "MM1R cells treated with RSL3, rep 2",           "MM1R",     "RSL3",      "2",       "P6",    "MM1R-T3-P6-RSL3",
  "GSM5534029", "Untreated MM1R cells, rep 2",                   "MM1R",     "Untreated", "2",       "P6",    "MM1R-T3-P6-UT",
  "GSM5534030", "MM1R cells treated with Fer-1 and RSL3, rep 3", "MM1R",     "RSL3_Fer1", "3",       "P7",    "MM1R-T3-P7-FR",
  "GSM5534031", "MM1R cells treated with RSL3, rep 3",           "MM1R",     "RSL3",      "3",       "P7",    "MM1R-T3-P7-RSL3",
  "GSM5534032", "Untreated MM1R cells, rep 3",                   "MM1R",     "Untreated", "3",       "P7",    "MM1R-T3-P7-UT",
  "GSM5534033", "MM1S cells treated with Fer-1 and RSL3, rep 1", "MM1S",     "RSL3_Fer1", "1",       "P16",   "MM1S-T3-P16-FR",
  "GSM5534034", "MM1S cells treated with RSL3, rep 1",           "MM1S",     "RSL3",      "1",       "P16",   "MM1S-T3-P16-RSL3",
  "GSM5534035", "Untreated MM1S cells, rep 1",                   "MM1S",     "Untreated", "1",       "P16",   "MM1S-T3-P16-UT",
  "GSM5534036", "MM1S cells treated with Fer-1 and RSL3, rep 2", "MM1S",     "RSL3_Fer1", "2",       "P22",   "MM1S-T3-P22-FR",
  "GSM5534037", "MM1S cells treated with RSL3, rep 2",           "MM1S",     "RSL3",      "2",       "P22",   "MM1S-T3-P22-RSL3",
  "GSM5534038", "Untreated MM1S cells, rep 2",                   "MM1S",     "Untreated", "2",       "P22",   "MM1S-T3-P22-UT",
  "GSM5534039", "MM1S cells treated with Fer-1 and RSL3, rep 3", "MM1S",     "RSL3_Fer1", "3",       "P23",   "MM1S-T3-P23-FR",
  "GSM5534040", "MM1S cells treated with RSL3, rep 3",           "MM1S",     "RSL3",      "3",       "P23",   "MM1S-T3-P23-RSL3",
  "GSM5534041", "Untreated MM1S cells, rep 3",                   "MM1S",     "Untreated", "3",       "P23",   "MM1S-T3-P23-UT"
)

expected_design <- meta %>%
  count(cell_line, condition, name = "n")

if (
  nrow(expected_design) != 6 ||
  any(expected_design$n != 3) ||
  !setequal(expected_design$cell_line, c("MM1R", "MM1S")) ||
  !setequal(expected_design$condition, c("Untreated", "RSL3", "RSL3_Fer1"))
) {
  stop("Locked metadata failed internal design validation.")
}

write_csv(meta, file.path(result_dir, "sample_metadata_locked.csv"))

# -------------------------------------------------------------------------
# 2. Download and read raw count matrix
# -------------------------------------------------------------------------

message("[2/9] Downloading raw supplementary counts directly from NCBI...")

count_filename <- paste0(accession, "_raw_readcounts_allsamples.txt.gz")
count_file <- file.path(raw_dir, count_filename)

count_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/",
  accession,
  "/suppl/",
  count_filename
)

options(timeout = max(1200, getOption("timeout")))

if (!file.exists(count_file) || file.info(count_file)$size == 0) {
  message("Downloading: ", count_filename)

  download_ok <- tryCatch(
    {
      utils::download.file(
        url = count_url,
        destfile = count_file,
        mode = "wb",
        method = "libcurl",
        quiet = FALSE
      )
      TRUE
    },
    error = function(e) {
      message("libcurl download failed: ", conditionMessage(e))
      FALSE
    }
  )

  if (!download_ok && .Platform$OS.type == "windows") {
    message("Retrying with Windows download method...")
    utils::download.file(
      url = count_url,
      destfile = count_file,
      mode = "wb",
      method = "wininet",
      quiet = FALSE
    )
  }
} else {
  message("Raw count file already exists; skipping download.")
}

if (!file.exists(count_file) || file.info(count_file)$size == 0) {
  stop("Raw count file was not downloaded successfully.")
}

message(
  "Raw count file ready: ",
  basename(count_file),
  " (",
  round(file.info(count_file)$size / 1024, 1),
  " KB)"
)

count_files <- count_file

raw_tbl <- data.table::fread(
  count_files[1],
  data.table = FALSE,
  check.names = FALSE
)

if (nrow(raw_tbl) == 0 || ncol(raw_tbl) < 19) {
  stop("Downloaded count table has unexpected dimensions.")
}

# -------------------------------------------------------------------------
# 3. Match count columns to GEO samples
# -------------------------------------------------------------------------

message("[3/9] Matching count columns to GEO samples...")

clean_key <- function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}

sample_fields <- bind_rows(
  meta %>% transmute(sample, key = clean_key(sample), key_type = "GSM"),
  meta %>% transmute(sample, key = clean_key(title), key_type = "title"),
  meta %>% transmute(
    sample,
    key = clean_key(verified_short_name),
    key_type = "verified_short_name"
  )
) %>%
  filter(nzchar(key)) %>%
  distinct()

count_names <- names(raw_tbl)
count_keys <- clean_key(count_names)

match_sample_column <- function(col_name, col_key) {
  # 1) Exact match against GSM/title/verified short name.
  exact <- sample_fields %>% filter(key == col_key)
  if (nrow(exact) == 1) return(exact$sample)

  # 2) Robust substring match against verified short names.
  #    This handles count columns that append/prepend BAM/file-processing text.
  short_keys <- clean_key(meta$verified_short_name)
  short_hits <- meta$sample[
    vapply(short_keys, function(k) {
      grepl(k, col_key, fixed = TRUE) || grepl(col_key, k, fixed = TRUE)
    }, logical(1))
  ]
  if (length(short_hits) == 1) return(short_hits)

  # 3) Fallback: GEO accession embedded in column name.
  gsm_hits <- meta$sample[
    vapply(meta$sample, function(g) {
      grepl(clean_key(g), col_key, fixed = TRUE)
    }, logical(1))
  ]
  if (length(gsm_hits) == 1) return(gsm_hits)

  NA_character_
}

matched_samples <- vapply(
  seq_along(count_names),
  function(i) match_sample_column(count_names[i], count_keys[i]),
  character(1)
)

sample_idx <- which(!is.na(matched_samples))

mapping_debug <- tibble(
  count_column = count_names,
  cleaned_count_column = count_keys,
  matched_sample = matched_samples
)

write_csv(
  mapping_debug,
  file.path(result_dir, "count_column_mapping_debug.csv")
)

if (length(sample_idx) != 18 || length(unique(matched_samples[sample_idx])) != 18) {
  message("Detected count columns:")
  message(paste(count_names, collapse = " | "))
  stop(
    "Could not uniquely map all 18 count columns to GEO samples. ",
    "See count_column_mapping_debug.csv."
  )
}

non_sample_idx <- setdiff(seq_along(count_names), sample_idx)

# Prefer an explicit symbol/name column, then gene ID, then first non-sample column.
preferred_symbol <- non_sample_idx[
  grepl("symbol|gene.?name", count_names[non_sample_idx], ignore.case = TRUE)
]
preferred_id <- non_sample_idx[
  grepl("gene.?id|ensembl|entrez", count_names[non_sample_idx], ignore.case = TRUE)
]

gene_idx <- if (length(preferred_symbol) > 0) {
  preferred_symbol[1]
} else if (length(preferred_id) > 0) {
  preferred_id[1]
} else {
  non_sample_idx[1]
}

gene_raw <- as.character(raw_tbl[[gene_idx]])

counts_df <- raw_tbl[, sample_idx, drop = FALSE]
names(counts_df) <- matched_samples[sample_idx]

# Reorder to GEO metadata.
counts_df <- counts_df[, meta$sample, drop = FALSE]

# -------------------------------------------------------------------------
# 4. Standardize gene identifiers
# -------------------------------------------------------------------------

message("[4/9] Standardizing gene identifiers...")

standardize_gene_id <- function(ids) {
  ids <- as.character(ids)
  ids_no_ver <- sub("\\.[0-9]+$", "", ids)

  ens_frac <- mean(grepl("^ENSG[0-9]+", ids_no_ver), na.rm = TRUE)
  entrez_frac <- mean(grepl("^[0-9]+$", ids_no_ver), na.rm = TRUE)

  if (ens_frac > 0.50) {
    symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = unique(ids_no_ver),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    mapped <- unname(symbols[ids_no_ver])
    id_type <- "ENSEMBL"
  } else if (entrez_frac > 0.50) {
    symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = unique(ids_no_ver),
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
    mapped <- unname(symbols[ids_no_ver])
    id_type <- "ENTREZID"
  } else {
    mapped <- ids_no_ver
    id_type <- "SYMBOL_OR_OTHER"
  }

  tibble(
    original_id = ids,
    normalized_id = ids_no_ver,
    gene = mapped,
    detected_id_type = id_type
  )
}

id_map <- standardize_gene_id(gene_raw)
write_csv(id_map, file.path(result_dir, "gene_identifier_mapping.csv"))

keep_mapped <- !is.na(id_map$gene) & nzchar(id_map$gene)

if (sum(keep_mapped) < 1000) {
  stop("Too few genes remained after identifier standardization.")
}

tmp <- tibble(gene = id_map$gene[keep_mapped]) %>%
  bind_cols(as_tibble(counts_df[keep_mapped, , drop = FALSE])) %>%
  mutate(across(-gene, as.numeric)) %>%
  group_by(gene) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

count_mat <- as.matrix(tmp[, -1])
rownames(count_mat) <- tmp$gene
storage.mode(count_mat) <- "numeric"

if (any(count_mat < 0, na.rm = TRUE)) {
  stop("Negative values detected in raw count matrix.")
}

integer_fraction <- mean(abs(count_mat - round(count_mat)) < 1e-8, na.rm = TRUE)

if (integer_fraction < 0.999) {
  stop(
    "Matrix is not sufficiently integer-like for DESeq2. Integer fraction = ",
    round(integer_fraction, 5)
  )
}

count_mat <- round(count_mat)

# -------------------------------------------------------------------------
# 5. QC
# -------------------------------------------------------------------------

message("[5/9] Running basic QC...")

qc <- tibble(
  sample = colnames(count_mat),
  library_size = colSums(count_mat),
  detected_genes = colSums(count_mat > 0)
) %>%
  left_join(meta, by = "sample")

write_csv(qc, file.path(result_dir, "library_qc.csv"))

# Global VST PCA for diagnostic purposes only.
qc_keep <- rowSums(count_mat >= 10) >= 3
qc_counts <- count_mat[qc_keep, , drop = FALSE]

qc_coldata <- meta %>%
  mutate(
    cell_line = factor(cell_line),
    condition = factor(condition)
  ) %>%
  as.data.frame()
rownames(qc_coldata) <- qc_coldata$sample

qc_dds <- DESeqDataSetFromMatrix(
  countData = qc_counts,
  colData = qc_coldata,
  design = ~ cell_line + condition
)

vsd <- vst(qc_dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)))

pca_tbl <- as_tibble(pca$x[, 1:2], rownames = "sample") %>%
  left_join(meta, by = "sample")

write_csv(pca_tbl, file.path(result_dir, "vst_pca_coordinates.csv"))

p <- ggplot(pca_tbl, aes(PC1, PC2, shape = condition)) +
  geom_point(size = 3) +
  facet_wrap(~ cell_line) +
  theme_bw() +
  labs(title = "GSE182638 VST PCA")

ggsave(
  file.path(result_dir, "vst_pca.png"),
  p,
  width = 8,
  height = 5,
  dpi = 180
)

corr <- cor(assay(vsd), method = "pearson")
png(
  file.path(result_dir, "vst_sample_correlation.png"),
  width = 1800,
  height = 1600,
  res = 180
)
pheatmap::pheatmap(corr, main = "GSE182638 VST sample correlation")
dev.off()

# -------------------------------------------------------------------------
# 6. Mechanistic universe
# -------------------------------------------------------------------------

message("[6/9] Loading ferroptosis master gene universe...")

fpt_universe <- read_csv(
  file.path(
    "ferroptosis-only",
    "config",
    "ferroptosis_master_gene_universe.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(
    rna_score_candidate == "YES",
    direction_for_rna %in% c(-1, 1)
  ) %>%
  mutate(direction = as.numeric(direction_for_rna))

presence <- fpt_universe %>%
  transmute(
    gene,
    module,
    provisional_status,
    evidence_tier,
    mechanistic_prior_M,
    present_in_GSE182638 = gene %in% rownames(count_mat)
  )

write_csv(
  presence,
  file.path(result_dir, "mechanistic_gene_presence.csv")
)

# -------------------------------------------------------------------------
# 7. Differential effects within each cell line
# -------------------------------------------------------------------------

message("[7/9] Fitting rescue-aware DE models...")

all_effects <- list()
all_mech_e <- list()

for (cl in c("MM1R", "MM1S")) {

  m <- meta %>%
    filter(cell_line == cl) %>%
    mutate(
      replicate = factor(replicate),
      condition = factor(
        condition,
        levels = c("Untreated", "RSL3_Fer1", "RSL3")
      )
    )

  sub_counts <- count_mat[, m$sample, drop = FALSE]

  filtered <- prefilter_counts(
    sub_counts,
    group = m$condition,
    min_count = 10L
  )

  coldata <- as.data.frame(m)
  rownames(coldata) <- coldata$sample

  # replicate acts as passage-matched blocking factor.
  dds <- fit_deseq2(filtered, coldata, ~ replicate + condition)

  contrasts <- list(
    RSL3_vs_Untreated = c("condition", "RSL3", "Untreated"),
    RSL3_vs_RSL3_Fer1 = c("condition", "RSL3", "RSL3_Fer1")
  )

  for (contrast_id in names(contrasts)) {

    eff <- ashr_effect_from_contrast(dds, contrasts[[contrast_id]]) %>%
      mutate(
        accession = accession,
        cell_line = cl,
        contrast_id = contrast_id
      )

    all_effects[[paste(cl, contrast_id, sep = "__")]] <- eff

    scored <- score_dictionary_genes(
      eff,
      fpt_universe %>% select(gene, direction)
    ) %>%
      left_join(
        fpt_universe %>%
          select(
            gene, module, role, provisional_status,
            evidence_tier, mechanistic_prior_M
          ),
        by = "gene"
      ) %>%
      mutate(
        accession = accession,
        cell_line = cl,
        contrast_id = contrast_id
      )

    all_mech_e[[paste(cl, contrast_id, sep = "__")]] <- scored
  }
}

effects <- bind_rows(all_effects)
mech_e <- bind_rows(all_mech_e)

write_csv(
  effects,
  file.path(result_dir, "per_contrast_effects_all_genes.csv")
)

write_csv(
  mech_e,
  file.path(result_dir, "per_contrast_mechanistic_gene_evidence.csv")
)

# -------------------------------------------------------------------------
# 8. Rescue-validated mechanistic evidence + all-gene ESR candidates
# -------------------------------------------------------------------------

message("[8/9] Computing rescue-validated evidence...")

# Mechanistic, direction-aware evidence.
mech_wide <- mech_e %>%
  select(
    gene, cell_line, contrast_id, e,
    module, role, provisional_status,
    evidence_tier, mechanistic_prior_M, direction
  ) %>%
  pivot_wider(
    names_from = contrast_id,
    values_from = e
  ) %>%
  mutate(
    e_rescue = rescue_validated_evidence(
      RSL3_vs_Untreated,
      RSL3_vs_RSL3_Fer1
    )
  )

write_csv(
  mech_wide,
  file.path(result_dir, "per_cell_line_rescue_validated_mechanistic_evidence.csv")
)

study_mech <- mech_wide %>%
  group_by(
    gene, module, role, provisional_status,
    evidence_tier, mechanistic_prior_M, direction
  ) %>%
  summarise(
    e_study = median(e_rescue, na.rm = TRUE),
    MM1R_e = e_rescue[cell_line == "MM1R"][1],
    MM1S_e = e_rescue[cell_line == "MM1S"][1],
    supporting_cell_lines = sum(is.finite(e_rescue)),
    .groups = "drop"
  ) %>%
  arrange(desc(e_study))

write_csv(
  study_mech,
  file.path(result_dir, "GSE182638_mechanistic_study_evidence.csv")
)

# ALL-GENE rescue-concordance table for future ESR signature derivation.
# No literature direction is imposed here.
effect_strength <- effects %>%
  group_by(cell_line, contrast_id) %>%
  add_effect_magnitude() %>%
  ungroup() %>%
  mutate(
    unsigned_strength = magnitude * confidence
  ) %>%
  select(
    gene, cell_line, contrast_id,
    beta, lfsr, magnitude, confidence, unsigned_strength
  ) %>%
  pivot_wider(
    names_from = contrast_id,
    values_from = c(beta, lfsr, magnitude, confidence, unsigned_strength)
  ) %>%
  mutate(
    same_rescue_direction =
      sign(beta_RSL3_vs_Untreated) ==
      sign(beta_RSL3_vs_RSL3_Fer1) &
      is.finite(beta_RSL3_vs_Untreated) &
      is.finite(beta_RSL3_vs_RSL3_Fer1) &
      beta_RSL3_vs_Untreated != 0 &
      beta_RSL3_vs_RSL3_Fer1 != 0,
    rescue_strength = ifelse(
      same_rescue_direction,
      pmin(
        unsigned_strength_RSL3_vs_Untreated,
        unsigned_strength_RSL3_vs_RSL3_Fer1
      ),
      0
    ),
    learned_direction = ifelse(
      same_rescue_direction,
      sign(beta_RSL3_vs_Untreated),
      0
    ),
    signed_rescue_score = learned_direction * rescue_strength
  )

write_csv(
  effect_strength,
  file.path(result_dir, "per_cell_line_all_gene_rescue_concordance.csv")
)

study_esr_candidates <- effect_strength %>%
  group_by(gene) %>%
  summarise(
    study_signed_rescue_score = median(signed_rescue_score, na.rm = TRUE),
    study_rescue_strength = abs(study_signed_rescue_score),
    MM1R_signed = signed_rescue_score[cell_line == "MM1R"][1],
    MM1S_signed = signed_rescue_score[cell_line == "MM1S"][1],
    cross_cell_direction_agreement =
      is.finite(MM1R_signed) &
      is.finite(MM1S_signed) &
      sign(MM1R_signed) == sign(MM1S_signed) &
      MM1R_signed != 0 &
      MM1S_signed != 0,
    .groups = "drop"
  ) %>%
  arrange(desc(study_rescue_strength))

write_csv(
  study_esr_candidates,
  file.path(result_dir, "GSE182638_empirical_ESR_candidates.csv")
)

# -------------------------------------------------------------------------
# 9. Audit summary
# -------------------------------------------------------------------------

message("[9/9] Writing audit summary...")

core_genes <- read_csv(
  file.path("ferroptosis-only", "config", "core_gene_review.csv"),
  show_col_types = FALSE
) %>%
  filter(review_decision == "CORE") %>%
  pull(gene)

summary_lines <- c(
  paste0("Accession: ", accession),
  paste0("Count file: ", basename(count_files[1])),
  paste0("Genes after ID standardization: ", nrow(count_mat)),
  paste0("Samples: ", ncol(count_mat)),
  paste0("Integer fraction: ", round(integer_fraction, 5)),
  paste0(
    "RNA-eligible mechanistic genes present: ",
    sum(presence$present_in_GSE182638), "/", nrow(presence)
  ),
  paste0(
    "Provisional CORE genes present: ",
    paste(core_genes[core_genes %in% rownames(count_mat)], collapse = ", ")
  ),
  "",
  "Primary contrasts:",
  "  RSL3 - Untreated",
  "  RSL3 - (RSL3 + Ferrostatin-1)",
  "",
  "Model within each cell line: ~ replicate + condition",
  "",
  "Important interpretation:",
  "  e_rescue > 0 supports the literature-expected RNA direction AND Fer-1 reversal.",
  "  e_rescue < 0 is contradictory evidence.",
  "  The all-gene ESR candidate file is exploratory and must not be frozen from this single study."
)

writeLines(
  summary_lines,
  file.path(result_dir, "analysis_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo.txt")
)

message("GSE182638 analysis completed.")
message("Please send back the result files in: ", result_dir)
