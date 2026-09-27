#!/usr/bin/env Rscript

# GSE287284 — direct neuronal ferroptosis validation
#
# Design:
#   RA-differentiated mouse N2a neuron-like cells
#   DMSO n=3 vs RSL3 n=3
#
# Input:
#   GEO per-sample processed mRNA FPKM files (~1.2 MB total series supplement)
#
# Statistical strategy:
#   log2(FPKM + 1) -> limma -> posterior effect with ashr
#
# Important:
#   This is a neuronal induction validation dataset, NOT a rescue dataset.
#   Frozen ESR candidates come only from GSE182638.
#
# Run from repository root:
#   source("ferroptosis-only/scripts/03_GSE287284_neuronal_validation.R")

required_pkgs <- c(
  "dplyr", "tidyr", "readr", "tibble", "ggplot2", "pheatmap",
  "limma", "ashr", "AnnotationDbi", "org.Mm.eg.db", "biomaRt"
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
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(limma)
  library(ashr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(biomaRt)
})

source("R/evidence_functions.R")

accession <- "GSE287284"
sample_file <- file.path(
  "ferroptosis-only", "config", "GSE287284_samples.csv"
)
raw_dir <- file.path(
  "ferroptosis-only", "data", "raw", accession
)
result_dir <- file.path(
  "ferroptosis-only", "results", accession
)
primary_candidate_file <- file.path(
  "ferroptosis-only", "results", "GSE182638",
  "GSE182638_empirical_ESR_candidates.csv"
)
ortholog_cache <- file.path(
  "ferroptosis-only", "data", "reference",
  "human_mouse_orthologs_ensembl101.csv"
)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ortholog_cache), recursive = TRUE, showWarnings = FALSE)

message("[1/9] Loading locked metadata...")

meta <- read_csv(sample_file, show_col_types = FALSE) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "RSL3")),
    replicate = factor(replicate)
  )

if (
  nrow(meta) != 6 ||
  sum(meta$condition == "DMSO") != 3 ||
  sum(meta$condition == "RSL3") != 3
) {
  stop("GSE287284 metadata does not match the locked 3-vs-3 design.")
}

write_csv(meta, file.path(result_dir, "sample_metadata_locked.csv"))

message("[2/9] Downloading processed mRNA FPKM files directly from GEO...")

geo_sample_bucket <- function(gsm) {
  paste0(substr(gsm, 1, 7), "nnn")
}

for (i in seq_len(nrow(meta))) {
  gsm <- meta$gsm[i]
  fn <- meta$processed_filename[i]
  dest <- file.path(raw_dir, fn)

  if (!file.exists(dest) || file.info(dest)$size == 0) {
    url <- paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/samples/",
      geo_sample_bucket(gsm), "/",
      gsm, "/suppl/", fn
    )

    message("Downloading ", fn, " ...")
    utils::download.file(
      url = url,
      destfile = dest,
      mode = "wb",
      method = "libcurl",
      quiet = FALSE
    )
  }
}

message("[3/9] Reading and standardizing gene-level FPKM...")

pick_gene_column <- function(tbl) {
  nms <- names(tbl)
  low <- tolower(nms)

  preferred <- which(
    grepl("gene.?symbol|symbol|gene.?name", low)
  )
  if (length(preferred) > 0) return(nms[preferred[1]])

  geneish <- which(grepl("gene", low))
  if (length(geneish) > 0) return(nms[geneish[1]])

  nms[1]
}

pick_fpkm_column <- function(tbl, gene_col) {
  nms <- names(tbl)
  low <- tolower(nms)

  f <- which(grepl("fpkm", low))
  if (length(f) > 0) return(nms[f[1]])

  numeric_cols <- nms[
    vapply(tbl, is.numeric, logical(1))
  ]
  numeric_cols <- setdiff(numeric_cols, gene_col)

  if (length(numeric_cols) == 0) {
    stop("Could not identify an FPKM/numeric expression column.")
  }

  numeric_cols[length(numeric_cols)]
}

standardize_mouse_gene <- function(ids) {
  ids <- as.character(ids)
  ids_no_ver <- sub("\\.[0-9]+$", "", ids)

  ens_frac <- mean(
    grepl("^ENSMUSG[0-9]+$", ids_no_ver),
    na.rm = TRUE
  )

  if (is.finite(ens_frac) && ens_frac > 0.5) {
    mapped <- AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = unique(ids_no_ver),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    return(unname(mapped[ids_no_ver]))
  }

  ids_no_ver
}

sample_tables <- lapply(seq_len(nrow(meta)), function(i) {
  gsm <- meta$gsm[i]
  path <- file.path(raw_dir, meta$processed_filename[i])

  x <- read_csv(path, show_col_types = FALSE)
  gene_col <- pick_gene_column(x)
  fpkm_col <- pick_fpkm_column(x, gene_col)

  out <- tibble(
    gene_raw = as.character(x[[gene_col]]),
    value = suppressWarnings(as.numeric(x[[fpkm_col]]))
  ) %>%
    mutate(
      gene = standardize_mouse_gene(gene_raw)
    ) %>%
    filter(
      !is.na(gene), nzchar(gene),
      is.finite(value), value >= 0
    ) %>%
    group_by(gene) %>%
    summarise(
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )

  names(out)[2] <- gsm
  out
})

expr_tbl <- Reduce(
  function(x, y) full_join(x, y, by = "gene"),
  sample_tables
)

expr_mat <- as.matrix(expr_tbl[, meta$gsm, drop = FALSE])
rownames(expr_mat) <- expr_tbl$gene
storage.mode(expr_mat) <- "numeric"
expr_mat[!is.finite(expr_mat)] <- 0

if (nrow(expr_mat) < 5000) {
  stop("Unexpectedly few genes after FPKM import.")
}

write_csv(
  as_tibble(expr_mat, rownames = "mouse_gene"),
  file.path(result_dir, "GSE287284_FPKM_matrix.csv")
)

message("[4/9] Running expression QC...")

qc <- tibble(
  gsm = colnames(expr_mat),
  total_FPKM = colSums(expr_mat),
  detected_genes = colSums(expr_mat > 0),
  genes_FPKM_ge_1 = colSums(expr_mat >= 1)
) %>%
  left_join(meta, by = "gsm")

write_csv(qc, file.path(result_dir, "expression_qc.csv"))

keep <- rowSums(expr_mat >= 0.5) >= 3
analysis_mat <- expr_mat[keep, , drop = FALSE]
log_expr <- log2(analysis_mat + 1)

vars <- apply(log_expr, 1, var)
top <- order(vars, decreasing = TRUE)[
  seq_len(min(5000, length(vars)))
]

pca <- prcomp(t(log_expr[top, , drop = FALSE]))

pca_tbl <- as_tibble(
  pca$x[, 1:2, drop = FALSE],
  rownames = "gsm"
) %>%
  left_join(meta, by = "gsm")

write_csv(
  pca_tbl,
  file.path(result_dir, "logFPKM_pca_coordinates.csv")
)

p <- ggplot(pca_tbl, aes(PC1, PC2, shape = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "GSE287284 log2(FPKM+1) PCA")

ggsave(
  file.path(result_dir, "logFPKM_pca.png"),
  p,
  width = 7,
  height = 5,
  dpi = 180
)

corr <- cor(log_expr, method = "pearson")
png(
  file.path(result_dir, "logFPKM_sample_correlation.png"),
  width = 1400,
  height = 1200,
  res = 180
)
pheatmap::pheatmap(
  corr,
  main = "GSE287284 sample correlation"
)
dev.off()

message("[5/9] Fitting limma model and posterior RSL3 effect...")

design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- sub("^condition", "", colnames(design))

fit <- limma::lmFit(log_expr, design)
cont <- limma::makeContrasts(
  RSL3_vs_DMSO = RSL3 - DMSO,
  levels = design
)
fit2 <- limma::contrasts.fit(fit, cont)
fit2 <- limma::eBayes(fit2)

beta_mle <- fit2$coefficients[, 1]
se_mle <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1]

ok <- is.finite(beta_mle) &
  is.finite(se_mle) &
  se_mle > 0

effects <- tibble(
  gene = rownames(log_expr),
  beta_mle = beta_mle,
  se_mle = se_mle,
  beta = NA_real_,
  posterior_sd = NA_real_,
  lfsr = NA_real_,
  p_value = fit2$p.value[, 1],
  adj_p_value = p.adjust(fit2$p.value[, 1], method = "BH")
)

ash_fit <- ashr::ash(
  betahat = beta_mle[ok],
  sebetahat = se_mle[ok],
  mixcompdist = "normal",
  method = "shrink"
)

effects$beta[ok] <- ash_fit$result$PosteriorMean
effects$posterior_sd[ok] <- ash_fit$result$PosteriorSD
effects$lfsr[ok] <- ash_fit$result$lfsr

write_csv(
  effects,
  file.path(result_dir, "RSL3_vs_DMSO_mouse_gene_effects.csv")
)

message("[6/9] Freezing GSE182638 ESR discovery candidates...")

if (!file.exists(primary_candidate_file)) {
  stop(
    "Missing GSE182638 discovery result: ",
    primary_candidate_file
  )
}

primary <- read_csv(
  primary_candidate_file,
  show_col_types = FALSE
)

frozen <- primary %>%
  filter(
    cross_cell_direction_agreement == TRUE,
    is.finite(study_signed_rescue_score),
    study_signed_rescue_score != 0
  ) %>%
  mutate(
    discovery_direction = sign(study_signed_rescue_score)
  ) %>%
  dplyr::select(
    gene,
    discovery_direction,
    study_signed_rescue_score,
    study_rescue_strength,
    MM1R_signed,
    MM1S_signed
  )

write_csv(
  frozen,
  file.path(
    result_dir,
    "GSE182638_ESR_discovery_snapshot.csv"
  )
)

message("[7/9] Building/caching one-to-one human-mouse ortholog map...")

all_human_genes <- unique(c(
  frozen$gene,
  read_csv(
    file.path(
      "ferroptosis-only",
      "config",
      "ferroptosis_master_gene_universe.csv"
    ),
    show_col_types = FALSE
  )$gene
))

if (!file.exists(ortholog_cache)) {
  human_mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl",
    version = 101
  )

  mouse_mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "mmusculus_gene_ensembl",
    version = 101
  )

  ortho_raw <- biomaRt::getLDS(
    attributes = "hgnc_symbol",
    filters = "hgnc_symbol",
    values = all_human_genes,
    mart = human_mart,
    attributesL = "mgi_symbol",
    martL = mouse_mart,
    uniqueRows = TRUE
  )

  ortho_raw <- as_tibble(ortho_raw)
  names(ortho_raw)[1:2] <- c(
    "human_gene",
    "mouse_gene"
  )

  ortho_one_to_one <- ortho_raw %>%
    filter(
      nzchar(human_gene),
      nzchar(mouse_gene)
    ) %>%
    distinct() %>%
    group_by(human_gene) %>%
    filter(n_distinct(mouse_gene) == 1) %>%
    ungroup() %>%
    group_by(mouse_gene) %>%
    filter(n_distinct(human_gene) == 1) %>%
    ungroup()

  write_csv(
    ortho_one_to_one,
    ortholog_cache
  )
}

orthologs <- read_csv(
  ortholog_cache,
  show_col_types = FALSE
)

message("[8/9] Testing frozen ESR candidates in neuronal RSL3 response...")

effects_human <- effects %>%
  rename(mouse_gene = gene) %>%
  inner_join(orthologs, by = "mouse_gene") %>%
  transmute(
    gene = human_gene,
    mouse_gene,
    beta_mle,
    se_mle,
    beta,
    posterior_sd,
    lfsr,
    p_value,
    adj_p_value
  )

dictionary <- frozen %>%
  transmute(
    gene,
    direction = discovery_direction
  )

neural_scored <- score_dictionary_genes(
  effects_human,
  dictionary
) %>%
  dplyr::select(
    gene,
    mouse_gene,
    direction,
    beta_mle,
    beta,
    lfsr,
    e
  ) %>%
  rename(
    discovery_direction = direction,
    neuronal_induction_e = e
  )

validation <- frozen %>%
  left_join(
    orthologs,
    by = c("gene" = "human_gene")
  ) %>%
  left_join(
    neural_scored %>%
      dplyr::select(
        gene,
        mouse_gene,
        beta_mle,
        beta,
        lfsr,
        neuronal_induction_e
      ),
    by = c("gene", "mouse_gene")
  ) %>%
  mutate(
    neural_status = case_when(
      is.na(mouse_gene) ~ "NO_ONE_TO_ONE_ORTHOLOG",
      !is.finite(neuronal_induction_e) ~ "NOT_TESTABLE",
      neuronal_induction_e > 0 ~ "SUPPORTS_DISCOVERY_DIRECTION",
      neuronal_induction_e < 0 ~ "CONTRADICTS_DISCOVERY_DIRECTION",
      TRUE ~ "NO_DIRECTIONAL_EVIDENCE"
    )
  ) %>%
  arrange(desc(neuronal_induction_e))

write_csv(
  validation,
  file.path(
    result_dir,
    "GSE287284_frozen_ESR_neuronal_validation.csv"
  )
)

# Mechanistic genes are reported descriptively, without assuming causal
# promoter/suppressor direction equals acute RSL3 transcript direction.
master <- read_csv(
  file.path(
    "ferroptosis-only",
    "config",
    "ferroptosis_master_gene_universe.csv"
  ),
  show_col_types = FALSE
)

mechanistic <- master %>%
  dplyr::select(
    gene,
    module,
    role,
    provisional_status,
    evidence_tier,
    mechanistic_prior_M
  ) %>%
  left_join(
    orthologs,
    by = c("gene" = "human_gene")
  ) %>%
  left_join(
    effects %>%
      rename(
        mouse_gene = gene,
        neuronal_beta = beta,
        neuronal_lfsr = lfsr,
        neuronal_beta_mle = beta_mle
      ) %>%
      dplyr::select(
        mouse_gene,
        neuronal_beta_mle,
        neuronal_beta,
        neuronal_lfsr
      ),
    by = "mouse_gene"
  )

write_csv(
  mechanistic,
  file.path(
    result_dir,
    "mechanistic_genes_neuronal_RSL3_effects.csv"
  )
)

message("[9/9] Writing audit summary...")

n_support <- sum(
  validation$neural_status ==
    "SUPPORTS_DISCOVERY_DIRECTION",
  na.rm = TRUE
)
n_contradict <- sum(
  validation$neural_status ==
    "CONTRADICTS_DISCOVERY_DIRECTION",
  na.rm = TRUE
)
n_no_ortho <- sum(
  validation$neural_status ==
    "NO_ONE_TO_ONE_ORTHOLOG",
  na.rm = TRUE
)
n_testable <- sum(
  validation$neural_status %in% c(
    "SUPPORTS_DISCOVERY_DIRECTION",
    "CONTRADICTS_DISCOVERY_DIRECTION",
    "NO_DIRECTIONAL_EVIDENCE"
  ),
  na.rm = TRUE
)

summary_lines <- c(
  paste0("Accession: ", accession),
  "Role: direct neuronal ferroptosis validation",
  "Species: Mus musculus",
  "Model: retinoic-acid differentiated N2a neuron-like cells",
  "Conditions: DMSO n=3; RSL3 n=3",
  "RSL3: 5 uM for 24 h",
  "",
  "Input: GEO processed gene-level mRNA FPKM",
  "Transformation: log2(FPKM + 1)",
  "Model: limma RSL3 - DMSO",
  "Posterior effect/confidence: ashr",
  "",
  paste0("Genes retained after expression filter: ", nrow(log_expr)),
  paste0("Frozen GSE182638 ESR candidates: ", nrow(frozen)),
  paste0("Testable one-to-one ortholog candidates: ", n_testable),
  paste0("Supports discovery direction: ", n_support),
  paste0("Contradicts discovery direction: ", n_contradict),
  paste0("No one-to-one ortholog: ", n_no_ortho),
  "",
  "Important interpretation:",
  "  This dataset has no Ferrostatin-1 rescue arm.",
  "  neuronal_induction_e tests cross-species neuronal transfer of the pre-specified GSE182638 direction.",
  "  It is not equivalent to rescue-validated evidence.",
  "  Mechanistic-gene acute RNA effects are descriptive only."
)

writeLines(
  summary_lines,
  file.path(result_dir, "analysis_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo.txt")
)

message("GSE287284 neuronal validation completed.")
message("Results: ", result_dir)
