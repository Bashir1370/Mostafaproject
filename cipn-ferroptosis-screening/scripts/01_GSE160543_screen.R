#!/usr/bin/env Rscript

# CIPN ferroptosis screening — GSE160543
#
# Practical question:
#   Does rat DRG from oxaliplatin- or paclitaxel-induced peripheral neuropathy
#   show transcriptomic evidence consistent with ferroptosis-associated activation?
#
# Primary evidence:
#   Vinik et al. 2024 validated 24-gene ferroptosis-vs-apoptosis biomarker panel
#
# Supporting evidence:
#   Gene Ontology Biological Process: Ferroptosis
#   WikiPathways Ferroptosis
#
# Context / specificity controls:
#   HALLMARK_APOPTOSIS
#   HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
#
# Run from repository root:
#   source("cipn-ferroptosis-screening/scripts/01_GSE160543_screen.R")

required_pkgs <- c(
  "DESeq2", "dplyr", "tidyr", "readr", "tibble", "data.table",
  "ggplot2", "pheatmap", "fgsea",
  "AnnotationDbi", "org.Rn.eg.db"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them and rerun the script."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
})

accession <- "GSE160543"

project_dir <- "cipn-ferroptosis-screening"
config_dir <- file.path(project_dir, "config")
raw_dir <- file.path(project_dir, "data", "raw", accession)
extract_dir <- file.path(raw_dir, "extracted")
result_dir <- file.path(project_dir, "results", accession)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(config_dir, "GSE160543_samples.csv")
vinik_file <- file.path(config_dir, "vinik_2024_24_biomarkers.csv")
vinik_mapping_file <- file.path(config_dir, "vinik_2024_24_human_to_rat_locked.csv")
msigdb_locked_file <- file.path(config_dir, "msigdb_2026.1_Hs_rat_gene_sets_locked.csv")
mechanistic_mapping_file <- file.path(config_dir, "mechanistic_panel_human_to_rat_locked.csv")
gene_set_manifest_file <- file.path(config_dir, "gene_set_manifest.csv")

# -------------------------------------------------------------------------
# 1. Locked metadata
# -------------------------------------------------------------------------

message("[1/10] Loading locked GSE160543 metadata...")

meta <- readr::read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    group = factor(
      group,
      levels = c("Vehicle", "Paclitaxel", "Oxaliplatin")
    ),
    replicate = factor(replicate)
  )

design_check <- meta %>% count(group, name = "n")

if (
  nrow(meta) != 12 ||
  nrow(design_check) != 3 ||
  any(design_check$n != 4)
) {
  stop("Locked GSE160543 metadata does not match the expected 3 groups x 4 replicates.")
}

write_csv(meta, file.path(result_dir, "sample_metadata.csv"))

# -------------------------------------------------------------------------
# 2. Download the small GEO processed-count archive
# -------------------------------------------------------------------------

message("[2/10] Downloading/extracting GEO processed count files...")

tar_name <- paste0(accession, "_RAW.tar")
tar_path <- file.path(raw_dir, tar_name)
tar_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/",
  accession,
  "/suppl/",
  tar_name
)

options(timeout = max(1200, getOption("timeout")))

if (!file.exists(tar_path) || file.info(tar_path)$size == 0) {
  utils::download.file(
    url = tar_url,
    destfile = tar_path,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )
} else {
  message("Processed GEO archive already exists; skipping download.")
}

if (length(list.files(extract_dir, recursive = TRUE)) == 0) {
  utils::untar(tar_path, exdir = extract_dir)
} else {
  message("Archive already extracted; skipping extraction.")
}

all_files <- list.files(
  extract_dir,
  recursive = TRUE,
  full.names = TRUE
)

if (length(all_files) < 12) {
  stop("Too few files were extracted from GSE160543_RAW.tar.")
}

# -------------------------------------------------------------------------
# 3. Reconstruct gene-level raw-count matrix
# -------------------------------------------------------------------------

message("[3/10] Reconstructing the 12-sample raw-count matrix...")

clean_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}

integer_fraction <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  if (!any(ok)) return(0)
  mean(abs(x[ok] - round(x[ok])) < 1e-8)
}

nonnegative_fraction <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  if (!any(ok)) return(0)
  mean(x[ok] >= 0)
}

choose_gene_column <- function(tbl) {
  nms <- names(tbl)
  low <- tolower(nms)

  preferred <- grep(
    "gene.?symbol|symbol|gene.?id|geneid|gene.?name",
    low
  )

  if (length(preferred) > 0) return(nms[preferred[1]])

  non_numeric <- nms[
    !vapply(tbl, is.numeric, logical(1))
  ]

  if (length(non_numeric) > 0) return(non_numeric[1])

  nms[1]
}

choose_count_column <- function(tbl, gene_col) {
  nms <- setdiff(names(tbl), gene_col)
  low <- tolower(nms)

  forbidden <- grepl(
    "fpkm|rpkm|tpm|cpm|length|len|ratio|fold|log",
    low
  )

  candidate_names <- nms[!forbidden]

  if (length(candidate_names) == 0) {
    stop("No candidate count columns remain after excluding normalized-expression columns.")
  }

  explicit <- candidate_names[
    grepl(
      "raw.?count|read.?count|expected.?count|counts|count",
      tolower(candidate_names)
    )
  ]

  score_column <- function(nm) {
    x <- suppressWarnings(as.numeric(tbl[[nm]]))
    c(
      integer = integer_fraction(x),
      nonnegative = nonnegative_fraction(x),
      finite = mean(is.finite(x))
    )
  }

  pool <- if (length(explicit) > 0) explicit else candidate_names

  scores <- t(vapply(
    pool,
    score_column,
    numeric(3)
  ))

  score_tbl <- tibble(
    column = pool,
    integer_fraction = scores[, "integer"],
    nonnegative_fraction = scores[, "nonnegative"],
    finite_fraction = scores[, "finite"]
  ) %>%
    arrange(
      desc(integer_fraction),
      desc(nonnegative_fraction),
      desc(finite_fraction)
    )

  selected <- score_tbl$column[1]

  if (
    score_tbl$nonnegative_fraction[1] < 0.99 ||
    score_tbl$finite_fraction[1] < 0.99
  ) {
    stop(
      "Could not identify a valid nonnegative raw-count column. ",
      "Top candidate: ", selected,
      "; nonnegative fraction = ", round(score_tbl$nonnegative_fraction[1], 4),
      "; finite fraction = ", round(score_tbl$finite_fraction[1], 4)
    )
  }

  list(selected = selected, scores = score_tbl)
}

standardize_rat_gene_ids <- function(ids) {
  ids <- as.character(ids)
  ids_no_ver <- sub("\\.[0-9]+$", "", ids)

  ens_frac <- mean(
    grepl("^ENSRNOG[0-9]+$", ids_no_ver),
    na.rm = TRUE
  )
  entrez_frac <- mean(
    grepl("^[0-9]+$", ids_no_ver),
    na.rm = TRUE
  )

  if (is.finite(ens_frac) && ens_frac > 0.50) {
    mapped <- AnnotationDbi::mapIds(
      org.Rn.eg.db,
      keys = unique(ids_no_ver),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )

    gene <- unname(mapped[ids_no_ver])
    id_type <- "ENSEMBL"
  } else if (is.finite(entrez_frac) && entrez_frac > 0.50) {
    mapped <- AnnotationDbi::mapIds(
      org.Rn.eg.db,
      keys = unique(ids_no_ver),
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )

    gene <- unname(mapped[ids_no_ver])
    id_type <- "ENTREZID"
  } else {
    gene <- ids_no_ver
    id_type <- "SYMBOL_OR_OTHER"
  }

  tibble(
    original_id = ids,
    normalized_id = ids_no_ver,
    gene = gene,
    detected_id_type = id_type
  )
}

sample_tables <- list()
audit_rows <- list()
id_maps <- list()

for (i in seq_len(nrow(meta))) {
  gsm <- meta$gsm[i]

  hits <- all_files[
    grepl(gsm, basename(all_files), fixed = TRUE)
  ]

  if (length(hits) != 1) {
    stop(
      "Expected exactly one processed file for ", gsm,
      "; found ", length(hits), "."
    )
  }

  f <- hits[1]
  x <- data.table::fread(
    f,
    data.table = FALSE,
    check.names = FALSE
  )

  if (nrow(x) < 1000 || ncol(x) < 2) {
    stop("Unexpected dimensions in ", basename(f), ".")
  }

  gene_col <- choose_gene_column(x)
  count_choice <- choose_count_column(x, gene_col)
  count_col <- count_choice$selected

  ids <- as.character(x[[gene_col]])
  id_map <- standardize_rat_gene_ids(ids)

  counts <- suppressWarnings(as.numeric(x[[count_col]]))

  keep <- (
    !is.na(id_map$gene) &
    nzchar(id_map$gene) &
    is.finite(counts) &
    counts >= 0
  )

  one <- tibble(
    gene = id_map$gene[keep],
    count = round(counts[keep])
  ) %>%
    group_by(gene) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )

  names(one)[2] <- gsm
  sample_tables[[gsm]] <- one

  audit_rows[[gsm]] <- tibble(
    gsm = gsm,
    file = basename(f),
    rows = nrow(x),
    columns = ncol(x),
    gene_column = gene_col,
    count_column = count_col,
    count_integer_fraction = integer_fraction(counts),
    count_nonnegative_fraction = nonnegative_fraction(counts),
    detected_id_type = unique(id_map$detected_id_type)[1]
  )

  id_maps[[gsm]] <- id_map %>%
    mutate(gsm = gsm)
}

audit <- bind_rows(audit_rows)

write_csv(
  audit,
  file.path(result_dir, "raw_file_column_audit.csv")
)

write_csv(
  bind_rows(id_maps),
  file.path(result_dir, "gene_identifier_mapping.csv")
)

count_tbl <- Reduce(
  function(a, b) full_join(a, b, by = "gene"),
  sample_tables
) %>%
  mutate(across(-gene, ~ replace_na(.x, 0)))

count_tbl <- count_tbl %>%
  dplyr::select(gene, all_of(meta$gsm))

count_mat <- as.matrix(count_tbl[, -1, drop = FALSE])
rownames(count_mat) <- count_tbl$gene
storage.mode(count_mat) <- "numeric"
count_mat <- round(count_mat)

if (ncol(count_mat) != 12 || nrow(count_mat) < 5000) {
  stop("Reconstructed count matrix has unexpected dimensions.")
}

if (any(count_mat < 0) || mean(abs(count_mat - round(count_mat)) < 1e-8) < 0.999) {
  stop("Reconstructed matrix is not a valid nonnegative integer count matrix.")
}

write_csv(
  as_tibble(count_mat, rownames = "gene"),
  file.path(result_dir, "count_matrix_gene_symbols.csv")
)

# -------------------------------------------------------------------------
# 4. DESeq2 and QC
# -------------------------------------------------------------------------

message("[4/10] Running DESeq2 and sample-level QC...")

coldata <- as.data.frame(meta)
rownames(coldata) <- coldata$gsm

keep <- rowSums(count_mat >= 10) >= 4
filtered_counts <- count_mat[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = coldata,
  design = ~ group
)

dds <- DESeq(dds, quiet = TRUE)

library_qc <- tibble(
  gsm = colnames(count_mat),
  library_size = colSums(count_mat),
  detected_genes = colSums(count_mat > 0),
  genes_count_ge_10 = colSums(count_mat >= 10)
) %>%
  left_join(
    meta %>% mutate(gsm = as.character(gsm)),
    by = "gsm"
  )

write_csv(
  library_qc,
  file.path(result_dir, "library_qc.csv")
)

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

pca <- prcomp(t(vst_mat))

pca_tbl <- as_tibble(
  pca$x[, 1:2, drop = FALSE],
  rownames = "gsm"
) %>%
  left_join(
    meta %>% mutate(gsm = as.character(gsm)),
    by = "gsm"
  )

write_csv(
  pca_tbl,
  file.path(result_dir, "pca_coordinates.csv")
)

p_pca <- ggplot(
  pca_tbl,
  aes(PC1, PC2, shape = group)
) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "GSE160543 DRG — VST PCA")

ggsave(
  file.path(result_dir, "pca.png"),
  p_pca,
  width = 7,
  height = 5,
  dpi = 180
)

sample_cor <- cor(vst_mat, method = "pearson")

png(
  file.path(result_dir, "sample_correlation.png"),
  width = 1600,
  height = 1400,
  res = 180
)
pheatmap(
  sample_cor,
  main = "GSE160543 VST sample correlation"
)
dev.off()

extract_de <- function(dds, treated, control = "Vehicle") {
  r <- results(
    dds,
    contrast = c("group", treated, control),
    independentFiltering = FALSE
  )

  tibble(
    gene = rownames(r),
    baseMean = r$baseMean,
    log2FoldChange = r$log2FoldChange,
    lfcSE = r$lfcSE,
    stat = r$stat,
    pvalue = r$pvalue,
    padj = r$padj
  )
}

de_oxa <- extract_de(dds, "Oxaliplatin")
de_ptx <- extract_de(dds, "Paclitaxel")

write_csv(
  de_oxa,
  file.path(result_dir, "DE_Oxaliplatin_vs_Vehicle.csv")
)

write_csv(
  de_ptx,
  file.path(result_dir, "DE_Paclitaxel_vs_Vehicle.csv")
)

# -------------------------------------------------------------------------
# 5. Build independent ferroptosis/control signatures
# -------------------------------------------------------------------------

message("[5/10] Loading LOCKED, audited ferroptosis/control signatures...")

# IMPORTANT:
# Inferential analysis never retrieves live MSigDB memberships or live orthologs.
# Provenance/species/version are locked in config/gene_set_manifest.csv and
# docs/GENE_SET_AUDIT.md. Update those files explicitly before changing sets.

manifest <- read_csv(gene_set_manifest_file, show_col_types = FALSE)
locked_msig <- read_csv(msigdb_locked_file, show_col_types = FALSE)
vinik <- read_csv(vinik_file, show_col_types = FALSE)
vinik_map <- read_csv(vinik_mapping_file, show_col_types = FALSE)

required_manifest <- c(
  "VINIK_2024_24_FERROPTOSIS_BIOMARKERS",
  "GOBP_FERROPTOSIS",
  "WP_FERROPTOSIS",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)

if (!all(required_manifest %in% manifest$signature)) {
  stop("gene_set_manifest.csv is missing one or more required signatures.")
}

if (
  nrow(vinik) != 24 ||
  nrow(vinik_map) != 24 ||
  anyDuplicated(vinik_map$human_symbol) ||
  any(is.na(vinik_map$dataset_gene_symbol)) ||
  any(!nzchar(vinik_map$dataset_gene_symbol))
) {
  stop("Locked Vinik-24 mapping is incomplete or malformed.")
}

if (!setequal(vinik$human_symbol, vinik_map$human_symbol)) {
  stop("Vinik source list and locked human-to-rat mapping do not contain the same 24 human genes.")
}

vinik_rat <- unique(vinik_map$dataset_gene_symbol)

get_locked_set <- function(name) {
  unique(
    locked_msig$rat_gene_symbol[
      locked_msig$signature == name &
      !is.na(locked_msig$rat_gene_symbol) &
      nzchar(locked_msig$rat_gene_symbol)
    ]
  )
}

gobp_genes <- get_locked_set("GOBP_FERROPTOSIS")
wp_genes <- get_locked_set("WP_FERROPTOSIS")
apoptosis_genes <- get_locked_set("HALLMARK_APOPTOSIS")
ros_genes <- get_locked_set("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")

expected_sizes <- c(
  GOBP_FERROPTOSIS = 30L,
  WP_FERROPTOSIS = 67L,
  HALLMARK_APOPTOSIS = 160L,
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = 50L
)

observed_sizes <- c(
  GOBP_FERROPTOSIS = length(gobp_genes),
  WP_FERROPTOSIS = length(wp_genes),
  HALLMARK_APOPTOSIS = length(apoptosis_genes),
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = length(ros_genes)
)

if (!identical(as.integer(observed_sizes), as.integer(expected_sizes))) {
  stop(
    "Locked MSigDB gene-set sizes changed unexpectedly. Observed: ",
    paste(names(observed_sizes), observed_sizes, sep = "=", collapse = "; ")
  )
}

signatures <- list(
  VINIK_2024_24_FERROPTOSIS_BIOMARKERS = vinik_rat,
  GOBP_FERROPTOSIS = gobp_genes,
  WP_FERROPTOSIS = wp_genes,
  HALLMARK_APOPTOSIS = apoptosis_genes,
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = ros_genes
)

signature_audit <- bind_rows(
  lapply(names(signatures), function(nm) {
    tibble(
      signature = nm,
      gene = unique(signatures[[nm]])
    )
  })
) %>%
  mutate(
    present_in_GSE160543 = gene %in% rownames(dds)
  )

vinik_present <- sum(
  signature_audit$signature == "VINIK_2024_24_FERROPTOSIS_BIOMARKERS" &
  signature_audit$present_in_GSE160543
)

if (vinik_present != 24L) {
  warning(
    "Audited Vinik-24 mapping is complete, but only ",
    vinik_present,
    "/24 genes passed the current DESeq2 expression filter."
  )
}

write_csv(
  vinik_map,
  file.path(result_dir, "vinik24_ortholog_mapping.csv")
)

write_csv(
  manifest,
  file.path(result_dir, "gene_set_manifest_used.csv")
)

write_csv(
  signature_audit,
  file.path(result_dir, "signature_members_used.csv")
)

# -------------------------------------------------------------------------
# 6. Whole-transcriptome GSEA
# -------------------------------------------------------------------------

message("[6/10] Running GSEA on both CIPN contrasts...")

run_gsea <- function(de_tbl, contrast_name) {
  rank_tbl <- de_tbl %>%
    filter(
      !is.na(gene),
      nzchar(gene),
      is.finite(stat)
    ) %>%
    arrange(desc(abs(stat))) %>%
    distinct(gene, .keep_all = TRUE)

  ranks <- rank_tbl$stat
  names(ranks) <- rank_tbl$gene
  ranks <- sort(ranks, decreasing = TRUE)

  fg <- fgsea::fgseaMultilevel(
    pathways = signatures,
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    eps = 0
  )

  as_tibble(fg) %>%
    mutate(
      contrast = contrast_name,
      leadingEdge = vapply(
        leadingEdge,
        paste,
        collapse = ";",
        FUN.VALUE = character(1)
      )
    ) %>%
    dplyr::select(
      contrast,
      pathway,
      NES,
      pval,
      padj,
      size,
      leadingEdge
    )
}

gsea_oxa <- run_gsea(
  de_oxa,
  "Oxaliplatin_vs_Vehicle"
)

gsea_ptx <- run_gsea(
  de_ptx,
  "Paclitaxel_vs_Vehicle"
)

gsea_all <- bind_rows(gsea_oxa, gsea_ptx) %>%
  arrange(contrast, padj)

write_csv(
  gsea_all,
  file.path(result_dir, "gsea_all_signatures.csv")
)

# -------------------------------------------------------------------------
# 7. Sample-level signature scores
# -------------------------------------------------------------------------

message("[7/10] Computing sample-level signature scores...")

row_z <- t(scale(t(vst_mat)))
row_z[!is.finite(row_z)] <- NA_real_

score_one_signature <- function(genes, signature_name) {
  present <- intersect(genes, rownames(row_z))

  if (length(present) < 3) {
    return(tibble())
  }

  vals <- colMeans(
    row_z[present, , drop = FALSE],
    na.rm = TRUE
  )

  tibble(
    gsm = names(vals),
    signature = signature_name,
    score = as.numeric(vals),
    genes_used = length(present)
  )
}

sample_scores <- bind_rows(
  lapply(names(signatures), function(nm) {
    score_one_signature(
      signatures[[nm]],
      nm
    )
  })
) %>%
  left_join(
    meta %>%
      mutate(gsm = as.character(gsm)) %>%
      dplyr::select(gsm, group, replicate),
    by = "gsm"
  )

write_csv(
  sample_scores,
  file.path(result_dir, "sample_level_signature_scores.csv")
)

test_score <- function(tbl, treated) {
  x <- tbl$score[tbl$group == treated]
  y <- tbl$score[tbl$group == "Vehicle"]

  tibble(
    treated = treated,
    mean_treated = mean(x, na.rm = TRUE),
    mean_vehicle = mean(y, na.rm = TRUE),
    mean_difference = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    t_pvalue = tryCatch(
      t.test(x, y)$p.value,
      error = function(e) NA_real_
    ),
    wilcoxon_pvalue = tryCatch(
      wilcox.test(x, y, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  )
}

score_tests <- sample_scores %>%
  group_by(signature) %>%
  group_modify(~ bind_rows(
    test_score(.x, "Oxaliplatin"),
    test_score(.x, "Paclitaxel")
  )) %>%
  ungroup() %>%
  mutate(
    t_padj = p.adjust(t_pvalue, method = "BH"),
    wilcoxon_padj = p.adjust(wilcoxon_pvalue, method = "BH")
  )

write_csv(
  score_tests,
  file.path(result_dir, "sample_level_score_tests.csv")
)

p_scores <- ggplot(
  sample_scores,
  aes(x = group, y = score)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2) +
  facet_wrap(~ signature, scales = "free_y") +
  theme_bw() +
  labs(
    title = "GSE160543 sample-level transcriptomic signature scores",
    x = NULL,
    y = "Mean within-sample gene z-score"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

ggsave(
  file.path(result_dir, "sample_level_signature_scores.png"),
  p_scores,
  width = 12,
  height = 8,
  dpi = 180
)

# -------------------------------------------------------------------------
# 8. Supporting mechanistic ferroptosis panel
# -------------------------------------------------------------------------

message("[8/10] Generating supporting mechanistic-gene panel...")

mechanistic_map <- read_csv(
  mechanistic_mapping_file,
  show_col_types = FALSE
) %>%
  transmute(
    human_symbol = human_symbol,
    rat_symbol = rat_gene_symbol,
    mapping_status = status
  )

if (
  nrow(mechanistic_map) != 13 ||
  anyDuplicated(mechanistic_map$human_symbol) ||
  any(is.na(mechanistic_map$rat_symbol))
) {
  stop("Locked mechanistic-panel mapping is incomplete or malformed.")
}

mechanistic_effects <- mechanistic_map %>%
  left_join(
    de_oxa %>%
      transmute(
        rat_symbol = gene,
        OXA_log2FC = log2FoldChange,
        OXA_stat = stat,
        OXA_padj = padj
      ),
    by = "rat_symbol"
  ) %>%
  left_join(
    de_ptx %>%
      transmute(
        rat_symbol = gene,
        PTX_log2FC = log2FoldChange,
        PTX_stat = stat,
        PTX_padj = padj
      ),
    by = "rat_symbol"
  )

write_csv(
  mechanistic_effects,
  file.path(result_dir, "mechanistic_panel_effects.csv")
)

panel_genes <- intersect(
  mechanistic_effects$rat_symbol,
  rownames(vst_mat)
)

if (length(panel_genes) >= 3) {
  panel_mat <- vst_mat[panel_genes, , drop = FALSE]
  panel_z <- t(scale(t(panel_mat)))
  panel_z[!is.finite(panel_z)] <- 0

  ann <- data.frame(
    group = meta$group
  )
  rownames(ann) <- meta$gsm

  human_labels <- mechanistic_effects$human_symbol[
    match(
      rownames(panel_z),
      mechanistic_effects$rat_symbol
    )
  ]

  rownames(panel_z) <- paste0(
    human_labels,
    " (",
    rownames(panel_z),
    ")"
  )

  png(
    file.path(result_dir, "mechanistic_panel_heatmap.png"),
    width = 1800,
    height = 1500,
    res = 180
  )
  pheatmap(
    panel_z,
    annotation_col = ann,
    main = "Supporting ferroptosis mechanistic genes — row-scaled VST"
  )
  dev.off()
}

# -------------------------------------------------------------------------
# 9. Compact evidence classification
# -------------------------------------------------------------------------

message("[9/10] Summarizing ferroptosis evidence...")

get_gsea_row <- function(tbl, pathway_name) {
  z <- tbl[tbl$pathway == pathway_name, , drop = FALSE]
  if (nrow(z) == 0) {
    return(tibble(
      NES = NA_real_,
      pval = NA_real_,
      padj = NA_real_
    ))
  }
  z[1, c("NES", "pval", "padj")]
}

summarize_contrast <- function(gsea_tbl, contrast_name) {
  v <- get_gsea_row(
    gsea_tbl,
    "VINIK_2024_24_FERROPTOSIS_BIOMARKERS"
  )
  k <- get_gsea_row(
    gsea_tbl,
    "GOBP_FERROPTOSIS"
  )
  w <- get_gsea_row(
    gsea_tbl,
    "WP_FERROPTOSIS"
  )
  a <- get_gsea_row(
    gsea_tbl,
    "HALLMARK_APOPTOSIS"
  )
  r <- get_gsea_row(
    gsea_tbl,
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
  )

  # Directional evidence is intentionally restricted to:
  #   (1) Vinik-24, a ferroptosis-vs-apoptosis biomarker panel, and
  #   (2) WP_FERROPTOSIS, used as supporting mechanistic pathway evidence.
  #
  # GOBP_FERROPTOSIS is reported but NOT counted as a directional activation
  # vote because GO process membership contains both ferroptosis-promoting and
  # ferroptosis-limiting genes (for example GPX4 and NFE2L2).

  strong <- (
    is.finite(v$NES) &&
    v$NES > 0 &&
    is.finite(v$padj) &&
    v$padj < 0.05 &&
    is.finite(w$NES) &&
    w$NES > 0 &&
    is.finite(w$padj) &&
    w$padj < 0.05
  )

  suggestive <- (
    (
      is.finite(v$NES) &&
      v$NES > 0 &&
      is.finite(v$pval) &&
      v$pval < 0.05
    ) ||
    (
      is.finite(w$NES) &&
      w$NES > 0 &&
      is.finite(w$padj) &&
      w$padj < 0.10
    )
  )

  evidence <- if (strong) {
    "STRONG_SUPPORT"
  } else if (suggestive) {
    "SUGGESTIVE_SUPPORT"
  } else {
    "NO_TRANSCRIPTOMIC_SUPPORT"
  }

  apoptosis_flag <- if (
    is.finite(a$NES) &&
    a$NES > 0 &&
    is.finite(a$padj) &&
    a$padj < 0.05
  ) {
    "APOPTOSIS_COINCIDENT"
  } else {
    "NO_SIGNIFICANT_APOPTOSIS_ENRICHMENT"
  }

  ros_context <- if (
    is.finite(r$NES) &&
    r$NES > 0 &&
    is.finite(r$padj) &&
    r$padj < 0.05
  ) {
    "ROS_ENRICHED"
  } else {
    "NO_SIGNIFICANT_ROS_ENRICHMENT"
  }

  tibble(
    contrast = contrast_name,
    ferroptosis_evidence = evidence,
    specificity_flag = apoptosis_flag,
    ros_context = ros_context,
    VINIK_NES = v$NES,
    VINIK_FDR = v$padj,
    GOBP_NES = k$NES,
    GOBP_FDR = k$padj,
    WP_NES = w$NES,
    WP_FDR = w$padj,
    APOPTOSIS_NES = a$NES,
    APOPTOSIS_FDR = a$padj,
    ROS_NES = r$NES,
    ROS_FDR = r$padj
  )
}

summary_tbl <- bind_rows(
  summarize_contrast(
    gsea_oxa,
    "Oxaliplatin_vs_Vehicle"
  ),
  summarize_contrast(
    gsea_ptx,
    "Paclitaxel_vs_Vehicle"
  )
)

write_csv(
  summary_tbl,
  file.path(result_dir, "evidence_summary.csv")
)

# -------------------------------------------------------------------------
# 10. Human-readable audit summary
# -------------------------------------------------------------------------

message("[10/10] Writing analysis summary...")

summary_lines <- c(
  paste0("Dataset: ", accession),
  "Tissue: rat dorsal root ganglion (bulk RNA-seq)",
  "Groups: Vehicle n=4; Paclitaxel n=4; Oxaliplatin n=4",
  paste0(
    "Gene-level count matrix after ID standardization: ",
    nrow(count_mat), " genes x ", ncol(count_mat), " samples"
  ),
  paste0(
    "Genes retained for DESeq2: ",
    nrow(dds)
  ),
  "",
  "Primary transcriptomic ferroptosis signature:",
  paste0(
    "  Vinik 2024 human 24-gene ferroptosis-vs-apoptosis panel; locked rat dataset symbols present after filtering: ",
    vinik_present, "/24"
  ),
  "",
  "Gene-set provenance:",
  "  Inferential memberships are loaded from locked config CSVs, not live databases.",
  "  MSigDB-derived sets are frozen from MSigDB 2026.1.Hs / msigdbr 26.1.1 human-to-rat ortholog output.",
  "  See config/gene_set_manifest.csv and docs/GENE_SET_AUDIT.md.",
  "",
  "Supporting signatures:",
  "  Gene Ontology GOBP_FERROPTOSIS (supporting membership set; NES is not a direct activation score)",
  "  WikiPathways WP_FERROPTOSIS",
  "",
  "Context controls:",
  "  HALLMARK_APOPTOSIS",
  "  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "",
  "Evidence classifications:",
  paste0(
    "  Oxaliplatin vs Vehicle: ",
    summary_tbl$ferroptosis_evidence[
      summary_tbl$contrast == "Oxaliplatin_vs_Vehicle"
    ],
    " | ",
    summary_tbl$specificity_flag[
      summary_tbl$contrast == "Oxaliplatin_vs_Vehicle"
    ],
    " | ",
    summary_tbl$ros_context[
      summary_tbl$contrast == "Oxaliplatin_vs_Vehicle"
    ]
  ),
  paste0(
    "  Paclitaxel vs Vehicle: ",
    summary_tbl$ferroptosis_evidence[
      summary_tbl$contrast == "Paclitaxel_vs_Vehicle"
    ],
    " | ",
    summary_tbl$specificity_flag[
      summary_tbl$contrast == "Paclitaxel_vs_Vehicle"
    ],
    " | ",
    summary_tbl$ros_context[
      summary_tbl$contrast == "Paclitaxel_vs_Vehicle"
    ]
  ),
  "",
  "Interpretation boundary:",
  "  These labels summarize transcriptomic evidence only.",
  "  Bulk DRG RNA-seq cannot by itself prove ferroptotic cell death or identify the exact responding cell type.",
  "",
  "GOBP_FERROPTOSIS is reported descriptively and is not used as a directional activation vote in the evidence classifier.",
  "The STRONG/SUGGESTIVE/NO_SUPPORT rules are pre-specified heuristic interpretation rules, not a clinically validated classifier."
)

writeLines(
  summary_lines,
  file.path(result_dir, "analysis_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo.txt")
)

message("GSE160543 CIPN ferroptosis screening completed.")
message("Results: ", result_dir)
