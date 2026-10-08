#!/usr/bin/env Rscript
# Step 02: all tested genes contribute to BH; significant protein-coding genes proceed.
local({
  explicit <- Sys.getenv("OIPN_PROJECT_DIR", unset = "")
  sourced <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  args <- grep("^--file=", commandArgs(), value = TRUE)
  if (nzchar(explicit)) {
    project <- normalizePath(explicit, mustWork = TRUE)
  } else if (length(sourced)) {
    project <- dirname(dirname(normalizePath(tail(sourced, 1), mustWork = TRUE)))
  } else if (length(args)) {
    project <- dirname(dirname(normalizePath(sub("^--file=", "", args[1]), mustWork = TRUE)))
  } else stop("Use source(full_script_path) or Rscript full_script_path.")
  lib <- file.path(project, ".r-library")
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  pkgs <- c("DESeq2", "jsonlite")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  check <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
  sha256 <- function(path) {
    check(nzchar(Sys.which("sha256sum")), "Ubuntu sha256sum is required.")
    result <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
    check(is.null(attr(result, "status")), paste("Cannot checksum", path))
    strsplit(result[1], "[[:space:]]+")[[1]][1]
  }
  input <- file.path(project, "results/01_dataset_audit")
  qc <- file.path(project, "results/01b_sample_qc")
  output <- file.path(project, "results/02_differential_expression")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(output, "SUCCESS.txt"))
  stage <- tempfile("de_stage_", tmpdir = output)
  dir.create(stage)
  csv <- function(x, file) write.csv(x, file.path(stage, file), row.names = FALSE, na = "NA")
  tryCatch({
    message("[1/5] Checking frozen sample design and count/QC provenance...")
    cfg_path <- file.path(project, "config/analysis_parameters.json")
    design_path <- file.path(project, "config/GSE286387_design.json")
    manifest_path <- file.path(project, "config/GSE286387_samples.csv")
    cfg <- jsonlite::fromJSON(cfg_path)
    frozen <- jsonlite::fromJSON(design_path)
    check(frozen$status == "accepted_for_step02_with_documented_metadata_limitations" &&
            frozen$accession == "GSE286387" && frozen$protocol_version == cfg$protocol_version,
          "Step 01 design is not accepted or version-matched.")
    check(frozen$formula == "~ condition" && cfg$differential_expression$model == frozen$formula &&
            frozen$reference_level == "control" && frozen$numerator == "oxaliplatin" &&
            frozen$denominator == "control" && !frozen$paired, "Unsupported frozen design.")
    de <- cfg$differential_expression
    check(de$engine == "DESeq2" && de$multiple_testing == "BH" &&
            identical(de$independent_filtering, FALSE) &&
            is.null(de$absolute_log2fc_cutoff) && is.null(de$baseMean_cutoff_after_DE) &&
            de$minimum_samples_rule == "smallest_group_size_across_two_groups" &&
            de$count_minimum == 10 && de$fdr_strictly_less_than == 0.05,
          "Unsupported thresholds; do not silently change the registered analysis.")
    check(!file.exists(file.path(input, "FAILURE.txt")) && !file.exists(file.path(qc, "FAILURE.txt")) &&
            file.exists(file.path(qc, "SUCCESS.txt")), "Input audit or local QC is incomplete/failed.")
    check(sha256(manifest_path) == frozen$sample_manifest_sha256, "Frozen sample manifest was modified.")
    counts_path <- file.path(input, "count_matrix.tsv")
    check(sha256(counts_path) == frozen$count_matrix_sha256, "Count matrix changed since design acceptance.")
    qc_hash <- read.csv(file.path(qc, "input_checksums.csv"), stringsAsFactors = FALSE)
    check(!anyDuplicated(qc_hash$input) &&
            identical(qc_hash$sha256[qc_hash$input == "count_matrix.tsv"], frozen$count_matrix_sha256),
          "QC did not use the accepted count matrix.")
    meta <- read.csv(manifest_path, stringsAsFactors = FALSE)
    original_meta <- read.csv(file.path(input, "sample_manifest.csv"), stringsAsFactors = FALSE)
    check(!anyDuplicated(meta$sample_id) && !anyDuplicated(meta$GEO_sample_id) &&
            all(meta$inclusion_status == "pass"), "Invalid or unaccepted samples.")
    check(setequal(meta$sample_id, original_meta$sample_id), "Original sample identities differ.")
    original_meta <- original_meta[match(meta$sample_id, original_meta$sample_id), ]
    check(identical(meta$GEO_sample_id, original_meta$GEO_sample_id) &&
            identical(meta$condition, original_meta$condition), "Original conditions or GSM identities differ.")
    check(identical(qc_hash$sha256[qc_hash$input == "sample_manifest.csv"],
                    sha256(file.path(input, "sample_manifest.csv"))), "QC sample metadata are stale.")
    tab <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
    check(names(tab)[1] == "stable_gene_id" && !anyDuplicated(tab$stable_gene_id) &&
            setequal(names(tab)[-1], meta$sample_id) && ncol(tab) == nrow(meta) + 1L,
          "Invalid count columns or duplicate stable IDs.")
    counts <- as.matrix(tab[, meta$sample_id, drop = FALSE])
    rownames(counts) <- tab$stable_gene_id
    check(is.numeric(counts) && all(is.finite(counts)) && all(counts >= 0) &&
            all(counts == trunc(counts)) && max(counts) <= .Machine$integer.max, "Invalid raw counts; no rounding allowed.")
    storage.mode(counts) <- "integer"
    ann_path <- file.path(input, "gene_annotations.csv")
    ann <- read.csv(ann_path, stringsAsFactors = FALSE, na.strings = "NA")
    check(!anyDuplicated(ann$stable_gene_id) && setequal(ann$stable_gene_id, rownames(counts)), "Annotation universe mismatch.")
    # Preserve the deposited biotype and all symbol candidates; verify against original export.
    raw_counts_path <- file.path(project, "data/raw/GSE286387/GSE286387_DRG_raw_counts.txt.gz")
    source <- read.csv(file.path(project, "docs/audits/source_manifest.csv"), stringsAsFactors = FALSE)
    check(identical(sha256(raw_counts_path), source$sha256[source$role == "counts"]), "Raw counts source changed.")
    raw <- read.delim(gzfile(raw_counts_path), check.names = FALSE, stringsAsFactors = FALSE)
    ann <- ann[match(rownames(counts), ann$stable_gene_id), ]
    check(all(ann$biotype[match(raw$ensembl_gene_id, ann$stable_gene_id)] == raw$gene_biotype),
          "Deposited biotype annotations were modified.")
    coldata <- data.frame(condition = factor(meta$condition, levels = c("control", "oxaliplatin")),
                          row.names = meta$sample_id)
    groups <- table(coldata$condition)
    check(!anyNA(coldata$condition) && groups["control"] == frozen$n_control &&
            groups["oxaliplatin"] == frozen$n_oxaliplatin &&
            min(groups) >= cfg$datasets$minimum_biological_units_per_group, "Invalid group sizes.")
    mm <- model.matrix(~ condition, coldata)
    check(qr(mm)$rank == ncol(mm), "Non-estimable design.")
    keep <- rowSums(counts >= de$count_minimum) >= min(groups)
    check(any(keep), "No genes pass the registered prefilter.")
    message("[2/5] Fitting DESeq2: oxaliplatin versus control, all retained biotypes...")
    dds <- DESeq2::DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], coldata, design = ~ condition)
    dds <- DESeq2::DESeq(dds, minReplicatesForReplace = 7, quiet = TRUE)
    result <- DESeq2::results(dds, contrast = c("condition", "oxaliplatin", "control"),
                             alpha = de$fdr_strictly_less_than, independentFiltering = FALSE, pAdjustMethod = "BH")
    message("[3/5] Auditing every gene and selecting significant protein-coding genes...")
    all <- data.frame(stable_gene_id = rownames(counts), symbol = ann$symbol_candidates,
                      biotype = ann$biotype, annotation_reference = "deposited_Ensembl_release_76",
                      prefilter_status = ifelse(keep, "pass", "excluded"),
                      baseMean = NA_real_, log2FC = NA_real_, lfcSE = NA_real_, statistic = NA_real_,
                      pvalue = NA_real_, padj = NA_real_, maxCooks = NA_real_,
                      na_reason = ifelse(keep, "", "count_prefilter_failed"))
    columns <- c(baseMean = "baseMean", log2FC = "log2FoldChange", lfcSE = "lfcSE", statistic = "stat",
                 pvalue = "pvalue", padj = "padj")
    for (name in names(columns)) all[keep, name] <- result[[columns[[name]]]]
    max_cooks <- S4Vectors::mcols(dds)$maxCooks
    if (!is.null(max_cooks)) all$maxCooks[keep] <- max_cooks
    cutoff <- qf(0.99, ncol(mm), nrow(mm) - ncol(mm))
    no_p <- keep & !is.finite(all$pvalue)
    all$na_reason[no_p] <- ifelse(is.finite(all$maxCooks[no_p]) & all$maxCooks[no_p] > cutoff,
                                  "cooks_outlier_default_Deseq2", "nonfinite_pvalue_other_or_unresolved")
    all$na_reason[keep & is.finite(all$pvalue) & !is.finite(all$padj)] <- "missing_adjusted_pvalue"
    all$valid_p_for_BH <- keep & is.finite(all$pvalue)
    all$significant <- keep & is.finite(all$padj) & all$padj < de$fdr_strictly_less_than
    all$passes_to_step03 <- all$significant & !is.na(all$biotype) & all$biotype == "protein_coding"
    all$direction <- ifelse(is.finite(all$log2FC), ifelse(all$log2FC > 0, "up", ifelse(all$log2FC < 0, "down", "zero")), NA_character_)
    all$protocol_version <- cfg$protocol_version
    valid <- all$valid_p_for_BH
    check(isTRUE(all.equal(all$padj[valid], p.adjust(all$pvalue[valid], method = "BH"), tolerance = 1e-12)),
          "BH universe/check failed.")
    discovery <- all[all$passes_to_step03, , drop = FALSE]
    discovery <- discovery[order(discovery$padj, discovery$stable_gene_id), , drop = FALSE]
    csv(all, "all_genes_de.csv")
    csv(all[keep, ], "tested_gene_universe.csv")
    csv(discovery, "significant_protein_coding_degs.csv")
    reasons <- ifelse(!keep, "COUNT_FILTER_FAILED", ifelse(!is.finite(all$padj), "DE_UNASSESSABLE",
                      ifelse(!all$significant, "NOT_SIGNIFICANT", ifelse(!all$passes_to_step03, "NONCODING", "PASS_TO_STEP03"))))
    csv(data.frame(step_id = "02", object_level = "gene", object_id = all$stable_gene_id,
                   status = ifelse(all$passes_to_step03, "pass", ifelse(keep & !is.finite(all$padj), "unassessable", "excluded")),
                   reason_code = reasons, reason_text = ifelse(nzchar(all$na_reason), all$na_reason, reasons),
                   input_source = "count_matrix.tsv", input_checksum = frozen$count_matrix_sha256,
                   protocol_version = cfg$protocol_version), "step_audit.csv")
    message("[4/5] Saving design, environment and fitted model...")
    saveRDS(dds, file.path(stage, "dds_fitted.rds"))
    csv(meta, "sample_manifest_used.csv")
    csv(data.frame(sample_id = meta$sample_id, mm, check.names = FALSE), "design_matrix.csv")
    inputs <- c(cfg_path, design_path, manifest_path, counts_path, ann_path)
    csv(data.frame(input = basename(inputs), sha256 = vapply(inputs, sha256, character(1))), "input_checksums.csv")
    writeLines(capture.output(sessionInfo()), file.path(stage, "sessionInfo.txt"))
    summary <- list(status = "STEP02_COMPLETE", accession = "GSE286387", protocol_version = cfg$protocol_version,
                    formula = frozen$formula, numerator = "oxaliplatin", denominator = "control",
                    input_genes = nrow(all), retained_genes = sum(keep), excluded_by_count = sum(!keep),
                    valid_pvalues_for_BH = sum(valid), nonfinite_pvalues_after_filter = sum(no_p),
                    significant_all_biotypes = sum(all$significant), protein_coding_DEGs = nrow(discovery),
                    protein_coding_up = sum(discovery$direction == "up"), protein_coding_down = sum(discovery$direction == "down"),
                    independentFiltering = FALSE, fold_change_cutoff = NULL, fdr = de$fdr_strictly_less_than,
                    DESeq2 = as.character(packageVersion("DESeq2")), R = as.character(getRversion()),
                    count_matrix_sha256 = frozen$count_matrix_sha256, metadata_limitations = frozen$remaining_limitations)
    jsonlite::write_json(summary, file.path(stage, "summary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    report <- c("# Step 02 — differential expression", "", "Status: STEP02_COMPLETE", "",
                paste0("- Design: ", frozen$formula, "; oxaliplatin / control; 5 RNA libraries per group; independent mice supported by public study design"),
                paste0("- Tested genes, all biotypes: ", sum(keep)), paste0("- Valid p-values for BH: ", sum(valid)),
                paste0("- Significant genes, all biotypes (BH FDR <0.05): ", sum(all$significant)),
                paste0("- Significant protein-coding genes: ", nrow(discovery)),
                paste0("- Protein-coding up / down: ", summary$protein_coding_up, " / ", summary$protein_coding_down),
                "- No fold-change cutoff; no post-DE baseMean threshold; independentFiltering=FALSE",
                paste0("- DESeq2 ", packageVersion("DESeq2"), "; R ", getRversion()), "",
                "All 10 samples retained; default DESeq2 Cook's handling kept. Missing results are not treated as nonsignificant evidence of resistance.",
                "Biotype is from the deposited historical annotation. Stable gene IDs remain distinct; symbols are not merged.",
                "Unreported batch and animal identifiers remain limitations; the reported 10-mouse study supports an unpaired design.",
                "DEG numbers need not equal the paper: our prefilter/BH universe differs and we impose no fold-change cutoff.",
                "Expression differences do not establish protein oxidation. Step 03 mapping follows this discovery list.")
    writeLines(report, file.path(stage, "de_report.md"))
    check(all(file.copy(list.files(stage, full.names = TRUE), output, overwrite = TRUE)), "Cannot publish DE outputs.")
    unlink(file.path(output, "FAILURE.txt"))
    writeLines("STEP02_COMPLETE", file.path(output, "SUCCESS.txt"))
    message("[5/5] Outputs saved: ", output)
    cat(paste(report[1:12], collapse = "\n"), "\n")
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(output, "FAILURE.txt"))
    stop("DE FAILED: ", conditionMessage(e), call. = FALSE)
  }, finally = unlink(stage, recursive = TRUE))
})
