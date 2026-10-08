#!/usr/bin/env Rscript
# Step 01b: exploratory sample QC, not differential-expression testing.
local({
  locate_project <- function() {
    explicit <- Sys.getenv("OIPN_PROJECT_DIR", unset = "")
    if (nzchar(explicit)) return(normalizePath(explicit, mustWork = TRUE))
    frames <- sys.frames()
    sourced <- unlist(lapply(frames, function(x) x$ofile), use.names = FALSE)
    if (length(sourced)) return(dirname(dirname(normalizePath(tail(sourced, 1), mustWork = TRUE))))
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(arg)) return(dirname(dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE))))
    stop("Run via source(full_script_path) or Rscript full_script_path.")
  }
  project <- locate_project()
  lib <- file.path(project, ".r-library")
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  needed <- c("DESeq2", "jsonlite")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages in project library: ", paste(missing, collapse = ", "))
  input <- file.path(project, "results", "01_dataset_audit")
  output <- file.path(project, "results", "01b_sample_qc")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(output, "SUCCESS.txt"))
  staging <- tempfile("qc_stage_", tmpdir = output)
  dir.create(staging)
  csv <- function(x, name, row.names = FALSE) {
    write.csv(x, file.path(staging, name), row.names = row.names, na = "NA")
  }
  check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  sha256 <- function(path) {
    check(nzchar(Sys.which("sha256sum")), "Ubuntu sha256sum command is required.")
    value <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
    check(is.null(attr(value, "status")), paste("Cannot checksum", path))
    strsplit(value[1], "[[:space:]]+")[[1]][1]
  }
  tryCatch({
    message("[1/5] Checking audited counts, metadata and source checksums...")
    check(!file.exists(file.path(input, "FAILURE.txt")), "Step 01 failed; rerun the Python audit successfully first.")
    config_path <- file.path(project, "config", "analysis_parameters.json")
    cfg <- jsonlite::fromJSON(config_path)
    check(identical(cfg$sample_qc$blind, TRUE) &&
            identical(cfg$sample_qc$initial_fit_type, "parametric") &&
            identical(cfg$sample_qc$pca_scale_features, FALSE) &&
            identical(cfg$sample_qc$automatic_sample_exclusion, FALSE),
          "Unsupported QC configuration; review before changing the registered method.")
    summary <- jsonlite::fromJSON(file.path(input, "audit_summary.json"))
    check(summary$status == "DATA_READY_QC_PENDING" && summary$accession == "GSE286387",
          "Unexpected input audit status/accession.")
    check(cfg$protocol_version == summary$protocol_version && cfg$datasets$primary_candidate == "GSE286387",
          "Configuration and input audit do not match.")
    source <- read.csv(file.path(input, "source_manifest.csv"), stringsAsFactors = FALSE)
    frozen_source <- read.csv(file.path(project, "docs", "audits", "source_manifest.csv"), stringsAsFactors = FALSE)
    check(!anyDuplicated(source$role) && setequal(source$role, c("soft", "counts")), "Unexpected source manifest.")
    for (i in seq_len(nrow(source))) {
      raw_path <- file.path(project, "data", "raw", "GSE286387", source$local_file[i])
      expected <- frozen_source$sha256[match(source$role[i], frozen_source$role)]
      check(length(expected) == 1L && !is.na(expected) && source$sha256[i] == expected,
            "Source changed since the committed audit; review before continuing.")
      check(sha256(raw_path) == expected, "Raw file checksum mismatch; rerun/review the audit.")
    }
    meta <- read.csv(file.path(input, "sample_manifest.csv"), stringsAsFactors = FALSE)
    frozen_meta <- read.csv(file.path(project, "docs", "audits", "sample_manifest.csv"), stringsAsFactors = FALSE)
    check(!anyDuplicated(meta$sample_id) && !anyDuplicated(meta$GEO_sample_id), "Duplicate sample identity.")
    check(setequal(meta$sample_id, frozen_meta$sample_id), "Unexpected samples.")
    frozen_meta <- frozen_meta[match(meta$sample_id, frozen_meta$sample_id), ]
    check(identical(meta$GEO_sample_id, frozen_meta$GEO_sample_id) && identical(meta$condition, frozen_meta$condition),
          "Sample identity or condition differs from audited mapping.")
    tab <- read.delim(file.path(input, "count_matrix.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
    check(names(tab)[1] == "stable_gene_id" && !anyDuplicated(tab$stable_gene_id), "Invalid/duplicate stable gene IDs.")
    check(setequal(names(tab)[-1], meta$sample_id) && ncol(tab) == nrow(meta) + 1L, "Count/sample columns mismatch.")
    counts <- as.matrix(tab[, meta$sample_id, drop = FALSE])
    rownames(counts) <- tab$stable_gene_id
    check(is.numeric(counts) && all(is.finite(counts)) && all(counts >= 0) && all(counts == trunc(counts)),
          "Counts must be finite nonnegative integers; no rounding is permitted.")
    check(max(counts) <= .Machine$integer.max && all(colSums(counts) > 0), "Invalid integer range or zero library.")
    # Verify every exported count against the deposited table, including annotation duplicate copies.
    raw <- read.delim(gzfile(file.path(project, "data", "raw", "GSE286387", source$local_file[source$role == "counts"])),
                      check.names = FALSE, stringsAsFactors = FALSE)
    check(setequal(unique(raw$ensembl_gene_id), rownames(counts)), "Exported gene universe differs from source.")
    check(all(counts[match(raw$ensembl_gene_id, rownames(counts)), , drop = FALSE] ==
                as.matrix(raw[, meta$sample_id, drop = FALSE])), "Processed count matrix differs from raw input.")
    check(nrow(counts) == summary$unique_genes && ncol(counts) == summary$n_samples, "Input summary dimensions mismatch.")
    group <- factor(meta$condition, levels = c("control", "oxaliplatin"))
    check(!anyNA(group) && all(table(group) == 5L), "Expected five samples per group.")
    check(cfg$differential_expression$minimum_samples_rule == "smallest_group_size_across_two_groups",
          "Unsupported prefilter rule; review implementation.")
    minimum_samples <- min(table(group))
    threshold <- cfg$differential_expression$count_minimum
    check(length(threshold) == 1L && is.finite(threshold) && threshold >= 1, "Invalid count minimum.")
    keep <- rowSums(counts >= threshold) >= minimum_samples
    check(sum(keep) >= 2L, "Too few genes pass the planned count filter.")
    csv(data.frame(stable_gene_id = rownames(counts), passes_planned_filter = keep), "gene_prefilter_audit.csv")
    storage.mode(counts) <- "integer"
    message("[2/5] Estimating size factors and blind variance-stabilizing transformation...")
    coldata <- data.frame(condition = group, row.names = meta$sample_id)
    # Intercept-only design is deliberate: group information does not determine QC transformation.
    dds <- DESeq2::DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], coldata, design = ~ 1)
    dds <- DESeq2::estimateSizeFactors(dds)
    vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
    vm <- SummarizedExperiment::assay(vsd)
    check(all(is.finite(vm)), "VST returned nonfinite values.")
    size <- DESeq2::sizeFactors(dds)
    metrics <- data.frame(sample_id = meta$sample_id, GEO_sample_id = meta$GEO_sample_id,
                          condition = meta$condition, library_size = colSums(counts),
                          detected_genes_count_gt0 = colSums(counts > 0),
                          genes_count_ge10 = colSums(counts >= 10), size_factor = unname(size),
                          inclusion_status = "held", reason = "manual_QC_and_metadata_review_pending")
    csv(metrics, "sample_metrics.csv")
    message("[3/5] Computing PCA, correlations and sample distances...")
    variance <- apply(vm, 1, var)
    variable <- which(is.finite(variance) & variance > 0)
    check(length(variable) >= 2L, "Insufficient variable genes for PCA.")
    ordered <- variable[order(-variance[variable], rownames(vm)[variable])]
    top_n <- cfg$sample_qc$pca_top_variable_genes
    check(length(top_n) == 1L && is.finite(top_n) && top_n >= 2L && top_n == trunc(top_n), "Invalid PCA gene limit.")
    selected <- head(ordered, as.integer(top_n))
    pc <- prcomp(t(vm[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
    pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    coords <- data.frame(sample_id = meta$sample_id, condition = meta$condition, pc$x[, 1:2, drop = FALSE])
    csv(coords, "pca_coordinates.csv")
    csv(data.frame(stable_gene_id = rownames(vm)[selected], VST_variance = variance[selected]), "pca_selected_genes.csv")
    distances <- as.matrix(dist(t(vm)))
    correlations <- cor(vm, method = "pearson")
    csv(distances, "sample_distances.csv", TRUE)
    csv(correlations, "sample_correlations.csv", TRUE)
    saveRDS(vsd, file.path(staging, "vst_qc.rds"))
    message("[4/5] Saving labeled plots and reproducibility records...")
    colors <- c(control = "#2673B8", oxaliplatin = "#C34A36")
    labels <- meta$sample_id
    plots <- list(
      library_size = function() {
        par(mar = c(10, 5, 4, 1))
        barplot(metrics$library_size / 1e6, names.arg = labels, las = 2,
                col = colors[meta$condition], ylab = "Deposited gene counts (millions)",
                main = "GSE286387: library totals (unique genes)")
        legend("topright", names(colors), fill = colors, bty = "n")
      },
      pca_top500 = function() {
        par(mar = c(5, 5, 4, 2))
        xr <- range(coords$PC1); yr <- range(coords$PC2)
        padding <- function(r) r + c(-1, 1) * max(diff(r), 1) * 0.2
        plot(coords$PC1, coords$PC2, pch = 19, col = colors[meta$condition],
             xlim = padding(xr), ylim = padding(yr),
             xlab = sprintf("PC1 (%.1f%%)", pct[1]), ylab = sprintf("PC2 (%.1f%%)", pct[2]),
             main = paste0("Blind VST PCA: top ", length(selected), " variable genes"))
        text(coords$PC1, coords$PC2, labels, pos = ifelse(grepl("_rep1$", labels), 1, 3), cex = 0.65)
        legend("topright", names(colors), col = colors, pch = 19, bty = "n")
      },
      sample_distance = function() {
        par(oma = c(2, 0, 3, 0))
        heatmap(distances, Rowv = as.dendrogram(hclust(as.dist(distances))), Colv = "Rowv",
                scale = "none", symm = TRUE, margins = c(11, 11), cexRow = 0.8, cexCol = 0.8,
                col = hcl.colors(100, "YlOrRd", rev = TRUE))
        mtext("VST sample distances (Euclidean)", side = 3, outer = TRUE, line = 1)
        mtext("Pale = closer; dark = more distant", side = 1, outer = TRUE)
      },
      sample_correlation = function() {
        par(oma = c(2, 0, 3, 0))
        heatmap(correlations, Rowv = as.dendrogram(hclust(as.dist(distances))), Colv = "Rowv",
                scale = "none", symm = TRUE, margins = c(11, 11), cexRow = 0.8, cexCol = 0.8,
                col = hcl.colors(100, "YlOrRd", rev = TRUE))
        mtext("VST sample correlations (Pearson)", side = 3, outer = TRUE, line = 1)
        mtext("Pale = lower correlation; dark = higher correlation", side = 1, outer = TRUE)
      }
    )
    for (name in names(plots)) {
      png(file.path(staging, paste0(name, ".png")), width = 2000, height = 1700, res = 180)
      tryCatch(plots[[name]](), finally = dev.off())
    }
    pdf(file.path(staging, "QC_plots.pdf"), width = 11, height = 9)
    tryCatch(for (draw in plots) draw(), finally = dev.off())
    files <- c(config_path, file.path(input, c("count_matrix.tsv", "sample_manifest.csv", "audit_summary.json")))
    csv(data.frame(input = basename(files), sha256 = vapply(files, sha256, character(1))), "input_checksums.csv")
    writeLines(capture.output(sessionInfo()), file.path(staging, "sessionInfo.txt"))
    report <- c("# Step 01b — sample QC", "", "Status: QC_GENERATED_REVIEW_PENDING", "",
                paste0("- Samples: ", ncol(counts), " (5 control, 5 oxaliplatin)"),
                paste0("- Unique input genes: ", nrow(counts)),
                paste0("- Genes for QC: ", sum(keep), " (count >= ", threshold, " in >= ", minimum_samples, " samples)"),
                paste0("- DESeq2: ", packageVersion("DESeq2"), "; R: ", getRversion()),
                sprintf("- PCA variance: PC1 %.2f%%; PC2 %.2f%%", pct[1], pct[2]),
                paste0("- Library size range: ", min(metrics$library_size), " to ", max(metrics$library_size)),
                "", "Transformation: blind VST, initial parametric fit; DESeq2 may automatically use a local dispersion trend if needed.",
                "PCA uses the top 500 variable filtered genes (or fewer if unavailable); distances and correlations use all filtered genes.",
                "No hypothesis tests, DEGs, sample exclusions or final inclusion approvals were generated.",
                "", "## Review before Step 02", "",
                "Inspect library totals, PCA, within/between-group distances and correlations together.",
                "A distant point or lack of group separation alone does not justify sample removal.",
                "Animal/pool IDs, age, batch and exact final-dose-to-collection interval remain unresolved in GEO metadata.",
                "Confirm independent biological units and review full study methods before freezing the design.",
                "Count-level QC cannot assess mapping rates, RNA degradation or read quality without additional evidence.",
                "", "Successful execution requires SUCCESS.txt and absence of FAILURE.txt in this output directory.")
    writeLines(report, file.path(staging, "qc_report.md"))
    csv(data.frame(step_id = "01", object_level = "sample", object_id = meta$GEO_sample_id,
                   status = "held", reason_code = "QC_REVIEW_PENDING", reason_text = "Plots generated; metadata/inclusion review pending",
                   input_source = "sample_manifest.csv", input_checksum = sha256(file.path(input, "sample_manifest.csv")),
                   protocol_version = cfg$protocol_version), "step_audit.csv")
    check(all(file.copy(list.files(staging, full.names = TRUE), output, overwrite = TRUE)), "Failed to publish QC outputs.")
    unlink(file.path(output, "FAILURE.txt"))
    writeLines("QC_GENERATED_REVIEW_PENDING", file.path(output, "SUCCESS.txt"))
    message("[5/5] QC saved to: ", output)
    cat(paste(report[1:12], collapse = "\n"), "\n")
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(output, "FAILURE.txt"))
    stop("QC FAILED: ", conditionMessage(e), call. = FALSE)
  }, finally = unlink(staging, recursive = TRUE))
})
