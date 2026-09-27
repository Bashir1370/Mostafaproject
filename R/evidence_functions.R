# Core evidence functions for the NEC–FPT classifier
# Version: v0.2.0
#
# Statistical principle:
#   raw counts -> DESeq2 effect + SE -> ashr posterior effect + lfsr
#   -> direction-aware per-contrast evidence e in [-1,1]
#
# No hard p-value threshold is used.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ashr)
  library(dplyr)
  library(tibble)
})

prefilter_counts <- function(count_mat, group, min_count = 10L) {
  stopifnot(ncol(count_mat) == length(group))
  smallest_group <- min(table(group))
  keep <- rowSums(count_mat >= min_count) >= smallest_group
  count_mat[keep, , drop = FALSE]
}

fit_deseq2 <- function(count_mat, coldata, design_formula) {
  stopifnot(all(colnames(count_mat) == rownames(coldata)))
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_mat)),
    colData = coldata,
    design = design_formula
  )
  dds <- DESeq(dds, quiet = TRUE)
  dds
}

ashr_effect_from_coef <- function(dds, coef_name) {
  res_mle <- results(dds, name = coef_name, independentFiltering = FALSE)

  ok <- is.finite(res_mle$log2FoldChange) &
        is.finite(res_mle$lfcSE) &
        res_mle$lfcSE > 0

  out <- tibble(
    gene = rownames(res_mle),
    beta_mle = res_mle$log2FoldChange,
    se_mle = res_mle$lfcSE,
    beta = NA_real_,
    posterior_sd = NA_real_,
    lfsr = NA_real_
  )

  fit <- ashr::ash(
    betahat = res_mle$log2FoldChange[ok],
    sebetahat = res_mle$lfcSE[ok],
    mixcompdist = "normal",
    method = "shrink"
  )

  out$beta[ok] <- fit$result$PosteriorMean
  out$posterior_sd[ok] <- fit$result$PosteriorSD
  out$lfsr[ok] <- fit$result$lfsr
  out
}

ashr_effect_from_contrast <- function(dds, contrast) {
  res_mle <- results(
    dds,
    contrast = contrast,
    independentFiltering = FALSE
  )

  ok <- is.finite(res_mle$log2FoldChange) &
        is.finite(res_mle$lfcSE) &
        res_mle$lfcSE > 0

  out <- tibble(
    gene = rownames(res_mle),
    beta_mle = res_mle$log2FoldChange,
    se_mle = res_mle$lfcSE,
    beta = NA_real_,
    posterior_sd = NA_real_,
    lfsr = NA_real_
  )

  fit <- ashr::ash(
    betahat = res_mle$log2FoldChange[ok],
    sebetahat = res_mle$lfcSE[ok],
    mixcompdist = "normal",
    method = "shrink"
  )

  out$beta[ok] <- fit$result$PosteriorMean
  out$posterior_sd[ok] <- fit$result$PosteriorSD
  out$lfsr[ok] <- fit$result$lfsr
  out
}

add_effect_magnitude <- function(effect_tbl) {
  effect_tbl %>%
    mutate(
      magnitude = ifelse(
        is.finite(beta),
        rank(abs(beta), ties.method = "average", na.last = "keep") /
          sum(is.finite(beta)),
        NA_real_
      ),
      confidence = ifelse(is.finite(lfsr), 1 - lfsr, NA_real_)
    )
}

directional_evidence <- function(effect_tbl, direction) {
  stopifnot(direction %in% c(-1, 1))

  add_effect_magnitude(effect_tbl) %>%
    mutate(
      signed_alignment = sign(direction * beta),
      e = signed_alignment * magnitude * confidence
    )
}

aggregate_signed_e <- function(x, method = "median") {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)

  if (method == "median") return(median(x))
  if (method == "mean") return(mean(x))
  stop("Unsupported aggregation method")
}

aggregate_empirical_E <- function(study_e) {
  # study_e must contain one value per independent study.
  x <- study_e[is.finite(study_e)]
  if (length(x) == 0) return(NA_real_)

  positive <- mean(pmax(x, 0))
  contradiction <- mean(pmax(-x, 0))

  tibble(
    E = max(0, positive - contradiction),
    positive_support = positive,
    contradiction = contradiction,
    independent_studies = length(x)
  )
}

rescue_validated_evidence <- function(e_induction, e_reversal) {
  # With contrast definitions:
  #   induction:  inducer - control
  #   reversal:   inducer - (inducer + rescue)
  # both should align to the target-state direction when rescue works.
  support <- pmin(pmax(e_induction, 0), pmax(e_reversal, 0))
  contradiction <- pmax(pmax(-e_induction, 0), pmax(-e_reversal, 0))
  support - contradiction
}

same_direction_mimicry <- function(effect_tbl, target_direction) {
  # IMPORTANT:
  # competitor/stress evidence is evaluated using the TARGET program's
  # expected direction, not the competitor program's own dictionary sign.
  directional_evidence(effect_tbl, target_direction)
}

compute_specificity_S <- function(E_target, competitor_same_direction, stress_same_direction,
                                  epsilon = 1e-8) {
  comps <- c(competitor_same_direction, stress_same_direction)
  comps <- comps[is.finite(comps)]
  C <- if (length(comps) == 0) 0 else max(comps, 0)

  if (!is.finite(E_target)) return(NA_real_)

  S <- max(0, (E_target - C) / (E_target + C + epsilon))

  tibble(
    S = S,
    competitor_mimicry_C = C
  )
}


score_dictionary_genes <- function(effect_tbl, dictionary_tbl) {
  # dictionary_tbl requires: gene, direction
  effect_tbl %>%
    inner_join(
      dictionary_tbl %>% dplyr::select(gene, direction),
      by = "gene"
    ) %>%
    add_effect_magnitude() %>%
    mutate(
      signed_alignment = sign(direction * beta),
      e = signed_alignment * magnitude * confidence
    )
}


three_way_rescue_validated_evidence <- function(e_induction, e_specific_vs_context, e_reversal) {
  # Example NEC GSE108621:
  # induction: TSZ - DMSO
  # specific_vs_context: TSZ - TNF
  # reversal: TSZ - (TSZ + Nec-1s)
  support <- pmin(
    pmax(e_induction, 0),
    pmax(e_specific_vs_context, 0),
    pmax(e_reversal, 0)
  )
  contradiction <- pmax(
    pmax(-e_induction, 0),
    pmax(-e_specific_vs_context, 0),
    pmax(-e_reversal, 0)
  )
  support - contradiction
}
