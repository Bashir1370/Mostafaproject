# Roadmap

## Milestone 1 — Reference-data ingestion
- download verified processed/raw reference datasets
- record platform, genome build, sample metadata and contrasts
- establish reproducible preprocessing per study
- never merge raw expression across studies

## Milestone 2 — Differential analysis
- estimate per-study effect sizes
- apply shrinkage where appropriate
- quantify sign confidence
- build standardized per-study effect tables

## Milestone 3 — Empirical reproducibility E
- compute direction-aware per-study evidence
- aggregate multiple contrasts within a study
- calculate positive support and contradiction penalty
- report independent-study depth separately

## Milestone 4 — Specificity S
- compare DPT evidence against same-direction FPT and generic-stress behavior
- mirror the procedure for FPT
- identify shared/non-discriminating genes

## Milestone 5 — Final gene weights
- W_MPS = M × E
- W_DSI = M × E × S
- document genes removed from RNA scoring despite high mechanistic relevance

## Milestone 6 — ESR signatures
- derive empirical DPT signature
- derive rescue-validated and cross-inducer FPT signatures
- build rank-based single-sample ESR scoring

## Milestone 7 — Calibration
- calculate MPS/MCI/ESR on all reference samples
- derive candidate thresholds
- bootstrap threshold stability
- leave-one-study-out validation
- retain continuous scores if thresholds are unstable

## Milestone 8 — Freeze v1.0
Before external application:
- freeze gene dictionaries
- freeze equations
- freeze weights
- freeze thresholds
- freeze software version
- archive validation results

## Success criterion

The framework should distinguish mechanism-specific DPT/FPT states while:
- not classifying generic oxidative stress as either state by default
- permitting Mixed and Indeterminate outcomes
- maintaining performance under held-out-study validation
- remaining interpretable at the stage and gene levels
