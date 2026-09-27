# Roadmap — NEC–FPT Framework

## Milestone 0 — Redesign completed
- replace DPT with NEC
- create necroptosis gene dictionary
- redefine specificity index as NFSI
- rebuild reference dataset library
- rebuild contrast registry
- retain FPT architecture with improved rescue/orthogonal references

## Milestone 1 — Reference-data inspection
- run `scripts/00_inspect_reference_datasets.R`
- record public file types and dimensions
- verify sample labels and replicate structure
- classify studies as READY_COUNTS / RAW_REPROCESS_REQUIRED / validation-only
- never merge raw expression across studies

## Milestone 2 — NEC primary calibration
Start with:
- GSE108621

Required evidence:
- TSZ vs DMSO
- TSZ vs TNF
- TSZ vs TSZ+Nec-1s

Then replicate with:
- GSE172027
- GSE154230

## Milestone 3 — FPT primary calibration
Start with:
- GSE182638

Then:
- GSE247883
- GSE255459

Hold:
- GSE319384 for orthogonal multi-inducer validation

## Milestone 4 — Empirical reproducibility E
- DESeq2 effect estimation
- ashr shrinkage / lfsr
- direction-aware evidence e
- hierarchical within-study aggregation
- cross-study positive support and contradiction penalty

## Milestone 5 — Specificity S
- evaluate NEC genes in FPT/stress references using NEC direction
- evaluate FPT genes in NEC/stress references using FPT direction
- calculate same-direction mimicry penalty
- expand negative comparator panel if needed

## Milestone 6 — ESR
- derive rescue-validated NEC empirical signature
- derive rescue/cross-inducer FPT empirical signature
- test signature stability across held-out studies
- build rank-based single-sample ESR

## Milestone 7 — Neural robustness
NEC:
- GSE287439 cell-type pseudobulk/scRNA validation

FPT:
- GSE287284
- GSE152988

## Milestone 8 — Score calibration
- NEC-MPS / NEC-MCI / NEC-ESR
- FPT-MPS / FPT-MCI / FPT-ESR
- NFSI
- bootstrap threshold stability
- leave-one-study-out validation
- report continuous evidence if thresholds are unstable

## Milestone 9 — Freeze v1.0
Freeze:
- dictionaries
- evidence tiers
- equations
- reference studies
- empirical weights
- ESR signatures
- thresholds
- software environment

## Success criteria

The framework should:
- distinguish NEC from FPT in held-out perturbation studies;
- not reduce to generic TNF inflammation or oxidative stress;
- preserve performance across different induction mechanisms;
- support Mixed and Indeterminate states;
- remain interpretable at gene, stage and study levels.
