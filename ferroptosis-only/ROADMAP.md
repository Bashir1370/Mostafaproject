# Roadmap — Ferroptosis-Only Model

## Phase 1 — Systematic gene curation — IN PROGRESS
Build a master ferroptosis gene universe from:
- primary mechanistic literature
- NCCD / curated pathway resources
- high-quality genetic perturbation studies
- validated ferroptosis suppressors/promoters

For every gene record:
- mechanistic module
- promoter / suppressor / context
- evidence tier
- direct perturbation evidence
- rescue evidence
- post-translational dependence
- RNA interpretability
- planned scoring role

## Phase 2 — Freeze mechanistic architecture
Define:
- Core
- Extended
- Mechanistic annotation
- Excluded from direct RNA scoring

## Phase 3 — Reference dataset QC — IN PROGRESS
Completed primary-rescue execution for GSE182638. GSE247883 independent rescue replication is now locked and prepared from raw reads.

Prioritize:
- GSE182638
- GSE247883
- GSE255459
- GSE319384
- GSE287284
- GSE152988

## Phase 4 — Empirical calibration
For each study:
- differential expression
- effect shrinkage
- sign confidence
- direction-aware evidence
- within-study aggregation
- cross-study empirical reproducibility E

## Phase 5 — Build empirical ferroptosis signature
Prefer:
- induction + Ferrostatin-1 reversal
- cross-inducer agreement
- cross-cell-line reproducibility
- neural validation

## Phase 6 — Final scoring
Calculate:
- stage scores
- FPT-MPS
- FPT-MCI
- FPT-ESR

## Phase 7 — Validation
Use:
- held-out multi-inducer dataset
- neuronal ferroptosis dataset
- generic oxidative-stress negative controls

## Phase 8 — Freeze v1.0
Freeze:
- gene universe
- core dictionary
- empirical signature
- formulas
- thresholds
- software environment


## Current empirical-calibration checkpoint
- GSE182638 primary rescue analysis completed locally.
- GSE182638 exposed the need to separate mechanistic susceptibility direction from acute transcriptional response.
- GSE247883 is locked as an independent replication dataset; no discovery re-selection is permitted.
- Raw SRA/ENA -> Salmon -> tximport -> DESeq2/ashr workflow prepared for GSE247883.

## Phase 1 current checkpoint
- Master gene universe expanded and referenced.
- Primary mechanistic bibliography created.
- No gene weight is final yet.
- Provisional CORE genes must pass a dedicated scientific review before empirical calibration.
