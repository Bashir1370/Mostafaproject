# Reference Dataset QC — NEC–FPT v0.2

This file records **pre-analysis data readiness**, not biological performance.

## Inclusion criteria for primary calibration

A study should ideally have:
1. interpretable induction of the target death mechanism;
2. biological replication;
3. a pathway-specific rescue/inhibitor or an orthogonal control;
4. count-level data or raw reads that can be reproducibly converted to counts;
5. metadata sufficient to reconstruct contrasts.

## Current readiness

| Accession | Program | Design strength | Public expression object | Readiness | Role |
|---|---|---|---|---|---|
| GSE108621 | NEC | TSZ + TNF-only + Nec-1s rescue | DESeq result tables; SRA | RAW_REPROCESS_REQUIRED | primary NEC calibration |
| GSE172027 | NEC | TSZ + Nec-1s rescue, human astrocytes | normalized log2; SRA | RAW_REPROCESS_REQUIRED | neural rescue replication |
| GSE154230 | NEC | two RIPK1-activating routes + Nec-1s, glia | normalized log; SRA | RAW_REPROCESS_REQUIRED | neuroimmune replication |
| GSE134234 | NEC | direct RIPK3 dimerization + TNF control | RPKM workbook; SRA | ORTHOGONAL_VALIDATION | direct-RIPK3 validation |
| GSE268650 | NEC | optogenetic RIPK3 + inhibitor/light controls | normalized log2 + DE tables; SRA | ORTHOGONAL_VALIDATION | human orthogonal validation |
| GSE287439 | NEC | TSZ +/- RIPK1 inhibitor in human iPSC tri-culture | 10x H5 | SCRNA_VALIDATION | cell-type validation |
| GSE182638 | FPT | RSL3 + Fer-1 rescue in 2 human lines | **raw read counts** | READY_COUNTS | primary FPT calibration |
| GSE247883 | FPT | RSL3 + Fer-1 rescue | FPKM; SRA | RAW_REPROCESS_REQUIRED | rescue replication |
| GSE255459 | FPT | Erastin + RSL3 in 3 human lines | **raw counts** | READY_COUNTS | cross-inducer calibration |
| GSE319384 | FPT | Erastin/RSL3/Ferroptocide, 2 genotypes | **raw counts** | READY_COUNTS | held-out multi-inducer validation |
| GSE287284 | FPT | neuronal RSL3 | FPKM; SRA | RAW_REPROCESS_REQUIRED | neural validation |
| GSE152988 | FPT | human iPSC-neuron genetic ferroptosis model | processed table; SRA | HUMAN_NEURAL_VALIDATION | genetic validation |
| GSE104664 | stress | H2O2 time course | **raw counts** | READY_COUNTS | oxidative-stress penalty |

## Priority order

### First executable count-level analyses
1. GSE182638
2. GSE255459
3. GSE319384 — keep held out until initial FPT weights/signature are defined
4. GSE104664

### NEC
GSE108621 has the best experimental design but requires raw-read reprocessing for the final DESeq2/ashr pipeline.

The already-deposited DESeq result tables can be used only for:
- metadata confirmation;
- contrast sanity checks;
- preliminary direction checks.

They should not silently replace final count-level inference when lfsr-based uncertainty is required.

## Outlier policy

Potential outliers are:
- reported;
- visualized by PCA/correlation;
- never silently removed.

Removal requires a documented technical reason and versioned decision.

## Cross-species policy

Mouse studies are analyzed in their native annotation first.
Human ortholog mapping occurs **after** per-study differential analysis.

Raw mouse/human matrices are never merged.
