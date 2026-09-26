# Reference Dataset QC — v0.1

This document records the initial data-readiness audit for the reference perturbation library.

## Inclusion criteria

A dataset is considered suitable for **primary empirical calibration** when it has:
1. a biologically interpretable perturbation/contrast;
2. adequate biological replication;
3. gene-level count data or raw reads that can be reproducibly converted to counts;
4. sample metadata sufficient to reconstruct the design matrix;
5. no unavoidable circularity for the gene being evaluated.

## Current audit

| Accession | Program | Replication | Public matrix | Final readiness | Planned contrast |
|---|---|---:|---|---|---|
| GSE282334 | DPT | 2/cell | normalized counts | Raw reprocessing required | (KO_Glu- − KO_Reg) − (EV_Glu- − EV_Reg) |
| GSE247883 | FPT | 3/group | FPKM | Raw reprocessing required | RSL3−DMSO plus Fer-1 reversal |
| GSE255459 | FPT | 2/condition/cell line | raw counts | Ready | Erastin−DMSO and RSL3−DMSO within each cell line |
| GSE131444 | FPT | 3/group | raw counts | Ready | Erastin−DMSO |
| GSE317656 | FPT | 3/group/cell line | gene counts | Ready | Erastin−Control within each cell line |
| GSE104664 | Generic stress | 3/group | raw counts | Ready | H2O2_16h−Control and H2O2_36h−Control |
| GSE55169 | Generic stress | no replicate per cell/fraction/time state | fractionated RNA profiles | Secondary only | Supporting annotation only |

## Critical decisions

### GSE282334
This is currently the strongest DPT transcriptomic reference in the seed library, but its supplementary matrix is normalized rather than raw counts. The final calibration should therefore reconstruct gene-level counts from SRA. The public normalized matrix may be used for exploratory PCA, sample-label checking and preliminary direction checks only.

### GSE247883
The rescue design is highly valuable mechanistically. However, the deposited supplementary matrix is FPKM. Final effect estimation and false-sign confidence should use reprocessed counts from SRA.

### GSE255459
This is immediately usable and particularly valuable because two canonical ferroptosis inducers are tested across three cell lines. Contrasts are first estimated within each cell line. Erastin and RSL3 evidence are then combined within this study before the study contributes one unit to cross-study E.

### GSE131444
Immediately usable count data. Because this dataset is mouse, ortholog mapping is performed only after per-study differential analysis. Raw mouse gene identifiers must not be mixed directly with human matrices.

### GSE317656
Immediately usable for validation after confirming processed columns are integer gene counts. The two cell lines are analyzed separately and aggregated within study.

### GSE104664
Selected as the primary generic oxidative-stress comparator because it has a clean untreated/H2O2 design, three biological replicates per group, and a public count matrix. The 16 h and 36 h contrasts are treated as two within-study contrasts and aggregated before specificity calculations.

### GSE55169
Retained only as supporting evidence. Its fractionated time-course design is scientifically interesting but inappropriate as the primary negative-control differential dataset because cell-line/fraction/time states lack biological replication.

## Next QC tasks

1. Import each READY_COUNTS matrix.
2. Confirm identifiers, integer/count nature, library sizes and zero inflation.
3. Reconstruct metadata from GEO sample labels.
4. PCA / sample-correlation QC.
5. Flag outliers before differential analysis using pre-specified rules.
6. Store standardized per-study count and metadata objects without merging studies.
