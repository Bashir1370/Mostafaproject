# Primary References and Reference Datasets

This file records the primary literature and public perturbation datasets currently anchoring the framework.

## Disulfidptosis / disulfide-stress foundation

1. **Liu X, Olszewski K, Zhang Y, et al.** Cystine transporter regulation of pentose phosphate pathway dependency and disulfide stress exposes a targetable metabolic vulnerability in cancer. *Nature Cell Biology*. 2020;22:476–486.  
   DOI: https://doi.org/10.1038/s41556-020-0496-x  
   Role in framework: SLC7A11–cystine loading, NADPH drain, PPP dependency, disulfide stress.

2. **Liu X, Nie L, Zhang Y, et al.** Actin cytoskeleton vulnerability to disulfide stress mediates disulfidptosis. *Nature Cell Biology*. 2023;25:404–414.  
   DOI: https://doi.org/10.1038/s41556-023-01091-2  
   Role in framework: definition of disulfidptosis, cytoskeletal disulfide bonding, F-actin collapse, WRC/Rac execution machinery.

3. **Tang M, et al.** Inhibition of thioredoxin reductase 1 sensitizes glucose-starved glioblastoma cells to disulfidptosis. *Cell Death & Differentiation*. 2025;32(4):598–612.  
   DOI: https://doi.org/10.1038/s41418-024-01440-0  
   Role in framework: TXNRD1 / thioredoxin defense and an empirical DPT transcriptomic reference.

## Ferroptosis foundation

4. **Doll S, Proneth B, Tyurina YY, et al.** ACSL4 dictates ferroptosis sensitivity by shaping cellular lipid composition. *Nature Chemical Biology*. 2017;13:91–98.  
   DOI: https://doi.org/10.1038/nchembio.2239

5. **Gao M, Monian P, Pan Q, et al.** Ferroptosis is an autophagic cell death process. *Cell Research*. 2016;26:1021–1032.  
   DOI: https://doi.org/10.1038/cr.2016.95

6. **Zou Y, Li H, Graham ET, et al.** Cytochrome P450 oxidoreductase contributes to phospholipid peroxidation in ferroptosis. *Nature Chemical Biology*. 2020;16:302–309.  
   DOI: https://doi.org/10.1038/s41589-020-0472-6

7. **Bersuker K, Hendricks JM, Li Z, et al.** The CoQ oxidoreductase FSP1 acts parallel to GPX4 to inhibit ferroptosis. *Nature*. 2019;575:688–692.  
   DOI: https://doi.org/10.1038/s41586-019-1705-2

8. **Doll S, Freitas FP, Shah R, et al.** FSP1 is a glutathione-independent ferroptosis suppressor. *Nature*. 2019;575:693–698.  
   DOI: https://doi.org/10.1038/s41586-019-1707-0

9. **Mao C, Liu X, Zhang Y, et al.** DHODH-mediated ferroptosis defence is a targetable vulnerability in cancer. *Nature*. 2021;593:586–590.  
   DOI: https://doi.org/10.1038/s41586-021-03539-7

10. **Soula M, Weber RA, Zilka O, et al.** Metabolic determinants of cancer cell sensitivity to canonical ferroptosis inducers. *Nature Chemical Biology*. 2020;16:1351–1360.  
    DOI: https://doi.org/10.1038/s41589-020-0613-y

## Verified public reference datasets

### DPT

- **GSE282334** — Homo sapiens; LN229-TAZ(4SA); EV/KO × regular/glucose-deprived; 8 samples.  
  Primary use: DPT interaction calibration.  
  Final-analysis note: GEO provides normalized counts as the supplementary matrix; raw SRA reads are available. Final effect estimation should therefore be performed from reprocessed count-level data.  
  Anti-circularity: TXNRD1 cannot use its own knockout contrast as empirical validation of TXNRD1 RNA behavior.

### FPT

- **GSE247883** — Homo sapiens A549; DMSO, RSL3, RSL3+Ferrostatin-1; 3 replicates/group.  
  Primary use: rescue-validated FPT response.  
  Final-analysis note: supplementary expression is FPKM; raw SRA should be reprocessed to gene counts for final calibration.

- **GSE255459** — Homo sapiens; BT549, HS578, SUM159; DMSO, Erastin, RSL3; duplicates.  
  Raw count matrices are provided.  
  Primary use: cross-inducer / cross-cell-line FPT consensus.

- **GSE131444** — Mus musculus MEF; DMSO vs Erastin; triplicates.  
  Raw count matrix is provided.  
  Primary use: independent cross-species FPT validation; human ortholog mapping is required before human-centered aggregation.

- **GSE317656** — Homo sapiens; SK-N-AS and KELLY neuroblastoma; control vs Erastin; triplicates.  
  Processed gene-count CSVs are provided.  
  Primary use: independent FPT validation in a neural-lineage context.

### Generic oxidative-stress comparators

- **GSE104664 — PRIMARY** — Homo sapiens HUVEC; untreated vs 200 µM H2O2 for 16 h and 36 h; triplicates.  
  Raw gene-count matrix is provided.  
  Primary use: same-direction generic oxidative-stress penalty used in specificity S.

- **GSE55169 — SECONDARY ONLY** — Homo sapiens MRC5/BJ fibroblasts; H2O2 time course; nuclear/cytosolic fractions.  
  Because individual cell-line/fraction/time comparisons do not contain adequate biological replication, this dataset is retained as supporting annotation rather than a primary specificity-calibration dataset.

## Reference policy

- Prefer primary mechanistic studies over review-derived gene lists.
- Public dataset metadata must be re-verified before analysis.
- Dataset inclusion does not imply equal evidence quality.
- Raw-count or reprocessed count-level data are preferred for effect estimation.
- Normalized expression may be used for exploratory checks but cannot silently substitute for final count-level calibration.
- The reference library is versioned; all additions/removals require a documented rationale.
