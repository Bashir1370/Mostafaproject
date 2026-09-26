# Primary References and Reference Datasets

This file records the primary literature and public perturbation datasets currently anchoring the framework.

## Disulfidptosis / disulfide-stress foundation

1. **Liu X, Olszewski K, Zhang Y, et al.** Cystine transporter regulation of pentose phosphate pathway dependency and disulfide stress exposes a targetable metabolic vulnerability in cancer. *Nature Cell Biology*. 2020;22:476–486.  
   DOI: https://doi.org/10.1038/s41556-020-0496-x  
   Role in framework: SLC7A11–cystine loading, NADPH drain, PPP dependency, disulfide stress.

2. **Liu X, Nie L, Zhang Y, et al.** Actin cytoskeleton vulnerability to disulfide stress mediates disulfidptosis. *Nature Cell Biology*. 2023;25.  
   DOI: https://doi.org/10.1038/s41556-023-01091-2  
   Role in framework: definition of disulfidptosis, cytoskeletal disulfide bonding, F-actin collapse, WRC/Rac execution machinery.

3. **Tang M, et al.** Inhibition of thioredoxin reductase 1 sensitizes glucose-starved glioblastoma cells to disulfidptosis. *Cell Death & Differentiation*. 2025;32(4):598–612.  
   PubMed linked from GEO GSE282334.  
   Role in framework: TXNRD1 / thioredoxin defense and an empirical DPT transcriptomic reference.

## Ferroptosis foundation

4. **Doll S, Proneth B, Tyurina YY, et al.** ACSL4 dictates ferroptosis sensitivity by shaping cellular lipid composition. *Nature Chemical Biology*. 2017;13:91–98.  
   DOI: https://doi.org/10.1038/nchembio.2239  
   Role: ACSL4, PUFA-phospholipid susceptibility.

5. **Gao M, Monian P, Pan Q, et al.** Ferroptosis is an autophagic cell death process. *Cell Research*. 2016;26:1021–1032.  
   DOI: https://doi.org/10.1038/cr.2016.95  
   Role: NCOA4-mediated ferritinophagy and labile iron.

6. **Zou Y, Li H, Graham ET, et al.** Cytochrome P450 oxidoreductase contributes to phospholipid peroxidation in ferroptosis. *Nature Chemical Biology*. 2020;16:302–309.  
   DOI: https://doi.org/10.1038/s41589-020-0472-6  
   Role: POR and lipid-peroxidation execution.

7. **Bersuker K, Hendricks JM, Li Z, et al.** The CoQ oxidoreductase FSP1 acts parallel to GPX4 to inhibit ferroptosis. *Nature*. 2019;575:688–692.  
   DOI: https://doi.org/10.1038/s41586-019-1705-2  
   Role: FSP1/AIFM2 anti-ferroptotic defense.

8. **Doll S, Freitas FP, Shah R, et al.** FSP1 is a glutathione-independent ferroptosis suppressor. *Nature*. 2019;575:693–698.  
   DOI: https://doi.org/10.1038/s41586-019-1707-0

9. **Mao C, Liu X, Zhang Y, et al.** DHODH-mediated ferroptosis defence is a targetable vulnerability in cancer. *Nature*. 2021;593:586–590.  
   DOI: https://doi.org/10.1038/s41586-021-03539-7  
   Role: mitochondrial ferroptosis-defense axis.

10. **Soula M, Weber RA, Zilka O, et al.** Metabolic determinants of cancer cell sensitivity to canonical ferroptosis inducers. *Nature Chemical Biology*. 2020;16:1351–1360.  
    DOI: https://doi.org/10.1038/s41589-020-0613-y  
    Role: ferroptosis metabolic requirements and GCH1/BH4-related defense context.

## Verified public reference datasets

### DPT

- **GSE282334** — Homo sapiens; LN229-TAZ(4SA); EV/KO × regular/glucose-deprived; 8 samples; RNA-seq.  
  GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282334  
  Proposed use: primary DPT empirical calibration via interaction term.  
  Important anti-circularity note: TXNRD1 itself cannot use the KO contrast to validate its own RNA diagnostic value.

### FPT

- **GSE247883** — Homo sapiens A549; DMSO, RSL3, RSL3+Ferrostatin-1; 3 replicates/group.  
  GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247883  
  Proposed use: rescue-validated FPT response.

- **GSE255459** — Homo sapiens; BT549, HS578, SUM159; DMSO, Erastin, RSL3; duplicates.  
  GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255459  
  Proposed use: cross-inducer / cross-cell-line FPT consensus.

- **GSE131444** — Mus musculus MEF; DMSO vs Erastin; triplicates.  
  GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131444  
  Proposed use: independent cross-species FPT validation.

- **GSE317656** — Homo sapiens neuroblastoma; SK-N-AS and KELLY; control vs Erastin; triplicates.  
  GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317656  
  Proposed use: independent FPT validation in a neural-lineage context.

### Generic oxidative-stress comparator

- **GSE55169** — Homo sapiens fibroblasts; H2O2 time-course; nuclear/cytosolic RNA.  
  GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55169  
  Proposed use: same-direction generic-stress penalty.

## Reference policy

- Prefer primary mechanistic studies over review-derived gene lists.
- Public dataset metadata must be re-verified before analysis.
- Dataset inclusion does not imply equal evidence quality.
- The reference library is versioned and may expand, but all additions require documented inclusion rationale.
