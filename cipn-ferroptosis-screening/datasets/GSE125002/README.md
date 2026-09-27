# GSE125002 analysis unit

**Dataset identity:** Mouse DRG / Oxaliplatin vs Vehicle / GEO GSE125002.

Run only after the shared mouse gene-set audit has been frozen:

```r
source("cipn-ferroptosis-screening/scripts/00_prepare_mouse_gene_sets.R") # once only
source("cipn-ferroptosis-screening/scripts/02_GSE125002_mouse_screen.R")
```

Outputs are written exclusively to `cipn-ferroptosis-screening/results/GSE125002/`.
The script downloads `GSE125002_RAW.tar`, uses only GSM3560840-845, excludes cisplatin, reconstructs gene-level mouse-symbol counts, runs DESeq2/QC, the five audited signatures, sample-level scores, and the mechanistic panel.
