# GSE286387 analysis unit

**Dataset identity:** Mouse DRG / chronic systemic Oxaliplatin vs Vehicle / GEO GSE286387.

Run only after the shared mouse gene-set audit has been frozen:

```r
source("cipn-ferroptosis-screening/scripts/00_prepare_mouse_gene_sets.R") # once only
source("cipn-ferroptosis-screening/scripts/03_GSE286387_mouse_screen.R")
```

Outputs are written exclusively to `cipn-ferroptosis-screening/results/GSE286387/`.
The loader retrieves GEO metadata/supplementary files, requires a genome-scale integer raw-count matrix, audits the expected 5 Vehicle + 5 Oxaliplatin RNA-seq samples, and refuses to substitute normalized or DEG-only tables for raw counts. Runtime sample mapping is exported for inspection.
