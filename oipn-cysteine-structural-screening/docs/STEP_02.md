# Step 02 — differential expression

## Entry and fixed method

Run after the original input audit and local sample QC succeed, using the accepted manifest and [design audit](DATASET_DESIGN_AUDIT.md). The script rejects modified frozen counts or sample manifest, stale QC provenance, noninteger counts and modified source biotypes. All ten libraries remain included.

Fit DESeq2 with `~ condition`; a positive log2FC means higher RNA expression after oxaliplatin. Prefilter: count >=10 in >=5 of all ten samples. Fit and apply BH across all retained biotypes with finite p-values, `independentFiltering=FALSE`, and default Cook's handling. No automatic count replacement occurs with five samples per group and the default replacement minimum of seven.

Entry to Step 03 requires **deposited protein-coding biotype and padj <0.05**. Both directions enter, with no absolute log2FC or post-analysis baseMean cutoff. Historical Ensembl annotations remain explicit; subsequent mapping must audit them against the protein reference. RNA expression does not demonstrate protein oxidation.

## Ubuntu execution

From the repository root, after pulling the updated files:

```bash
Rscript oipn-cysteine-structural-screening/scripts/02_DE_GSE286387.R
cat oipn-cysteine-structural-screening/results/02_differential_expression/de_report.md
```

Alternatively, in the RStudio Console:

```r
source("/home/bashir/Desktop/Mostafaproject/oipn-cysteine-structural-screening/scripts/02_DE_GSE286387.R")
```

The script loads the project's R library and requires DESeq2 and jsonlite. It does not install packages. It stages output and creates `SUCCESS.txt` only after publication finishes; `FAILURE.txt` or a missing success marker means outputs must not be used.

## Outputs and exit gate

`results/02_differential_expression/` contains:

| Output | Meaning |
|---|---|
| `all_genes_de.csv` | Every input stable gene ID, including prefilter exclusions and missing-result reasons |
| `tested_gene_universe.csv` | All count-filter-passing genes; retains noncoding and Cook's-unassessable genes |
| `significant_protein_coding_degs.csv` | Exactly the genes entering Step 03, sorted by padj and stable ID |
| `step_audit.csv` | One terminal status and reason for every input gene |
| `dds_fitted.rds` | Fitted DESeq2 object for reproducibility |
| `sample_manifest_used.csv`, `design_matrix.csv` | Actual sample inclusion and model matrix |
| `input_checksums.csv`, `sessionInfo.txt` | Input and environment provenance |
| `summary.json`, `de_report.md`, `SUCCESS.txt` | Machine/human summary and successful execution marker |

Review the summary, missing-result reasons, sample/model identity and exact discovery contract before mapping proteins. A missing p-value is unassessable, not evidence of biological resistance. No downstream structural filtering occurs in Step 02. DEG totals need not reproduce the manuscript: our registered count filter, BH universe and lack of fold-change cutoff differ.

## Implementation validation

A real-data implementation run with R 4.3.3 / DESeq2 1.42.0 tested 18,538 genes; 18,503 had valid p-values. It selected 893 significant genes across biotypes and 854 significant coding genes (477 up, 377 down). These are implementation-validation results; the user's DESeq2 1.42.1 run remains the authoritative workstation reproduction pending review. The small snapshot is under `docs/audits/DE_validation/`.

Tests independently recompute BH over all valid biotypes, verify complete gene auditing and exact significant-coding selection including small fold changes, and reject tampered frozen counts/sample manifests. Run after local Step 02:

```bash
python3 -m unittest discover -s oipn-cysteine-structural-screening/tests -v
```

Integration tests require Rscript with DESeq2/jsonlite; result-contract tests skip explicitly if Step 02 results are unavailable.
