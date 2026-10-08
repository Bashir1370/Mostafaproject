# Step 01b — sample QC

Status: QC_GENERATED_REVIEW_PENDING

- Samples: 10 (5 control, 5 oxaliplatin)
- Unique input genes: 40481
- Genes for QC: 18538 (count >= 10 in >= 5 samples)
- DESeq2: 1.42.0; R: 4.3.3
- PCA variance: PC1 35.23%; PC2 24.73%
- Library size range: 52963470 to 60984426

Transformation: blind VST, initial parametric fit; DESeq2 may automatically use a local dispersion trend if needed.
PCA uses the top 500 variable filtered genes (or fewer if unavailable); distances and correlations use all filtered genes.
No hypothesis tests, DEGs, sample exclusions or final inclusion approvals were generated.

## Review before Step 02

Inspect library totals, PCA, within/between-group distances and correlations together.
A distant point or lack of group separation alone does not justify sample removal.
Animal/pool IDs, age, batch and exact final-dose-to-collection interval remain unresolved in GEO metadata.
Confirm independent biological units and review full study methods before freezing the design.
Count-level QC cannot assess mapping rates, RNA degradation or read quality without additional evidence.

Successful execution requires SUCCESS.txt and absence of FAILURE.txt in this output directory.
