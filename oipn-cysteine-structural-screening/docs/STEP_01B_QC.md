# Step 01b — sample QC in R

This is the remaining computational QC substep of Step 01, not a new scientific step. Input is the actual Step 01 raw-count audit output. It uses no old classifier dictionaries, ferroptosis signatures, symbol collapsing or existing ferroptosis helper functions.

## Entry checks

Requires DESeq2, jsonlite, Ubuntu sha256sum, the original downloaded raw files and Step 01 outputs. Root resolves from Rscript/source path, or explicit OIPN_PROJECT_DIR for test fixtures. Both project library and output paths resolve independently of the current working directory.

The script checks protocol/accession status, source SHA256 against the committed snapshot, sample IDs/conditions against the committed audited mapping, unique stable gene IDs, exact sample-column membership, count integer/range/zero-library validity and all exported count values against all deposited count rows. Changed sources require review, not silent acceptance. Missing/failed Python audits terminate QC.

## Computation

1. Preserve the complete 40,481-gene matrix as input; apply the declared count >=10 in at least five samples rule for the QC transform (18,538 genes). Write a gene prefilter audit. The unfiltered source is not overwritten.
2. Estimate DESeq2 size factors on retained genes and use varianceStabilizingTransformation with blind=TRUE and an intercept-only dataset. Initial fitType is parametric; DESeq2 can report/use automatic local trend fallback. The full transform is used rather than fast vst subsampling, avoiding an implicit nsub size requirement. Warnings remain visible to the user.
3. PCA on up to 500 highest-variance retained VST genes, centered but not feature-scaled. Ties resolve by stable gene ID. Sample distances (Euclidean) and sample correlations (Pearson) use all retained VST genes.
4. Record raw unique-gene library totals, detected features, features >=10 and size factors. Group labels are used to color plots, not to select PCA features or estimate QC dispersion.
5. Generate four labeled PNGs and a four-page PDF; save numerical coordinates/matrices, selected PCA genes, VST object, SHA256 of processed inputs/configuration and sessionInfo.

No Wald tests, likelihood-ratio tests, adjusted p-values or DEGs are computed. Normalized/VST data are for QC only and will never become DESeq2 count inputs.

## Output and review gate

`results/01b_sample_qc/`: library_size.png, pca_top500.png, sample_distance.png, sample_correlation.png, QC_plots.pdf, sample_metrics.csv, pca_coordinates.csv, pca_selected_genes.csv, sample_distances.csv, sample_correlations.csv, gene_prefilter_audit.csv, vst_qc.rds, input_checksums.csv, sessionInfo.txt, qc_report.md, step_audit.csv and SUCCESS.txt.

Status: QC_GENERATED_REVIEW_PENDING. Every sample remains held. Missing animal/pool identity, age and batch are not converted into invented labels. Sample retention/exclusion and final DE design require review; PCA distance alone is insufficient technical evidence for deletion. No criterion requires control and treatment to separate visibly in PCA.

Files are generated in staging and copied only on computational success. On rerun SUCCESS.txt is removed first. Failure writes FAILURE.txt; old outputs without SUCCESS.txt must not be interpreted. A successful run removes FAILURE.txt. A completed computational QC run does not certify biological independence or absence of batch confounding.

## Validation and architecture

All 83 current repository files were read and the independent root classifier, ferroptosis-only, CIPN ferroptosis and new structural screening architecture cross-checked. The new script imports only this project's count audit/configuration; other workspaces are unchanged.

Actual official inputs tested with R 4.3.3, DESeq2 1.42.0 and jsonlite 1.8.8 on Linux. Both Rscript and source() execution paths were exercised. Seven tests pass: four existing count-integrity tests plus real-data R QC, altered count rejection and altered condition rejection. Full run produces all ten held samples and 18,538 genes, with PC1 35.23% and PC2 24.73% in that tested environment. These numbers are computational observations, not evidence that any particular sample should be removed. User reports DESeq2 1.42.1; exact workstation reproduction remains pending.

```bash
python3 -m unittest discover -s oipn-cysteine-structural-screening/tests -v
```

The R integration tests skip explicitly if Rscript/DESeq2/jsonlite or real audit inputs are unavailable; they do not pretend to have executed. No unnecessary tests run scientific inference.

Method reference: [official DESeq2 vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), especially transformation, PCA and sample-distance sections.
