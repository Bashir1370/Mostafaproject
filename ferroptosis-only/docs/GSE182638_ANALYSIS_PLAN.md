# GSE182638 Primary Ferroptosis Calibration Plan

## Why this dataset is first

GSE182638 contains 18 human RNA-seq samples from two multiple-myeloma cell lines:
- MM1R
- MM1S

Each cell line contains:
- untreated control
- RSL3
- Ferrostatin-1 pre-treatment + RSL3

with three biological replicates per condition.

The deposited GEO supplementary data include a raw read-count matrix, so the final DESeq2 pipeline can start directly from public counts.

## Locked contrasts

Within each cell line:

1. **Induction**
   ```text
   RSL3 - Untreated
   ```

2. **Rescue reversal**
   ```text
   RSL3 - (RSL3 + Ferrostatin-1)
   ```

Because both contrasts point toward the RSL3 state, a genuine rescue-consistent transcript should show the same effect direction in both contrasts.

## Blocking factor

Biological replicates correspond to matched cell-passage blocks within each cell line.

Therefore the differential-expression model is:

```text
~ replicate + condition
```

rather than a simple unblocked condition-only model.

## Mechanistic evidence

For each literature-curated RNA-eligible gene:

```text
e_induction
e_reversal
```

are combined using:

```text
support =
  min(max(e_induction,0),
      max(e_reversal,0))

contradiction =
  max(max(-e_induction,0),
      max(-e_reversal,0))

e_rescue = support - contradiction
```

The two cell lines are then aggregated within this study so GSE182638 contributes one independent study-level evidence value per gene.

## Empirical ESR discovery

A separate all-gene analysis is retained for later empirical-signature construction.

A gene is rescue-concordant within a cell line when:
- RSL3 changes it relative to untreated;
- RSL3 differs from the Fer-1 rescued state in the same direction.

No final ESR gene set will be selected from GSE182638 alone.

Cross-study replication with independent ferroptosis datasets remains mandatory.

## Primary outputs

- sample_metadata_geo.csv
- library_qc.csv
- vst_pca_coordinates.csv
- vst_pca.png
- vst_sample_correlation.png
- gene_identifier_mapping.csv
- mechanistic_gene_presence.csv
- per_contrast_effects_all_genes.csv
- per_contrast_mechanistic_gene_evidence.csv
- per_cell_line_rescue_validated_mechanistic_evidence.csv
- GSE182638_mechanistic_study_evidence.csv
- per_cell_line_all_gene_rescue_concordance.csv
- GSE182638_empirical_ESR_candidates.csv
- analysis_summary.txt
- sessionInfo.txt

## Interpretation rule

This study calibrates **transcriptomic behavior under experimentally rescued ferroptosis**.

It does not prove that a transcript alone is a ferroptosis execution marker.
