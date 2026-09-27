# GSE287284 Direct Neuronal Validation Plan

## Role

GSE287284 is a **direct neuronal ferroptosis validation dataset**.

It is not used to discover a replacement signature and it is not treated as rescue evidence because there is no Ferrostatin-1 rescue arm.

## Verified experimental design

GEO reports:
- organism: Mus musculus
- model: N2a cells differentiated toward a neuronal phenotype with retinoic acid
- DMSO: n=3
- RSL3: n=3
- differentiation: 10 uM retinoic acid for 96 h
- ferroptosis induction: 5 uM RSL3 for 24 h
- platform: Illumina NovaSeq 6000
- processed mRNA files: per-sample gene-level FPKM
- genome assembly reported by GEO: GRCm38

## Scientific question

Do the directions discovered in the human GSE182638 Fer-1-rescue experiment transfer to an independent neuronal ferroptosis model?

The test is:

```text
GSE182638 frozen human ESR direction
            ↓
human-to-mouse one-to-one ortholog
            ↓
GSE287284 RSL3 vs DMSO
            ↓
same direction = neuronal support
opposite direction = neuronal contradiction
```

## No data leakage

The ESR candidate universe is frozen from GSE182638 before GSE287284 is inspected.

GSE287284 can:
- support a frozen candidate
- contradict a frozen candidate
- return no directional evidence
- return not testable/no one-to-one ortholog

It cannot add a new candidate to the primary ESR signature.

## Cross-species mapping

Human discovery genes are translated to mouse using Ensembl release 101 orthology.

Only **one-to-one human-mouse mappings** are used for the primary validation table.

This avoids treating simple capitalization equivalence as proof of orthology.

The ortholog table is cached at:

```text
ferroptosis-only/data/reference/human_mouse_orthologs_ensembl101.csv
```

## Processed-data analysis

GEO provides FPKM, not raw integer counts.

Therefore DESeq2 is not used.

Primary analysis:
1. download the six per-sample mRNA FPKM files
2. merge gene-level expression
3. filter genes expressed at FPKM >= 0.5 in at least 3 samples
4. transform as log2(FPKM + 1)
5. fit limma model
6. contrast RSL3 - DMSO
7. use limma posterior variance to obtain standard errors
8. apply ashr to obtain posterior effect and lfsr
9. score only frozen GSE182638 directions

## Why ashr is retained

Using limma coefficient standard errors followed by ashr keeps the validation evidence conceptually comparable to the discovery pipeline:

```text
effect size
+
uncertainty
↓
posterior beta + lfsr
↓
direction-aware evidence
```

The result is still labeled **neuronal induction evidence**, not rescue evidence.

## Mechanistic genes

ACSL4, LPCAT3, POR, GPX4 and the extended mechanistic universe are reported separately.

Their causal promoter/suppressor roles are not imposed as expected acute transcriptional directions.

This preserves the distinction learned from GSE182638:

```text
causal susceptibility role != acute transcriptional response
```

## Main outputs

- sample_metadata_locked.csv
- GSE287284_FPKM_matrix.csv
- expression_qc.csv
- logFPKM_pca_coordinates.csv
- logFPKM_pca.png
- logFPKM_sample_correlation.png
- RSL3_vs_DMSO_mouse_gene_effects.csv
- GSE182638_ESR_discovery_snapshot.csv
- GSE287284_frozen_ESR_neuronal_validation.csv
- mechanistic_genes_neuronal_RSL3_effects.csv
- analysis_summary.txt
- sessionInfo.txt

## Interpretation

A positive neuronal_induction_e means the neuronal RSL3 response agrees with the direction pre-specified by GSE182638.

A negative value means cross-species neuronal contradiction.

This dataset does not establish ferroptosis specificity by itself because it lacks a rescue arm. Its value is **neural transfer/generalization**.
