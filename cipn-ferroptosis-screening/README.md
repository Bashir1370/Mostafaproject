# CIPN Ferroptosis Screening

A deliberately simple, transcriptome-based workflow to answer one practical question:

> Does dorsal-root-ganglion (DRG) RNA-seq from chemotherapy-induced peripheral neuropathy (CIPN) show transcriptomic evidence consistent with ferroptosis?

This workspace is intentionally separate from the more complex ferroptosis-modeling experiments in `ferroptosis-only/`.

## First target dataset

**GSE160543** (rat DRG, bulk RNA-seq):
- Vehicle: n = 4
- Paclitaxel: n = 4
- Oxaliplatin: n = 4
- processed counts are available from GEO in a small supplementary archive

Primary contrasts:
1. Oxaliplatin vs Vehicle
2. Paclitaxel vs Vehicle

## Evidence strategy

We do **not** infer ferroptosis from one mechanistic gene such as GPX4, ACSL4, or SLC7A11.

Instead we ask whether several independent transcriptomic resources converge.

### Primary directional state signature
**Vinik 2024 — 24 validated ferroptosis-vs-apoptosis biomarkers**

This set was derived from transcriptomic ferroptosis/apoptosis datasets and experimentally validated with multiple ferroptosis inducers, apoptosis inducers, qRT-PCR, and in vivo experiments.

For the rat DRG analysis, human biomarkers are mapped to rat orthologs.

### Supporting ferroptosis pathways
- KEGG Ferroptosis
- WikiPathways Ferroptosis (WP_FERROPTOSIS)

These contain mechanistic pathway genes and are treated as supporting pathway-level evidence, not as a directional execution signature.

### Specificity/context controls
- HALLMARK_APOPTOSIS
- HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY

These help distinguish a ferroptosis-associated signal from a broad cell-death or oxidative-stress response.

## Analysis

The main script performs:

```text
GSE160543 raw processed counts
        ↓
gene-ID standardization
        ↓
DESeq2
        ↓
Oxaliplatin vs Vehicle
Paclitaxel vs Vehicle
        ↓
whole-transcriptome ranking by DESeq2 Wald statistic
        ↓
fgsea
 ├─ VINIK_2024_24_FERROPTOSIS_BIOMARKERS
 ├─ KEGG_FERROPTOSIS
 ├─ WP_FERROPTOSIS
 ├─ HALLMARK_APOPTOSIS
 └─ HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
        ↓
sample-level module scores from VST expression
        ↓
supporting mechanistic-gene heatmap
        ↓
compact evidence summary
```

## Interpretation

The workflow uses four practical labels:

- **STRONG_SUPPORT**  
  Vinik-24 has positive FDR-significant enrichment and at least one mechanistic ferroptosis pathway is also positively FDR-significant.

- **SUGGESTIVE_SUPPORT**  
  Ferroptosis signatures trend positively but do not meet the stronger convergence rule.

- **MIXED_NON_SPECIFIC**  
  Ferroptosis signatures are positive, but apoptosis/general ROS programs are comparably strong and interpretation requires caution.

- **NO_TRANSCRIPTOMIC_SUPPORT**  
  No reproducible positive enrichment of the ferroptosis resources is observed.

These labels refer only to **transcriptomic evidence**. Bulk DRG RNA-seq cannot by itself prove ferroptotic cell death or identify the exact cell type producing the signal.

## Run

From the repository root:

```r
source("cipn-ferroptosis-screening/scripts/01_GSE160543_screen.R")
```

Results are written to:

```text
cipn-ferroptosis-screening/results/GSE160543/
```

## Main outputs

- `sample_metadata.csv`
- `library_qc.csv`
- `pca.png`
- `DE_Oxaliplatin_vs_Vehicle.csv`
- `DE_Paclitaxel_vs_Vehicle.csv`
- `gsea_all_signatures.csv`
- `signature_members_used.csv`
- `vinik24_ortholog_mapping.csv`
- `sample_level_signature_scores.csv`
- `sample_level_score_tests.csv`
- `mechanistic_panel_heatmap.png`
- `evidence_summary.csv`
- `analysis_summary.txt`
- `sessionInfo.txt`

## Key references

- GEO GSE160543: rat DRG vehicle/paclitaxel/oxaliplatin RNA-seq.
- KEGG Ferroptosis pathway: hsa/rno04216.
- WikiPathways Ferroptosis: WP4313 / MSigDB WP_FERROPTOSIS.
- Vinik Y et al. Advanced Science 2024. DOI: 10.1002/advs.202307263.
