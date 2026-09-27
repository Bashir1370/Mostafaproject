# GSE247883 Independent Rescue Replication Plan

## Role

GSE247883 is an **independent rescue-replication dataset**, not a discovery dataset.

Verified design:
- Homo sapiens
- A549 cells
- DMSO, n=3
- RSL3 2 uM for 24 h, n=3
- RSL3 2 uM + Ferrostatin-1 5 uM for 24 h, n=3
- paired-end RNA-seq
- Illumina NovaSeq 6000
- GEO reports GRCh38 / Ensembl release 101 processing
- raw reads are available through SRA/ENA
- GEO processed supplement is FPKM and is not used as DESeq2 input

## Primary scientific question

Do transcriptional features discovered in the GSE182638 Fer-1 rescue experiment reproduce, with the same pre-specified direction, in an independent A549 RSL3/Fer-1 experiment?

## No data leakage rule

GSE247883 must not select a replacement discovery signature.

The discovery universe is generated only from GSE182638:
- require cross-cell-line rescue-direction agreement in MM1R and MM1S
- keep the GSE182638 signed rescue direction
- freeze that table before evaluating GSE247883
- save an MD5 checksum of the frozen snapshot

GSE247883 can:
- support a frozen candidate
- contradict a frozen candidate
- leave a candidate not testable

It cannot create a new ESR candidate for the primary signature.

## Mechanistic genes

Mechanistic causal role and acute RNA response are treated separately.

For ACSL4, LPCAT3, POR, GPX4 and the extended mechanistic universe:
- report acute RSL3 effect
- report Fer-1 reversal effect
- report rescue-concordance descriptively
- do **not** assume the causal promoter/suppressor sign is the expected acute transcript direction

This prevents conflating baseline susceptibility with compensatory transcriptional response.

## Raw-data workflow

1. Resolve the locked SRX accessions against ENA.
2. Record exact run accession, FASTQ URLs, MD5 and compressed size.
3. Download paired FASTQ files only after reviewing total size.
4. Quantify against Ensembl release 101 transcriptome with Salmon.
5. Import Salmon estimates with tximport.
6. Fit DESeq2 model:
   `~ condition`
7. Estimate posterior effects with ashr.
8. Evaluate two contrasts:
   - RSL3 - DMSO
   - RSL3 - (RSL3 + Fer-1)
9. Test only the frozen GSE182638 ESR candidate universe.
10. Save QC, replication evidence and audit metadata.

## Why the model is ~ condition

Unlike GSE182638, the replicate labels in GSE247883 are independent biological replicates and are not passage-matched blocks across conditions. Therefore replicate is not included as a blocking factor.

## Execution order on Windows

From repository root:

```powershell
powershell -ExecutionPolicy Bypass -File "ferroptosis-only\scripts\02a_GSE247883_prepare_raw.ps1"
```

This is a preflight only. It writes an ENA manifest and reports total compressed FASTQ size without downloading reads.

After reviewing the size:

```powershell
powershell -ExecutionPolicy Bypass -File "ferroptosis-only\scripts\02a_GSE247883_prepare_raw.ps1" -DownloadReads
```

Quantification uses the official Salmon Docker image:

```powershell
powershell -ExecutionPolicy Bypass -File "ferroptosis-only\scripts\02b_GSE247883_salmon_quant.ps1" -Threads 8
```

Then in R, from repository root:

```r
source("ferroptosis-only/scripts/02_GSE247883_rescue_replication.R")
```

## Main outputs

- sample_metadata_locked.csv
- salmon_mapping_qc.csv
- vst_pca.png
- vst_sample_correlation.png
- per_contrast_effects_all_genes.csv
- mechanistic_genes_acute_RNA_effects.csv
- mechanistic_genes_rescue_concordance_descriptive.csv
- GSE182638_ESR_discovery_snapshot.csv
- GSE182638_ESR_discovery_snapshot.md5.txt
- GSE247883_frozen_ESR_candidate_replication.csv
- analysis_summary.txt
- sessionInfo.txt

## Interpretation

A positive `replication_e` means that both:
- RSL3 induction aligns with the direction pre-specified by GSE182638, and
- the RSL3 vs RSL3+Fer-1 contrast also aligns with that direction.

A negative value is contradictory evidence.

No hard replication threshold is frozen at this stage; cross-study calibration will be defined only after the rescue and cross-inducer studies are assembled.
