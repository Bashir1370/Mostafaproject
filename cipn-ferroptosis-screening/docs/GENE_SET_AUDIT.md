# GSE160543 gene-set provenance and species audit

## Scope

This document locks the provenance, species handling, and interpretation boundaries for every signature used by `01_GSE160543_screen.R`.

Dataset: GSE160543, rat dorsal root ganglion bulk RNA-seq (Rattus norvegicus).

## Primary signature: Vinik-24

Source: Vinik Y et al. *Programming a Ferroptosis-to-Apoptosis Transition Landscape Revealed Ferroptosis Biomarkers and Repressors for Cancer Therapy*. Advanced Science. 2024;11:2307263. DOI: 10.1002/advs.202307263. PMID: 38441406.

The source study selected 24 biomarkers from dataset-derived and gradient-derived candidates and validated them as ferroptosis-versus-apoptosis transcriptomic biomarkers. They are human biomarkers; they are **not** a rat-validated signature.

For this project, the human panel is mapped to rat before analysis. The locked mapping is in:
`config/vinik_2024_24_human_to_rat_locked.csv`.

Two mappings require explicit reconciliation:
- Human GARS (current HGNC symbol GARS1) maps to current rat Gars1, but the GSE160543 annotation contains the legacy symbol `Gars`; the analysis therefore uses `Gars`.
- Human CALM2 was previously mapped by `babelgene(top=TRUE)` to rat `Calm1`. NCBI/Alliance identifies rat `Calm2` as the ortholog of human CALM2, so the locked analysis uses `Calm2`.

The previous exploratory run therefore used 23/24 Vinik members in the dataset and used Calm1 for CALM2. The audited rerun must use all 24 locked dataset symbols and should be treated as the authoritative run.

## MSigDB signatures

The exploratory run used `msigdbr 26.1.1`, corresponding to MSigDB 2026.1. The call `msigdbr(species = "Rattus norvegicus")` starts from the human MSigDB database by default and computationally maps genes to rat orthologs.

To prevent future database/package updates from silently changing the analysis, the exact rat memberships recovered in the audited exploratory run are frozen in:
`config/msigdb_2026.1_Hs_rat_gene_sets_locked.csv`.

The source/version/count metadata are frozen in:
`config/gene_set_manifest.csv`.

Included sets:
- GOBP_FERROPTOSIS — GO:0097707 / MSigDB M46862; supporting evidence.
- WP_FERROPTOSIS — WikiPathways WP4313 / MSigDB M39768; supporting evidence.
- HALLMARK_APOPTOSIS — MSigDB M5902; context control.
- HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY — MSigDB M5938; context control.

## Interpretation boundaries

1. Ortholog mapping is computational. A pathway mapped from human to rat is not equivalent to a pathway experimentally validated in rat DRG.
2. GOBP_FERROPTOSIS is a process-membership set and includes genes that can promote or limit ferroptosis. Its GSEA NES must not be read as a simple ferroptosis activation/deactivation score.
3. WP_FERROPTOSIS is a mechanistic pathway membership set, not a directional transcriptional signature.
4. HALLMARK_APOPTOSIS and HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY are context controls and must not contribute positive ferroptosis evidence.
5. Bulk DRG RNA-seq can support a ferroptosis-associated transcriptional program but cannot prove ferroptotic cell death or identify the responding cell type.
6. The STRONG/SUGGESTIVE/NO_SUPPORT labels are project-specific heuristic summaries, not clinically validated classifiers.

## Reproducibility rule

The screening script must load the locked CSVs above for the inferential analysis. Live `msigdbr` or live ortholog mapping may be used only for an explicit audit/update step, never silently during the main analysis.
