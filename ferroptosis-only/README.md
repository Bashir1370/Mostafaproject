# Ferroptosis-Only Modeling Workspace

This folder contains an independent, disease-agnostic systems-biology framework for modeling **ferroptosis alone**, without requiring comparison against another regulated cell-death program.

The parent NEC–FPT framework remains unchanged.

## Objective

Build a reproducible transcriptome-based model that answers:

1. Is the molecular environment permissive for ferroptosis?
2. Are the major mechanistic modules of ferroptosis jointly supported?
3. Does the sample resemble experimentally validated ferroptotic transcriptomic states?
4. How strong and reproducible is the evidence across independent perturbation datasets?

## Core design

The model separates:

- **Mechanistic Core Score**
- **Regulatory/Context Score**
- **Empirical Ferroptosis Resemblance (FPT-ESR)**
- **Mechanistic Completeness (FPT-MCI)**

Unlike the NEC–FPT comparison model, this project does **not** require a cross-program specificity index.

## Initial mechanistic modules

- F1 — Iron availability
- F2 — PUFA-phospholipid susceptibility
- F3 — Lipid-peroxidation machinery
- F4 — Anti-ferroptotic defense
- F5 — Context/regulatory modifiers
- FPT-ESR — empirical ferroptosis state resemblance

## Repository structure

```text
ferroptosis-only/
├── README.md
├── ROADMAP.md
├── docs/
│   ├── PROTOCOL.md
│   ├── GENE_CURATION.md
│   └── REFERENCES.md
└── config/
    ├── ferroptosis_master_gene_universe.csv
    └── reference_datasets.csv
```

## Current status

**v0.2.0 — empirical rescue calibration in progress**

Completed:
- expanded ferroptosis Master Gene Universe
- evidence-tier framework
- RNA-interpretability rules
- primary mechanistic reference library
- CORE / EXTENDED / MECHANISTIC / EMPIRICAL_ONLY / PENDING classes

Completed additionally:
- critical review of provisional CORE genes
- locked GSE182638 primary rescue-analysis design
- executable GSE182638 DESeq2 + ashr + Fer-1 rescue calibration script

Next:
- GSE182638 primary rescue analysis completed
- GSE247883 locked as independent rescue replication
- preflight ENA manifest, raw-read download, Docker/Salmon quantification, and DESeq2/ashr replication scripts added
- next: run GSE247883 preflight and inspect raw-read download size before acquisition.

See:
- [GENE_CURATION.md](docs/GENE_CURATION.md)
- [REFERENCES.md](docs/REFERENCES.md)
