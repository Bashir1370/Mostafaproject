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
│   └── PROTOCOL.md
└── config/
    ├── ferroptosis_master_gene_universe.csv
    └── reference_datasets.csv
```

## Current status

**v0.1.0 — workspace initialized**

Next task: systematic curation of the ferroptosis gene universe before any scoring weights are frozen.
