# Mostafaproject — Necroptosis–Ferroptosis Classification Framework

A **disease-agnostic, mechanism-first, evidence-weighted and empirically calibrated** framework for distinguishing **necroptosis (NEC)** from **ferroptosis (FPT)** using transcriptomic data.

## Scope

This project is intentionally independent of any specific disease, tissue, drug, or clinical application.

The goal is to build a reproducible molecular classification framework that separately measures:

1. **NEC-MPS** — necroptosis molecular permissiveness.
2. **FPT-MPS** — ferroptosis molecular permissiveness.
3. **MCI** — mechanistic completeness of each program.
4. **ESR** — empirical state resemblance to experimentally confirmed perturbation states.
5. **NFSI** — Necroptosis–Ferroptosis Specificity Index.

## Why these two programs?

They have:
- biologically distinct execution mechanisms;
- strong causal experimental foundations;
- multiple independent perturbation transcriptomes;
- inhibitor/rescue designs;
- orthogonal induction strategies;
- neural or neuroimmune validation datasets.

Mechanistically:

```text
Necroptosis:
RIPK1/RIPK3 necrosome competence
        ->
MLKL activation / membrane disruption

Ferroptosis:
Iron availability
        ->
PUFA-phospholipid susceptibility
        ->
phospholipid peroxidation
        ->
failure of anti-ferroptotic defenses
```

## Core principles

- Mechanism before association.
- RNA expression is interpreted as **molecular competence/permissiveness**, not direct proof of a post-translational execution event.
- Genes whose function depends primarily on phosphorylation, oligomerization, catalytic activity or isoform balance are not naively interpreted from mRNA.
- Mechanistic priors are recalibrated using independent perturbation transcriptomes.
- Directly manipulated genes cannot validate themselves in the same experiment.
- Rescue designs are prioritized.
- Independent studies, not individual samples or contrasts, are the replication unit.
- Generic oxidative/inflammatory responses are explicitly penalized during specificity calibration.
- The classifier allows **NEC-dominant, FPT-dominant, Mixed, and Indeterminate** outcomes.
- Thresholds are learned from reference data and frozen before any future application dataset is evaluated.

## Repository structure

```text
Mostafaproject/
├── README.md
├── CHANGELOG.md
├── ROADMAP.md
├── docs/
│   ├── PROTOCOL.md
│   ├── SCORING_SPEC.md
│   ├── DECISION_LOG.md
│   ├── REFERENCES.md
│   ├── REFERENCE_DATA_QC.md
│   └── EMPIRICAL_CALIBRATION_WORKFLOW.md
├── config/
│   ├── necroptosis_gene_dictionary.csv
│   ├── fpt_gene_dictionary.csv
│   ├── reference_datasets.csv
│   ├── reference_contrasts.csv
│   └── analysis_parameters.csv
├── R/
│   └── evidence_functions.R
└── scripts/
    ├── 00_inspect_reference_datasets.R
    └── 01_GSE255459_fpt.R
```

## Current status

**v0.2.0 — NEC–FPT redesign**

Completed:
- NEC/FPT mechanistic architecture
- initial mechanistic evidence dictionaries
- rescue/orthogonal reference-data library
- direction-aware empirical evidence logic
- anti-circularity policy
- MPS/MCI/ESR/NFSI definitions
- pre-specified contrast registry

Next:
- inspect all reference matrices
- reprocess SRA-only studies
- estimate per-study empirical evidence
- calculate E and S
- derive empirical NEC/FPT signatures
- calibrate thresholds
- leave-one-study-out validation
- freeze v1.0
