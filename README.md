# Mostafaproject — DPT–FPT Classification Framework

A **disease-agnostic, mechanism-first, evidence-weighted and empirically calibrated** framework for distinguishing **disulfidptosis (DPT)** from **ferroptosis (FPT)** using transcriptomic data.

## Scope

This project is intentionally independent of any specific disease, tissue, drug, or clinical application.

The goal is to build a **reproducible classifier** that can answer four distinct questions:

1. **DPT-MPS** — Is the molecular environment permissive for disulfidptosis?
2. **FPT-MPS** — Is the molecular environment permissive for ferroptosis?
3. **MCI** — Is the required mechanism complete, rather than only partially present?
4. **DSI** — When mechanism-specific evidence exists, is it more consistent with DPT or FPT?

A separate **ESR (Empirical State Resemblance)** layer measures similarity to experimentally confirmed reference states.

## Core principles

- Mechanism before association.
- Causal regulators, protective systems, modifiers, and damaged substrates are not treated as equivalent.
- Direction of effect is explicit.
- Protein-level damage is not automatically converted into an RNA marker.
- Shared redox biology is retained for permissiveness but down-weighted for DPT-vs-FPT specificity.
- Mechanistic priors are recalibrated using independent perturbation transcriptomes.
- A directly manipulated gene cannot validate itself in the same experiment.
- Independent studies, not individual samples or contrasts, are the replication unit.
- The classifier supports **DPT-dominant, FPT-dominant, Mixed, and Indeterminate** outputs.
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
│   └── REFERENCES.md
└── config/
    ├── dpt_gene_dictionary.csv
    ├── fpt_gene_dictionary.csv
    └── reference_datasets.csv
```

## Current status

**v0.1.0 — architecture freeze candidate**

Completed:
- mechanistic stage architecture
- initial DPT/FPT evidence matrices
- direction-aware scoring logic
- mechanistic prior definition
- empirical reproducibility concept
- specificity concept
- MPS/MCI/ESR/DSI definitions
- anti-circularity rules
- classification logic
- initial reference perturbation library

Not yet completed:
- re-analysis of reference datasets
- empirical score `E`
- specificity score `S`
- final per-gene weights
- threshold calibration
- leave-one-study-out validation
- production implementation and tests

See [docs/PROTOCOL.md](docs/PROTOCOL.md) for the full workflow.
