# Decision Log — NEC–FPT Framework

## DL-001 — Disease-agnostic development
The classifier is developed independently of any future disease, tissue, drug or clinical application.

## DL-002 — Necroptosis replaces disulfidptosis
The active framework was redesigned from DPT–FPT to **NEC–FPT** because necroptosis combines:
- strong mechanistic separation from ferroptosis;
- direct causal evidence in neuropathic-injury biology;
- multiple independent GEO perturbation datasets;
- rescue/inhibitor controls;
- orthogonal RIPK3 activation strategies;
- human neural/neuroimmune validation datasets.

The earlier DPT design remains recoverable through Git history but is not part of the active model.

## DL-003 — Mechanism-first design
The framework is not a flat death-gene score. NEC and FPT are decomposed into mechanistic components.

## DL-004 — NEC core is deliberately small
Direct NEC RNA scoring begins with:
- RIPK1/RIPK3 necrosome competence;
- MLKL execution competence.

Upstream triggers and checkpoints are not allowed to dominate the core score.

## DL-005 — Post-translational execution is not inferred from RNA
pRIPK3, pMLKL and MLKL oligomerization define execution more directly than transcript abundance. RNA scores therefore represent competence/permissiveness and state resemblance, not proof of execution.

## DL-006 — CASP8/FADD/CFLAR are annotations initially
Their effects depend on catalytic state, complex composition and/or isoforms. Bulk RNA direction is not safely converted into a simple positive/negative NEC score.

## DL-007 — RIPK1 is important but not universally required
RIPK1 receives lower prior/specificity status than RIPK3/MLKL because RIPK3 can be activated through RIPK1-independent routes such as ZBP1.

## DL-008 — MPS and MCI are distinct
MPS measures overall molecular permissiveness.
MCI measures completeness and penalizes a missing required module.

## DL-009 — ESR remains independent
Empirical transcriptomic resemblance is not merged into the literature-defined mechanistic score.

## DL-010 — Three-component gene weighting
Final weights combine:
- mechanistic prior M;
- empirical reproducibility E;
- specificity S.

## DL-011 — No p-value-only weighting
Effect direction, effect magnitude and sign confidence are used. Contradictions are retained.

## DL-012 — Anti-circularity
A directly manipulated gene cannot validate itself in the same experiment.

## DL-013 — Independent study is the replication unit
Multiple contrasts, cell lines or time points from one study are aggregated before cross-study E.

## DL-014 — GSE108621 is primary NEC calibration
Its TNF-only arm helps remove generic inflammatory transcription and its Nec-1s arm provides pharmacological rescue.

## DL-015 — Orthogonal NEC validation is mandatory
Direct/optogenetic RIPK3 datasets are held for robustness testing rather than allowing all NEC references to share the same TSZ induction chemistry.

## DL-016 — GSE182638 is primary FPT rescue calibration
It provides raw counts, two independent human cell lines, RSL3 induction and Ferrostatin-1 reversal.

## DL-017 — Cross-inducer FPT evidence is required
Erastin, RSL3 and additional inducers are used to reduce drug-specific transcriptomic artifacts.

## DL-018 — Generic stress is a competitor
Oxidative stress datasets are used to down-weight nonspecific stress genes. Additional inflammatory/apoptotic comparators may be added later under a versioned protocol.

## DL-019 — Specificity index renamed NFSI
`NFSI` = Necroptosis–Ferroptosis Specificity Index.
Positive values favor NEC-specific evidence; negative values favor FPT-specific evidence.

## DL-020 — Four top-level classes
- NEC-dominant
- FPT-dominant
- Mixed
- Indeterminate / Neither

## DL-021 — Thresholds are learned
No final biological cutoff is hand-picked. Thresholds require bootstrap and whole-study holdout validation.

## DL-022 — Freeze before external application
Dictionaries, equations, empirical signatures and thresholds must be versioned and frozen before external use.
