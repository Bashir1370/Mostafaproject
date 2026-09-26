# Decision Log

## DL-001 — Disease-agnostic development
The classifier is developed independently of any disease, tissue, drug, or clinical application.

## DL-002 — Mechanism-first design
The framework is not a conventional flat gene-set score. DPT and FPT are decomposed into mechanistic stages.

## DL-003 — Separate evidence classes
Causal regulators, protective systems, metabolic modifiers, and proteins damaged post-translationally are not treated as equivalent.

## DL-004 — DPT cytoskeletal substrates are not direct RNA activity markers
ACTB, MYH9, FLNA and related proteins remain mechanistic annotations unless empirical RNA evidence demonstrates diagnostic utility.

## DL-005 — Direction-aware scoring
Promoters and protective genes have opposite expected directions. Protective genes are scored inversely.

## DL-006 — SLC7A11 is shared/contextual
SLC7A11 is central to DPT permissiveness but also central to ferroptosis resistance. It contributes strongly to DPT context but must not directly dominate DSI.

## DL-007 — DPT initiation uses a biological gate
DPT requires both cystine-loading capacity and reducing-capacity vulnerability. D1 and D2 are coupled using a geometric interaction.

## DL-008 — MPS and MCI measure different concepts
MPS measures overall permissiveness/intensity. MCI measures mechanistic completeness and is weakest-link sensitive.

## DL-009 — ESR remains separate
Empirical transcriptomic resemblance is kept separate from mechanistic permissiveness.

## DL-010 — Three-component weighting
Final gene weights combine:
- mechanistic prior M
- empirical reproducibility E
- specificity S

## DL-011 — No p-value-only weighting
Effect direction, effect magnitude and uncertainty are prioritized. Contradictions are explicitly penalized.

## DL-012 — Anti-circularity
A directly manipulated gene cannot use the same experiment to validate its own transcriptomic behavior.

## DL-013 — Independent study is the replication unit
Multiple contrasts or cell lines from one study are aggregated before cross-study evidence combination.

## DL-014 — Rescue designs receive special value
FPT responses reversed by Ferrostatin-1 are more specific than inducer-vs-control behavior alone.

## DL-015 — Generic oxidative stress is a competitor
H2O2 and similar stress datasets are used to down-weight nonspecific stress-response genes.

## DL-016 — Four top-level classes
Possible outputs:
- DPT-dominant
- FPT-dominant
- Mixed
- Indeterminate / Neither

The classifier is not forced to assign every sample to one mechanism.

## DL-017 — Thresholds are learned
Cutoffs must be learned from reference data and evaluated using bootstrap and leave-one-study-out validation.

## DL-018 — Freeze before external application
Gene dictionaries, equations, calibration rules and thresholds are versioned and frozen before future application.
