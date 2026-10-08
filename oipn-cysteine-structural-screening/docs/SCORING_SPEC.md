# Scoring specification — v0.1.0

## Primary target

Relative structural prioritization of eligible free-thiol-compatible cysteine sites. This is not calibrated for all ROS chemistry, oxidation kinetics, S-glutathionylation specifically, or ROS exposure in a particular compartment. No assumed DRG pH is needed for this rank-based index.

## Reference and feature definitions

Reference = all unique quality/state-eligible discovery sites with complete SG SASA and pKa, after selecting one primary structure context per site. Freeze site membership and raw values. Do not duplicate a site across equivalent assembly chains to weight it multiple times.

- A_j: SASA of SG, in Å², probe radius 1.4 Å; same implementation/radii across structures.
- K_j: predicted Cys pKa from the recorded primary tool/preparation.
- D_j: distance to nearest other represented cysteine SG; supplemental only.

For N >= 2, rank x ascending with average ranks for exact ties (smallest rank=1), and define R(x)=(rank(x)-1)/(N-1). All equal values give R=0.5. N=1 gives R=0.5 and an insufficient_reference flag. N=0 yields empty outputs and no score. No early rounding; sort/order on full precision.

S_j = 0.5 R(A_j) + 0.5 R(-K_j).

Equal weights are an a priori exploratory choice, not experimentally established weights. Higher accessibility and lower pKa are working ranking directions to be evaluated independently. Values lie in [0,1] but are not probabilities. Scores depend on the reference population and cannot be directly compared across independently normalized datasets.

## Protein score and ties

S_p = max(S_j) over complete eligible sites of protein p.

All Cys attaining the maximum are reported as best sites. Display ties by stable protein ID; do not invent biological order within ties. Canonical proteins map back to stable mouse gene IDs; enrichment counts a gene only once.

Report n_cys_total, n_cys_structurally_covered, n_cys_quality_eligible, n_cys_state_eligible and n_cys_scored. Scored-Cys coverage = n_cys_scored / n_cys_total. Also report observed sequence coverage when available; it is a different metric.

## Shortlist

Threshold = empirical 0.75 quantile of protein scores using linear interpolation (R type 7 / NumPy linear convention). Include S_p >= threshold, preserving all ties. With extensive ties this can exceed 25%; report actual membership rather than force an arbitrary count. No absolute score cutoff labels a protein biologically sensitive/resistant.

## Missing and special cases

- Missing pKa or SASA: no primary site score; retain failure reason.
- Zero valid SASA: scoreable, not automatically resistant.
- No eligible Cys: protein score=NA and explicit status.
- No other represented Cys: D=NA; no penalty.
- Structural disulfide/metal/incompatible PTM: separate branch, not forced into free-thiol calculations.
- Partial coverage: maximum refers only to assessed sites. Missing sites can hide the true most-sensitive site.

## Fixed sensitivity scenarios

- Accessibility/pKa weights: 0.25/0.75 and 0.75/0.25; main remains 0.5/0.5.
- Predicted neighbors: pLDDT >=90 instead of >=70; target stays >=90.
- Experimental-only and predicted-only subsets; additionally compare on a common protein/site set to distinguish selection effects.
- Suitable alternative conformation/assembly where available, selected independently of score.
- Protein mean of up to three highest-scoring eligible sites, as a secondary aggregation only.
- Examine length, n_cys_total, n_cys_scored and coverage associations with max score.

Freeze the eligible reference membership for weight-only scenarios. Quality/source scenarios change eligibility and must explicitly export their reference sets; do not interpret normalization changes as physical changes in proteins.
