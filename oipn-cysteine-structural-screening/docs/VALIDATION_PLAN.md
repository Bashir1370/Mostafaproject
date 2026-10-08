# Independent evaluation and robustness

## What is evaluated

The frozen equal-weight index is evaluated as a ranking of cysteine thiol oxidation susceptibility. H2O2 chemoproteomic measurements are a candidate benchmark, not direct proof of oxidation in mouse DRG/OIPN. Keep oxidation chemistry, cell type, species, oxidant concentration, duration, abundance correction and measurement coverage visible.

## Reference acceptance

Accept a measured site only with traceable accession/sequence, unambiguous residue mapping and a documented experimental response. Negatives are measured low/nonresponding sites under the tested conditions, not every unreported Cys. Preserve continuous responses if available. Declare label rules from experimental methodology before examining predicted scores. Ambiguous/unmeasured sites are excluded from performance calculations but counted in the audit.

Candidate source: Fu et al. 2017, PMID 28827280, DOI 10.1074/mcp.RA117.000108. Access and supplementary-table audit are pending. Its human-cell measurements must not be described as mouse-DRG validation. Benchmark structures are acquired/prepared under the same stated rules for the benchmark species; do not force human experimental labels onto mouse sites without a separate explicit mapping analysis.

## Comparisons and uncertainty

Compare the composite with SASA-only and -pKa-only on identical evaluated sites. When labels permit, report AUROC and AUPRC, class prevalence, precision at the predeclared shortlist fraction and protein-block bootstrap confidence intervals (seed 20261008; 2000 resamples). Resample proteins rather than treating their sites as independent. Skip a metric with a documented reason if a resample has only one class. For continuous outcomes, report a rank correlation and compatible uncertainty.

Report protein-level dependence on length, Cys number and coverage; use matched/permutation controls if sample sizes permit. Known redox-site annotations cannot simultaneously supply positive score points and serve as independent validation.

## Robustness outputs

Run the scenarios in SCORING_SPEC.md. Report eligible counts, common-set Spearman rank correlation, shortlist Jaccard overlap, per-protein rank shifts and per-candidate shortlist retention. Missingness due to changed eligibility is separate from falling below a threshold. Retention is descriptive; no post-hoc cutoff is used to claim robustness.

## Interpretation gate

If benchmark data are inaccessible, label the index unvalidated. If performance does not justify the composite or is highly unstable, publish/report it only as an exploratory feature-based index. If it has independent support, limit the claim to the tested chemistry/context. No achieved metric by itself establishes OIPN-specific oxidation.

## Model amendments

v0.1.0 weights remain fixed. A later trained/calibrated model is a new version with protein/family-separated training and test sets and preprocessing fitted on training only. Do not tune on the OIPN shortlist or reuse the tuning data as the independent test.

## Independent OIPN expression evidence

GSE125002 is a complementary dataset candidate, not a pooled extension of GSE286387. Analyze its own matched contrast and report effect direction/FDR for discovery candidates. It does not rediscover a replacement shortlist or redefine the structural score. Differences in route/time/pooling limit interpretation of absent replication.
