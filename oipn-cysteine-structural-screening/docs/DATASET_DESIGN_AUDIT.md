# GSE286387 design acceptance — 2026-10-08

## Primary source and reported facts

[Full study methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC13356913/), DOI [10.1007/s12035-025-05463-7](https://doi.org/10.1007/s12035-025-05463-7), and [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286387) were reviewed with the deposited counts and sample metadata.

The study reports male C57BL/6 mice aged 8–10 weeks at purchase, weekly intraperitoneal oxaliplatin 10 mg/kg for eight weeks, and 5% dextrose controls. Bilateral lumbar L1–L6 DRGs were collected; euthanasia occurred within a week after the experiments. RNA integrity exceeded 8. Indexed libraries were pooled for sequencing on NovaSeq 6000. Alignment used STAR and Ensembl release 76. Figure 4 reports five RNA-seq samples per group; the data-availability statement describes ten mice. The manuscript's staged DEG selection uses additional fold-change thresholds.

## Interpretation and frozen decision

The ten-mouse description, five-per-group figure and ten deposited samples support treating the RNA libraries as ten independent mouse units. This is an interpretation of the public study design, not a recovered animal-ID crosswalk. Library pooling occurs after indexing and does not establish cross-animal tissue pooling. Original animal identifiers and sequencing batch remain unknown. `analyst_unit_id` labels libraries for audit only and must not be presented as original animal IDs.

Accept all ten samples for Step 02 with the documented limitations. Use the unpaired model `~ condition`, reference `control`, contrast `oxaliplatin / control`. Unknown batch cannot be included as a covariate and cannot be assumed absent. Exact age at collection and time since final injection remain unresolved. This acceptance does not establish that every mouse had a separately documented behavioral phenotype.

## Reviewed user QC

The user's committed [Ubuntu snapshot](audits/ubuntu_2026-10-08/) contains R 4.3.3 / DESeq2 1.42.1 diagnostics. Counts and source checksums agree with the official-input audit. Ten samples passed execution integrity; library totals are 52,963,470–60,984,426 and pairwise VST correlations are approximately 0.9866–0.9962. PCA explains 35.23% and 24.73% on its first two axes. `DRG_10mg_rep1` differs from the other treatment libraries on PCA, but the available count-level metrics provide no independent technical or identity justification to exclude it. It is retained. Count QC cannot replace read-level mapping or degradation metrics.

## Reproducible acceptance records

- `config/GSE286387_samples.csv`: accepted manifest; unchanged GSM identities and conditions, explicit unknowns, no exclusions.
- `config/GSE286387_design.json`: model, contrast, sample-manifest and count-matrix SHA256, evidence and limitations.
- Original `results/01_dataset_audit/sample_manifest.csv` and the user's QC snapshot remain historical inputs; their `held` status reflects the earlier pre-review stage.

This resolves the Step 01 design gate under public-metadata limitations. No scientific threshold, protein score or protocol version was changed. A later discovery of technical failure or a different biological-unit design requires a logged revision and rerunning downstream analysis.
