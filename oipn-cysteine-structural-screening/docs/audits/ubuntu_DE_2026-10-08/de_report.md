# Step 02 — differential expression

Status: STEP02_COMPLETE

- Design: ~ condition; oxaliplatin / control; 5 RNA libraries per group; independent mice supported by public study design
- Tested genes, all biotypes: 18538
- Valid p-values for BH: 18503
- Significant genes, all biotypes (BH FDR <0.05): 893
- Significant protein-coding genes: 854
- Protein-coding up / down: 477 / 377
- No fold-change cutoff; no post-DE baseMean threshold; independentFiltering=FALSE
- DESeq2 1.42.1; R 4.3.3

All 10 samples retained; default DESeq2 Cook's handling kept. Missing results are not treated as nonsignificant evidence of resistance.
Biotype is from the deposited historical annotation. Stable gene IDs remain distinct; symbols are not merged.
Unreported batch and animal identifiers remain limitations; the reported 10-mouse study supports an unpaired design.
DEG numbers need not equal the paper: our prefilter/BH universe differs and we impose no fold-change cutoff.
Expression differences do not establish protein oxidation. Step 03 mapping follows this discovery list.
