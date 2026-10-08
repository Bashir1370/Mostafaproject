# Step 01 — official dataset audit and local reproduction

## Scope

This substep downloads official GEO SOFT metadata and the submitter-designated raw gene count table for GSE286387. It checks species, tissue, exact treatment labels, ten sample titles, sample-column identity, integer/nonnegative values and stable-ID duplication. It does not select DEGs, fit DESeq2, exclude samples or declare the complete Step 01 passed.

## Run

From the repository root on Ubuntu:

```bash
python3 oipn-cysteine-structural-screening/scripts/01_audit_GSE286387.py
```

Requires Python >=3.8, standard library only. Downloads go to `data/raw/GSE286387`; outputs go to `results/01_dataset_audit`, both inside this project directory. Paths resolve from the script, independently of the terminal working directory. Existing local files can be re-audited with `--offline`. Online runs fetch fresh files and recompute checksums. A changed checksum requires review against the committed snapshot before continuing. A failed audit returns nonzero and writes FAILURE.txt; successful files from a failed attempt must not be used. The script removes stale matrix/summary/report on failure. An interrupted process must likewise be rerun successfully before outputs are used.

## Actual results, 2026-10-08

| Check | Result |
|---|---|
| Official input rows | 40,757 |
| Unique mouse Ensembl IDs | 40,481 |
| Unique protein-coding gene IDs in deposited annotation | 21,534 |
| Samples | 5 vehicle + 5 oxaliplatin |
| Duplicate IDs / redundant rows | 167 / 276 |
| Duplicate count conflicts | 0 |
| Fractional, negative or nonfinite counts | 0 |
| All-zero unique genes | 11,082 |
| Planned >=10 counts in >=5 samples filter | 18,538 unique genes pass; filter not yet applied to exported matrix |

The supplementary table explicitly labels raw counts in GEO processing metadata. It includes all-zero genes and many genes with maximum count below 100, rather than only the paper's 91 selected DEGs. The paper's normalized 91-gene matrix is not used. The export is genome-scale; completeness of every reference feature is not independently established.

All duplicated IDs have exactly identical ten-sample counts and identical biotype. Repeated annotation rows are consistent with one-to-many annotation export. Keep one count vector per stable ID, never sum these copies; retain all source annotations and the deduplication audit. Conflicting duplicates terminate the script. Mapping later uses stable Ensembl IDs rather than blindly selecting an Entrez alternative. This is an input-integrity implementation decision, not a new DEG threshold.

## Metadata and remaining review

- GEO sample metadata explicitly provides DRG, Mus musculus, C57BL/6 and treatment. The series summary/design supplies male mice, intraperitoneal 10 mg/kg weekly for eight weeks and collection after that treatment period.
- GSM8726281–GSM8726285 are treatment replicates 1–5; GSM8726286–GSM8726290 are control replicates 1–5. Count-column order differs from GEO order and is mapped by exact titles, never position or regex guesses.
- GEO reports five RNA samples per group. It does not expose per-animal IDs or explicit pooling details. `animal_or_pool_id=NA` is intentional. The study-level design suggests independent samples, but animal-level independence/pooling needs full-methods review before freezing inclusion. Neither matching replicate suffixes nor absence of animal IDs implies pairing.
- Age and sequencing batch are not specified in GEO sample fields. No batch or pairing term is invented. The exact interval between final injection and collection is unresolved.
- Alignment metadata specifies STAR 2.5.1a and Ensembl release 76. Later current sequence mapping must record release/ID changes; deposited biotype is a historical annotation.
- Original processing says a max-count <100 filter was used; many such genes remain in the deposited matrix. This suggests filtering occurred downstream of this export, but its exact stage is unresolved. Our declared >=10 in >=5 samples filter remains unchanged.
- R sample QC must inspect library size, expressed-feature counts, transformed PCA and sample distances. Outliers trigger review, not automatic deletion. The design and final sample manifest are frozen before DE fitting.

## Outputs and verification

`sample_manifest.csv` (all samples held for QC), `source_manifest.csv` (URL/checksum/time), `count_matrix.tsv` (all unique gene count vectors), `gene_annotations.csv`, `source_gene_annotations.csv`, `duplicate_gene_audit.csv`, `step_audit.csv`, `audit_summary.json`, and `qc_report.md`.

A small real-run snapshot is committed under `docs/audits/`; large count inputs and regenerated outputs are not committed. Tests:

```bash
python3 -m unittest discover -s oipn-cysteine-structural-screening/tests -v
```

Four tests pass: identical annotation duplicates do not inflate counts; conflicting duplicates fail; fractional/negative/missing/nonfinite values fail; mismatched sample identities fail. Actual official files were audited successfully using Python 3.12.14. No R analysis has yet been executed.

Sources: [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286387), exact official file URLs/checksums in `docs/audits/source_manifest.csv`, [study DOI](https://doi.org/10.1007/s12035-025-05463-7).
