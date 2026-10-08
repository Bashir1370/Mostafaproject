# Step 03b — resolving ambiguous gene-to-protein mappings

This is a substep of scientific Step 03, not an eleventh study step. It preserves the original first-pass outputs and refines held multi-entry mappings using independent representative evidence. Mapping-resolution implementation version: 1.0.0; protocol: 0.1.0.

## Coverage objective

The user requires mapping most discovery genes before cysteine inventory. Coverage is reported explicitly; there is no retrospectively chosen threshold or permission to force an uncertain mapping. The combined output must be reproduced and reviewed before Step 04. Unresolved genes remain in the audit, and conclusions must reflect mapping coverage.

## Entry

The checksum-bound 854-gene Step 02 snapshot and a successful local Step 03 run, including its cached UniProt JSON responses, are required. The resolver verifies the Step 03 input/helper hashes, raw response hashes, gene universe, original mapping/candidate statuses and accepted sequence hashes. A modified baseline or failed/incomplete input stops the run.

## Representative selection: all conditions required

1. Query every candidate accession of all multi-entry held genes using the EMBL-EBI Proteins API gene-centric service. Preserve pagination, raw responses and source metadata.
2. Identify gene groups through explicit **Ensembl stable gene identifiers** in returned group members. Gene names and MGI labels do not substitute for the exact source ID.
3. Require exactly one representative identifier across all relevant groups. Every original candidate accession must be accounted for by that gene's group. Duplicate identical group evidence is harmless; conflicting parents remain held.
4. The representative must be an original, exact-gene mouse UniProt candidate, not an unverified outside accession. Apply the original displayed-isoform, taxonomy and full-sequence checks again. A representative isoform identifier must match the current displayed isoform when supplied.

Reviewed status, protein length, function, gene significance and eventual structural score are not selection criteria. This permits a valid unreviewed representative and refuses an unsupported reviewed one. A service-proposed representative lacking the original exact gene link remains held for a separate bridge audit; it is not silently substituted.

## What remains outside this resolver

Six no-current-link genes require a separate historical/current Ensembl identity review. Missing cross-references do not by themselves establish that a gene was retired. Fragment sequences are not repaired or inferred. Representative conflicts and unresolved isoform associations remain held. These cases can be revisited with documented gene/transcript/translation evidence; no symbol-only rescue is allowed.

## Ubuntu execution

Python 3.8+ standard library only; no new package installation. After the completed first pass, run:

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/03b_resolve_GSE286387.py
cat oipn-cysteine-structural-screening/results/03b_mapping_resolution/resolution_report.md
```

The script resolves the project path automatically and supports `--project-dir`. The first live run retrieves 58 accession batches, with up to three simultaneous requests, and may take several minutes. Repeat runs validate and reuse cached evidence under `data/raw/gene_centric/step03b/<discovery_sha256>/`. Transient errors are retried; failed requests are not converted into biological missingness. Corrupt or incomplete cache pairs stop execution.

## Combined outputs and gate

`results/03b_mapping_resolution/` contains:

| File | Meaning |
|---|---|
| `gene_protein_mapping.csv` | All 854 genes, original status, resolution method/version, updated accepted mappings or explicit held reasons |
| `canonical_sequences.fasta` | All unique accepted proteins from baseline plus resolved mappings, with gene associations |
| `gene_centric_evidence.csv` | Candidate group representatives, exact-gene members, candidate coverage and resolution reason |
| `step_audit.csv` | One terminal status per input gene |
| `source_manifest.csv` | Every gene-centric page URL, response SHA256, access time and available release header |
| `input_checksums.csv` | Discovery/configuration, baseline mapping/source manifest, original helper and resolver hashes |
| `summary.json`, `resolution_report.md` | Coverage, gene/protein counts, reasons and reference limitations |
| `SUCCESS.txt`, `FAILURE.txt` | Generation success versus failure; outputs are invalid without success or when failure exists |

The status is **STEP03B_GENERATED_REVIEW_PENDING**. Baseline accepted mappings remain unchanged; they are not dropped when gene-centric information is unavailable. The original `results/03_protein_mapping/` remains intact. After review, this combined mapping/FASTA is the prospective Step 04 input; the baseline FASTA alone is not the expanded set. No cysteine count filter or structural score is calculated here.

## Reference limitation

The original UniProt sequences/links remain bound to their recorded release. If the gene-centric service supplies a release header, it must match that release; otherwise execution stops. In validation the service supplied no release header. Its response content hashes and retrieval times are recorded, and its relationship evidence is checked against the original UniProt entries. **Do not label the gene-centric snapshot as release 2026_03 merely because the sequence snapshot uses that release.** A live relationship change may affect workstation reproduction and requires review.

## Real-input validation

The accepted 854-gene input and UniProt sequence release 2026_03 produced 504 additional resolved genes. Combined mapping: **789 genes (92.4%), 787 unique proteins**, 64 held and one unassessable fragment. Held reasons: 49 service representatives lack an original exact gene link; seven have only noncanonical/unresolved isoform links; six still have no current exact link; one has multiple displayed isoform identifiers; one has conflicting gene-centric representatives. These are implementation-validation results, not a workstation Step 03b snapshot.

Tests cover explicit representative evidence, repeated group rows, absence of reviewed-status/symbol shortcuts, complete candidate accounting, parent conflicts, absent original gene links, noncanonical representatives, sequence quality, baseline input tampering and combined gene/FASTA integrity. Earlier Steps 01–03 tests remain in the suite. Cache replay reproduces mapping totals.

## Primary source

[EMBL-EBI Proteins API documentation](https://www.ebi.ac.uk/proteins/api/doc/) provides the gene-centric service alongside protein/genome information. API relationships are identity evidence for selecting a representative; they do not establish the isoform expressed in DRG or actual oxidation.
