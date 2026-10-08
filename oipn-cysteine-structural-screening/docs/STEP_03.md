# Step 03 — exact mouse gene mapping and canonical sequences

## Input and objective

The accepted Step 02 snapshot supplies 854 unique mouse stable gene IDs. Its discovery-file SHA256 is bound in `config/GSE286387_step03_input.json`. Step 03 checks this binding before requesting data. It preserves every input gene and does not require a cysteine; Cys inventory is Step 04.

The deposited annotation is Ensembl release 76. This implementation retrieves **current UniProtKB entries with exact Ensembl GeneId cross-references**, preserving the versioned gene, transcript and protein identifiers found there. Removing the numeric version suffix permits comparison of stable IDs; it does not demonstrate unchanged annotation or reconstruct release-76 proteins. Missing current links remain held for a later version-aware rescue audit. No symbol, ortholog or longest-sequence fallback is used.

## Fixed first-pass criteria

| Check | Pass | Otherwise |
|---|---|---|
| Gene association | Exact stable Ensembl GeneId in returned entry | `held`: no current exact link |
| Species | Taxonomy 10090 | `held`: wrong species |
| Accessions | Exactly one mouse UniProtKB primary accession for the gene | `held`: multiple candidates; keep all options |
| Canonical association | Explicit displayed-isoform cross-reference, or entry-level cross-reference without a specified isoform | `held` if only noncanonical/unresolved isoforms are linked |
| Sequence | Nonempty, matching declared length and MD5 when provided, standard amino acids plus U/O, not flagged Fragment | `unassessable`: missing, ambiguous, inconsistent or fragment sequence |

`canonical_status=UniProtKB_displayed_sequence` means the representative sequence supplied in that UniProtKB entry. It does not mean the dominant DRG transcript, an experimentally measured sequence, or necessarily a mature protein. Explicit isoform association and entry-level unspecified association are distinguished in the output. The canonical isoform is not assumed to be accession-1. Reviewed/unreviewed status is recorded; neither is an automatic tie breaker. Multi-entry genes are held even when one candidate is reviewed. Fragment holding avoids treating an incomplete sequence as a complete protein in later Cys inventory.

These conservative implementation choices operationalize the registered ambiguous-mapping hold policy. They define a first pass, not permission to discard held genes from the research question. The held cohort must be inspected and any rescue strategy fixed and logged before implementation; candidate identities, oxidation score and functional familiarity must not decide a rescue.

## Ubuntu execution

No extra Python package or R installation is needed. Python 3.8+ standard library and HTTPS access to `rest.uniprot.org` suffice. Run in the Ubuntu terminal:

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/03_map_GSE286387.py
cat oipn-cysteine-structural-screening/results/03_protein_mapping/mapping_report.md
```

The script resolves its own project root; `--project-dir` permits an explicit project path. It retrieves 35 batches of 25 genes with up to three concurrent GET requests and follows all result pages. Typical downloads are tens of MB and can take several minutes. It retries transient HTTP/network failures; failed requests stop execution rather than becoming unmapped genes.

Responses and release headers are cached under `data/raw/UniProt/step03/<discovery_sha256>/`. A repeat run verifies cached content checksums before reuse. Every page must identify the same UniProt release. A corrupt/incomplete cache or mixed-release download stops the run; preserve the failed evidence before correcting a cache. No cache is silently refreshed. Existing cached data may belong to an earlier release than a fresh download, which is recorded explicitly.

## Outputs and review gate

Outputs are under `results/03_protein_mapping/`:

| File | Content |
|---|---|
| `gene_protein_mapping.csv` | One terminal row per input gene; candidates, mapping status/reason, accepted identity, canonical evidence, sequence/version provenance |
| `mapping_candidates.csv` | Every exact gene–entry candidate and its cross-references, including held alternatives |
| `canonical_sequences.fasta` | One displayed sequence per unique accepted accession; gene associations in header |
| `step_audit.csv` | One pass/held/unassessable status for every input gene |
| `source_manifest.csv` | Query/page URLs, UTC access times, response SHA256 and UniProt release/date |
| `input_checksums.csv` | Discovery, accepted binding, protocol config and executable script hashes |
| `summary.json`, `mapping_report.md` | Counts, reasons and limitations |
| `SUCCESS.txt` / `FAILURE.txt` | Successful generation versus failed/incomplete run |

The generated status is **STEP03_GENERATED_REVIEW_PENDING**. Success confirms complete execution, not final biological approval. Verify sequence hashes/lengths, exact IDs, gene–protein multiplicity, held reasons and release consistency before Step 04. Multiple genes may map to one accession; genes remain distinct and FASTA proteins are deduplicated by accession, never by sequence or symbol. Zero passed genes is a documented outcome, not a reason to relax identity criteria.

## Real-input implementation validation

Validation used the uploaded workstation 854-gene input and UniProt release 2026_03 (release date 02-September-2026). It found 285 passing genes corresponding to 283 unique proteins, 562 genes held for multiple mouse entries, six held for no current exact link, and one unassessable fragment. This is an implementation run; workstation reproduction remains pending. No Cys filtering or structural scoring has occurred.

The 562 ambiguous genes are a major coverage limitation of the strict first pass. Their candidate list is retained for a separate identity-resolution review; they are not biologically negative controls. The accepted subset alone does not represent every DEG-associated protein.

Tests cover exact versioned IDs, explicit canonical isoforms other than -1, reviewed/unreviewed ambiguity, wrong species, noncanonical-only associations, sequence problems, corrupt cache, missing pagination, modified discovery input, and complete real-output gene/FASTA integrity. Existing Steps 01–02 regression tests remain in the suite.

## Primary documentation

- [EMBL-EBI: Sequence and isoforms](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/exploring-a-uniprotkb-entry/the-entry-view/sequence-isoforms/): representative sequence and alternative isoforms.
- [EMBL-EBI: Programmatic UniProt access](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/getting-data-from-uniprot/accessing-uniprot-data-programmatically/): API entry/query retrieval and release provenance.
- [UniProt REST endpoint](https://rest.uniprot.org/uniprotkb/search): live queries and JSON cross-reference/sequence fields were exercised on the actual mouse input during validation.

## Implemented ambiguity-resolution substep

The user reports successful first-pass workstation execution with the same 285-gene / 283-protein totals. [Step 03b](STEP_03B.md) now resolves ambiguous entries using explicit gene-centric representative evidence and all-candidate accounting, preserving baseline outputs. Its validation expands coverage to 789 genes (92.4%); workstation reproduction is pending. Do not use the baseline-only FASTA as the expanded Step 04 input. The remaining six no-current-link genes, fragment and unresolved representative cases stay auditable.
