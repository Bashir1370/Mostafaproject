# Step 05a — structural candidate catalogue

Protocol 0.1.0; implementation 1.0.1 (diagnostic revision). This is the metadata-catalogue substep of scientific Step 05. It does not add an eleventh scientific step or replace site-level structure review.

## Reviewed input

The user's Ubuntu Step 04 snapshot is reviewed at source commit e419a4e98e8b51368a0d5123b08fc128415eee1c. All nine snapshot SHA256 records and all 27 upstream/script checksum records pass. Complete site/protein/FASTA/audit tables reproduce independent enumeration. See docs/audits/ubuntu_cysteine_2026-10-09/REVIEW.md. config/GSE286387_step05_input.json binds the accepted 760 C-positive proteins and 10,799 canonical Cys positions. The script also validates prior sequence provenance and compares every site with the frozen canonical sequence.

No-Cys accessions and unresolved mappings remain in the prior full audits; they are not fabricated structural negatives. The historically frozen analysis_parameters.json, earlier scripts and scientific thresholds remain unchanged.

## Queries and fixed criteria

For every accepted accession, query both services, even if an experimental candidate exists:

1. PDBe SIFTS best_structures API: retain every reported PDB/chain/interval record, independent of API ordering, coverage or resolution ranking. Mouse taxonomy 10090 is required for a primary candidate. Missing taxonomy is held; nonmouse chains are excluded from the primary candidate set. Intervals outside the frozen sequence are held. The coarse start/end coordinates are not residue alignment or proof of observed Cys/SG atoms. This API catalogue is not an exhaustive search of unannotated PDB structures.
2. AlphaFold DB prediction API: retain every returned model, version, chain, sequence interval and service-provided file URL. Require exact accession, mouse taxonomy and agreement with the frozen full canonical sequence or exact represented subsequence. Explicit alternate isoforms, sequence differences, out-of-range intervals and complexes are held for separate mapping review. Partial-model boundaries are flagged. No average pLDDT, entry review status or eventual score selects a model.

Use current fields modelEntityId, sequence, sequenceStart/End, with explicit compatible legacy fallback. If both forms disagree, stop on a schema conflict. Download URLs come from the API; no hard-coded model-v4/v6 URL is constructed. A metadata sequence match still requires coordinate validation later. Reported sequence checksums are retained as metadata; internal sequence SHA256 is computed independently.

Candidate status is candidate/held/excluded, while structural_eligibility is always not_assessed. Protein audit status is held (candidate(s) require local review, or offered records have unresolved identity) or unassessable (no eligible record reported and no held record to resolve). No protein/site receives a structural pass in this substep. Experimental structures retain priority only after local suitability is established in Step 05b.

## Runtime and provenance

Python >=3.8 standard library; no R or new package required. Default eight concurrent protein queries, configurable 1–16 via --workers. HTTPS access to www.ebi.ac.uk and alphafold.ebi.ac.uk is required on first run. Each service response is cached under data/raw/structure_catalogue/step05a/<accepted-FASTA-SHA256> with raw response SHA256, HTTP status, URL, access UTC and service release header (explicitly absent where not supplied). A completed cache is replayed without refreshing; --offline requires the full cache. Partial or corrupt cache pairs stop execution rather than silently refreshing.

A 404 or empty valid response means that service did not report data at retrieval time; it is not biological absence. HTTP 429/5xx and network timeouts are retried up to four attempts; other errors/schema failures stop the run, with no success marker. Network failures never become no-structure counts. Retries reuse previously completed cache entries. Success requires both service queries for every input accession; 1,520 response provenance rows are expected. Failure cancels queued jobs and preserves completed cache data.

## Ubuntu execution

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/05a_catalogue_GSE286387.py
cat oipn-cysteine-structural-screening/results/05a_structure_catalogue/catalogue_report.md
```

Later cache-only replay:

```bash
python3 oipn-cysteine-structural-screening/scripts/05a_catalogue_GSE286387.py --offline
```

## Output contract

| File | Content |
|---|---|
| structure_candidates.csv | Every deduplicated service candidate record: accession, source, structure/chain/version, taxonomy, reported interval, sequence-check result, source URLs/hash, metadata decision/reason, structural_eligibility=not_assessed. |
| protein_structure_inventory.csv | One record per 760 accessions, with both service record counts, mouse/sequence-matched candidate counts and pending/unassessable status. |
| source_manifest.csv | Both response provenance records for every accession, including cached 404s. |
| step_audit.csv | One provisional Step 05 protein status per input accession; no site inclusion approval. |
| input_checksums.csv | Reviewed binding, snapshots, upstream sequence inputs and executed scripts. |
| summary.json; catalogue_report.md | Separate source/candidate/protein counts, model versions and environment. |
| SUCCESS.txt / FAILURE.txt | STEP05A_GENERATED_REVIEW_PENDING requires no FAILURE.txt. |

All outputs are under results/05a_structure_catalogue, now trackable in Git. Add/commit/push actual Ubuntu output after reproduction; validation outputs are recorded separately. Tables remaining after failure must not be consumed.

## Next within Step 05

Step 05b remains to be implemented: coordinate/assembly and local-validation retrieval; exact canonical-to-structure mapping; mutation/engineered-site review; target SG and required atom completeness; alternate occupancy handling; local pLDDT >=90 at target and >=70 at all other modeled residues within 6 A of SG; boundary and PAE/domain/assembly context review. Incomplete experimental evidence remains held. Catalogue counts must not be presented as quality-eligible structures, covered Cys sites or scored proteins.

## Primary documentation

- PDBe APIs: https://www.ebi.ac.uk/pdbe/pdbe-rest-api
- SIFTS API documentation: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
- AlphaFold API field changes: https://www.ebi.ac.uk/pdbe/news/breaking-changes-afdb-predictions-api
- AlphaFold versions/download-URL guidance: https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes

## Real-input implementation validation

The complete 1,520-query run and final offline replay cover all 760 proteins. Results: 94 with mouse PDBe candidates, 731 with exact-full-sequence AlphaFold candidates, 94 with both, 29 without a metadata-eligible candidate. Candidate rows total 1,872 (1,533 metadata candidates, 324 held, 15 nonmouse-chain exclusions). All 10,799 Cys sites remain prospective input; none is structurally approved here. The response/model versions and input/script hashes are recorded under docs/audits/Structure_catalogue_validation. Actual Ubuntu Step 05a reproduction remains pending.

Validation: all 42 regression tests passed, including eleven structure-catalogue identity/cache/error/real-output checks and all earlier input/QC/DE/mapping/inventory tests.

## HTTP access diagnostics — revision 1.0.1

The user reported HTTP 403 on Ubuntu before a catalogue was generated. The initial exception text omitted the failing endpoint; the report does not establish whether PDBe, AlphaFold, a proxy or another network component rejected the request. The same initial service queries succeed (expected PDBe 404 / AlphaFold 200) from the validation environment, so workstation access failure is not independently reproduced.

Fatal HTTP errors now include the exact requested URL and write failure_context.json beside FAILURE.txt, with status, response URL, allowlisted diagnostic headers and a bounded 1,024-byte response excerpt. Cookies and unrelated headers are not collected. Error responses are never added to the accepted metadata cache or converted into missing-structure results. The request endpoints, headers, identity gates and success-output tables remain unchanged. Existing historical validation files retain their original executable checksum/version.

After pulling the revision, rerun the normal command and provide failure_context.json if the request fails. A missing catalogue_report.md after failure is a consequence of the incomplete run. No network workaround or root cause is claimed until the endpoint/response evidence is available.

Diagnostic revision validation: all 44 regression tests pass; final cache replay reproduces all four complete catalogue/query/protein-audit tables byte-for-byte. Workstation 403 cause and resolution remain pending.

## Frozen reference mode — revision 1.0.2

The second Ubuntu attempt failed during HTTPS transport with SSL UNEXPECTED_EOF_WHILE_READING. The failing endpoint/root cause is still unknown; TLS verification is unchanged. To reproduce the completed validation catalogue without depending on workstation API connectivity, revision 1.0.2 supports a bundled frozen public response reference:

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/05a_catalogue_GSE286387.py --reference-cache
cat oipn-cysteine-structural-screening/results/05a_structure_catalogue/catalogue_report.md
```

The 629,151-byte archive contains exactly 1,520 query response/metadata pairs captured on 2026-10-09 UTC. Its SHA256 is f486decef5424cfe215a523d0df3e64dfe0171646401d573f484985f7743a2e3. The descriptor binds the archive, original validation source ledger, and accepted canonical FASTA. The loader verifies the complete URL universe, member names/types, response checksums, HTTP statuses and original retrieval times before importing. It never uses unrestricted tar extraction. An existing altered reference cache is rejected.

Reference bytes live in a separate data/raw/structure_catalogue/reference/<archive SHA256> namespace; partial workstation/network caches remain separate. This option forces offline execution, including when --offline is omitted. summary.json and the report explicitly identify frozen_reference metadata and the archive checksum; input_checksums.csv includes the descriptor/archive/source ledger. source_manifest.csv retains the original upstream acquisition timestamps. These are reused public reference responses, not new downloads on the workstation. A later database update requires a separately documented reference revision.

This resolves reproducibility of Step 05a only. It does not fix the workstation HTTPS connection, download coordinates, approve structures, or score oxidative susceptibility. The existing structural/identity gates and frozen scientific protocol are unchanged. Successful output still requires SUCCESS.txt and review before Step 05b.

## Actual Ubuntu reproduction accepted

The uploaded Python 3.13.13 execution at f8a1a4e8c60593101425ebc09ba55c2081393bba passes the full input and output review. See audits/ubuntu_catalogue_2026-10-09/REVIEW.md. Its complete immutable snapshot and checksum manifest are bound in ../config/GSE286387_step05b_input.json. Metadata-only acceptance does not change structural inclusion gates or claim fresh API acquisition.
