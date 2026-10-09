# Step 05b1 — raw coordinate acquisition and archive

Protocol 0.1.0 is unchanged. This is the acquisition substep within Step 05b; residue mapping, assembly/context review and local quality remain future implementations. Downloading a file is not a structural/site pass.

## Input and scope

Verify config/GSE286387_step05b_input.json, all accepted snapshot hashes, the recorded 42 provenance inputs and the upstream canonical binding before network access. Use the accepted 760-protein / 10,799-Cys universe. Of 1,872 offered metadata records, download only the 1,533 mouse/identity metadata candidates. The 324 held and 15 excluded records remain in candidate_file_links.csv; they are not deleted or silently resolved.

Deduplicate files by kind and URL: 445 full deposited PDB asymmetric-unit mmCIF files, 731 AlphaFold mmCIF files, 731 confidence JSON files and 731 PAE JSON files: 1,176 required coordinate files and 2,638 total requested files. All PDB chains, atoms, ligands and assembly-operation annotations present in each deposited file are retained. No biological assembly is generated/selected in this substep. Deposited local validation reports and appropriate generated/downloaded assemblies remain later acquisition/context work before local eligibility can be approved. Legacy PDB-format files are not duplicated alongside mmCIF. AlphaFold URLs/versions come from the frozen candidate table; they are not guessed or upgraded.

## Ubuntu commands

Uses Python standard library only; no additional R/Python packages or administrative installation is required.

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/05b1_download_GSE286387.py --probe &&
python3 oipn-cysteine-structural-screening/scripts/05b1_download_GSE286387.py --workers 4 --archive
cat oipn-cysteine-structural-screening/results/05b1_structure_download/download_report.md
```

Probe checks one file of every planned kind (four files). It produces PROBE_SUCCESS.txt only if all probe transports work and both coordinate samples are available. It never produces full-download SUCCESS.txt. A failing probe prevents bulk requests. Read failure_context.json for exact affected URLs/statuses. The optional --plan-only writes the complete plan without acquisition; --offline requires existing verified caches. Do not execute multiple instances concurrently against the same cache/output directory.

Full acquisition checks all files, preserves the complete plan and reports failures without treating them as biological absence. HTTP 404 is explicit unavailable evidence; missing required coordinates prevent completion. Missing optional confidence/PAE remains a recorded gap even when all coordinates succeed; no confidence/context gate is passed on that basis. Network/SSL/HTTP 403 errors are failed, not unavailable. Retry transient errors at most four times, 120-second request timeout (revision 1.0.1 defaults); TLS verification remains enabled. Unexpected redirected hosts, payload identity, coordinate-tag/content-length failures or corrupted caches fail acquisition. mmCIF checks here are transport/content sanity checks, not full parser/sequence/atom validation.

Successful completion requires all 1,176 coordinate files plus no failed requests and SUCCESS.txt. Status is STEP05B1_GENERATED_REVIEW_PENDING. The report retains auxiliary 404 gaps. Incomplete runs use FAILURE.txt; rerun the same full command to reuse byte-verified successful files and retry failed/unavailable endpoints. An orphan body without metadata is never trusted; confirmed corrupted cache pairs are reported rather than silently replaced. Original retrieval dates/ETag/Last-Modified and raw SHA256 are preserved on reuse.

## Outputs and storage

- results/05b1_structure_download/: full plan, all candidate/file associations, download manifest, input checksums, status summary/report, failure URL context and success/failure marker.
- data/raw/structure_coordinates/step05b1/<plan SHA256>/: unchanged file bytes and acquisition metadata. Keys derive from kind/URL, avoiding unsafe or colliding server filenames.
- data/archives/step05b1_<plan SHA256 prefix>.tar.gz: raw files with sidecars, all frozen input files/executable and output audit tables. archive_manifest.csv binds every archived source file; the script reads back every tar member and verifies SHA256 before replacing the archive. archive_receipt.json records the final archive path/hash/size and remains outside that archive to avoid self-reference.

An archive can capture a partial run: its included status/manifest says incomplete; archive existence does not mean completion. Rebuilding after resumed acquisition replaces that plan's archive only after verification. Keep enough local disk space for raw files plus their compressed archive; the exact total depends on deposited structures/PAE sizes. Back up the final archive separately. data/ remains Git-ignored; commit compact results/audits to GitHub, not bulk coordinate files or the tar.gz. The actual workstation receipt is the archive's hash reference.

## Official format source

https://www.rcsb.org/docs/programmatic-access/file-download-services documents mmCIF entry downloads at https://files.rcsb.org/download/<PDB ID>.cif and distinguishes biological-assembly files. Predicted coordinate/auxiliary URLs are the exact AlphaFold API URLs preserved in the accepted Step 05a reference.

## Validation and limits

The complete frozen plan contains exactly 2,638 unique files and links all 1,872 offered records. Real online probe downloads succeeded for one file of each kind (four), then cached replay was tested. All 56 then-existing regression tests passed; the additional full acquisition/archive/offline-resume test also passed (nine acquisition tests total). This is implementation validation, not the user's bulk download. Bulk live acquisition has not been performed in the validation environment. Workstation HTTPS transport remains independently unverified; use the probe first.

## Interrupted workstation acquisition — revision 1.0.1

Actual first bulk run: 2,121 successful files / 517 failed. All AlphaFold coordinates and confidence files succeeded; seven PDB coordinates and 510 PAE files remain. See audits/ubuntu_download_incomplete_2026-10-09/REVIEW.md. HTTP 400 cause is unresolved because the old executable did not retain response-body evidence. Two unchanged failed PAE URLs work from the validation environment; do not alter their model versions or substitute another endpoint without provenance.

The revised executable defaults to timeout=120 seconds and four attempts for transient errors, including IncompleteRead and remote disconnection. --timeout (1–600) and --attempts (1–8) are explicit CLI controls. --probe --retry-probe reads the prior download manifest and selects one previous failed/unavailable URL per file kind, verifying its ID/URL against the frozen plan; otherwise it uses the first planned URL. HTTP errors now include status, response URL, an allowlisted header JSON string and at most 2,048 response bytes in the failure context/CSV. Cookies are omitted. HTTP 400 is not blindly retried or relabeled unavailable.

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/05b1_download_GSE286387.py --probe --retry-probe --timeout 120 --attempts 4 &&
python3 oipn-cysteine-structural-screening/scripts/05b1_download_GSE286387.py --workers 1 --timeout 120 --attempts 4 --archive
```

If the failed-URL probe still fails, inspect failure_context.json; the bulk command is not run by this conditional shell sequence. Original successful caches/metadata are preserved; a completed archive replaces the previous plan archive only after read-back checksum verification. Partial archive reports remain explicitly incomplete. Validation: all 60 tests passed plus full JSON downloads of two previously failed PAE URLs from the validation environment. Workstation recovery is not claimed.

## Actual Ubuntu completion accepted

Source b2205ac65ebf0881653dd985ee279aef2f9fa994 completes all 2,638 planned files and all 1,176 coordinates. Full input, frozen-plan, acquisition-history and archive-manifest reconciliation passes; see audits/ubuntu_download_2026-10-09/REVIEW.md. The receipt binds the local 447,690,451-byte tar.gz. config/GSE286387_step05b2_input.json freezes complete actual manifests and requires local byte-hash verification before mapping. The archive/raw bodies were not remotely inspected; residue mapping/context/quality implementations remain pending.
