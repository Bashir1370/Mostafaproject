# Ubuntu Step 05a snapshot review

Status: STEP05A_SNAPSHOT_ACCEPTED_FOR_STEP05B

Source commit: f8a1a4e8c60593101425ebc09ba55c2081393bba. Actual workstation: Python 3.13.13 / Ubuntu Linux. Implementation 1.0.2; scientific protocol 0.1.0 unchanged.

All eight uploaded files are preserved byte-for-byte. All 42 recorded input checksums match the accepted inputs, executable and frozen public reference. The four complete candidate/protein/query/audit tables match the independent offline reference reproduction byte-for-byte. Every accepted accession and both service queries are accounted for; no structural pass was issued. A checksum manifest and config/GSE286387_step05b_input.json bind this reviewed snapshot separately from mutable results/.

Counts: 760 proteins / 10,799 prospective Cys sites; 1,872 offered records (1,533 metadata candidates, 324 held, 15 nonmouse exclusions). 731 proteins have sequence-matching AlphaFold options; 94 also have mouse PDBe options. 29 lack a metadata-eligible candidate: 27 are unassessable and two retain unresolved offered records and are held. Overall 733 protein records are held pending structural work; 27 are unassessable at the metadata stage.

Provenance: frozen_reference metadata, original acquisition on 2026-10-09 UTC, archive SHA256 f486decef5424cfe215a523d0df3e64dfe0171646401d573f484985f7743a2e3. This workstation execution made no fresh API queries. Source timestamps/releases are upstream reference provenance. HTTP/SSL transport issues remain unresolved and may affect future coordinate retrieval.

Acceptance is limited to using the audited catalogue as Step 05b input. It is not approval of structures, local Cys quality, chemical states or oxidative susceptibility. Keep all 760 proteins and 10,799 sites in downstream accounting; do not turn missing metadata into zero susceptibility. Retain every offered candidate and resolve structural/context eligibility before choosing a conformation. No coordinates, SASA, pKa or structural scores have been computed.
