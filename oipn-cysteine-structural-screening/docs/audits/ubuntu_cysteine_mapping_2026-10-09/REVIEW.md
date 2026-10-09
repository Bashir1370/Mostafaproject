# Actual Ubuntu Step 05b2 review

Source commit: 7dd9d5dde4e2a887a3af46f5c6e0c7b90889641c. Review date: 2026-10-09.

All nine uploaded files were fetched at this immutable commit and their Git blob hashes verified. Every executable/input checksum matches the frozen chain. The 25,026 mapping rows reproduce the summary status/reason counts; all 1,872 offered candidates and all 760 proteins / 10,799 sites remain accounted for. Reconstructed site audits exactly match the uploaded CSV. All 2,638 raw/sidecar ledger hashes agree with the accepted acquisition/archive manifests.

9,741 unique sites have mapped SG options, all with AlphaFold options. 394 of those sites additionally have experimental mapped options. Record counts are 9,741 AlphaFold mapped and 4,569 experimental mapped; experimental options also have 10,028 unassessable and 688 held rows. Among 1,058 sites without any mapped option, 265 remain mapping-evidence unresolved and 793 lack metadata-eligible structures. All 14,310 mapped records report no absent CYS heavy atom names; this alone does not establish a coherent conformer.

The actual Ubuntu run reports Python 3.13.13 / Gemmi 0.7.5 and successful parsing of all 1,176 coordinate files after the local raw hash gate. Raw coordinates and the archive remain on the workstation and were not independently parsed remotely in full. This review accepts mapping/provenance for local-confidence assessment, not local quality, chemical state, biological assembly or final structural inclusion.
