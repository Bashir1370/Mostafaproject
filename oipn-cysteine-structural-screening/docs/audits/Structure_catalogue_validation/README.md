# Step 05a implementation validation

The reviewed user's Ubuntu Step 04 Cys-positive snapshot is the input. These outputs are implementation validation, not the user's Step 05a workstation run. All 760 proteins and both services (1,520 completed queries) were queried live, followed by checksum-verified offline replay with the final script. Service response SHA256/access times are preserved in source_manifest.csv; raw metadata cache remains under ignored data/raw. No coordinate retrieval or structural eligibility pass occurs here. Canonical sequence release is independent of reported AlphaFold model versions.

Catalogue: 94 proteins with mouse PDBe/SIFTS candidates, 731 with exact-full-sequence AlphaFold candidates, 94 with both, 29 without a metadata-eligible candidate. The 1,872 offered records include held alternate isoforms/taxonomy/sequence cases and excluded nonmouse chains; raw record counts are not counts of eligible proteins or sites. The all-protein audit retains unresolved metadata cases as held separately from absent eligible records.

Validation: all 42 regression tests passed, including eleven structure-catalogue identity/cache/error/real-output checks and all earlier input/QC/DE/mapping/inventory tests.
