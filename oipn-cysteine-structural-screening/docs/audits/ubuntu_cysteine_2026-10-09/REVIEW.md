# Ubuntu Step 04 review — accepted for Step 05

Source commit: e419a4e98e8b51368a0d5123b08fc128415eee1c; recorded 2026-10-09 (Asia/Tehran). This snapshot is the user's actual workstation run: Python 3.13.13, Linux-7.0.0-34-generic-x86_64-with-glibc2.39, generated UTC 2026-10-09T07:02:11.198151+00:00.

All nine snapshot SHA256 records and all 27 upstream/executing-script checksum records pass. The five complete protein/site/FASTA/protein-audit/gene-audit files are byte-identical to independent enumeration. Remote tracked results/04_cysteine_inventory files have the same Git blob identities as the corresponding uploaded snapshot, including environment/summary provenance. No validation-environment file is relabeled as a workstation run.

Accepted sequence inventory: 787 input accessions, 760 C-positive and 27 no-C accessions; 31 have one C, 729 have multiple C; 10,799 unique canonical sites. Full gene audit: 854 genes = 760 C-positive gene mappings + 29 no-C gene mappings + 64 held mappings + one unassessable mapping. Gene/accession counts differ because some genes share a no-C accession; no merging of stable genes or duplication of sites occurs.

Cys positions retain full canonical precursor numbering, including U/O without counting U as C. The 27 excluded proteins are excluded only from the Cys-based method. The 65 upstream mapping cases remain unresolved and are not treated as Cys-free.

Freeze the C-positive sequence/site inventory in config/GSE286387_step05_input.json. This acceptance does not establish free thiol state, structure availability, local quality, oxidation or susceptibility. Step 05a catalogues metadata candidates; coordinate retrieval, sequence-to-residue mapping, assembly/mutation review and local QC remain Step 05b work. Experimental structures retain priority only when suitable at the individual site. Historical configuration and earlier scripts remain unchanged.
