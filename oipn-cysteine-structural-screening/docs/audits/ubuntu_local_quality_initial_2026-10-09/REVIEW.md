# Actual Ubuntu Step 05b3 initial run — correction required

Source commit: 7bb88bf0b17542d6218454c0c098f1b75fd476db. Review date: 2026-10-09. Implementation 1.0.0.

All nine uploaded Git blob hashes and all 80 input/executable hashes were verified before the executable amendment. The frozen input chain, 14,310 unique mapped-option associations, all 731 per-model summaries, complete 760-protein/10,799-site accounting and all 2,638 raw/sidecar ledger hashes pass. The uploaded summary and per-option gate decisions agree with the frozen thresholds. The raw coordinate corpus remains local to Ubuntu, not independently parsed remotely in full.

Valid option results: 4,580 predicted local-pass sites in 667 proteins; 5,036 options fail target pLDDT, 125 fail neighboring pLDDT; 4,569 experimental mapped alternatives await manual validation/context. Every predicted PAE matrix parsed descriptively; all 9,741 predicted mapped-site options are full canonical models. No final structural approvals or scores exist.

Defect: the site aggregation received only mapped options and a prior summary reason prioritizing mapped evidence. It consequently lost pending original mapping/metadata evidence when a predicted mapped option failed. 44 initially excluded sites have held upstream residue-mapping options; 2,034 have unresolved offered metadata, with union 2,042 sites. These must remain held, not terminally excluded, while evidence is unresolved.

Revision 1.0.1 passes explicit per-site original mapping-review counts and metadata-review flags into the terminal site audit. Expected corrected site counts are 4,580 local-pass/context-pending sites, 2,384 other held sites, 3,042 excluded and 793 unassessable (6,964 held total). The 4,580/667 biological candidate counts and all per-option confidence decisions are unchanged. Workstation rerun/review is required before a downstream input binding. This initial snapshot is preserved as historical provenance and is not accepted as final Step 05b3 site audit.
