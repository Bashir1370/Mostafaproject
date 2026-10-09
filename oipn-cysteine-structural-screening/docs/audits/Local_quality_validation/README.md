# Local-confidence implementation validation

90 full regression tests passed in 40.775 seconds; the final 16 local-quality tests passed again in 2.568 seconds after numeric-schema and audit refinements. Real accepted-hash AlphaFold model/confidence/PAE sample AF-A0A087WRH0-F1 yields 21 local-pass / 14 target-confidence-fail options among 35 sites. All neighbor sets and local gate outcomes match independent brute-force all-atom distances. Parsed PAE is descriptive and has no threshold.

The end-to-end representative fixture includes its three raw files plus 11gl: 43 mapped options (35 predicted / eight experimental); all 760 proteins and 10,799 sites remain audited, with experimental options held. A corrupted raw file removes SUCCESS and creates FAILURE. Full 731-model execution is pending on Ubuntu; local-pass samples do not imply final structural/chemical approval. See STEP_05B3.md.
