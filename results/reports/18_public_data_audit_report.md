# Stage 18: The second phase of public data and coordinate audit

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: 2026-09-11
- Verification Status: VERIFIED
- Version Label: phase2_stage18_v1

This stage performs local container, row and column number, and coordinate contract audits on the second phase candidate sources. 3 local files have been fully verified; the number of data sets that can directly enter strict pair analysis is 1.

- GSE74353 HEK293T in-vitro icSHAPE: admission stage 19, paired with existing in-vivo files in the same GRCh37.74 transcript coordinate system.
- GSE162356 HeLa SAC-seq: Native seven-column BED with coordinate base contract passed; final semantics for numeric columns continue to retain source fields until stage 20 and are not renamed to stoichiometry without confirmation.
- GSE162356 HEK293: Only close-match sensitivity analysis of HEK293T and cannot be marked as strict cell line replication.
- Human RNA MaP with GSE50676: currently missing quantitative m6A matching cell conditions, remains in audit status, and does not splice raw scores.

Machine table:`results/tables/18_public_data_audit.csv`;Independent source list:`metadata/manifests/source/phase2_source_manifest.csv`。
