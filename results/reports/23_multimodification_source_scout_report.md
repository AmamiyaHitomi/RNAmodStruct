# Stage 23: The third phase of multi-modification public data reconnaissance

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: deep-research/source-scout
- Origin Date: 2026-09-12
- Verification Status: SOURCE_METADATA_VERIFIED
- Version Label: phase3_stage23_v1

In this phase, 7 candidate modifications were checked, forming a download queue of 5 first batches of small processed files. The modifications that meet the download preparation conditions are Nm, m5C, m7G, and the stage threshold is`GO_TO_DOWNLOAD_AUDIT`. This only authorizes access to post-download field and coordinate auditing, which does not mean that cross-modification modeling can already be carried out.

## First batch priority

1. **m5C/GSE140995**: HeLa bsRNA-seq, post-processing workbook containing candidate sites and non-conversion ratios, prioritized as continuous or semi-continuous track audit.
2. **m7G / GSE112276**: Both HeLa and HEK293T have HighFC files, but it is enrichment/fold change semantics and cannot be used as a stoichiometric ratio.
3. **Nm/GSE90164**: Revised compact BED site table for HeLa and HEK, treated by binary site orbitals.

## Suspension source

- m1A: GSE97909 proposes a cap structure cross-reactivity explanation for a broad range of antibody peaks, not suitable as a general mRNA m1A layer.
- Ψ: GSE255287 matches HEK293T, but the three processed matrices total about 12.5 GB, and the site extraction process needs to be designed first.
- ac4C: GSE162043 matches HeLa, but the mRNA ac4C detection remains methodologically controversial and must be reanalyzed using untreated versus NAT10-knockout controls.
- A-to-I: REDIportal v3 is primarily from GTEx/TCGA and does not strictly match existing HeLa/HEK293T culture conditions; engineered HEK293T edited data is more suitable for the Phase 4 perturbation route.

The next mandatory threshold is to check the reference version file by file after downloading, starting from 0/1, link direction, field semantics, repeat consistency and actual overlap rate with the existing structure table. Only after at least two modifications pass this threshold can the unified feature space be launched.
