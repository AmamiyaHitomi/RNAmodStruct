# Stage 24: Multi-modified files, coordinates and structure coverage audit

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: 2026-09-12
- Verification Status: VERIFIED
- Version Label: phase3_stage24_v1

5/5 downloaded files passed the fixed SHA-256 check and were parsed successfully. The three modifications m5C, m7G and Nm have enough records to fall into the structural exons of the corresponding cell lines and the transcript lengths are exactly the same. The stage gate is`GO_TO_STRUCTURE_FEATURE_EXTRACTION`。

## Coordinates and measurement conclusions

- **m5C/HeLa**: 1,034 unique hg38 single-base sites; 100.00% reference base coincidence; 688 sites fall into transcript exons with HeLa icSHAPE data. The three non-transformed scales B/C/E were retained as independent replicates.
- **m7G/HeLa, HEK293T**: 6,022 and 6,290 hg38 enrichment intervals, respectively; they are not single base sites, nor stoichiometric, and are only allowed for interval-level sensitivity analysis. Structural exon overlap was 549 and 132 intervals, respectively.
- **Nm/HeLa, HEK**: 2,103 and 699 sites of the current GEO revision file are used directly by hg38, with structural exon overlap of 1,692 and 522 respectively. The original paper method describes hg19, but the revised file has hg38 characteristics: some coordinates exceed the corresponding hg19 chromosome length, and the HeLa direct hg38 overlap is much higher than the overlap after the wrong liftOver. Since the file was submitted with`start=end`Indicates the point position, and the ±1 nt coordinate sensitivity analysis must still be retained in the next stage; the HEK tag can only be used as an approximate match for HEK293T.

## Explain boundaries

“Structural overlap” here means that the modification record falls into a transcript exon with corresponding cell line structural data, and does not mean that there is already a qualified ±10 nt icSHAPE window at that position. Window coverage, representative transcript selection, and final sample size will be determined after site-by-site extraction at stage 25. Raw values ​​for the three techniques were not combined.
