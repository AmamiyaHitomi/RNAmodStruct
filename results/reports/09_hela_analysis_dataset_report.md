# Stage 09 HeLa analysis dataset report

## Inputs and reference

- HeLa m6A: GLORI `GSM6432595_Hela-1` (137,662 sites) and `GSM6432596_Hela-2` (140,705 sites);
  exact `(Chr, Sites, Strand)` union of 164,761 sites.
- HeLa structure: `GSE145805_RAW.tar` member `GSM4333258_HeLa.out.txt.gz` (64,616 versioned Ensembl
  transcripts, hg38), read directly from the archive without materializing the immutable download area.
- Reference: Ensembl release 88 / GRCh38 GTF and cDNA, matched by `transcript_id.transcript_version`.
- 60,350 transcripts are eligible with exact declared-length agreement across icSHAPE, GTF exon sums,
  and cDNA. Excluded: 3,040 without a GTF transcript and 1,226 without a cDNA sequence in release 88.
- HeLa is on GRCh38, so GLORI sites map directly onto Ensembl 88 exons (no hg19→hg38 liftOver).

## Attrition

| Step | Stage | Sites | Retained |
|---|---:|---:|
| 1 | GLORI union | 164,761 | — |
| 2 | Matched to HeLa icSHAPE transcript | 86,932 | 52.8% |
| 3 | Detected in both GLORI replicates | 65,050 | 74.8% |
| 4 | Main structure dataset (≥70% valid ±10 nt) | 25,996 | 40.0% |
| 5 | Model-ready full 201-nt sequence | 24,960 | 96.0% |

## QC

- Main dataset: 25,996 sites from 4,829 analysis genes; 57,826 of 86,932 mapped sites have multiple
  isoforms (maximum 59); 253,714 total isoform mappings retained.
- GRCh38 center-reference base failures: 0 (expected A on `+`, T on `-`).
- Replicate agreement (main dataset): Pearson 0.9836, Spearman 0.9707, median absolute ratio difference
  0.0228.
- The HeLa icSHAPE abundance field is the literal `*`, so `icshape_abundance_rpkm` is missing for every
  site (recorded deviation; the association stage excludes that covariate, the transfer stage imputes it
  with the frozen HEK293T median).

## Outputs

`data/final/09_hela_site_level_dataset.csv.gz`, `09_hela_main_analysis_dataset.csv.gz`,
`09_hela_all_isoform_mappings.csv.gz`; summaries under `results/09_hela_*`.
