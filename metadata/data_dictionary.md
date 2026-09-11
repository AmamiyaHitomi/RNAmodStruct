# Data dictionary (W1 draft)

Checked: 2026-09-10. This dictionary describes observed file structure. Coordinate origins and biological
interpretation that cannot be proven from the file itself remain explicitly unresolved.

## HEK293T icSHAPE processed reactivity

File: `GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz`

- Format: gzip-compressed, headerless, tab-delimited text.
- Rows: 10,164 transcripts.
- Reference stated by GEO: Ensembl GRCh37 release 74 / hg19 transcriptome.
- Transcript IDs are unversioned in this file, for example `ENST00000539264`.
- Every row passed the invariant `number of score fields == declared transcript length`.
- Declared bases: 13,671,469; valid scores: 5,101,164; `NULL`: 8,570,305.
- Observed valid-score fraction: 37.31%; observed score range: 0 to 1.

| Position | Observed meaning |
|---:|---|
| 1 | Ensembl transcript ID |
| 2 | Transcript length in nucleotides |
| 3 | RPKM |
| 4 onward | One icSHAPE reactivity value per transcript nucleotide; `NULL` means unmeasured |

## HeLa icSHAPE processed reactivity

Container: `GSE145805_RAW.tar`  
Required member: `GSM4333258_HeLa.out.txt.gz`

- The TAR contains six final `.out.txt.gz` members: HEK293, HeLa, K562, HepG2, H9, and mES.
- Format of the HeLa member: gzip-compressed, headerless, tab-delimited text.
- Rows: 64,616 transcripts.
- Reference stated by GEO: hg38.
- All transcript IDs are versioned in the observed file, for example `ENST00000606214.5`.
- Every row passed the invariant `number of score fields == declared transcript length`.
- Declared bases: 111,656,194; valid scores: 42,998,767; `NULL`: 68,657,427.
- Observed valid-score fraction: 38.51%; observed score range: 0 to 1.

| Position | Observed meaning |
|---:|---|
| 1 | Versioned Ensembl transcript ID |
| 2 | Transcript length in nucleotides |
| 3 | Literal `*` in observed rows; biological meaning not applicable/unspecified in submitted format |
| 4 onward | One icSHAPE reactivity value per transcript nucleotide; `NULL` means no confident score |

The submitted DMSO sample `GSM4333257` and NAI-N3 sample `GSM4333258` are processing inputs; the
archive exposes the final HeLa output under `GSM4333258_HeLa.out.txt.gz`. They must not be averaged as
two independent structure replicates.

## GLORI processed site tables

Files:

- `GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv.gz`: 214,715 rows.
- `GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv.gz`: 215,021 rows.
- `GSM6432595_Hela-1_35bp_m2.totalm6A.FDR.csv.gz`: 137,662 rows.
- `GSM6432596_Hela-2_35bp_m2.totalm6A.FDR.csv.gz`: 140,705 rows.

Despite the `.csv.gz` suffix, all four files are gzip-compressed **tab-delimited** tables. They have the
same 14 columns, both strands are present, and the initial full streaming parse found no malformed rows.
The human reference stated by GEO is GRCh38.

| Column | Observed or provisional meaning | Status |
|---|---|---|
| `Chr` | Chromosome label such as `chr10` | observed |
| `Sites` | 1-based genomic coordinate of the candidate modified A itself; genomic base is A on `+` and T on `-` | confirmed from GLORI-tools coordinate and strand/base code |
| `Strand` | `+` or `-` | observed |
| `Gene` | Submitted gene symbol/label; `ELSE` occurs | observed; mapping rules to audit |
| `Transcript` | Submitted transcript accession/label; RefSeq-style IDs and `ELSE` occur | observed; version mapping to audit |
| `NonCR` | Gene-level non-conversion background: aggregated retained A counts divided by aggregated A+G counts across qualifying A loci (`AG >= 15`) | proportion; confirmed from code |
| `AGcov` | Retained A-or-G read count at the site after pipeline filters | integer effective read/molecule count |
| `Acov` | Retained A read count at the site | integer effective read/molecule count |
| `Genecov` | Mean A+G depth across the gene's qualifying A loci: aggregated A+G count divided by number of those loci | reads per qualifying A locus |
| `Ratio` | Site A rate, `round(Acov / AGcov, 5)`; primary downstream m6A proportion | confirmed in all 708,103 rows |
| `NormeRatio` | Conversion-scaled site ratio, `Ratio * (1 - NonCR)`; sensitivity variable only | exact identity in all 708,103 rows |
| `Pvalue` | Submitted site-level P value | observed |
| `P_adjust` | Submitted adjusted P value | adjustment family/method to confirm |
| `Sample` | Constant `m6A.sites` in all four audited files | observed |

## Coordinate/reference decisions for W2

1. `Sites` is confirmed as 1-based and denotes the transcript-oriented modified A itself.
2. `Ratio` is the primary m6A response; `NormeRatio` is retained for sensitivity analysis.
3. HEK293T is pinned to Ensembl release 74 / GRCh37 as declared by GEO.
4. HeLa is pinned operationally to Ensembl release 88 / GRCh38 because releases 88 and 89 tie for the
   best tested ID-plus-exon-length match. GEO does not state the exact Ensembl release, so this remains a
   documented inference rather than an author-declared fact.
5. GLORI reproduction assets are pinned to RefSeq release 109.20190905 / GRCh38.p13 as specified by
   GLORI-tools. Transcript-aware GRCh37 mapping followed by recorded liftOver is still required for HEK293T.

See `metadata/01_glori_field_semantics.md`, `metadata/02_reference_manifest.csv`, and
`metadata/02_reference_compatibility.csv` for evidence and hashes.
