# Candidate data source verification

Checked: 2026-09-10 (Asia/Shanghai)  
Mode: source fact-check against authoritative repository records and the associated primary papers.

## Overall assessment

Six candidate input records were reviewed: one HEK293T icSHAPE processed file, two HEK293T GLORI
replicates, one HeLa icSHAPE archive/pair, and two HeLa GLORI replicates. Five are confirmed at the
individual processed-file level. The HeLa icSHAPE source is confirmed at the GEO Series/archive level,
but its internal member filenames and final-score readiness remain to be audited after obtaining the TAR.

## Verified sources

### HEK293T experimental structure

- GEO Series: [GSE74353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74353).
- Candidate processed file: `GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz`.
- Relevant inputs listed by GEO: DMSO GSM2075714/GSM2075715 and in-vivo NAI-N3
  GSM2075718/GSM2075719.
- GEO sample metadata states that reads were mapped to the human transcriptome `GRCh37.74 (hg19)`.
  The processed table is transcript-based: transcript ID, length, RPKM, then one score per base;
  `NULL` denotes an unmeasured base.
- Verdict: **confirmed**. It is not natively coordinate-compatible with the GRCh38 GLORI files.

### HEK293T m6A

- GEO Series: [GSE210563](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210563).
- Replicate 1: [GSM6432590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6432590),
  `GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv.gz`.
- Replicate 2: [GSM6432591](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6432591),
  `GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv.gz`.
- Both are WT HEK293T, GLORI+, Day 0, poly(A) RNA and are labelled rep1/rep2. GEO specifies GRCh38
  for the human GLORI processing.
- The associated primary paper is Liu et al., “Absolute quantification of single-base m6A methylation
  in the mammalian transcriptome using GLORI,” *Nature Biotechnology* 41, 355–366 (2023),
  [DOI 10.1038/s41587-022-01487-9](https://doi.org/10.1038/s41587-022-01487-9).
- Verdict: **confirmed**.

### HeLa experimental structure

- GEO Series: [GSE145805](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145805).
- Samples: [GSM4333257](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4333257) is HeLaD
  (DMSO control); GSM4333258 is HeLaN (in-vivo NAI-N3).
- GEO explicitly says the in-vivo structure was calculated from D and N with `icSHAPE-pipe` and that
  human reads were processed against hg38. Thus D/N are the paired control/treatment inputs to the
  chemical-probing calculation, not biological replicates to average as independent measurements.
- Available Series bundle: `GSE145805_RAW.tar` (TAR of TXT, 424.4 MB on the checked GEO record).
- Associated primary paper: Sun et al., “Predicting dynamic cellular protein–RNA interactions by deep
  learning using in vivo RNA structures,” *Cell Research* 31, 495–516 (2021),
  [DOI 10.1038/s41422-021-00476-y](https://doi.org/10.1038/s41422-021-00476-y).
- Verdict: **confirmed with caveat**. The archive must be listed before assuming which member is the
  final HeLa per-nucleotide reactivity table.

### HeLa m6A control set

- Replicate 1: [GSM6432595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6432595),
  `GSM6432595_Hela-1_35bp_m2.totalm6A.FDR.csv.gz`.
- Replicate 2: [GSM6432596](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6432596),
  `GSM6432596_Hela-2_35bp_m2.totalm6A.FDR.csv.gz`.
- Both are WT HeLa, `hypoxia-`, GLORI+, Day 0, poly(A) RNA. They are the intended no-hypoxia
  control pair; the hypoxia-treated samples are GSM6432598/GSM6432599 and are not part of this candidate set.
- Verdict: **confirmed**.

## Reference-version consequence

| Input | Submitted reference |
|---|---|
| HEK293T GSE74353 icSHAPE | GRCh37.74 / hg19 transcriptome |
| HEK293T and HeLa GSE210563 GLORI | GRCh38 |
| HeLa GSE145805 icSHAPE | hg38 |

The discovery pair has a real reference-version mismatch. A transcript-aware mapping or a documented
cross-version conversion is therefore mandatory before joining HEK293T structure and GLORI sites.

## Download size and access

| Candidate input | GEO-reported download size |
|---|---:|
| HEK293T in-vivo icSHAPE processed reactivity | 8.8 MB |
| HEK293T GLORI rep1 + rep2 | 7.7 + 7.7 = 15.4 MB |
| HeLa icSHAPE Series TAR bundle | 424.4 MB |
| HeLa GLORI control rep1 + rep2 | 5.0 + 5.1 = 10.1 MB |
| **Candidate download total** | **approximately 458.7 MB** |

These are public GEO records and the linked processed files are downloadable without a paid subscription.
The total above is for the processed candidate inputs in the current roadmap. It does not include raw SRA
FASTQ data, reference genome/annotation assets, extracted copies, caches, intermediates, or results. Reserve
at least 2–3 GB for initial download/extraction work and retain the roadmap's 10–20 GB whole-project budget.

## Remaining limitations

- All six candidate downloads are present locally and passed container/full-stream integrity checks.
  SHA-256 values are recorded in `source_manifest.csv` and `download_checksums.csv`.
- GLORI field semantics and coordinate convention have now been checked against the paper, the
  GLORI-tools implementation, and all 708,103 locally downloaded result rows. The evidence and working
  decisions are recorded in `01_glori_field_semantics.md`.
- `GSE145805_RAW.tar` was listed successfully. The required final HeLa member is
  `GSM4333258_HeLa.out.txt.gz`; its nested gzip stream also passed integrity checking.
- The mapping references are now downloaded and version-pinned in `data/reference/`; provenance,
  versions, sizes, SHA-256 values, and integrity results are recorded in `02_reference_manifest.csv`.
  Ensembl release 88 is the selected operational GRCh38 transcript annotation because it gives the best
  tested agreement with the versioned HeLa transcript IDs, but GEO specifies only “hg38”, not an exact
  Ensembl release. This selection is therefore a documented compatibility inference rather than an
  author-declared version.

## AI-assistance disclosure

This verification record was prepared with AI-assisted web retrieval and local environment inspection.
Repository accessions, filenames, conditions, and reference builds were checked against linked primary
GEO records; unresolved points are marked rather than inferred.
