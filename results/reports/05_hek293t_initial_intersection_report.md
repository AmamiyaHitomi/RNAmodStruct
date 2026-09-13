# HEK293T preliminary GLORI–icSHAPE intersection

Run date: 2026-09-10  
Status: **PASS**

## Fixed interpretation and coordinate path

- GLORI sites use `(Chr, Sites, Strand)`, where `Sites` is the 1-based modified genomic A position.
- The icSHAPE score vector is materialized as a 1-based, transcript-oriented coordinate: the first score
  corresponds to cDNA base 1. This is checked by exact agreement among the declared vector length,
  GRCh37.74 cDNA length, and summed GTF exon length.
- HEK293T mapping follows `transcript position → GRCh37.74 exon → hg19 genomic base →
  hg19ToHg38 chain → GRCh38 genomic base`.
- The main structure window is ±10 nt with the center excluded. Both sides must independently contain at
  least 70% valid icSHAPE values. A separate ±50 coverage value is retained for display/QC.
- Transcript isoforms are retained as separate mapping rows. Genomic-site summaries count a site once.

## Audit and coordinate gate

- All 6 source data assets and 9 reference assets passed size/SHA-256 and full-stream integrity checks.
- All six datasets had zero malformed rows; both icSHAPE datasets had zero vector-length failures.
- The GLORI `Ratio` and `NormeRatio` identities had zero failures across 708,103 rows.
- The coordinate prototype contains 30 stratified cases: both strands, exon-internal and splice-boundary,
  and GLORI-matched/non-matched positions. Twelve are GLORI matches.
- All 30 have transcript base A, the correct GRCh38 strand-oriented center base, and a successful
  liftOver round trip. There were zero failed assertions.

## Preliminary intersection

| Quantity | Result |
|---|---:|
| GLORI replicate 1 unique sites | 214,715 |
| GLORI replicate 2 unique sites | 215,021 |
| Exact replicate-common sites | 176,649 |
| Replicate-union sites | 253,087 |
| Union sites with ≥1 HEK293T icSHAPE transcript position | 20,559 (8.12%) |
| Replicate-1 sites with an icSHAPE position | 18,276 (8.51%) |
| Replicate-2 sites with an icSHAPE position | 18,444 (8.58%) |
| Replicate-common sites with an icSHAPE position | 16,161 (9.15%) |
| Matched transcript/isoform rows | 33,962 |
| Non-`ELSE` GLORI gene labels among matched sites | 2,593 |
| Matched sites with ≥1 valid ±10 structure window | 5,415 (26.34%) |
| Valid transcript-level ±10 windows | 9,547 / 33,962 (28.11%) |
| Transcript rows with a non-missing center reactivity | 10,959 / 33,962 (32.27%) |
| Mean / median complete-window ±10 coverage | 0.333 / 0.000 |
| GRCh38 center-base failures among matched sites | 0 |

The zero median is not a coding of missingness as zero: missing values remain NaN. It reflects the
block-like sparsity of the submitted icSHAPE tracks; the mean and the explicit 70%-per-side rule are more
informative for eligibility.

## Why union sites did not match

| Reason | Sites | Fraction of GLORI union |
|---|---:|---:|
| No exon belonging to an eligible icSHAPE transcript | 230,866 | 91.220% |
| hg19→hg38 chain position unmapped | 939 | 0.371% |
| Chain mapping ambiguous | 722 | 0.285% |
| Transcript-oriented center was not A | 1 | <0.001% |

The single non-A case is mitochondrial `chrM:13710:+` (`ND5`) and is retained in the status table rather
than silently discarded. It does not undermine the nuclear-coordinate checks but should remain excluded
from the joined structure analysis.

## Output interpretation

- `05_hek293t_glori_site_mapping_status.csv.gz` has one row per GLORI union site and is the denominator and
  attrition record.
- `05_hek293t_glori_icshape_intersection.csv.gz` has one row per genomic-site/transcript pair and preserves
  isoform multiplicity, replicate-specific `Ratio`, center reactivity, and window coverage.
- These results are a coordinate/coverage feasibility analysis, not yet a biological association test.
  The 5,415 sites with at least one valid main window clear the roadmap's provisional ≥500-site gate,
  but split construction, isoform resolution, covariates, and analysis thresholds still need to be
  frozen before testing biological hypotheses.

## Reproduction

Run from the repository root with the project Python interpreter:

```powershell
& 'E:\ancd\envs\my_pytorch\python.exe' 'src/03_audit_inputs.py'
& 'E:\ancd\envs\my_pytorch\python.exe' 'src/04_validate_coordinate_mapping.py'
& 'E:\ancd\envs\my_pytorch\python.exe' 'src/05_compute_initial_intersection.py'
& 'E:\ancd\envs\my_pytorch\python.exe' -m unittest discover -s tests -p 'test_*.py' -v
```

Detailed machine-readable counts are in `05_hek293t_initial_intersection_summary.csv`; run completion is
recorded in `05_hek293t_initial_intersection_run_status.csv`, and each stage has a uniquely prefixed log
under `results/logs/`.
