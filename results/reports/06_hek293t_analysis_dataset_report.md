# Stage 06: HEK293T analysis-dataset build

Run date: 2026-09-10  
Status: **PASS**

## Frozen primary rules

- Primary site set: sites detected in both GLORI replicates.
- Primary m6A response: count-weighted pooled A rate,
  `(Acov_rep1 + Acov_rep2) / (AGcov_rep1 + AGcov_rep2)`.
- No new `AGcov` threshold is imposed after the deposited GLORI FDR filtering.
- Primary experimental-structure feature uses ±10 nt and excludes the center. Each side separately needs
  at least 7/10 valid icSHAPE measurements.
- A full, unpadded 201-nt transcript sequence is required only for the later model-ready subset.
- Missing icSHAPE measurements remain null and are never converted to zero.

All decisions are recorded in `config/06_analysis_config.yaml` before the association-analysis stage.

## Deterministic isoform resolution

For each genomic site, transcript mappings are ranked by:

1. valid primary ±10 window;
2. larger valid-value count in the ±50 display window;
3. agreement between GLORI and Ensembl gene names;
4. ascending transcript ID as the deterministic final tie-break.

The selected representative produces the site-level table. All 33,962 mappings remain available in the
all-isoform table for sensitivity analyses, with exactly one rank-1 representative per genomic site.

## Dataset sizes

| Dataset stage | Sites |
|---|---:|
| GLORI replicate union | 253,087 |
| Mapped to at least one icSHAPE transcript | 20,559 |
| Detected in both GLORI replicates | 16,161 |
| Primary structure-analysis dataset | 4,409 |
| Model-ready subset with full 201-nt sequence | 4,096 |

The primary structure dataset contains 1,324 analysis genes. The largest gene contributes 38 sites.
Among the 20,559 mapped genomic sites, 7,266 have more than one compatible transcript isoform; the
maximum is 19 isoforms at one site.

## Replicate QC in the primary dataset

| Metric | Result |
|---|---:|
| Pearson correlation of replicate `Ratio` | 0.9879 |
| Spearman correlation of replicate `Ratio` | 0.9728 |
| Median absolute replicate difference | 0.01549 |

These values describe replicate agreement only. They are not tests of the relationship between m6A and
RNA structure.

## Exclusions after transcript mapping

| Site-level outcome | Sites |
|---|---:|
| Included in primary structure dataset | 4,409 |
| Insufficient ±10 icSHAPE coverage | 11,752 |
| Not detected in both GLORI replicates | 4,398 |

## Outputs

- `data/final/06_hek293t_site_level_dataset.csv.gz`: one deterministic representative per mapped site,
  retaining inclusion flags and exclusion reasons.
- `data/final/06_hek293t_main_analysis_dataset.csv.gz`: the 4,409-row primary structure dataset.
- `data/final/06_hek293t_all_isoform_mappings.csv.gz`: all transcript mappings for sensitivity analysis.
- `results/tables/06_hek293t_attrition_summary.csv`: stage-by-stage site attrition.
- `results/tables/06_hek293t_dataset_qc_summary.csv`: machine-readable QC metrics.
- `results/tables/06_hek293t_gene_site_counts.csv`: gene-cluster sizes for split and dependence checks.

Automated tests verify site uniqueness, pooled-ratio arithmetic, main-window eligibility, fixed JSON
window lengths, 201-nt center A, and exactly one representative isoform per site.
