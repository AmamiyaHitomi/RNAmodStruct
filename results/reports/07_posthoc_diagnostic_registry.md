# Stage 07 post hoc diagnostic registry

## Transcript-region stratification

- Status: post hoc; not part of the frozen primary or sensitivity family.
- Trigger: the required statistical-fallacy scan called for an explicit check of whether the pooled direction could conceal region-specific reversal (Simpson's paradox).
- Method: repeat the adjusted gene-clustered model within each transcript region having at least 200 sites and 100 genes; apply Benjamini–Hochberg correction within this diagnostic family.
- Permitted use: direction-consistency and heterogeneity diagnosis only.
- Prohibited use: confirmatory support for the primary claim or replacement of the prespecified sensitivity family.
- Output: `results/tables/07_hek293t_posthoc_region_diagnostic.csv`.
