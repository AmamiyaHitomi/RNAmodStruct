# Stage 08 HEK293T prediction-model report

## Scope and freeze

The model population contains 4,096 sites with complete 201-nt sequence windows and valid
experimental flank reactivity. The deterministic joint grouping joins every site from the same gene and
also joins genes connected by an identical 201-nt sequence. It produced 1,205
independent groups. A single frozen split assigned 3,271 sites to development and
825 sites to the held-out test set. Gene, identical-sequence, and joint-group overlap were all zero.

All imputation, scaling, category discovery, and Ridge alpha selection were fitted within development
folds. The candidate alphas were [0.01, 0.1, 1.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0]; selected values were M0=0.1, M1=300,
M2=300, M3=300, and M4=300. Predictions were clipped to [0, 1]
for every model.

### Protocol deviation

An initial implementation was evaluated on this same test split before the continuous 3-mer and
experimental-reactivity features were correctly standardized within each training fold. That first run
gave M1 MAE 0.2068117, M3 MAE 0.2068004, and ΔMAE 0.00113 percentage points. The preprocessing defect
was then corrected without changing the population, split, candidate alphas, feature definitions, or
model family, and the results below were generated. Because test outcomes had already been viewed, the
corrected result is transparently classified as **test-reused after a preprocessing correction**, not as
a fully untouched confirmatory test. A later acceptance audit found that the original alpha grid ended
at the development-selected boundary (100); the grid was widened using development CV only, with a
hard failure now preventing boundary optima from being frozen. No alternative split was searched.
Independent confirmation comes from the external HeLa stage.

## Held-out result

| Model | Inputs | MAE | RMSE | Spearman | R2 |
|---|---|---:|---:|---:|---:|
| M0 | common background | 0.20820 | 0.25053 | 0.3284 | 0.1247 |
| M1 | background + sequence | 0.20546 | 0.24890 | 0.3719 | 0.1361 |
| M2 | background + sequence + predicted structure | 0.21193 | 0.25697 | 0.3218 | 0.0792 |
| M3 | background + sequence + experimental R_flank10 | 0.20549 | 0.24875 | 0.3710 | 0.1371 |
| M4 | background + sequence + predicted + experimental structure | 0.21188 | 0.25680 | 0.3220 | 0.0804 |

The primary complete-matrix comparison gives ΔMAE = MAE(M2) − MAE(M4) =
0.005 percentage points (paired 286-group
bootstrap 95% CI -0.018 to
0.027; 2,000 replicates). Positive values
mean lower error after adding experimental structure.

Predicted structure did not improve this Ridge sequence baseline: M2 MAE exceeded M1 MAE by
0.647 percentage points. This comparison is
descriptive within the reused test set and does not justify searching for a favorable alternative split.

## Interpretation boundary

The full Ridge matrix M0--M4 is complete. The M4-vs-M2 comparison estimates the incremental predictive
value of experimental R_flank10 beyond the frozen sequence and partition-function-derived structure
baseline. Predictive gain from experimental structure is not evidence of a causal mechanism, and the held-out interval covers
test-group sampling uncertainty for fixed fitted models rather than batch or training uncertainty.
