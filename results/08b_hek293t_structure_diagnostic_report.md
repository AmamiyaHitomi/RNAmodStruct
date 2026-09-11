# Stage 08B development-only structure diagnostic report

## Scope

All feature-block comparisons use only the frozen HEK293T development set (3,271 sites,
922 genes, 919 joint groups).
The reused internal test set was not consulted for block selection. Every imputation and scaling step was
fitted within its GroupKFold training fold.

## Expanded experimental-structure blocks

| Block | Selected alpha | Mean CV MAE | ΔMAE vs M1 (percentage points) | Role |
|---|---:|---:|---:|---|
| R_none | 300 | 0.20451 | 0.0000 | candidate |
| R_core | 300 | 0.20454 | -0.0029 | candidate |
| R_directional | 300 | 0.20461 | -0.0092 | candidate |
| R_profile_no_center | 300 | 0.20629 | -0.1779 | candidate |
| R_profile_all | 300 | 0.20633 | -0.1812 | sensitivity_only |

Among the prespecified selectable blocks, `R_none` had the lowest development CV MAE
(0.20451). No experimental-structure block improved on the sequence baseline, so the
development-stage decision is to retain no expanded R block. `R_core` was the least harmful
non-null block (ΔMAE -0.0029 percentage points). This is a
development-stage selection result, not held-out evidence. The
all-position profile includes the modified center and is sensitivity-only because that position may
contain local probing or modification-linked signal.

## GLORI replicate measurement reference

Across the development set, replicate-1 versus replicate-2 MAE was 0.02651
(2.65 percentage points), RMSE was 0.04220,
Spearman was 0.9744, and Pearson was 0.9875.
Depth-stratified results are retained in the accompanying table. These discrepancies quantify observed
repeat variability; they are not a formal, identifiable upper bound on model performance because the
combined outcome and each replicate have different measurement-error properties.

## Decision boundary

No expanded-R block is selected. `R_none` is frozen as the development decision for external HeLa
evaluation; `R_core` may be retained only as a prespecified descriptive sensitivity analysis. Association,
predictive increment, and repeat agreement remain distinct estimands.
