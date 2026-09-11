# Stage 08 prediction-modeling decision and deviation log

## Frozen minimum-version decision

- Population: the 4,096 stage-06 records with complete 201-nt sequence and valid experimental
  `R_flank10`.
- Independence unit: connected components formed from `analysis_gene` and identical 201-nt sequences.
- Split: deterministic 80% development / 20% test assignment with seed 20260910.
- Development selection: five-fold `GroupKFold`; final Ridge alpha candidates 0.01, 0.1, 1, 10, 30,
  100, 300, 1000, and 3000. Boundary optima are rejected rather than frozen.
- Available model matrix: M0 (common background), M1 (background + sequence), and M3 (background +
  sequence + experimental `R_flank10`). M2 and M4 remain pending because predicted-structure features
  are absent.
- Primary available comparison: ΔMAE = MAE(M1) − MAE(M3), with a paired test-group bootstrap.

## 2026-09-10 preprocessing correction after first test view

The first executable version standardized the common continuous covariates, but left normalized 3-mer
frequencies and `R_flank10` on their raw scales. Ridge penalization is scale-dependent, so this was an
implementation defect. The first test output had already been inspected:

| Version | M1 MAE | M3 MAE | ΔMAE (percentage points) | 95% CI (percentage points) |
|---|---:|---:|---:|---:|
| Initial, defective scaling | 0.2068117 | 0.2068004 | 0.00113 | −0.00965 to 0.01065 |
| Corrected training-fold scaling | 0.2053234 | 0.2052843 | 0.00391 | −0.03431 to 0.03873 |

The correction was fitted within each development fold and did not change the population, split,
candidate alphas, feature definitions, or learner family. Nevertheless, reusing the same test set after
viewing its labels means the corrected HEK293T result is not described as fully untouched. Creating a
new favorable HEK293T split would compound the problem, so the assignment is retained and the planned
HeLa analysis is designated as the next independent predictive validation.

## 2026-09-10 stage 08A completion

ViennaRNA 2.7.2 partition-function features were generated at 37°C with dangles=2 for all 4,092 unique
201-nt sequences and mapped to 4,096 sites. This completed M2 and M4. On the reused HEK293T test split,
After the final development-only alpha-grid audit, M2 was worse than M1 by 0.6466 MAE percentage points
(paired group-bootstrap 95% CI 0.2635 to 1.0328 points worse). M4 improved on M2 by only 0.0052 points
(95% CI −0.0176 to 0.0267). These results are descriptive because the HEK293T test outcomes had already
been viewed before 08A was added.

## 2026-09-10 stage 08B development decision

Five experimental-structure blocks were compared using only the frozen development groups. None
improved mean GroupKFold MAE over the M1 sequence baseline. `R_none` was therefore selected; `R_core`
was merely the least harmful non-null block. The GLORI replicate-to-replicate MAE was 2.65 percentage
points overall and decreased from 5.44 points in the lowest depth quartile to 0.98 points in the highest.
Observed replicate disagreement is much smaller than the approximately 20.5-point prediction MAE, so
replicate noise alone does not explain the model error.

## 2026-09-11 alpha-grid acceptance correction

An acceptance audit detected that M1--M4 selected alpha=100 at the upper edge of the original grid and
that development CV MAE was still decreasing. Without consulting HEK293T test or HeLa labels, the grid
was expanded to 0.01--3000. M0 selected 0.1 and M1--M4 selected 300, all strictly inside the grid. The
pipeline now raises an error whenever an optimum lies on either boundary. Stages 08, 08B, and 09C were
regenerated. The substantive conclusion was unchanged: no experimental-structure predictive increment
was detected internally, and the frozen increment remained negative on HeLa.
