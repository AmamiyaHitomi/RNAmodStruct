# Stage 08 / 08A / 08B validation report

## Verdict

**PASS with interpretive caution.** The full M0--M4 Ridge matrix, ViennaRNA predicted-structure feature
generation, development-only expanded experimental-structure diagnostic, leakage audit, saved models,
and repeat-agreement analysis are complete.

## Stage 08A evidence

- ViennaRNA 2.7.2, 37°C, dangles=2, partition function enabled.
- 4,092 unique complete 201-nt sequences; 4,096 mapped model sites.
- Every unpaired-probability vector has 201 finite values bounded by [0,1].
- Saved P features include position-wise unpaired probabilities, MFE/nt, ensemble free energy/nt,
  ensemble diversity, pairing-state entropy, and prespecified window summaries.
- M2 test MAE was 0.21193 versus 0.20546 for M1: predicted structure worsened MAE by 0.6466 percentage
  points (paired group-bootstrap 95% CI 0.2635 to 1.0328 points worse).
- M4 test MAE was 0.21188. Its improvement over M2 was 0.0052 percentage points, with 95% CI −0.0176
  to 0.0267; no reliable incremental experimental-structure benefit was detected conditional on P.

Because predicted structure is a deterministic nonlinear transformation of the same sequence, these
results assess whether this particular representation helps a Ridge model. They do not establish that
computational RNA structure is biologically irrelevant.

## Stage 08B evidence

Only the frozen development set was used: 3,271 sites, 922 genes, and 919 joint groups. The compared
blocks were no R, flank10, directional up/down, the −50..+50 profile without center, and an all-position
sensitivity block. All preprocessing was fitted inside five-fold GroupKFold.

No R block improved development MAE over M1. `R_none` was selected. The least harmful non-null block,
`R_core`, worsened MAE by 0.0029 percentage points; the no-center full profile worsened it by 0.1779
points. Expanded experimental profiles therefore should not be promoted to external confirmatory input.

GLORI replicate agreement was high (Spearman 0.9744; Pearson 0.9875), with replicate MAE 2.65
percentage points. Replicate MAE declined from 5.44 points in the lowest coverage quartile to 0.98 in
the highest. This is a measurement-variability reference, not a formal noise ceiling, but it shows that
repeat disagreement is much smaller than model prediction error.

## Validation boundary

Gene, identical-sequence, and joint-group overlap across the original HEK293T split are all zero.
However, HEK293T test outcomes were viewed before a preprocessing correction and before 08A was added.
All current test comparisons are therefore descriptive/reused rather than untouched confirmation. The
alpha grid was later widened using development CV only; all selected optima are now interior and guarded
by a boundary check. Confirmatory predictive evidence comes from the HeLa frozen-model transfer, which
did not use HeLa labels for feature or parameter selection.
