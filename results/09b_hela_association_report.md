# Stage 09B HeLa association replication report

## Verdict

**Association replicates.** Re-estimating the frozen stage-07 model on HeLa (25,996 sites, 4,829 genes)
gives a positive and statistically distinguishable association between local experimental icSHAPE
reactivity and pooled GLORI m6A proportion, in the same direction as HEK293T.

## Primary result

| Quantity | HeLa | HEK293T (frozen stage 07) |
|---|---:|---:|
| Sites / genes | 25,996 / 4,829 | 4,409 / 1,324 |
| β per SD reactivity | +0.00947 | +0.01430 |
| Cluster-robust 95% CI (per SD) | +0.00633 to +0.01261 | +0.00595 to +0.02264 |
| Cluster-bootstrap 95% CI (per SD, 1,000) | +0.00643 to +0.01271 | +0.00650 to +0.02276 |
| p value | 3.53×10⁻⁹ | 0.000787 |
| Approximate percentage points per SD | ≈0.95 | ≈1.43 |
| R² | 0.157 | 0.174 |

Both intervals exclude zero with the same (positive) direction; HeLa shows a slightly smaller magnitude.
One predictor SD is 0.13538 icSHAPE reactivity units in HeLa.

## Sensitivity family (7 prespecified, BH adjusted)

All seven estimates are positive (range 0.0058–0.0115 per SD), all adjusted q ≤ 0.00048, including
excluding −2..+2, complete flank coverage only, each GLORI replicate as outcome, per-gene equal weight,
single-isoform sites only, and canonical DRACH sites only.

## Diagnostics

- Heteroscedasticity (Breusch–Pagan p ≈ 0) → inference uses CR1 gene-cluster covariance and cluster
  bootstrap.
- Functional form: RESET p = 1.1×10⁻¹⁹ and spline-minus-linear AIC = −2.70 → modest descriptive
  preference for a nonlinear form; the linear coefficient is an average local trend.
- Post hoc region stratification (BH adjusted): 3′UTR +0.0124/SD (q = 2.5×10⁻⁹), CDS +0.0040/SD
  (q = 0.37), 5′UTR +0.0054/SD (q = 0.43), noncoding/other +0.0083/SD (q = 0.43). All regions are
  positive; the effect is concentrated in 3′UTR.

## Deviation

The HeLa icSHAPE file carries no abundance/RPKM field (literal `*`), so the `log1p_icshape_abundance_rpkm`
covariate is unavailable and excluded from the HeLa association design (recorded in the run status and
decision log). This is a technical deviation, not a change to the frozen predictor or outcome definition.

## Interpretation boundary

Observational, cross-study integration; adjustment does not identify a causal effect. The estimand is
site-level, conditional on GLORI detection in both replicates and adequate icSHAPE flank coverage.
