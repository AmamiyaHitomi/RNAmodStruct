## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: validate
- Origin Date: 2026-09-10
- Verification Status: VERIFIED
- Version Label: validation_v1
- Upstream Dependencies: stage-06 HEK293T analysis dataset and frozen stage-07 configuration

## Validation Report

- **Source**: `stage07_hek293t_association`
- **Overall Confidence**: **CAUTION**
- **Verification basis**: deterministic full re-run, automated invariants, rendered-figure QA, and statistical-fallacy scan

The selected-site analysis supports a small positive association between local experimental RNA
reactivity and pooled GLORI m6A proportion. It does not establish that RNA structure causes m6A, that
m6A causes structural change, or that the result generalizes to all GLORI sites.

### Statistical findings

| Metric | Test | Value | Effect size | Confidence |
|---|---|---:|---:|---|
| Primary association | Adjusted OLS with CR1 gene-cluster covariance | p = 0.000787 | +0.01430 m6A proportion per 1 SD reactivity; 95% CI 0.00595–0.02264 | SOLID within the selected analysis population |
| Cluster-bootstrap interval | 1,000 gene-cluster resamples | 1,000 successful | 95% CI 0.00650–0.02276 per SD | SOLID |
| Sensitivity family | Seven prespecified analyses, BH adjusted | all q ≤ 0.0343 | all estimates positive, range 0.01020–0.02051 per SD | SOLID for direction under the prespecified variants |
| Linear-model fit | Adjusted OLS | R² = 0.1742 | descriptive only | CAUTION |
| Functional form | RESET and 3-df spline sensitivity | RESET p = 1.88×10⁻15; spline AIC − linear AIC = −1.72 | modest descriptive preference for a nonlinear form | CAUTION |
| Region diagnostic | Post hoc region-stratified adjusted models, BH adjusted | 3′UTR q = 0.00588; CDS q = 0.01782; 5′UTR q = 0.3667; other q = 0.6392 | positive in 3′UTR/CDS; negative but imprecise in 5′UTR/other | CAUTION |

The primary population contains 4,409 sites from 1,324 genes. One predictor SD is 0.12309 icSHAPE
reactivity units. The fitted +0.01430 effect therefore corresponds to about **1.43 percentage points**
in pooled m6A proportion, conditional on the recorded covariates. This is statistically distinguishable
from zero but small in absolute magnitude.

### Model and data warnings

| Type | Detail | Affected |
|---|---|---|
| Heteroscedasticity | Breusch–Pagan p = 5.10×10⁻66. CR1 gene-cluster covariance and cluster bootstrap were used for inference. | Standard OLS standard errors must not be used. |
| Functional-form mismatch | RESET strongly rejects the simple linear specification; the spline improves AIC only slightly and preserves evidence of association. | Interpret the linear coefficient as an average local trend, not a universal slope. |
| Structure-availability selection | Valid-window versus insufficient-window SMDs were −0.619 for local GC, +0.583 for transcript position, +0.370 for GLORI coverage, +0.208 for abundance, and +0.362 for isoform count. | Generalization beyond the 4,409 selected sites. |
| Region heterogeneity | 3′UTR and CDS estimates are positive, whereas smaller 5′UTR and other-region estimates are negative with intervals crossing zero. | Claims of a uniform transcript-wide direction. |
| Cross-study integration | GLORI and icSHAPE measurements come from separate HEK293T experiments and are linked through reference annotation. | Batch, condition, and unmeasured biological confounding. |
| Bounded response | The response lies in [0,1]; OLS produced no out-of-range fitted values here, but distributional assumptions remain approximate. | Parametric model interpretation. |

### Fallacy scan

- **Coverage**: 11/11 fallacy types checked

| Fallacy | Severity | Detail | Recommendation |
|---|---|---|---|
| Simpson's paradox | CAUTION | The pooled positive direction is reproduced in 3′UTR and CDS, but not in the smaller 5′UTR and other-region strata. | Report region estimates as post hoc diagnostics and avoid claiming complete subgroup consistency. |
| Ecological fallacy | NOTE | The unit is an individual mapped site, with dependence clustered by gene; no gene-average effect is projected onto every site. | Keep the estimand at site level. |
| Berkson's bias | CAUTION | Entry requires GLORI detection and adequate icSHAPE coverage; jointly selected measurements can induce association. | Treat the finding as conditional on ascertainment and coverage. |
| Collider bias | CAUTION | Conditioning on successful cross-platform mapping and structure-window validity may condition on common consequences of biological/technical variables. | Do not interpret the adjusted coefficient causally without a defensible causal graph and independent validation. |
| Base-rate neglect | NOTE | No classifier, diagnostic threshold, PPV, or prevalence inference is used. | Not applicable to the reported association. |
| Regression to the mean | NOTE | There is no pre/post selection on an extreme baseline value. | No material issue identified. |
| Survivorship bias | CAUTION | Only 20,559 of 253,087 union GLORI sites map to icSHAPE transcripts, and only 4,409 enter the primary structure analysis. | Carry the full attrition table and unmatched reasons into any manuscript or downstream model. |
| Multiple comparisons / look-elsewhere effect | NOTE | One frozen primary model was used; seven prespecified sensitivity p-values were BH adjusted. The region analysis is separately labelled post hoc and BH adjusted. | Do not promote exploratory profile features or subgroup results to new confirmatory claims. |
| Researcher degrees of freedom | NOTE | Stage-06 and stage-07 YAML configurations were frozen before fitting; deterministic rules and logs are retained. There was no external preregistration. | Preserve configs and clearly distinguish prespecified from post hoc analyses. |
| Causal over-interpretation | CAUTION | The design is observational and cross-study. Adjustment does not identify an intervention effect. | Use “associated with,” never “causes,” “promotes,” or “inhibits.” |
| Reverse causality | CAUTION | Concurrent data cannot determine whether structure affects modification, modification affects structure, or both reflect another process. | Directionality requires perturbation or temporally resolved data. |

### Reproducibility

- **Method**: full deterministic re-run with identical configuration and seed, followed by SHA-256 comparison
- **Verdict**: **REPRODUCIBLE**

| Metric | Original | Re-run | Diff | Status |
|---|---:|---:|---:|---|
| Primary coefficient per SD | 0.01429675084 | 0.01429675084 | 0 | MATCH |
| Bootstrap interval per SD | 0.00649952053–0.02276398591 | 0.00649952053–0.02276398591 | 0 | MATCH |
| Deterministic result tables | 9 files | 9 files | all SHA-256 identical | MATCH |
| Automated tests | 14 | 14 | 0 failures | MATCH |

File-level hashes are recorded in `results/07_hek293t_reproducibility_check.csv`.

### Figure validation

The 180 × 105 mm overview figure passed source preflight (21/21 checks), panel-alignment audit,
rendered collision audit (0 failures, 0 warnings), PDF font audit (minimum 6 pt), and final-size visual
inspection. Panel-c profile intervals are pointwise gene-cluster bootstrap intervals and must not be read
as simultaneous bands.

### Final interpretation boundary

Within HEK293T sites detected in both GLORI replicates and having sufficient icSHAPE coverage, higher
mean ±10 nt icSHAPE reactivity is associated with a modestly higher pooled m6A proportion after the
specified adjustment. The result is robust across the seven prespecified sensitivity analyses, but
selection into the integrated dataset, residual nonlinearity, regional heterogeneity, cross-study
confounding, and unresolved directionality limit generalization and causal interpretation.

