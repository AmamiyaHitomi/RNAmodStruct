# ModStruct first phase method expansion summary

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run + validate
- Origin Date: 2026-09-13
- Verification Status: VERIFIED
- Version Label: phase1_methods_extension_v1
- Evidence boundary: exploratory expansion after the main results are known; HeLa only conducts external reviews after freezing the process and does not interpret predictions or significance as causal mechanisms.

## Completion status

12–17 Each of the six implementation packages features frozen YAML, executable scripts, machine-readable tables, technical reports, and proprietary tests. The deterministic full-process rerun compared 16 machine products, and the result was `REPRODUCIBLE`. Original 06–10 mainline not written back.

## Main results

1. **Dual cell line summary**: Fixed effect β/SD=0.01007 (95% CI 0.00713–0.01301); DL random effect β/SD=0.01027 (95% CI 0.00675–0.01379), I²=11.2%. With only two cell lines, heterogeneity estimates are unstable.
2. **Conditional dependence**: HEK293T residual correlation=0.0222, group sign flip p=0.169; HeLa residual correlation=0.0141, p=0.02639. This is not a high-dimensional conditional mutual information estimate.
3. **Structural entropy**: After adjusting for covariates and R_flank10, HEK293T β/SD=-0.00328 (BH q=0.6873); HeLa β/SD=-0.00460 (BH q=0.01167). Stable independent effects were not supported in the development set, and the small negative result for HeLa was therefore only an exploratory idiosyncratic signal.
4. **Isomer sensitivity**: The β/SD range of four calibers of HEK293T is 0.01430–0.02051; HeLa is 0.00943–0.01151, and the direction is stable. Weight does not represent heterogeneitybody abundance.
5. **401 nt window**: The full window retains positions 2,143 for HEK293T and 16,818 for HeLa. HeLa M1 MAE=0.22758, M4 MAE=0.22875, structure expansion without external gain.
6. **Nonlinear model**: HeLa M4 MAE: Ridge=0.21463, HGB=0.18942. HGB is better than Ridge with the same characteristics, but the M4 of HGB is ΔMAE=-0.001687 relative to M2, and the experimental structure still does not provide stable external gain.

## Comprehensive conclusion

The first period strengthened the robustness boundary: the small positive correlation of the main line maintained direction at meta and isomeric calibers, but neither conditional dependence, 401 nt nor nonlinear analysis showed stable additional predictive value from the experimental structure. Structural entropy results are inconsistent between development and external sets and cannot be upgraded to mechanistic conclusions. Both positive and negative results are retained.

## Statistical fallacy scan

Checked for type 11/11 risks; full machine table in `results/tables/phase1_statistical_validation.csv`. The main reservation caveats are post hoc expansion, unstable heterogeneity between the two studies, selectivity of the full 401 nt, and that all observational results should not be interpreted causally.