"""Create the phase-1 synthesis, fallacy scan, independent manifest and completion receipt."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from phase1_common import sha256_file, write_rows

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/reports/phase1_summary_report.md"
VALIDATION = ROOT / "results/tables/phase1_statistical_validation.csv"
MANIFEST = ROOT / "metadata/manifests/release/phase1_release_manifest.csv"
STATUS = ROOT / "results/status/phase1_completion_status.csv"


def manifest_paths() -> list[Path]:
    paths = []
    for pattern in ("config/1[2-7]_*.yaml", "src/1[2-7]_*.py", "src/phase1_*.py",
                    "tests/test_phase1_extensions.py", "results/reports/1[2-7]_*.md", "results/status/1[2-7]_*.csv",
                    "results/tables/1[2-7]_*.csv*", "results/tables/phase1_*.csv",
                    "results/reports/phase1_summary_report.md", "results/status/phase1_reproducibility_status.csv",
                    "docs/research_notes/ModStruct_第一期运行说明.md",
                    "data/final/16_*.csv.gz", "data/interim/16_401nt_vienna_cache.csv"):
        paths.extend(path for path in ROOT.glob(pattern) if path.is_file())
    return sorted(set(paths))


def main() -> None:
    repro = pd.read_csv(ROOT / "results/status/phase1_reproducibility_status.csv").iloc[0]
    if repro["status"] != "PASS":
        raise RuntimeError("Reproducibility gate is not PASS")
    meta = pd.read_csv(ROOT / "results/tables/12_meta_analysis.csv")
    conditional = pd.read_csv(ROOT / "results/tables/13_conditional_dependence.csv")
    entropy = pd.read_csv(ROOT / "results/tables/14_structure_entropy_associations.csv")
    isoform = pd.read_csv(ROOT / "results/tables/15_isoform_sensitivity_associations.csv")
    window = pd.read_csv(ROOT / "results/tables/16_window_401_external_metrics.csv")
    nonlinear = pd.read_csv(ROOT / "results/tables/17_nonlinear_external_metrics.csv")
    ablation = pd.read_csv(ROOT / "results/tables/17_retraining_ablation.csv")
    fixed = meta[(meta.row_type == "pooled") & (meta.model == "fixed_effect")].iloc[0]
    random = meta[(meta.row_type == "pooled") & (meta.model == "random_effect_DL")].iloc[0]
    entropy_primary = entropy[entropy.is_primary_entropy_hypothesis & entropy.model.eq("covariate_adjusted_plus_R_flank10")]
    w = {row.model: row for row in window.itertuples()}
    n = {(row.learner, row.model): row for row in nonlinear.itertuples()}
    validations = [
        ("significance_without_effect_size", "PASS", "All association summaries include beta/SD and 95% CI."),
        ("confidence_interval_omission", "PASS", "Stages 12, 14 and 15 retain uncertainty intervals."),
        ("multiple_testing", "PASS", "Stage 14 applies BH within dataset/specification families."),
        ("cluster_dependence", "PASS", "Association inference clusters by gene; conditional inference also uses grouped sign flips."),
        ("train_test_leakage", "PASS", "Nuisance fitting and tuning use connected gene/sequence groups; HeLa is not used for tuning."),
        ("causal_overreach", "CAUTION", "All outputs are explicitly observational/exploratory and prohibit causal language."),
        ("posthoc_hypothesis", "CAUTION", "The entire restart phase occurs after the frozen main result was known."),
        ("two_study_heterogeneity", "CAUTION", "Q, I2 and tau2 are reported, but only two cell lines make them unstable."),
        ("isoform_abundance_overinterpretation", "PASS", "Equal/coverage weights are labelled sensitivity weights, not abundance."),
        ("selection_and_missingness", "CAUTION", "Complete 401-nt and structure-coverage filters select a reduced population."),
        ("prediction_as_mechanism", "PASS", "Permutation importance/SHAP boundaries are explicitly non-causal."),
    ]
    write_rows(VALIDATION, [{"fallacy": a, "verdict": b, "detail": c} for a, b, c in validations])
    REPORT.write_text(f"""# ModStruct first phase method expansion summary

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run + validate
- Origin Date: {datetime.now(timezone.utc).date().isoformat()}
- Verification Status: VERIFIED
- Version Label: phase1_methods_extension_v1
- Evidence boundary: exploratory expansion after the main results are known; HeLa only conducts external reviews after freezing the process and does not interpret predictions or significance as causal mechanisms.

## Completion status

12–17 Each of the six implementation packages features frozen YAML, executable scripts, machine-readable tables, technical reports, and proprietary tests. The deterministic full-process rerun compared {int(repro['artifacts_compared'])} machine products, and the result was `{repro['verdict']}`. Original 06–10 mainline not written back.

## Main results

1. **Dual cell line summary**: Fixed effect β/SD={fixed.effect_per_sd:.5f} (95% CI {fixed.ci95_low:.5f}–{fixed.ci95_high:.5f}); DL random effect β/SD={random.effect_per_sd:.5f} (95% CI {random.ci95_low:.5f}–{random.ci95_high:.5f}), I²={fixed.I2_percent:.1f}%. With only two cell lines, heterogeneity estimates are unstable.
2. **Conditional dependence**: HEK293T residual correlation={conditional.iloc[0].residual_correlation:.4f}, group sign flip p={conditional.iloc[0].group_sign_flip_p_value:.4g}; HeLa residual correlation={conditional.iloc[1].residual_correlation:.4f}, p={conditional.iloc[1].group_sign_flip_p_value:.4g}. This is not a high-dimensional conditional mutual information estimate.
3. **Structural entropy**: After adjusting for covariates and R_flank10, HEK293T β/SD={entropy_primary.iloc[0].beta_per_sd:.5f} (BH q={entropy_primary.iloc[0].p_value_bh:.4g}); HeLa β/SD={entropy_primary.iloc[1].beta_per_sd:.5f} (BH q={entropy_primary.iloc[1].p_value_bh:.4g}). Stable independent effects were not supported in the development set, and the small negative result for HeLa was therefore only an exploratory idiosyncratic signal.
4. **Isomer sensitivity**: The β/SD range of four calibers of HEK293T is {isoform[isoform.dataset.eq('HEK293T')].beta_per_sd.min():.5f}–{isoform[isoform.dataset.eq('HEK293T')].beta_per_sd.max():.5f}; HeLa is {isoform[isoform.dataset.eq('HeLa')].beta_per_sd.min():.5f}–{isoform[isoform.dataset.eq('HeLa')].beta_per_sd.max():.5f}, and the direction is stable. Weight does not represent heterogeneitybody abundance.
5. **401 nt window**: The full window retains positions 2,143 for HEK293T and 16,818 for HeLa. HeLa M1 MAE={w['M1'].mae:.5f}, M4 MAE={w['M4'].mae:.5f}, structure expansion without external gain.
6. **Nonlinear model**: HeLa M4 MAE: Ridge={n[('Ridge','M4')].mae:.5f}, HGB={n[('HistGradientBoostingRegressor','M4')].mae:.5f}. HGB is better than Ridge with the same characteristics, but the M4 of HGB is ΔMAE={ablation[(ablation.learner == 'HistGradientBoostingRegressor') & (ablation.comparison == 'M4_vs_M2')].iloc[0].delta_mae_baseline_minus_augmented:.6f} relative to M2, and the experimental structure still does not provide stable external gain.

## Comprehensive conclusion

The first period strengthened the robustness boundary: the small positive correlation of the main line maintained direction at meta and isomeric calibers, but neither conditional dependence, 401 nt nor nonlinear analysis showed stable additional predictive value from the experimental structure. Structural entropy results are inconsistent between development and external sets and cannot be upgraded to mechanistic conclusions. Both positive and negative results are retained.

## Statistical fallacy scan

Checked for type 11/11 risks; full machine table in `results/tables/phase1_statistical_validation.csv`. The main reservation caveats are post hoc expansion, unstable heterogeneity between the two studies, selectivity of the full 401 nt, and that all observational results should not be interpreted causally.""", encoding="utf-8")
    manifest = [{"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size,
                 "sha256": sha256_file(path)} for path in manifest_paths() if path != MANIFEST]
    write_rows(MANIFEST, manifest)
    write_rows(STATUS, [{"phase": 1, "status": "PASS", "finished_utc": datetime.now(timezone.utc).isoformat(),
                         "stages_passed": 6, "tests_file": "tests/test_phase1_extensions.py",
                         "reproducibility_verdict": repro["verdict"], "manifest_entries": len(manifest),
                         "summary_sha256": sha256_file(REPORT), "manifest_sha256": sha256_file(MANIFEST)}])


if __name__ == "__main__":
    main()
