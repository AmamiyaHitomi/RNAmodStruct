"""Independent HeLa re-estimate of the frozen stage-07 association model (stage 09B)."""

from __future__ import annotations

import csv
import json
import logging
import math
import statistics
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml
from patsy import dmatrix
from scipy.stats import kurtosis, skew
from statsmodels.stats.diagnostic import het_breuschpagan, linear_reset
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "09b_hela_association.yaml"
INPUT = ROOT / "data" / "final" / "09_hela_main_analysis_dataset.csv.gz"
SITE_INPUT = ROOT / "data" / "final" / "09_hela_site_level_dataset.csv.gz"
TABLES = ROOT / "results" / "tables"
RESULTS = ROOT / "results"
LOG = RESULTS / "logs" / "09b_run_hela_association.log"


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def load_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, compression="gzip", low_memory=False)


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        if not rows:
            return
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows(rows)


class DesignBuilder:
    """Frozen stage-07 design with the HeLa-unavailable abundance covariate excluded."""

    continuous_names = [
        "gc_fraction_21",
        "transcript_position_fraction",
        "log1p_combined_agcov",
        "coverage_flank10",
        "distance_to_stop_codon_tx",
        "distance_to_nearest_splice_edge_tx",
    ]

    def __init__(self, df: pd.DataFrame):
        self.continuous = {
            "gc_fraction_21": pd.to_numeric(df["gc_fraction_21"]),
            "transcript_position_fraction": pd.to_numeric(df["transcript_position_fraction"]),
            "log1p_combined_agcov": np.log1p(pd.to_numeric(df["combined_agcov"])),
            "coverage_flank10": pd.to_numeric(df["coverage_flank10"]),
            "distance_to_stop_codon_tx": pd.to_numeric(df["distance_to_stop_codon_tx"], errors="coerce"),
            "distance_to_nearest_splice_edge_tx": pd.to_numeric(df["distance_to_nearest_splice_edge_tx"], errors="coerce"),
        }
        self.stats = {}
        for name, series in self.continuous.items():
            median = float(series.median())
            filled = series.fillna(median)
            sd = float(filled.std(ddof=1))
            self.stats[name] = (median, sd if sd > 0 else 1.0, float(filled.mean()))
        categorical = pd.get_dummies(
            df[["drach_subtype", "transcript_region"]].astype(str),
            prefix=["drach", "region"],
            drop_first=True,
            dtype=float,
        )
        self.categories = list(categorical.columns)

    def build(self, df: pd.DataFrame, predictor: str) -> pd.DataFrame:
        output = pd.DataFrame(index=df.index)
        output["Intercept"] = 1.0
        output[predictor] = pd.to_numeric(df[predictor])
        raw_series = {
            "gc_fraction_21": pd.to_numeric(df["gc_fraction_21"]),
            "transcript_position_fraction": pd.to_numeric(df["transcript_position_fraction"]),
            "log1p_combined_agcov": np.log1p(pd.to_numeric(df["combined_agcov"])),
            "coverage_flank10": pd.to_numeric(df["coverage_flank10"]),
            "distance_to_stop_codon_tx": pd.to_numeric(df["distance_to_stop_codon_tx"], errors="coerce"),
            "distance_to_nearest_splice_edge_tx": pd.to_numeric(df["distance_to_nearest_splice_edge_tx"], errors="coerce"),
        }
        for name, series in raw_series.items():
            median, sd, mean = self.stats[name]
            output[f"z_{name}"] = (series.fillna(median) - mean) / sd
            if name in {"distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx"}:
                output[f"missing_{name}"] = series.isna().astype(float)
        categories = pd.get_dummies(
            df[["drach_subtype", "transcript_region"]].astype(str),
            prefix=["drach", "region"],
            dtype=float,
        )
        for name in self.categories:
            output[name] = categories[name] if name in categories else 0.0
        return output.astype(float)


def fit_model(df, builder, outcome, predictor, weights=None):
    design = builder.build(df, predictor)
    constant_columns = [
        name for name in design.columns
        if name != "Intercept" and design[name].nunique(dropna=False) <= 1
    ]
    if constant_columns:
        design = design.drop(columns=constant_columns)
    y = pd.to_numeric(df[outcome]).astype(float)
    model = sm.OLS(y, design) if weights is None else sm.WLS(y, design, weights=weights)
    ordinary = model.fit()
    clustered = model.fit(cov_type="cluster", cov_kwds={"groups": df["analysis_gene"], "use_correction": True})
    return ordinary, clustered, design, y


def cluster_bootstrap_beta(design, y, groups, predictor, replicates, seed):
    x = design.to_numpy(dtype=float)
    yv = y.to_numpy(dtype=float)
    genes, inverse = np.unique(groups.astype(str).to_numpy(), return_inverse=True)
    rng = np.random.default_rng(seed)
    betas = []
    predictor_index = list(design.columns).index(predictor)
    for _ in range(replicates):
        sampled = rng.integers(0, len(genes), size=len(genes))
        gene_weights = np.bincount(sampled, minlength=len(genes)).astype(float)
        weights = gene_weights[inverse]
        keep = weights > 0
        root_w = np.sqrt(weights[keep])
        try:
            beta = np.linalg.lstsq(x[keep] * root_w[:, None], yv[keep] * root_w, rcond=None)[0]
            betas.append(beta[predictor_index])
        except np.linalg.LinAlgError:
            continue
    return np.asarray(betas)


def association_row(name, df, builder, outcome, predictor, weights=None):
    ordinary, clustered, design, y = fit_model(df, builder, outcome, predictor, weights)
    beta = float(clustered.params[predictor])
    se = float(clustered.bse[predictor])
    ci = clustered.conf_int().loc[predictor]
    predictor_sd = float(pd.to_numeric(df[predictor]).std(ddof=1))
    row = {
        "analysis": name,
        "outcome": outcome,
        "predictor": predictor,
        "sites": len(df),
        "genes": df["analysis_gene"].nunique(),
        "beta_raw": beta,
        "se_cluster": se,
        "ci95_low_raw": float(ci.iloc[0]),
        "ci95_high_raw": float(ci.iloc[1]),
        "p_value": float(clustered.pvalues[predictor]),
        "predictor_sd": predictor_sd,
        "beta_per_sd": beta * predictor_sd,
        "ci95_low_per_sd": float(ci.iloc[0]) * predictor_sd,
        "ci95_high_per_sd": float(ci.iloc[1]) * predictor_sd,
        "r_squared": float(ordinary.rsquared),
    }
    return row, ordinary, clustered


def exclude_center_neighborhood_predictor(df, minimum_fraction):
    result = df.copy()
    means = []
    valid = []
    for encoded in result["reactivity_window_m50_p50_json"]:
        window = json.loads(encoded)
        left = [window[index] for index in range(40, 48) if window[index] is not None]
        right = [window[index] for index in range(53, 61) if window[index] is not None]
        okay = len(left) / 8 >= minimum_fraction and len(right) / 8 >= minimum_fraction
        valid.append(okay)
        means.append(float(np.mean(left + right)) if okay else np.nan)
    result["reactivity_mean_flank10_exclude_m2_p2"] = means
    return result[np.asarray(valid)].copy()


def missingness_rows(site_df):
    common = site_df[site_df["detected_both_replicates"].astype(str).eq("True")].copy()
    common["structure_window_group"] = np.where(
        common["main_window_valid"].astype(str).eq("True"), "valid_flank10", "insufficient_flank10"
    )
    metrics = ["combined_ratio", "combined_agcov", "gc_fraction_21", "transcript_position_fraction", "isoform_count"]
    rows = []
    groups = {name: part for name, part in common.groupby("structure_window_group")}
    for metric in metrics:
        numeric = {name: pd.to_numeric(part[metric], errors="coerce").dropna() for name, part in groups.items()}
        valid = numeric.get("valid_flank10", pd.Series(dtype=float))
        missing = numeric.get("insufficient_flank10", pd.Series(dtype=float))
        pooled = math.sqrt(((len(valid) - 1) * valid.var(ddof=1) + (len(missing) - 1) * missing.var(ddof=1)) / max(len(valid) + len(missing) - 2, 1))
        smd = (valid.mean() - missing.mean()) / pooled if pooled > 0 else math.nan
        for name, series in numeric.items():
            rows.append({"metric": metric, "group": name, "n_nonmissing": len(series), "mean": series.mean(), "median": series.median(), "sd": series.std(ddof=1), "standardized_mean_difference_valid_minus_insufficient": smd})
    return rows


def binned_trend_rows(df, bins=20):
    work = df.copy()
    work["reactivity_bin"] = pd.qcut(work["reactivity_mean_flank10"], q=bins, duplicates="drop")
    rows = []
    for index, (_, part) in enumerate(work.groupby("reactivity_bin", observed=True), start=1):
        rows.append({"bin": index, "sites": len(part), "genes": part["analysis_gene"].nunique(), "reactivity_mean": part["reactivity_mean_flank10"].mean(), "ratio_mean": part["combined_ratio"].mean(), "ratio_median": part["combined_ratio"].median(), "reactivity_min": part["reactivity_mean_flank10"].min(), "reactivity_max": part["reactivity_mean_flank10"].max()})
    return rows


def posthoc_region_diagnostic_rows(df, builder, outcome, predictor):
    rows = []
    for region, subset in df.groupby("transcript_region", observed=True):
        if len(subset) < 200 or subset["analysis_gene"].nunique() < 100:
            continue
        row, _, _ = association_row(f"transcript_region_{region}", subset.copy(), builder, outcome, predictor)
        row["analysis_role"] = "posthoc_fallacy_diagnostic"
        row["diagnostic_reason"] = "Simpson_paradox_direction_check"
        rows.append(row)
    if rows:
        adjusted = multipletests([row["p_value"] for row in rows], method="fdr_bh")[1]
        for row, qvalue in zip(rows, adjusted):
            row["p_adjust_bh_region_diagnostic"] = qvalue
    return rows


def positional_profile_rows(df, replicates, seed):
    work = df.copy()
    work["ratio_quartile"], boundaries = pd.qcut(
        work["combined_ratio"], q=4, labels=["Q1", "Q2", "Q3", "Q4"], retbins=True, duplicates="raise"
    )
    boundary_rows = [
        {"boundary": name, "combined_ratio": value}
        for name, value in zip(["minimum", "Q1_Q2", "Q2_Q3", "Q3_Q4", "maximum"], boundaries)
    ]
    output = []
    rng = np.random.default_rng(seed)
    for quartile, part in work.groupby("ratio_quartile", observed=True):
        matrix = np.asarray([[
            np.nan if value is None else float(value)
            for value in json.loads(encoded)
        ] for encoded in part["reactivity_window_m50_p50_json"]], dtype=float)
        genes, inverse = np.unique(part["analysis_gene"].astype(str).to_numpy(), return_inverse=True)
        gene_sums = np.zeros((len(genes), matrix.shape[1]), dtype=float)
        gene_counts = np.zeros((len(genes), matrix.shape[1]), dtype=float)
        for gene_index in range(len(genes)):
            selected = matrix[inverse == gene_index]
            gene_sums[gene_index] = np.nansum(selected, axis=0)
            gene_counts[gene_index] = np.sum(~np.isnan(selected), axis=0)
        boot = np.full((replicates, matrix.shape[1]), np.nan)
        for b in range(replicates):
            sampled = rng.integers(0, len(genes), size=len(genes))
            weights = np.bincount(sampled, minlength=len(genes))
            numerator = weights @ gene_sums
            denominator = weights @ gene_counts
            boot[b] = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)
        for column, relative_pos in enumerate(range(-50, 51)):
            observed = matrix[:, column]
            valid = observed[~np.isnan(observed)]
            output.append({
                "ratio_quartile": str(quartile),
                "relative_position": relative_pos,
                "mean_reactivity": np.mean(valid) if len(valid) else "",
                "ci95_low_pointwise": np.nanpercentile(boot[:, column], 2.5) if len(valid) else "",
                "ci95_high_pointwise": np.nanpercentile(boot[:, column], 97.5) if len(valid) else "",
                "valid_sites": len(valid),
                "valid_genes": part.loc[~np.isnan(observed), "analysis_gene"].nunique(),
                "quartile_sites": len(part),
                "quartile_genes": len(genes),
            })
    return output, boundary_rows


def main() -> None:
    setup_logging()
    started = datetime.now().isoformat(timespec="seconds")
    with CONFIG_PATH.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    df = load_csv(INPUT)
    site_df = load_csv(SITE_INPUT)
    if df["site_id"].duplicated().any():
        raise RuntimeError("Stage-09 primary population contains duplicate site IDs")
    builder = DesignBuilder(df)
    predictor = config["primary_model"]["predictor"]
    outcome = config["primary_model"]["outcome"]
    logging.info("fitting primary model: n=%d genes=%d", len(df), df["analysis_gene"].nunique())
    primary, ordinary, clustered = association_row("primary_adjusted_OLS", df, builder, outcome, predictor)
    design = builder.build(df, predictor)
    bootstrap = cluster_bootstrap_beta(
        design,
        pd.to_numeric(df[outcome]),
        df["analysis_gene"],
        predictor,
        int(config["primary_model"]["cluster_bootstrap_replicates"]),
        int(config["primary_model"]["bootstrap_seed"]),
    )
    sd_r = primary["predictor_sd"]
    primary.update({
        "bootstrap_successful_replicates": len(bootstrap),
        "bootstrap_ci95_low_raw": np.percentile(bootstrap, 2.5),
        "bootstrap_ci95_high_raw": np.percentile(bootstrap, 97.5),
        "bootstrap_ci95_low_per_sd": np.percentile(bootstrap * sd_r, 2.5),
        "bootstrap_ci95_high_per_sd": np.percentile(bootstrap * sd_r, 97.5),
    })
    coefficient_rows = []
    clustered_ci = clustered.conf_int()
    for name in design.columns:
        coefficient_rows.append({"term": name, "estimate": clustered.params[name], "cluster_se": clustered.bse[name], "ci95_low": clustered_ci.loc[name, 0], "ci95_high": clustered_ci.loc[name, 1], "p_value": clustered.pvalues[name]})

    sensitivity_specs = []
    excluded = exclude_center_neighborhood_predictor(
        df, float(config["sensitivity_family"]["exclude_minus2_plus2_minimum_valid_fraction_each_side"])
    )
    sensitivity_specs.append(("exclude_minus2_to_plus2", excluded, outcome, "reactivity_mean_flank10_exclude_m2_p2", None))
    sensitivity_specs.append(("complete_flank10_coverage", df[(df["coverage_up10"] == 1) & (df["coverage_down10"] == 1)].copy(), outcome, predictor, None))
    sensitivity_specs.append(("GLORI_replicate_1", df.copy(), "ratio_rep1", predictor, None))
    sensitivity_specs.append(("GLORI_replicate_2", df.copy(), "ratio_rep2", predictor, None))
    gene_sizes = df.groupby("analysis_gene")["site_id"].transform("count").to_numpy()
    sensitivity_specs.append(("equal_weight_per_gene", df.copy(), outcome, predictor, 1.0 / gene_sizes))
    sensitivity_specs.append(("single_isoform_sites", df[df["isoform_count"] == 1].copy(), outcome, predictor, None))
    sensitivity_specs.append(("canonical_DRACH_sites", df[df["drach_subtype"] != "non_DRACH"].copy(), outcome, predictor, None))
    sensitivity_rows = []
    for name, subset, sensitivity_outcome, sensitivity_predictor, weights in sensitivity_specs:
        logging.info("fitting sensitivity %s: n=%d", name, len(subset))
        row, _, _ = association_row(name, subset, builder, sensitivity_outcome, sensitivity_predictor, weights)
        sensitivity_rows.append(row)
    adjusted = multipletests([row["p_value"] for row in sensitivity_rows], method="fdr_bh")[1]
    for row, qvalue in zip(sensitivity_rows, adjusted):
        row["p_adjust_bh_sensitivity_family"] = qvalue

    residuals = np.asarray(ordinary.resid)
    fitted = np.asarray(ordinary.fittedvalues)
    bp = het_breuschpagan(residuals, design.to_numpy())
    reset = linear_reset(ordinary, power=2, use_f=True)
    spline_basis = dmatrix(
        "cr(x, df=3, constraints='center') - 1",
        {"x": pd.to_numeric(df[predictor]).to_numpy()},
        return_type="dataframe",
    )
    spline_design = design.drop(columns=[predictor]).copy()
    spline_names = []
    for index in range(spline_basis.shape[1]):
        name = f"spline_reactivity_{index + 1}"
        spline_design[name] = spline_basis.iloc[:, index].to_numpy()
        spline_names.append(name)
    spline_ordinary = sm.OLS(pd.to_numeric(df[outcome]), spline_design).fit()
    spline_clustered = sm.OLS(pd.to_numeric(df[outcome]), spline_design).fit(cov_type="cluster", cov_kwds={"groups": df["analysis_gene"], "use_correction": True})
    restriction = np.zeros((len(spline_names), spline_design.shape[1]))
    for index, name in enumerate(spline_names):
        restriction[index, list(spline_design.columns).index(name)] = 1
    spline_wald = spline_clustered.wald_test(restriction, scalar=True)
    diagnostics = [
        {"metric": "n_sites", "value": len(df), "interpretation": "Primary population"},
        {"metric": "n_genes", "value": df["analysis_gene"].nunique(), "interpretation": "Cluster count"},
        {"metric": "r_squared_linear", "value": ordinary.rsquared, "interpretation": "Descriptive fit, not causal variance explained"},
        {"metric": "adjusted_r_squared_linear", "value": ordinary.rsquared_adj, "interpretation": "Descriptive fit"},
        {"metric": "residual_skewness", "value": skew(residuals), "interpretation": "OLS residual diagnostic"},
        {"metric": "residual_excess_kurtosis", "value": kurtosis(residuals), "interpretation": "OLS residual diagnostic"},
        {"metric": "breusch_pagan_p", "value": bp[1], "interpretation": "Heteroscedasticity diagnostic; inference uses cluster-robust covariance"},
        {"metric": "ramsey_reset_p", "value": float(reset.pvalue), "interpretation": "Quadratic functional-form diagnostic"},
        {"metric": "predictions_below_zero", "value": int(np.sum(fitted < 0)), "interpretation": "Bounded-outcome limitation"},
        {"metric": "predictions_above_one", "value": int(np.sum(fitted > 1)), "interpretation": "Bounded-outcome limitation"},
        {"metric": "spline_r_squared", "value": spline_ordinary.rsquared, "interpretation": "Prespecified 3-df spline sensitivity"},
        {"metric": "spline_joint_wald_p", "value": float(spline_wald.pvalue), "interpretation": "Joint association of spline terms, not a linearity test"},
        {"metric": "spline_minus_linear_AIC", "value": spline_ordinary.aic - ordinary.aic, "interpretation": "Negative favors spline descriptively"},
    ]

    profile_rows, quartile_boundaries = positional_profile_rows(
        df,
        int(config["descriptive_profiles"]["bootstrap_replicates"]),
        int(config["descriptive_profiles"]["bootstrap_seed"]),
    )
    missing_rows = missingness_rows(site_df)
    trend_rows = binned_trend_rows(df)
    region_diagnostic_rows = posthoc_region_diagnostic_rows(df, builder, outcome, predictor)

    TABLES.mkdir(parents=True, exist_ok=True)
    write_rows(TABLES / "09b_hela_primary_association.csv", [primary])
    write_rows(TABLES / "09b_hela_primary_model_coefficients.csv", coefficient_rows)
    write_rows(TABLES / "09b_hela_sensitivity_associations.csv", sensitivity_rows)
    write_rows(TABLES / "09b_hela_model_diagnostics.csv", diagnostics)
    write_rows(TABLES / "09b_hela_missingness_bias_summary.csv", missing_rows)
    write_rows(TABLES / "09b_hela_binned_reactivity_ratio_trend.csv", trend_rows)
    write_rows(TABLES / "09b_hela_reactivity_profile_by_ratio_quartile.csv", profile_rows)
    write_rows(TABLES / "09b_hela_ratio_quartile_boundaries.csv", quartile_boundaries)
    write_rows(TABLES / "09b_hela_posthoc_region_diagnostic.csv", region_diagnostic_rows)
    status = "PASS" if len(bootstrap) >= 950 and len(sensitivity_rows) == 7 else "FAIL"
    write_rows(
        RESULTS / "09b_hela_association_run_status.csv",
        [{"status": status, "started_at": started, "finished_at": datetime.now().isoformat(timespec="seconds"), "sites": len(df), "genes": df["analysis_gene"].nunique(), "bootstrap_successful": len(bootstrap), "sensitivity_models": len(sensitivity_rows), "unavailable_covariates": json.dumps(config["primary_model"].get("unavailable_covariates", []), sort_keys=True), "config": str(CONFIG_PATH.relative_to(ROOT)).replace("\\", "/")}],
    )
    logging.info("stage 09B analysis finished: status=%s beta_per_sd=%.6g bootstrap_ci=[%.6g, %.6g]", status, primary["beta_per_sd"], primary["bootstrap_ci95_low_per_sd"], primary["bootstrap_ci95_high_per_sd"])
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
