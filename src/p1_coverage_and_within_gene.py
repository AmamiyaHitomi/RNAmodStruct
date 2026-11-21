"""Post hoc P1 coverage, weighting, within-gene, and spline analyses."""

from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy import build_design_matrices, dmatrix
from scipy.linalg import qr
from scipy.special import expit
from scipy.stats import norm

from p1_reactivity_partial_r2 import MODELS
from presubmission_stage3_association_context import categorical_balance, numerical_balance
from replay_presubmission_stage4_models import sha256


PREDICTOR = "reactivity_mean_flank10"
OUTCOME = "combined_ratio"
GENE = "analysis_gene"
ROOT = Path(__file__).resolve().parents[1]
Z95 = norm.ppf(.975)
NUMERIC_BALANCE = (
    OUTCOME, "combined_agcov", "gc_fraction_21", "transcript_position_fraction",
    "isoform_count", "distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx",
)
CATEGORICAL = ("drach_subtype", "transcript_region")
WINDOWS = ((10, 7), (10, 6), (10, 8), (10, 10), (5, 4), (15, 11), (20, 14))


def design_matrix(stage, frame: pd.DataFrame) -> pd.DataFrame:
    design = stage.DesignBuilder(frame).build(frame, PREDICTOR)
    keep = design.nunique().gt(1) | design.columns.isin(["Intercept", PREDICTOR])
    return design.loc[:, keep]


def coefficient_row(fit, frame, sd: float, **labels) -> dict:
    low, high = fit.conf_int().loc[PREDICTOR]
    return {
        **labels, "sites": len(frame), "genes": frame[GENE].nunique(), "reference_sd": sd,
        "beta_raw": fit.params[PREDICTOR], "beta_pp_per_reference_sd": fit.params[PREDICTOR] * sd * 100,
        "ci95_low_pp": low * sd * 100, "ci95_high_pp": high * sd * 100,
    }


def coverage_windows(stage, mapped, primary, window, sd):
    rows = []
    for flank, threshold in WINDOWS:
        left = window[:, 50-flank:50]
        right = window[:, 51:51+flank]
        left_n = np.isfinite(left).sum(axis=1)
        right_n = np.isfinite(right).sum(axis=1)
        physical = mapped.transcript_pos_1based.gt(flank) & (
            mapped.transcript_pos_1based + flank <= mapped.transcript_length
        )
        valid = physical & (left_n >= threshold) & (right_n >= threshold)
        trial = mapped.loc[valid].copy()
        x = (np.nansum(left[valid], axis=1) + np.nansum(right[valid], axis=1)) / (left_n[valid] + right_n[valid])
        trial[PREDICTOR] = x
        trial["coverage_flank10"] = (left_n[valid] + right_n[valid]) / (2 * flank)
        if (flank, threshold) == (10, 7):
            np.testing.assert_array_equal(np.sort(primary.site_id), np.sort(trial.site_id))
            ordered = trial.set_index("site_id").loc[primary.site_id]
            np.testing.assert_allclose(ordered[PREDICTOR], primary[PREDICTOR], atol=1e-7, rtol=0)
            np.testing.assert_allclose(ordered["coverage_flank10"], primary["coverage_flank10"], atol=1e-7, rtol=0)
            trial = primary
        _, fit, _, _ = stage.fit_model(trial, stage.DesignBuilder(trial), OUTCOME, PREDICTOR)
        rows.append(coefficient_row(fit, trial, sd, flank_nt=flank, min_valid_per_side=threshold, comparison="eligible_under_rule"))

        same = primary.loc[primary.site_id.isin(trial.site_id)].copy()
        alternative = same.copy()
        alternative[PREDICTOR] = trial.set_index("site_id").loc[same.site_id, PREDICTOR].to_numpy()
        builder = stage.DesignBuilder(same)
        for label, frame in (("common_sites_primary_window", same), ("common_sites_alternative_window", alternative)):
            _, fit, _, _ = stage.fit_model(frame, builder, OUTCOME, PREDICTOR)
            rows.append(coefficient_row(fit, frame, sd, flank_nt=flank, min_valid_per_side=threshold, comparison=label))
    return rows


def coverage_balance(mapped):
    data = mapped.assign(coverage_group=np.where(mapped.main_analysis_included, "eligible", "ineligible"))
    rows = []
    for column in NUMERIC_BALANCE:
        rows.extend(numerical_balance(data, column, np.log1p if column == "combined_agcov" else None))
    for column in CATEGORICAL:
        rows.extend(categorical_balance(data, column))
    return rows


def selection_design(target):
    features = pd.DataFrame(index=target.index)
    for column in NUMERIC_BALANCE:
        value = target[column].astype(float)
        if column == "combined_agcov":
            value = np.log1p(value)
        if value.isna().any():
            features[f"{column}_missing"] = value.isna().astype(float)
        filled = value.fillna(value.median())
        if filled.std() > 0:
            features[column] = (filled - filled.mean()) / filled.std()
    other = pd.get_dummies(target[list(CATEGORICAL)].fillna("<missing>"), dtype=float, drop_first=True)
    features = pd.concat([features, other], axis=1)
    features = features.loc[:, features.nunique().gt(1)]
    return sm.add_constant(features, has_constant="add")


def probability_model(x, s, frequency):
    """Ridge-stabilized logistic IRLS, with no penalty on the intercept."""
    fitted = np.zeros(x.shape[1])
    penalty = np.eye(x.shape[1])
    penalty[0, 0] = 0
    for _ in range(60):
        p = expit(x @ fitted)
        information = x.T @ (x * (frequency * p * (1 - p))[:, None]) + penalty
        step = np.linalg.solve(information, x.T @ (frequency * (s - p)) - penalty @ fitted)
        fitted += step
        if np.max(np.abs(step)) < 1e-7:
            return expit(x @ fitted)
    raise RuntimeError("Selection probability IRLS failed to converge")


def reweight(stage, mapped, primary, sd, replicates, seed):
    physically_possible = mapped.transcript_pos_1based.gt(10) & (
        mapped.transcript_pos_1based + 10 <= mapped.transcript_length
    )
    target = mapped.loc[physically_possible].copy()
    selected = target.main_analysis_included.to_numpy(dtype=bool)
    frame = primary.set_index("site_id").loc[target.loc[selected, "site_id"]].copy()
    groups, names = pd.factorize(target[GENE], sort=True)
    x_frame = selection_design(target)
    x = x_frame.to_numpy()
    d = design_matrix(stage, frame).drop(columns="z_coverage_flank10")
    # Retain the predictor and remove exact dependent nuisance columns.
    nuisance = d.drop(columns=PREDICTOR)
    _, r, order = qr(nuisance.to_numpy(), pivoting=True, mode="economic")
    rank = np.linalg.matrix_rank(r)
    d = d[[PREDICTOR, *nuisance.columns[order[:rank]]]]
    design = d.to_numpy()
    y = frame[OUTCOME].to_numpy()
    scheme = ("unweighted", "ipw_untrimmed", "ipw_p01_p99")

    def weights_for(frequency):
        probabilities = probability_model(x, selected, frequency)
        observed = frequency[selected]
        stabilized = (np.dot(frequency, selected) / frequency.sum()) / probabilities[selected]
        repeated = np.repeat(stabilized, observed.astype(int))
        q1, q99 = np.quantile(repeated, [.01, .99])
        weights = np.stack([np.ones_like(stabilized), stabilized, np.clip(stabilized, q1, q99)])
        return probabilities, weights, (q1, q99)

    def coefficients(weights, frequency):
        values = []
        for w in weights * frequency[selected]:
            root_w = np.sqrt(w)
            weighted_design = design * root_w[:, None]
            solution, _, rank, _ = np.linalg.lstsq(weighted_design, root_w * y, rcond=None)
            if rank < design.shape[1]:
                nuisance_rank = np.linalg.matrix_rank(weighted_design[:, 1:])
                if rank <= nuisance_rank:
                    raise RuntimeError("Reactivity coefficient is unidentified in a gene bootstrap sample")
            values.append(solution[0] * sd * 100)
        return values

    ones = np.ones(len(target), dtype=int)
    probabilities, weights, limits = weights_for(ones)
    effects = coefficients(weights, ones)
    rng = np.random.default_rng(seed)
    boot = np.empty((replicates, len(scheme)))
    for i in range(replicates):
        multiplicity = np.bincount(rng.integers(0, len(names), len(names)), minlength=len(names))
        frequency = multiplicity[groups]
        _, w, _ = weights_for(frequency)
        boot[i] = coefficients(w, frequency)
    bounds = np.quantile(boot, [.025, .975], axis=0)
    differences = np.asarray(effects) - effects[0]
    difference_ci = np.quantile(boot - boot[:, [0]], [.025, .975], axis=0)
    rows = []
    for j, (name, w) in enumerate(zip(scheme, weights)):
        rows.append(dict(method=name, target_sites=len(target), target_genes=len(names),
                         excluded_structurally_ineligible=int((~physically_possible).sum()),
                         sites=len(frame), genes=frame[GENE].nunique(), reference_sd=sd,
                         beta_pp_per_reference_sd=effects[j], ci95_low_pp=bounds[0, j], ci95_high_pp=bounds[1, j],
                         difference_from_unweighted_pp=differences[j], difference_ci95_low_pp=difference_ci[0, j],
                         difference_ci95_high_pp=difference_ci[1, j], bootstrap_replicates=replicates, bootstrap_seed=seed,
                         effective_sites=w.sum()**2 / (w @ w), min_weight=w.min(), max_weight=w.max(),
                         median_weight=np.median(w), trim_low=limits[0], trim_high=limits[1]))
    propensity = []
    for label, mask in (("retained", selected), ("excluded", ~selected)):
        p = probabilities[mask]
        quantiles = np.quantile(p, [0, .01, .05, .5, .95, .99, 1])
        propensity.append(dict(group=label, sites=mask.sum(), **dict(zip(("min", "p01", "p05", "median", "p95", "p99", "max"), quantiles))))
    balance = []
    for name, w in zip(scheme, weights):
        for column, values in x_frame.drop(columns="const").items():
            pooled = values.std(ddof=1)
            balance.append(dict(method=name, variable=column,
                                target_mean=values.mean(), weighted_retained_mean=np.average(values[selected], weights=w),
                                standardized_difference=(np.average(values[selected], weights=w) - values.mean()) / pooled))
    return rows, propensity, balance


def within_gene(stage, primary, sd):
    summaries = primary.groupby(GENE)[PREDICTOR].agg(["size", "min", "max"])
    rows = []
    for threshold in (2, 3, 5):
        eligible = summaries.index[summaries["size"].ge(threshold) & summaries["max"].gt(summaries["min"])]
        frame = primary.loc[primary[GENE].isin(eligible)]
        nominal = design_matrix(stage, frame)
        grouped = frame[GENE]
        adjusted = sm.OLS(frame[OUTCOME], nominal).fit(cov_type="cluster", cov_kwds={"groups": grouped})
        d = nominal - nominal.groupby(grouped).transform("mean")
        d = d.loc[:, d.abs().max().gt(1e-12)]
        y = frame[OUTCOME] - frame.groupby(GENE)[OUTCOME].transform("mean")
        fe = sm.OLS(y, d).fit(cov_type="cluster", cov_kwds={"groups": grouped, "use_correction": False})
        n, g = len(frame), grouped.nunique()
        rank = np.linalg.matrix_rank(d.to_numpy())
        correction = g / (g - 1) * (n - 1) / (n - g - rank)
        half = Z95 * fe.bse[PREDICTOR] * np.sqrt(correction)
        beta = fe.params[PREDICTOR]
        row = coefficient_row(adjusted, frame, sd, min_sites_per_gene=threshold, method="pooled_same_sites")
        rows.append(row)
        rows.append({**row, "method": "gene_fixed_effects", "beta_raw": beta,
                     "beta_pp_per_reference_sd": beta * sd * 100,
                     "ci95_low_pp": (beta-half) * sd * 100, "ci95_high_pp": (beta+half) * sd * 100})
    return rows


def splines(stage, primary):
    x = primary[PREDICTOR].to_numpy()
    design = design_matrix(stage, primary)
    basis = dmatrix("cr(x, df=3, constraints='center') - 1", {"x": x}, return_type="dataframe")
    model = np.column_stack([design.drop(columns=PREDICTOR), basis])
    fit = sm.OLS(primary[OUTCOME], model).fit(
        cov_type="cluster", cov_kwds={"groups": primary[GENE], "use_correction": True}
    )
    grid = np.linspace(*np.quantile(x, [.01, .99]), 140)
    extrapolate_matrix = np.asarray(build_design_matrices([basis.design_info], {"x": grid})[0])
    at_grid = np.column_stack([np.tile(model.mean(axis=0)[:-3], (len(grid), 1)), extrapolate_matrix])
    pred = at_grid @ fit.params.to_numpy()
    se = np.sqrt(np.einsum("ij,jk,ik->i", at_grid, fit.cov_params(), at_grid))
    counts, edges = np.histogram(x, bins=np.linspace(0, 1, 46))
    if counts.sum() != len(x):
        raise ValueError("Histogram does not account for every site")
    return (
        pd.DataFrame({"reactivity_flank10": grid, "adjusted_ratio": pred, "ci95_low": pred-Z95*se, "ci95_high": pred+Z95*se}),
        pd.DataFrame({"left": edges[:-1], "right": edges[1:], "sites": counts}),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New output directory")
    parser.add_argument("--bootstrap", type=int, default=500)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    parts = {}
    provenance = {}
    for index, (cell, module, _) in enumerate(MODELS):
        stage = importlib.import_module(module)
        primary = stage.load_csv(stage.INPUT)
        site = stage.load_csv(stage.SITE_INPUT)
        mapped = site.loc[site.detected_both_replicates].copy()
        np.testing.assert_array_equal(np.sort(mapped.loc[mapped.main_analysis_included, "site_id"]), np.sort(primary.site_id))
        windows = np.asarray([json.loads(encoded) for encoded in mapped.reactivity_window_m50_p50_json], dtype=float)
        sd = primary[PREDICTOR].std(ddof=1)
        print(f"{cell}: windows, balance, within-gene models and curve", flush=True)
        result = {
            "window_sensitivity": coverage_windows(stage, mapped, primary, windows, sd),
            "coverage_balance": coverage_balance(mapped),
            "gene_fixed_effects": within_gene(stage, primary, sd),
        }
        result["adjusted_spline"], result["reactivity_histogram"] = splines(stage, primary)
        attrition = pd.read_csv(ROOT / "results/tables" / ("06_hek293t_attrition_summary.csv" if index == 0 else "09_hela_attrition_summary.csv")).head(4)
        attrition["fraction_from_initial"] = attrition.sites / attrition.sites.iloc[0]
        attrition["fraction_from_previous"] = attrition.sites / attrition.sites.shift()
        result["attrition"] = attrition.drop(columns="retained_from_previous")
        print(f"{cell}: re-estimating weights in {args.bootstrap} gene bootstraps", flush=True)
        result["ipw_estimates"], result["ipw_overlap"], result["ipw_balance"] = reweight(stage, mapped, primary, sd, args.bootstrap, 2026112102 + index)
        for name, rows in result.items():
            parts.setdefault(name, []).append(pd.DataFrame(rows).assign(cell_line=cell))
        provenance[cell] = {path.name: sha256(path) for path in (stage.INPUT, stage.SITE_INPUT, stage.CONFIG_PATH)}
    for name, frames in parts.items():
        pd.concat(frames, ignore_index=True).to_csv(args.output / f"{name}.csv", index=False, mode="x")
    (args.output / "input_sha256.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
