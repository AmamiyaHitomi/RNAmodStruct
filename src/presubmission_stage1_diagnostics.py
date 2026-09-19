"""Read-only frozen-prediction diagnostics for main_v2 presubmission stage 1.

The script fits no model and creates a new result directory only.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "presubmission_stage1"
PHASE0 = ROOT / "metadata" / "manifests" / "presubmission_phase0"
HEK_PRED = ROOT / "results/tables/08_hek293t_test_predictions.csv"
HELA_PRED = ROOT / "results/tables/09c_hela_transfer_predictions.csv"
HEK_METRICS = ROOT / "results/tables/08_hek293t_test_metrics.csv"
HELA_METRICS = ROOT / "results/tables/09c_hela_transfer_metrics.csv"
HEK_ASSIGNMENTS = ROOT / "results/tables/08_hek293t_group_assignments.csv"
HEK_DATA = ROOT / "data/final/06_hek293t_main_analysis_dataset.csv.gz"
MODELS = ["M0", "M1", "M2", "M3", "M4", "C_dev_median"]
COHORTS = [
    "HEK_reused_test",
    "HeLa_all_qualifying",
    "HeLa_development_overlap_excluded",
    "HeLa_hek_main_4409_overlap_excluded",
]
COLORS = {
    "M0": "#2166ac", "M1": "#b2182b", "M2": "#ef8a62",
    "M3": "#4d9221", "M4": "#9970ab", "C_dev_median": "#333333",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def add_constant(frame: pd.DataFrame, value: float) -> pd.DataFrame:
    result = frame.copy()
    result["prediction_C_dev_median"] = value
    return result


def metrics(cohort: str, model: str, frame: pd.DataFrame) -> dict:
    y = frame["combined_ratio"].to_numpy(dtype=float)
    pred = frame[f"prediction_{model}"].to_numpy(dtype=float)
    if not np.isfinite(y).all() or not np.isfinite(pred).all():
        raise ValueError(f"Nonfinite label/prediction in {cohort}/{model}")
    residual = pred - y
    absolute = np.abs(residual)
    rho = None if np.ptp(pred) == 0 else float(spearmanr(y, pred).statistic)
    return {
        "cohort": cohort,
        "model": model,
        "sites": len(frame),
        "genes": int(frame["analysis_gene"].nunique()),
        "joint_groups": int(frame["joint_group"].nunique()),
        "mae": float(mean_absolute_error(y, pred)),
        "rmse": float(np.sqrt(mean_squared_error(y, pred))),
        "spearman": rho,
        "r2": float(r2_score(y, pred)),
        "mean_signed_error_prediction_minus_observed": float(residual.mean()),
        "median_absolute_error": float(np.median(absolute)),
        "p90_absolute_error": float(np.quantile(absolute, 0.9)),
        "observed_mean": float(y.mean()),
        "observed_sd": float(y.std(ddof=0)),
        "prediction_mean": float(pred.mean()),
        "prediction_sd": float(pred.std(ddof=0)),
        "prediction_min": float(pred.min()),
        "prediction_q05": float(np.quantile(pred, 0.05)),
        "prediction_median": float(np.median(pred)),
        "prediction_q95": float(np.quantile(pred, 0.95)),
        "prediction_max": float(pred.max()),
        "at_clip_0_count": int((pred == 0).sum()),
        "at_clip_1_count": int((pred == 1).sum()),
        "at_either_clip_fraction": float(((pred == 0) | (pred == 1)).mean()),
    }


def calibration_rows(cohort: str, model: str, frame: pd.DataFrame, edges: np.ndarray) -> list[dict]:
    pred = frame[f"prediction_{model}"].to_numpy(dtype=float)
    y = frame["combined_ratio"].to_numpy(dtype=float)
    # Internal cut points come only from HEK development labels. These are
    # fixed prediction-value bins, not evaluation-outcome quantiles.
    bins = np.concatenate(([-np.inf], edges[1:-1], [np.inf]))
    index = np.searchsorted(bins[1:-1], pred, side="right")
    rows = []
    for bin_index in range(len(bins) - 1):
        mask = index == bin_index
        count = int(mask.sum())
        rows.append({
            "cohort": cohort,
            "model": model,
            "bin_index": bin_index + 1,
            "development_label_quantile_low": float(edges[bin_index]),
            "development_label_quantile_high": float(edges[bin_index + 1]),
            "bin_count": count,
            "mean_prediction": float(pred[mask].mean()) if count else None,
            "mean_observed": float(y[mask].mean()) if count else None,
            "mean_signed_error_prediction_minus_observed": float((pred[mask] - y[mask]).mean()) if count else None,
        })
    if sum(row["bin_count"] for row in rows) != len(frame):
        raise AssertionError("Calibration bins did not cover all sites")
    return rows


def make_figures(metric_frame: pd.DataFrame, calibration_frame: pd.DataFrame, site_frame: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9), sharex=True, sharey=True)
    for ax, cohort in zip(axes.flat, COHORTS):
        ax.plot([0, 1], [0, 1], color="#aaaaaa", linestyle="--", linewidth=1, label="perfect")
        for model in MODELS:
            part = calibration_frame[
                calibration_frame["cohort"].eq(cohort) & calibration_frame["model"].eq(model)
                & calibration_frame["bin_count"].gt(0)
            ]
            ax.plot(part["mean_prediction"], part["mean_observed"], marker="o", markersize=3,
                    linewidth=1.3, color=COLORS[model], label=model)
        ax.set(title=cohort.replace("_", " "), xlim=(0, 1), ylim=(0, 1))
        ax.grid(alpha=0.2)
        ax.set_xlabel("Mean predicted modification proportion")
        ax.set_ylabel("Mean observed modification proportion")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, frameon=False)
    fig.suptitle("Frozen prediction calibration; bins fixed by HEK development labels")
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(OUT / "calibration.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), sharex=True, sharey=True)
    for ax, cohort in zip(axes.flat, COHORTS):
        part = site_frame[site_frame["cohort"].eq(cohort)]
        for model in MODELS[:-1]:
            ax.hist(part[f"prediction_{model}"], bins=np.linspace(0, 1, 41), histtype="step",
                    density=True, linewidth=1.4, color=COLORS[model], label=model)
        ax.axvline(float(part["prediction_C_dev_median"].iloc[0]), color=COLORS["C_dev_median"],
                   linestyle="--", linewidth=1.4, label="C dev median")
        ax.set(title=cohort.replace("_", " "), xlim=(0, 1))
        ax.grid(alpha=0.2)
        ax.set_xlabel("Predicted modification proportion")
        ax.set_ylabel("Density")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False)
    fig.suptitle("Frozen prediction distributions")
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(OUT / "prediction_distributions.png", dpi=180)
    plt.close(fig)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing result directory: {OUT}")
    assignments = pd.read_csv(HEK_ASSIGNMENTS)
    hek_labels = pd.read_csv(HEK_DATA, usecols=["site_id", "combined_ratio"])
    development = assignments.loc[assignments["split"].eq("development"), ["site_id"]].merge(
        hek_labels, on="site_id", validate="one_to_one")
    if len(development) != 3271:
        raise ValueError("Unexpected HEK development size")
    constant = float(development["combined_ratio"].median())
    edges = np.quantile(development["combined_ratio"].to_numpy(dtype=float), np.linspace(0, 1, 11))
    if len(np.unique(edges)) != len(edges):
        raise ValueError("Development label quantile bins are not unique")

    hek = add_constant(pd.read_csv(HEK_PRED), constant)
    hela = add_constant(pd.read_csv(HELA_PRED), constant)
    if len(hek) != 825 or len(hela) != 24960:
        raise ValueError("Frozen prediction table counts changed")
    test = assignments[assignments["split"].eq("test")]
    if set(hek["site_id"]) != set(test["site_id"]):
        raise ValueError("HEK frozen predictions do not match fixed test IDs")
    group_check = hek[["site_id", "joint_group"]].merge(
        test[["site_id", "joint_group"]], on="site_id", validate="one_to_one",
        suffixes=("_prediction", "_assignment"),
    )
    if not group_check.eval("joint_group_prediction == joint_group_assignment").all():
        raise ValueError("HEK test joint groups changed")

    freeze = json.loads((PHASE0 / "freeze_manifest.json").read_text(encoding="utf-8"))
    flags = pd.read_csv(PHASE0 / "hela_site_cohort_flags.csv")
    if set(hela["site_id"]) != set(flags["site_id"]):
        raise ValueError("Frozen HeLa IDs disagree with phase 0")
    hela = hela.merge(flags[["site_id", "joint_group", "development_overlap_excluded",
                              "hek_main_4409_overlap_excluded"]], on="site_id",
                      validate="one_to_one", suffixes=("", "_phase0"))
    if not hela["joint_group"].eq(hela["joint_group_phase0"]).all():
        raise ValueError("Frozen HeLa joint groups disagree with phase 0")
    hela.drop(columns="joint_group_phase0", inplace=True)
    sets = {
        "HEK_reused_test": hek,
        "HeLa_all_qualifying": hela,
        "HeLa_development_overlap_excluded": hela[hela["development_overlap_excluded"]].copy(),
        "HeLa_hek_main_4409_overlap_excluded": hela[hela["hek_main_4409_overlap_excluded"]].copy(),
    }
    for cohort, key in [("HeLa_all_qualifying", "all_qualifying"),
                        ("HeLa_development_overlap_excluded", "development_overlap_excluded"),
                        ("HeLa_hek_main_4409_overlap_excluded", "hek_main_4409_overlap_excluded")]:
        ids = sorted(sets[cohort]["site_id"].astype(str))
        digest = hashlib.sha256("".join(f"{item}\n" for item in ids).encode()).hexdigest()
        if digest != freeze["cohorts"][key]["site_ids_sha256"]:
            raise ValueError(f"Phase 0 ID hash changed for {cohort}")

    metric_rows, calibration, sites = [], [], []
    for cohort, frame in sets.items():
        for model in MODELS:
            metric_rows.append(metrics(cohort, model, frame))
            calibration.extend(calibration_rows(cohort, model, frame, edges))
        subset = frame[["site_id", "analysis_gene", "joint_group", "combined_ratio"] +
                       [f"prediction_{model}" for model in MODELS]].copy()
        subset.insert(0, "cohort", cohort)
        sites.append(subset)
    metric_frame = pd.DataFrame(metric_rows)
    calibration_frame = pd.DataFrame(calibration)
    site_frame = pd.concat(sites, ignore_index=True)

    checks = []
    published_hek = pd.read_csv(HEK_METRICS)
    published_hela = pd.read_csv(HELA_METRICS)
    for cohort, subset in [("HEK_reused_test", None), ("HeLa_all_qualifying", "all_qualifying"),
                           ("HeLa_development_overlap_excluded", "strict")]:
        published = published_hek if subset is None else published_hela[published_hela["subset"].eq(subset)]
        for model in MODELS[:-1]:
            actual = metric_frame[metric_frame["cohort"].eq(cohort) & metric_frame["model"].eq(model)].iloc[0]
            expected = published[published["model"].eq(model)].iloc[0]
            for field in ["mae", "rmse", "spearman", "r2"]:
                difference = float(actual[field] - expected[field])
                checks.append({"cohort": cohort, "model": model, "metric": field,
                               "recomputed": actual[field], "published": expected[field],
                               "absolute_difference": abs(difference), "pass_at_1e_minus_8": abs(difference) <= 1e-8})
            if int(actual["sites"]) != int(expected["test_sites"] if subset is None else expected["sites"]):
                raise ValueError(f"Published site count changed for {cohort}/{model}")
    checks_frame = pd.DataFrame(checks)
    if not checks_frame["pass_at_1e_minus_8"].all():
        raise ValueError("A frozen prediction metric did not reproduce")

    # HeLa medians use evaluation labels. Keep them separate from deployable baselines.
    oracle = []
    for cohort in COHORTS[1:]:
        frame = sets[cohort]
        value = float(frame["combined_ratio"].median())
        oracle.append({"cohort": cohort, "status": "non_deployable_oracle_uses_evaluation_labels",
                       "HeLa_evaluation_label_median": value,
                       "oracle_mae": float(np.abs(frame["combined_ratio"] - value).mean())})

    OUT.mkdir(parents=False, exist_ok=False)
    metric_frame.to_csv(OUT / "diagnostic_metrics.csv", index=False, lineterminator="\n")
    calibration_frame.to_csv(OUT / "calibration.csv", index=False, lineterminator="\n")
    site_frame.to_csv(OUT / "site_predictions_with_constant.csv", index=False, lineterminator="\n")
    checks_frame.to_csv(OUT / "published_metric_reproduction.csv", index=False, lineterminator="\n")
    pd.DataFrame(oracle).to_csv(OUT / "oracle_diagnostics.csv", index=False, lineterminator="\n")
    make_figures(metric_frame, calibration_frame, site_frame)
    provenance = {
        "fitted_models": False,
        "HEK_development_label_median": constant,
        "calibration_edges_from_HEK_development_labels": edges.tolist(),
        "calibration_rule": "bin prediction values by fixed development-label decile boundaries; report mean prediction versus mean observed within each bin",
        "signed_error_rule": "prediction_minus_observed",
        "frozen_input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in
                                [HEK_PRED, HELA_PRED, HEK_METRICS, HELA_METRICS, HEK_ASSIGNMENTS, HEK_DATA,
                                 PHASE0 / "freeze_manifest.json", PHASE0 / "hela_site_cohort_flags.csv"]},
        "published_metric_tolerance": 1e-8,
        "checks_count": len(checks_frame),
        "checks_passed": int(checks_frame["pass_at_1e_minus_8"].sum()),
        "output_sha256": {p.name: sha256(p) for p in OUT.iterdir() if p.is_file()},
    }
    with (OUT / "provenance.json").open("x", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"constant": constant, "checks": len(checks_frame),
                      "max_metric_abs_difference": float(checks_frame["absolute_difference"].max()),
                      "mae": metric_frame.pivot(index="cohort", columns="model", values="mae").to_dict("index")}, indent=2))


if __name__ == "__main__":
    main()
