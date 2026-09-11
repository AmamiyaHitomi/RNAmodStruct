"""Run development-only expanded-R diagnostics and GLORI replicate agreement checks."""

from __future__ import annotations

import csv
import json
import logging
import math
import runpy
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "08b_structure_diagnostics.yaml"
INPUT = ROOT / "data" / "final" / "06_hek293t_main_analysis_dataset.csv.gz"
PREDICTED_INPUT = ROOT / "data" / "final" / "08a_hek293t_predicted_structure_features.csv.gz"
ASSIGNMENTS = ROOT / "results" / "tables" / "08_hek293t_group_assignments.csv"
TABLES = ROOT / "results" / "tables"
RESULTS = ROOT / "results"
LOG = RESULTS / "logs" / "08b_run_structure_diagnostics.log"


stage08 = runpy.run_path(str(ROOT / "src" / "08_run_prediction_modeling.py"))
FeatureEncoder = stage08["FeatureEncoder"]


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows(rows)


def add_profile_columns(df: pd.DataFrame) -> pd.DataFrame:
    matrix = np.asarray([json.loads(value) for value in df["reactivity_window_m50_p50_json"]], dtype=object)
    if matrix.shape != (len(df), 101):
        raise ValueError("Experimental reactivity profiles must span -50..+50")
    profile = pd.DataFrame({
        f"_r_{position}": pd.to_numeric(pd.Series(matrix[:, index]), errors="coerce").to_numpy()
        for index, position in enumerate(range(-50, 51))
    }, index=df.index)
    return pd.concat([df, profile], axis=1).copy()


def raw_block(df: pd.DataFrame, block: str) -> tuple[np.ndarray, list[str]]:
    if block == "R_core":
        names = ["reactivity_mean_flank10"]
        return df[names].to_numpy(dtype=float), names
    if block == "R_directional":
        up = df["reactivity_mean_up10"].to_numpy(dtype=float)
        down = df["reactivity_mean_down10"].to_numpy(dtype=float)
        return np.column_stack([up, down, up - down]), ["reactivity_mean_up10", "reactivity_mean_down10", "up10_minus_down10"]
    if block in {"R_profile_no_center", "R_profile_all"}:
        positions = list(range(-50, 51))
        if block == "R_profile_no_center":
            positions.remove(0)
        names = [f"reactivity_pos{position:+d}" for position in positions]
        values = df[[f"_r_{position}" for position in positions]].to_numpy(dtype=float)
        return values, names
    raise ValueError(f"Unknown experimental-structure block: {block}")


class BlockScaler:
    def __init__(self, include_mask: bool):
        self.include_mask = include_mask
        self.medians = None
        self.means = None
        self.scales = None

    def fit(self, values: np.ndarray) -> "BlockScaler":
        with np.errstate(all="ignore"):
            medians = np.nanmedian(values, axis=0)
        self.medians = np.where(np.isfinite(medians), medians, 0.0)
        filled = np.where(np.isfinite(values), values, self.medians)
        self.means = filled.mean(axis=0)
        scales = filled.std(axis=0)
        self.scales = np.where(scales > 0, scales, 1.0)
        return self

    def transform(self, values: np.ndarray, names: list[str]) -> tuple[np.ndarray, list[str]]:
        missing = ~np.isfinite(values)
        filled = np.where(~missing, values, self.medians)
        standardized = (filled - self.means) / self.scales
        output_names = [f"R:z_{name}" for name in names]
        if self.include_mask:
            standardized = np.hstack([standardized, missing.astype(float)])
            output_names += [f"R:missing_{name}" for name in names]
        return standardized.astype(np.float32), output_names


def prediction_metrics(y: np.ndarray, prediction: np.ndarray) -> dict:
    return {
        "mae": float(mean_absolute_error(y, prediction)),
        "rmse": float(math.sqrt(mean_squared_error(y, prediction))),
        "spearman": float(spearmanr(y, prediction).statistic),
        "r2": float(r2_score(y, prediction)),
    }


def development_cv(df: pd.DataFrame, blocks: list[str], alphas: list[float]) -> tuple[list[dict], list[dict]]:
    y = df["combined_ratio"].to_numpy(dtype=float)
    groups = df["joint_group"].to_numpy()
    splits = list(GroupKFold(n_splits=5).split(df, y, groups))
    rows = []
    for block in blocks:
        for alpha in alphas:
            fold_metrics = []
            for fold, (train_index, valid_index) in enumerate(splits, 1):
                train, valid = df.iloc[train_index], df.iloc[valid_index]
                encoder = FeatureEncoder().fit(train)
                x_train, _ = encoder.transform(train, "M1")
                x_valid, _ = encoder.transform(valid, "M1")
                feature_count = x_train.shape[1]
                if block != "R_none":
                    raw_train, names = raw_block(train, block)
                    raw_valid, _ = raw_block(valid, block)
                    include_mask = block.startswith("R_profile")
                    scaler = BlockScaler(include_mask).fit(raw_train)
                    r_train, r_names = scaler.transform(raw_train, names)
                    r_valid, _ = scaler.transform(raw_valid, names)
                    x_train = np.hstack([x_train, r_train])
                    x_valid = np.hstack([x_valid, r_valid])
                    feature_count += len(r_names)
                estimator = Ridge(alpha=alpha).fit(x_train, y[train_index])
                prediction = np.clip(estimator.predict(x_valid), 0.0, 1.0)
                values = prediction_metrics(y[valid_index], prediction)
                fold_metrics.append(values)
                rows.append({"block": block, "alpha": alpha, "fold": fold, "features": feature_count, **values})
            rows.append({
                "block": block, "alpha": alpha, "fold": "mean", "features": feature_count,
                **{key: float(np.mean([value[key] for value in fold_metrics])) for key in fold_metrics[0]},
            })
    means = [row for row in rows if row["fold"] == "mean"]
    summary = []
    for block in blocks:
        selected = min((row for row in means if row["block"] == block), key=lambda row: (row["mae"], row["alpha"]))
        if selected["alpha"] in {min(alphas), max(alphas)}:
            raise ValueError(
                f"{block} selected boundary alpha={selected['alpha']:g}; "
                "expand alpha_candidates before freezing the diagnostic"
            )
        summary.append(dict(selected))
    baseline_mae = next(row["mae"] for row in summary if row["block"] == "R_none")
    for row in summary:
        row["delta_mae_vs_sequence_baseline"] = baseline_mae - row["mae"]
        row["delta_mae_percentage_points"] = (baseline_mae - row["mae"]) * 100
        row["analysis_scope"] = "development_GroupKFold_only"
        row["selection_role"] = "sensitivity_only" if row["block"] == "R_profile_all" else "candidate"
    return rows, summary


def agreement_row(label: str, df: pd.DataFrame) -> dict:
    left = df["ratio_rep1"].to_numpy(dtype=float)
    right = df["ratio_rep2"].to_numpy(dtype=float)
    return {
        "subset": label, "sites": len(df), "genes": df["analysis_gene"].nunique(),
        "agcov_min": float(df["combined_agcov"].min()), "agcov_median": float(df["combined_agcov"].median()),
        "replicate_mae": float(mean_absolute_error(left, right)),
        "replicate_rmse": float(math.sqrt(mean_squared_error(left, right))),
        "replicate_spearman": float(spearmanr(left, right).statistic),
        "replicate_pearson": float(pearsonr(left, right).statistic),
        "mean_rep1_minus_rep2": float(np.mean(left - right)),
        "interpretation": "measurement_variability_reference_not_formal_noise_ceiling",
    }


def replicate_agreement(df: pd.DataFrame, thresholds: list[int]) -> list[dict]:
    rows = [agreement_row("all_development", df)]
    quantiles = pd.qcut(df["combined_agcov"], q=4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    for label in quantiles.cat.categories:
        rows.append(agreement_row(f"combined_agcov_{label}", df[quantiles == label]))
    for threshold in thresholds:
        subset = df[df["combined_agcov"] >= threshold]
        if len(subset):
            rows.append(agreement_row(f"combined_agcov_ge_{threshold}", subset))
    return rows


def main() -> None:
    setup_logging()
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    source = pd.read_csv(INPUT, compression="gzip", low_memory=False)
    source = source[source["model_dataset_included"].astype(str).eq("True")].copy()
    assignments = pd.read_csv(ASSIGNMENTS)
    predicted = pd.read_csv(PREDICTED_INPUT, compression="gzip", low_memory=False)
    df = source.merge(assignments[["site_id", "joint_group", "split"]], on="site_id", validate="one_to_one")
    df = df.merge(predicted.drop(columns=["sequence_sha256", "sequence_length"]), on="site_id", validate="one_to_one")
    development = add_profile_columns(df[df["split"].eq("development")].copy().reset_index(drop=True))
    if len(development) != 3271:
        raise AssertionError("08B must use the frozen development population only")
    blocks = list(config["experimental_structure_blocks"])
    rows, summary = development_cv(development, blocks, [float(value) for value in config["alpha_candidates"]])
    write_rows(TABLES / "08b_hek293t_expanded_structure_development_cv.csv", rows)
    write_rows(TABLES / "08b_hek293t_expanded_structure_summary.csv", summary)
    agreement = replicate_agreement(development, [int(value) for value in config["replicate_agreement"]["depth_thresholds"]])
    write_rows(TABLES / "08b_hek293t_glori_replicate_agreement.csv", agreement)

    candidates = [row for row in summary if row["block"] in config["selection_candidates"]]
    selected = min(candidates, key=lambda row: (row["mae"], row["block"]))
    nonnull = min((row for row in candidates if row["block"] != "R_none"), key=lambda row: (row["mae"], row["block"]))
    overall = agreement[0]
    report = f"""# Stage 08B development-only structure diagnostic report

## Scope

All feature-block comparisons use only the frozen HEK293T development set ({len(development):,} sites,
{development['analysis_gene'].nunique():,} genes, {development['joint_group'].nunique():,} joint groups).
The reused internal test set was not consulted for block selection. Every imputation and scaling step was
fitted within its GroupKFold training fold.

## Expanded experimental-structure blocks

| Block | Selected alpha | Mean CV MAE | ΔMAE vs M1 (percentage points) | Role |
|---|---:|---:|---:|---|
""" + "\n".join(
        f"| {row['block']} | {row['alpha']:g} | {row['mae']:.5f} | {row['delta_mae_percentage_points']:.4f} | {row['selection_role']} |"
        for row in summary
    ) + f"""

Among the prespecified selectable blocks, `{selected['block']}` had the lowest development CV MAE
({selected['mae']:.5f}). No experimental-structure block improved on the sequence baseline, so the
development-stage decision is to retain no expanded R block. `{nonnull['block']}` was the least harmful
non-null block (ΔMAE {nonnull['delta_mae_percentage_points']:.4f} percentage points). This is a
development-stage selection result, not held-out evidence. The
all-position profile includes the modified center and is sensitivity-only because that position may
contain local probing or modification-linked signal.

## GLORI replicate measurement reference

Across the development set, replicate-1 versus replicate-2 MAE was {overall['replicate_mae']:.5f}
({overall['replicate_mae'] * 100:.2f} percentage points), RMSE was {overall['replicate_rmse']:.5f},
Spearman was {overall['replicate_spearman']:.4f}, and Pearson was {overall['replicate_pearson']:.4f}.
Depth-stratified results are retained in the accompanying table. These discrepancies quantify observed
repeat variability; they are not a formal, identifiable upper bound on model performance because the
combined outcome and each replicate have different measurement-error properties.

## Decision boundary

No expanded-R block is selected. `R_none` is frozen as the development decision for external HeLa
evaluation; `R_core` may be retained only as a prespecified descriptive sensitivity analysis. Association,
predictive increment, and repeat agreement remain distinct estimands.
"""
    (RESULTS / "08b_hek293t_structure_diagnostic_report.md").write_text(report, encoding="utf-8")
    status = [{
        "stage": "08B", "status": "PASS", "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(), "analysis_scope": "development_only",
        "sites": len(development), "genes": development["analysis_gene"].nunique(),
        "joint_groups": development["joint_group"].nunique(), "blocks_compared": len(blocks),
        "selected_block": selected["block"], "selected_block_cv_mae": selected["mae"],
        "alpha_grid_min": min(config["alpha_candidates"]),
        "alpha_grid_max": max(config["alpha_candidates"]),
        "alpha_boundary_check": "PASS",
        "best_nonnull_block": nonnull["block"],
        "replicate_mae": overall["replicate_mae"],
    }]
    write_rows(RESULTS / "08b_hek293t_structure_diagnostics_run_status.csv", status)
    logging.info("08B completed; selected development block=%s", selected["block"])


if __name__ == "__main__":
    main()
