"""Stage 13: cross-fitted residual-on-residual (GCM-style) conditional-dependence test."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

from phase1_common import TabularEncoder, make_joint_groups, sha256_file, write_rows

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "13_conditional_dependence.yaml"
TABLE = ROOT / "results" / "tables" / "13_conditional_dependence.csv"
RESIDUALS = ROOT / "results" / "tables" / "13_cross_fitted_residuals.csv.gz"
STATUS = ROOT / "results" / "status" / "13_conditional_dependence_run_status.csv"
REPORT = ROOT / "results" / "reports" / "13_conditional_dependence_report.md"


def cross_fitted_residuals(df: pd.DataFrame, folds: int, alpha: float) -> pd.DataFrame:
    work = df.reset_index(drop=True).copy()
    work["joint_group"] = make_joint_groups(work, "sequence_window_201")
    groups = work["joint_group"].to_numpy()
    y_m = work["combined_ratio"].to_numpy(float)
    y_r = work["reactivity_mean_flank10"].to_numpy(float)
    pred_m, pred_r = np.full(len(work), np.nan), np.full(len(work), np.nan)
    fold_id = np.zeros(len(work), dtype=int)
    splitter = GroupKFold(n_splits=folds)
    for fold, (train, valid) in enumerate(splitter.split(work, y_m, groups), 1):
        encoder = TabularEncoder().fit(work.iloc[train])
        x_train, _ = encoder.transform(work.iloc[train], "M1")
        x_valid, _ = encoder.transform(work.iloc[valid], "M1")
        pred_m[valid] = Ridge(alpha=alpha).fit(x_train, y_m[train]).predict(x_valid)
        pred_r[valid] = Ridge(alpha=alpha).fit(x_train, y_r[train]).predict(x_valid)
        fold_id[valid] = fold
    if np.isnan(pred_m).any() or np.isnan(pred_r).any():
        raise AssertionError("Cross-fitting did not cover every row")
    return pd.DataFrame({"site_id": work["site_id"], "analysis_gene": work["analysis_gene"],
                         "joint_group": groups, "fold": fold_id, "residual_m": y_m - pred_m,
                         "residual_r": y_r - pred_r, "prediction_m": pred_m, "prediction_r": pred_r})


def infer(residuals: pd.DataFrame, replicates: int, seed: int) -> dict:
    product = residuals["residual_m"].to_numpy() * residuals["residual_r"].to_numpy()
    design = np.ones((len(product), 1))
    fit = sm.OLS(product, design).fit(cov_type="cluster", cov_kwds={"groups": residuals["analysis_gene"]})
    estimate, se = float(fit.params[0]), float(fit.bse[0])
    grouped = pd.DataFrame({"joint_group": residuals["joint_group"], "product": product}).groupby("joint_group")["product"].sum().to_numpy()
    rng = np.random.default_rng(seed)
    null = np.empty(replicates)
    denominator = len(product)
    for index in range(replicates):
        null[index] = float(np.sum(grouped * rng.choice([-1.0, 1.0], len(grouped))) / denominator)
    p_perm = float((1 + np.sum(np.abs(null) >= abs(estimate))) / (replicates + 1))
    correlation = float(np.corrcoef(residuals["residual_m"], residuals["residual_r"])[0, 1])
    return {"sites": len(product), "genes": residuals["analysis_gene"].nunique(),
            "joint_groups": residuals["joint_group"].nunique(), "gcm_mean_residual_product": estimate,
            "cluster_se": se, "cluster_ci95_low": estimate - 1.96 * se,
            "cluster_ci95_high": estimate + 1.96 * se, "cluster_p_value": float(fit.pvalues[0]),
            "group_sign_flip_p_value": p_perm, "residual_correlation": correlation,
            "permutation_replicates": replicates}


def main() -> None:
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    all_residuals, rows = [], []
    folds = int(config["cross_fitting"]["folds"])
    alpha = float(config["cross_fitting"]["ridge_alpha"])
    for offset, (dataset, item) in enumerate(config["datasets"].items()):
        source = pd.read_csv(ROOT / item["input"], low_memory=False)
        source = source[source["model_dataset_included"].astype(str).eq("True")].copy()
        residuals = cross_fitted_residuals(source, folds, alpha)
        residuals.insert(0, "dataset", dataset)
        all_residuals.append(residuals)
        row = {"dataset": dataset, "role": item["role"], "nuisance_model": "Ridge",
               "ridge_alpha": alpha, "folds": folds,
               **infer(residuals, int(config["inference"]["group_sign_flip_replicates"]),
                       int(config["inference"]["seed"]) + offset)}
        rows.append(row)
    write_rows(TABLE, rows)
    pd.concat(all_residuals).to_csv(RESIDUALS, index=False, compression={"method": "gzip", "mtime": 0})
    write_rows(STATUS, [{"stage": 13, "status": "PASS", "started_utc": started.isoformat(),
                          "finished_utc": datetime.now(timezone.utc).isoformat(), "datasets": len(rows),
                          "config_sha256": sha256_file(CONFIG), "table_sha256": sha256_file(TABLE),
                          "residuals_sha256": sha256_file(RESIDUALS)}])
    lines = ["# Stage 13: Conditional dependency testing", "", "## Material Passport", "",
             "- Origin: phase-1 exploratory extension", "- Verification Status: ANALYZED", "",
             "A residual-on-residual/GCM-style test using gene-wise cross-fitting to the same 201 nt sequence connected group. The", ""]
    for row in rows:
        lines.append(f"- {row['dataset']} ({row['role']}): residual correlation {row['residual_correlation']:.4f}, clustered p={row['cluster_p_value']:.4g}, and grouped sign-flip p={row['group_sign_flip_p_value']:.4g}.")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
