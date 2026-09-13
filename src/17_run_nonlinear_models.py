"""Stage 17: Ridge/HGB comparison, nested ablations and grouped permutation explanation."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

from phase1_common import TabularEncoder, make_joint_groups, model_metrics, sha256_file, write_rows

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/17_nonlinear_models.yaml"
TABLES = ROOT / "results/tables"
STATUS = ROOT / "results/status/17_nonlinear_models_run_status.csv"
REPORT = ROOT / "results/reports/17_nonlinear_models_report.md"


def load_dataset(item: dict) -> pd.DataFrame:
    source = pd.read_csv(ROOT / item["input"], low_memory=False)
    source = source[source["model_dataset_included"].astype(str).eq("True")].copy()
    predicted = pd.read_csv(ROOT / item["predicted"], low_memory=False)
    return source.merge(predicted, on="site_id", validate="one_to_one").reset_index(drop=True)


def make_estimator(learner: str, parameters: dict):
    if learner == "Ridge":
        return Ridge(alpha=float(parameters["alpha"]), solver="lsqr")
    return HistGradientBoostingRegressor(random_state=20260917, early_stopping=False, **parameters)


def candidate_grid(config: dict):
    yield from (("Ridge", {"alpha": float(x)}) for x in config["learners"]["Ridge"]["alpha_candidates"])
    yield from (("HistGradientBoostingRegressor", dict(x)) for x in config["learners"]["HistGradientBoostingRegressor"]["candidates"])


def select_parameters(data: pd.DataFrame, groups: np.ndarray, config: dict) -> tuple[list[dict], dict]:
    rows, selected = [], {}
    folds = int(config["development_cv"]["folds"])
    for model_name in config["model_matrix"]:
        for learner, parameters in candidate_grid(config):
            fold_scores = []
            for fold, (train, valid) in enumerate(GroupKFold(folds).split(data, groups=groups), 1):
                encoder = TabularEncoder().fit(data.iloc[train])
                x_train, _ = encoder.transform(data.iloc[train], model_name)
                x_valid, _ = encoder.transform(data.iloc[valid], model_name)
                fitted = make_estimator(learner, parameters).fit(x_train, data.iloc[train]["combined_ratio"])
                prediction = np.clip(fitted.predict(x_valid), 0, 1)
                score = model_metrics(data.iloc[valid]["combined_ratio"].to_numpy(), prediction)
                fold_scores.append(score)
                rows.append({"model": model_name, "learner": learner,
                             "parameters": json.dumps(parameters, sort_keys=True), "fold": fold, **score})
            rows.append({"model": model_name, "learner": learner,
                         "parameters": json.dumps(parameters, sort_keys=True), "fold": "mean",
                         **{key: float(np.mean([x[key] for x in fold_scores])) for key in fold_scores[0]}})
        for learner in ("Ridge", "HistGradientBoostingRegressor"):
            best = min((row for row in rows if row["model"] == model_name and row["learner"] == learner
                        and row["fold"] == "mean"), key=lambda row: row["mae"])
            selected[(model_name, learner)] = json.loads(best["parameters"])
    return rows, selected


def external_evaluation(hek: pd.DataFrame, hela: pd.DataFrame, config: dict, selected: dict) -> list[dict]:
    rows = []
    for model_name in config["model_matrix"]:
        for learner in ("Ridge", "HistGradientBoostingRegressor"):
            encoder = TabularEncoder().fit(hek)
            x_hek, _ = encoder.transform(hek, model_name)
            x_hela, _ = encoder.transform(hela, model_name)
            parameters = selected[(model_name, learner)]
            fitted = make_estimator(learner, parameters).fit(x_hek, hek["combined_ratio"])
            prediction = np.clip(fitted.predict(x_hela), 0, 1)
            rows.append({"dataset": "HeLa", "model": model_name, "learner": learner,
                         "parameters_selected_in_HEK293T": json.dumps(parameters, sort_keys=True),
                         **model_metrics(hela["combined_ratio"].to_numpy(), prediction)})
    return rows


def grouped_permutation(data: pd.DataFrame, groups: np.ndarray, parameters: dict,
                        folds: int, repeats: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    deltas = {name: [] for name in ("common_covariates", "sequence_3mer", "predicted_structure", "experimental_structure")}
    for train, valid in GroupKFold(folds).split(data, groups=groups):
        encoder = TabularEncoder().fit(data.iloc[train])
        x_train, names = encoder.transform(data.iloc[train], "M4")
        x_valid, _ = encoder.transform(data.iloc[valid], "M4")
        fitted = make_estimator("HistGradientBoostingRegressor", parameters).fit(x_train, data.iloc[train]["combined_ratio"])
        y = data.iloc[valid]["combined_ratio"].to_numpy()
        baseline = model_metrics(y, np.clip(fitted.predict(x_valid), 0, 1))["mae"]
        blocks = {
            "common_covariates": [i for i, name in enumerate(names) if name.startswith("C:")],
            "sequence_3mer": [i for i, name in enumerate(names) if name.startswith("X:")],
            "predicted_structure": [i for i, name in enumerate(names) if name.startswith("P:")],
            "experimental_structure": [i for i, name in enumerate(names) if name.startswith("R:")],
        }
        for block, columns in blocks.items():
            for _ in range(repeats):
                permuted = x_valid.copy()
                order = rng.permutation(len(valid))
                permuted[:, columns] = permuted[order][:, columns]
                score = model_metrics(y, np.clip(fitted.predict(permuted), 0, 1))["mae"]
                deltas[block].append(score - baseline)
    return [{"learner": "HistGradientBoostingRegressor", "model": "M4", "feature_block": block,
             "delta_mae_mean": float(np.mean(values)), "delta_mae_sd": float(np.std(values, ddof=1)),
             "delta_mae_q025": float(np.quantile(values, .025)), "delta_mae_q975": float(np.quantile(values, .975)),
             "fold_repeat_estimates": len(values)} for block, values in deltas.items()]


def main() -> None:
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    hek, hela = load_dataset(config["datasets"]["HEK293T"]), load_dataset(config["datasets"]["HeLa"])
    hek["joint_group"] = make_joint_groups(hek, "sequence_window_201")
    groups = hek["joint_group"].to_numpy()
    cv_rows, selected = select_parameters(hek, groups, config)
    external = external_evaluation(hek, hela, config, selected)
    importance = grouped_permutation(
        hek, groups, selected[("M4", "HistGradientBoostingRegressor")],
        int(config["development_cv"]["folds"]), int(config["explanation"]["repeats"]), int(config["explanation"]["seed"]),
    )
    ablations = []
    for learner in ("Ridge", "HistGradientBoostingRegressor"):
        metrics = {row["model"]: row for row in external if row["learner"] == learner}
        for baseline, augmented in (("M0", "M1"), ("M1", "M2"), ("M1", "M3"), ("M2", "M4")):
            ablations.append({"dataset": "HeLa", "learner": learner, "comparison": f"{augmented}_vs_{baseline}",
                              "baseline": baseline, "augmented": augmented,
                              "delta_mae_baseline_minus_augmented": metrics[baseline]["mae"] - metrics[augmented]["mae"]})
    cv_path = TABLES / "17_nonlinear_development_cv.csv"
    external_path = TABLES / "17_nonlinear_external_metrics.csv"
    write_rows(cv_path, cv_rows); write_rows(external_path, external)
    write_rows(TABLES / "17_grouped_permutation_importance.csv", importance)
    write_rows(TABLES / "17_retraining_ablation.csv", ablations)
    write_rows(STATUS, [{"stage": 17, "status": "PASS", "started_utc": started.isoformat(),
                         "finished_utc": datetime.now(timezone.utc).isoformat(), "hek293t_sites": len(hek),
                         "hela_sites": len(hela), "external_fits": len(external),
                         "config_sha256": sha256_file(CONFIG), "cv_sha256": sha256_file(cv_path),
                         "external_sha256": sha256_file(external_path)}])
    ridge = next(row for row in external if row["learner"] == "Ridge" and row["model"] == "M4")
    hgb = next(row for row in external if row["learner"] == "HistGradientBoostingRegressor" and row["model"] == "M4")
    REPORT.write_text(f"""# Stage 17: Nonlinear Model and Interpretation

## Material Passport

- Origin: phase-1 exploratory extension
- Verification Status: ANALYZED

Both Ridge and HistGradientBoostingRegressor retain M0–M4 nested ablation, and only the gene-identical sequence connected group cross-validation of HEK293T is used for parameter adjustment; HeLa does not participate in parameter adjustment. M4 external MAE: Ridge={ridge['mae']:.5f}, HGB={hgb['mae']:.5f}, ΔMAE of HGB relative to Ridge={ridge['mae']-hgb['mae']:.6f}.

Interpretation of the main results as feature block permutation importance in grouped cross-validation on HEK293T, supplemented by retraining ablation on HeLa. Significance is not a causal explanation; SHAP does not serve as a basis for conclusions, nor does it expand the hyperparameter space based on HeLa performance.""", encoding="utf-8")


if __name__ == "__main__":
    main()
