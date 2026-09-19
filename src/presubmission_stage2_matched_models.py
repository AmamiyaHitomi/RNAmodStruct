"""Run the phase-0-fixed, post hoc matched model comparisons for main_v2.

Only the frozen HEK293T development partition is used for fitting and tuning.
HEK test and HeLa are evaluated after all candidate selection is complete.
The script refuses to overwrite a previous stage-2 result directory.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold

from phase1_common import TabularEncoder


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "presubmission_stage2"
PHASE0 = ROOT / "metadata" / "manifests" / "presubmission_phase0"
SPEC_PATH = PHASE0 / "analysis_spec.json"
FREEZE_PATH = PHASE0 / "freeze_manifest.json"
CONFIG17 = ROOT / "config/17_nonlinear_models.yaml"
INPUTS = {
    "hek": ROOT / "data/final/06_hek293t_main_analysis_dataset.csv.gz",
    "hela": ROOT / "data/final/09_hela_main_analysis_dataset.csv.gz",
    "hek_predicted": ROOT / "data/final/08a_hek293t_predicted_structure_features.csv.gz",
    "hela_predicted": ROOT / "data/final/09a_hela_predicted_structure_features.csv.gz",
    "assignments": ROOT / "results/tables/08_hek293t_group_assignments.csv",
    "hela_frozen_predictions": ROOT / "results/tables/09c_hela_transfer_predictions.csv",
    "phase0_flags": PHASE0 / "hela_site_cohort_flags.csv",
    "spec": SPEC_PATH,
    "freeze": FREEZE_PATH,
    "config17": CONFIG17,
}
SCENARIOS = [
    {"name": "C1_M0_plus_R", "encoding": "original_201_onehot_plus_3mer",
     "drop_abundance": False, "models": ["M0", "M0R"], "learners": ["Ridge"]},
    {"name": "C2_no_abundance", "encoding": "original_201_onehot_plus_3mer",
     "drop_abundance": True, "models": ["M0", "M0R", "M1", "M3"], "learners": ["Ridge"]},
    {"name": "C3_matched_3mer", "encoding": "normalized_3mer_only",
     "drop_abundance": False, "models": ["M0", "M1", "M2", "M3", "M4"], "learners": ["Ridge", "HGB"]},
]
COHORT_FLAGS = {
    "HeLa_all_qualifying": None,
    "HeLa_development_overlap_excluded": "development_overlap_excluded",
    "HeLa_hek_main_4409_overlap_excluded": "hek_main_4409_overlap_excluded",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def id_sha256(ids: list[str]) -> str:
    return hashlib.sha256("".join(f"{item}\n" for item in sorted(ids)).encode("utf-8")).hexdigest()


def stage08_feature_encoder_class():
    path = ROOT / "src/08_run_prediction_modeling.py"
    spec = importlib.util.spec_from_file_location("frozen_stage08", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the original FeatureEncoder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.FeatureEncoder


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, dict]:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
    freeze = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
    if spec["status"] != "post_hoc_plan_fixed_after_original_HEK_test_and_HeLa_evaluations":
        raise ValueError("Unexpected analysis-spec status")
    hek = pd.read_csv(INPUTS["hek"], low_memory=False)
    hek = hek[hek["model_dataset_included"].eq(True)].copy()
    hela = pd.read_csv(INPUTS["hela"], low_memory=False)
    hela = hela[hela["model_dataset_included"].eq(True)].copy()
    for frame, source in [(hek, INPUTS["hek_predicted"]), (hela, INPUTS["hela_predicted"])]:
        predicted = pd.read_csv(source, low_memory=False)
        predicted = predicted.drop(columns=["sequence_sha256", "sequence_length"])
        merged = frame.merge(predicted, on="site_id", how="left", validate="one_to_one", sort=False)
        if merged["predicted_unpaired_probabilities_json"].isna().any():
            raise ValueError(f"Missing predicted features from {source}")
        if frame is hek:
            hek = merged
        else:
            hela = merged
    assignments = pd.read_csv(INPUTS["assignments"])
    hek = hek.merge(assignments[["site_id", "split", "joint_group"]], on="site_id",
                    validate="one_to_one", sort=False)
    if len(hek) != 4096 or int(hek["split"].eq("development").sum()) != 3271:
        raise ValueError("Frozen HEK split changed")
    development = hek[hek["split"].eq("development")].copy().reset_index(drop=True)
    test = hek[hek["split"].eq("test")].copy().reset_index(drop=True)
    if len(test) != 825 or set(development["joint_group"]) & set(test["joint_group"]):
        raise ValueError("Frozen HEK test/group assignment changed")
    hela_frozen = pd.read_csv(INPUTS["hela_frozen_predictions"],
                              usecols=["site_id", "joint_group"])
    hela = hela.merge(hela_frozen, on="site_id", validate="one_to_one", sort=False)
    flags = pd.read_csv(INPUTS["phase0_flags"],
                        usecols=["site_id", "development_overlap_excluded",
                                 "hek_main_4409_overlap_excluded"])
    hela = hela.merge(flags, on="site_id", validate="one_to_one", sort=False)
    if len(hela) != 24960 or not hela["sequence_window_201_full"].eq(True).all():
        raise ValueError("Frozen HeLa model population changed")
    for name, flag in COHORT_FLAGS.items():
        part = hela if flag is None else hela[hela[flag]]
        key = name.removeprefix("HeLa_")
        if len(part) != freeze["cohorts"][key]["site_count"]:
            raise ValueError(f"Phase-0 cohort count changed: {name}")
        if id_sha256(part["site_id"].astype(str).tolist()) != freeze["cohorts"][key]["site_ids_sha256"]:
            raise ValueError(f"Phase-0 cohort ID hash changed: {name}")
    return development, test, hela, spec, freeze


def encoder_for(scenario: dict, original_class, train: pd.DataFrame):
    if scenario["encoding"] == "normalized_3mer_only":
        return TabularEncoder().fit(train)
    encoder = original_class()
    if scenario["drop_abundance"]:
        encoder.numeric = [name for name in encoder.numeric
                           if name != "log1p_icshape_abundance_rpkm"]
    return encoder.fit(train)


def matrix(encoder, scenario: dict, frame: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str]]:
    if model == "M0R":
        common, common_names = encoder.common(frame)
        structure, structure_names = encoder.structure(frame)
        values, names = np.hstack([common, structure]), common_names + structure_names
    else:
        values, names = encoder.transform(frame, model)
    if scenario["drop_abundance"] and any("abundance" in name for name in names):
        raise AssertionError("Abundance survived the feature deletion")
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite features in {scenario['name']}/{model}")
    return values, names


def candidates(spec: dict, config17: dict, learner: str) -> list[dict]:
    if learner == "Ridge":
        values = next(x for x in spec["comparisons"] if x["id"] == "C3_matched_Ridge_HGB_3mer")
        return [{"alpha": float(x)} for x in values["Ridge_candidates"]]
    return [dict(x) for x in config17["learners"]["HistGradientBoostingRegressor"]["candidates"]]


def estimator(scenario: dict, learner: str, params: dict):
    if learner == "HGB":
        return HistGradientBoostingRegressor(random_state=20260917,
                                             early_stopping=False, **params)
    solver = "lsqr" if scenario["encoding"] == "normalized_3mer_only" else "auto"
    return Ridge(alpha=float(params["alpha"]), solver=solver)


def evaluate(y: np.ndarray, pred: np.ndarray) -> dict:
    rho = spearmanr(y, pred).statistic
    return {
        "mae": float(mean_absolute_error(y, pred)),
        "rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "spearman": float(rho),
        "r2": float(r2_score(y, pred)),
        "mean_signed_error_prediction_minus_observed": float((pred - y).mean()),
        "at_clip_fraction": float(((pred == 0) | (pred == 1)).mean()),
    }


def model_key(scenario: str, learner: str, model: str) -> str:
    return f"{scenario}__{learner}__{model}"


def select_candidates(development: pd.DataFrame, spec: dict, config17: dict,
                      original_class) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    splits = list(GroupKFold(n_splits=5).split(
        development, groups=development["joint_group"].to_numpy()))
    score_by_candidate = defaultdict(list)
    fold_rows = []
    for scenario in SCENARIOS:
        print("CV scenario", scenario["name"], flush=True)
        for fold, (train_idx, valid_idx) in enumerate(splits, 1):
            train, valid = development.iloc[train_idx], development.iloc[valid_idx]
            encoder = encoder_for(scenario, original_class, train)
            y_train = train["combined_ratio"].to_numpy(dtype=float)
            y_valid = valid["combined_ratio"].to_numpy(dtype=float)
            for model in scenario["models"]:
                x_train, names_train = matrix(encoder, scenario, train, model)
                x_valid, names_valid = matrix(encoder, scenario, valid, model)
                if names_train != names_valid:
                    raise AssertionError("Fold feature order changed")
                for learner in scenario["learners"]:
                    for candidate_id, params in enumerate(candidates(spec, config17, learner)):
                        fitted = estimator(scenario, learner, params).fit(x_train, y_train)
                        pred = np.clip(fitted.predict(x_valid), 0, 1)
                        mae = float(mean_absolute_error(y_valid, pred))
                        key = (scenario["name"], learner, model, candidate_id)
                        score_by_candidate[key].append(mae)
                        fold_rows.append({"scenario": scenario["name"], "learner": learner,
                                          "model": model, "candidate_id": candidate_id,
                                          "parameters": json.dumps(params, sort_keys=True),
                                          "fold": fold, "validation_sites": len(valid), "mae": mae})
            print("  fold", fold, "complete", flush=True)
    mean_rows, selected = [], {}
    for scenario in SCENARIOS:
        for learner in scenario["learners"]:
            grid = candidates(spec, config17, learner)
            for model in scenario["models"]:
                for candidate_id, params in enumerate(grid):
                    key = (scenario["name"], learner, model, candidate_id)
                    values = score_by_candidate[key]
                    if len(values) != 5:
                        raise ValueError(f"Incomplete CV for {key}")
                    mean_rows.append({"scenario": scenario["name"], "learner": learner,
                                      "model": model, "candidate_id": candidate_id,
                                      "parameters": json.dumps(params, sort_keys=True),
                                      "fold": "mean", "validation_sites": None,
                                      "mae": float(np.mean(values))})
                options = [row for row in mean_rows if row["scenario"] == scenario["name"]
                           and row["learner"] == learner and row["model"] == model]
                best = min(options, key=lambda row: (row["mae"], row["candidate_id"]))
                selected[model_key(scenario["name"], learner, model)] = {
                    "scenario": scenario["name"], "learner": learner, "model": model,
                    "candidate_id": int(best["candidate_id"]),
                    "parameters": grid[int(best["candidate_id"])],
                    "cv_mean_mae": float(best["mae"]),
                    "candidate_boundary": best["candidate_id"] in {0, len(grid) - 1},
                }
    return pd.DataFrame(fold_rows + mean_rows), pd.DataFrame(selected.values()), selected


def fit_and_predict(development: pd.DataFrame, test: pd.DataFrame, hela: pd.DataFrame,
                    original_class, selected: dict) -> tuple[dict, pd.DataFrame, pd.DataFrame, list]:
    test_wide = test[["site_id", "analysis_gene", "joint_group", "combined_ratio"]].copy()
    hela_wide = hela[["site_id", "analysis_gene", "joint_group", "combined_ratio",
                      "development_overlap_excluded", "hek_main_4409_overlap_excluded"]].copy()
    bundles = {}
    selected_rows = []
    for scenario in SCENARIOS:
        print("Final fit scenario", scenario["name"], flush=True)
        encoder = encoder_for(scenario, original_class, development)
        for model in scenario["models"]:
            x_train, train_names = matrix(encoder, scenario, development, model)
            x_test, test_names = matrix(encoder, scenario, test, model)
            x_hela, hela_names = matrix(encoder, scenario, hela, model)
            if train_names != test_names or train_names != hela_names:
                raise AssertionError("Final train/test/HeLa feature order changed")
            for learner in scenario["learners"]:
                key = model_key(scenario["name"], learner, model)
                choice = selected[key]
                fitted = estimator(scenario, learner, choice["parameters"]).fit(
                    x_train, development["combined_ratio"].to_numpy(dtype=float))
                test_wide[key] = np.clip(fitted.predict(x_test), 0, 1)
                hela_wide[key] = np.clip(fitted.predict(x_hela), 0, 1)
                state = encoder.to_dict() if scenario["encoding"].startswith("original") else encoder.__dict__
                bundles[key] = {
                    "model": fitted, "encoder_state": state, "encoder_kind": scenario["encoding"],
                    "drop_abundance": scenario["drop_abundance"],
                    "feature_names": train_names, "selected_parameters": choice["parameters"],
                    "training_site_ids_sha256": id_sha256(development["site_id"].astype(str).tolist()),
                    "prediction_clip": [0.0, 1.0], "analysis_role": "post_hoc_sensitivity",
                }
                selected_rows.append({**choice, "key": key, "features": len(train_names),
                                      "HEK_reused_test_mae": float(mean_absolute_error(
                                          test["combined_ratio"], test_wide[key])),
                                      "HeLa_all_mae": float(mean_absolute_error(
                                          hela["combined_ratio"], hela_wide[key]))})
                print("  fitted", key, "CV MAE", round(choice["cv_mean_mae"], 6), flush=True)
    return bundles, test_wide, hela_wide, selected_rows


def evaluation_table(test_wide: pd.DataFrame, hela_wide: pd.DataFrame,
                     keys: list[str]) -> pd.DataFrame:
    rows = []
    cohorts = {"HEK_reused_test": test_wide}
    for name, flag in COHORT_FLAGS.items():
        cohorts[name] = hela_wide if flag is None else hela_wide[hela_wide[flag]]
    for name, frame in cohorts.items():
        y = frame["combined_ratio"].to_numpy(dtype=float)
        for key in keys:
            scenario, learner, model = key.split("__")
            rows.append({"cohort": name, "scenario": scenario, "learner": learner,
                         "model": model, "key": key, "sites": len(frame),
                         "genes": int(frame["analysis_gene"].nunique()),
                         "joint_groups": int(frame["joint_group"].nunique()),
                         "role": "reused_descriptive" if name == "HEK_reused_test" else "post_hoc_HeLa_sensitivity",
                         **evaluate(y, frame[key].to_numpy(dtype=float))})
    return pd.DataFrame(rows)


def comparisons() -> list[tuple[str, str, str, str]]:
    rows = []
    for scenario, learner, pairs in [
        ("C1_M0_plus_R", "Ridge", [("M0", "M0R")]),
        ("C2_no_abundance", "Ridge", [("M0", "M0R"), ("M1", "M3")]),
        ("C3_matched_3mer", "Ridge", [("M0", "M1"), ("M1", "M2"),
                                        ("M1", "M3"), ("M2", "M4")]),
        ("C3_matched_3mer", "HGB", [("M0", "M1"), ("M1", "M2"),
                                      ("M1", "M3"), ("M2", "M4")]),
    ]:
        rows.extend(("feature_increment", scenario,
                     model_key(scenario, learner, base), model_key(scenario, learner, aug))
                    for base, aug in pairs)
    for model in ["M0", "M1", "M2", "M3", "M4"]:
        rows.append(("matched_learner", "C3_matched_3mer",
                     model_key("C3_matched_3mer", "Ridge", model),
                     model_key("C3_matched_3mer", "HGB", model)))
    return rows


def paired_contrasts(test_wide: pd.DataFrame, hela_wide: pd.DataFrame,
                     keys: list[str], spec: dict) -> tuple[pd.DataFrame, dict]:
    definitions = comparisons()
    rows = []
    index_hashes = {}
    for cohort, flag in COHORT_FLAGS.items():
        frame = hela_wide if flag is None else hela_wide[hela_wide[flag]]
        groups, codes = np.unique(frame["joint_group"].astype(str).to_numpy(), return_inverse=True)
        group_n = len(groups)
        counts = np.bincount(codes, minlength=group_n).astype(float)
        draw = np.random.default_rng(int(spec["paired_bootstrap"]["seed"])).integers(
            0, group_n, size=(int(spec["paired_bootstrap"]["replicates"]), group_n), dtype=np.int32)
        index_hashes[cohort] = hashlib.sha256(draw.tobytes()).hexdigest()
        denominator = counts[draw].sum(axis=1)
        y = frame["combined_ratio"].to_numpy(dtype=float)
        mae, bootstrap = {}, {}
        for key in keys:
            errors = np.abs(y - frame[key].to_numpy(dtype=float))
            sums = np.bincount(codes, weights=errors, minlength=group_n)
            mae[key] = float(errors.mean())
            bootstrap[key] = sums[draw].sum(axis=1) / denominator
        for kind, scenario, baseline, augmented in definitions:
            values = bootstrap[baseline] - bootstrap[augmented]
            point = mae[baseline] - mae[augmented]
            rows.append({"cohort": cohort, "comparison_kind": kind, "scenario": scenario,
                         "baseline_key": baseline, "augmented_key": augmented,
                         "delta_mae_baseline_minus_augmented": point,
                         "delta_mae_percentage_points": 100 * point,
                         "ci95_low": float(np.quantile(values, 0.025)),
                         "ci95_high": float(np.quantile(values, 0.975)),
                         "ci95_low_percentage_points": float(100 * np.quantile(values, 0.025)),
                         "ci95_high_percentage_points": float(100 * np.quantile(values, 0.975)),
                         "bootstrap_groups": group_n,
                         "bootstrap_replicates": len(values),
                         "CI_scope": "fixed_fitted_models",
                         "role": "post_hoc_HeLa_sensitivity"})
        print("Bootstrap complete", cohort, "groups", group_n, flush=True)
    # HEK test contrasts are descriptive because that test split was reused.
    y = test_wide["combined_ratio"].to_numpy(dtype=float)
    for kind, scenario, baseline, augmented in definitions:
        point = float(np.abs(y - test_wide[baseline]).mean()
                      - np.abs(y - test_wide[augmented]).mean())
        rows.append({"cohort": "HEK_reused_test", "comparison_kind": kind,
                     "scenario": scenario, "baseline_key": baseline, "augmented_key": augmented,
                     "delta_mae_baseline_minus_augmented": point,
                     "delta_mae_percentage_points": 100 * point,
                     "ci95_low": None, "ci95_high": None,
                     "ci95_low_percentage_points": None, "ci95_high_percentage_points": None,
                     "bootstrap_groups": None, "bootstrap_replicates": 0,
                     "CI_scope": "not_reported_reused_test", "role": "reused_descriptive"})
    return pd.DataFrame(rows), index_hashes


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    started = datetime.now(timezone.utc)
    development, test, hela, spec, freeze = load_inputs()
    original_class = stage08_feature_encoder_class()
    config17 = yaml.safe_load(CONFIG17.read_text(encoding="utf-8"))
    cv, selected_cv, selected = select_candidates(development, spec, config17, original_class)
    bundles, test_wide, hela_wide, selected_rows = fit_and_predict(
        development, test, hela, original_class, selected)
    keys = list(bundles)
    # The C1 M0 result anchors the new fitting pipeline to the original frozen M0.
    original_test = pd.read_csv(ROOT / "results/tables/08_hek293t_test_predictions.csv",
                                usecols=["site_id", "prediction_M0"])
    anchor = test_wide[["site_id", "C1_M0_plus_R__Ridge__M0"]].merge(
        original_test, on="site_id", validate="one_to_one")
    anchor_max_abs = float(np.max(np.abs(anchor["C1_M0_plus_R__Ridge__M0"]
                                          - anchor["prediction_M0"])))
    if anchor_max_abs > 1e-6:
        raise ValueError(f"C1 M0 does not reproduce frozen M0: {anchor_max_abs}")
    metrics = evaluation_table(test_wide, hela_wide, keys)
    contrasts, index_hashes = paired_contrasts(test_wide, hela_wide, keys, spec)
    if len(metrics) != 64 or len(contrasts) != 64:
        raise AssertionError("Unexpected comparison matrix size")
    OUT.mkdir(parents=False, exist_ok=False)
    (OUT / "models").mkdir(exist_ok=False)
    cv.to_csv(OUT / "development_cv_candidates.csv", index=False, lineterminator="\n")
    pd.DataFrame(selected_rows).to_csv(OUT / "selected_models.csv", index=False, lineterminator="\n")
    metrics.to_csv(OUT / "evaluation_metrics.csv", index=False, lineterminator="\n")
    contrasts.to_csv(OUT / "paired_contrasts.csv", index=False, lineterminator="\n")
    test_wide.to_csv(OUT / "HEK_reused_test_predictions.csv", index=False, lineterminator="\n")
    hela_wide.to_csv(OUT / "HeLa_predictions.csv", index=False, lineterminator="\n")
    for key, bundle in bundles.items():
        joblib.dump(bundle, OUT / "models" / f"{key}.joblib")
    provenance = {
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS_post_hoc_sensitivity",
        "frozen_development_sites": len(development),
        "reused_test_sites": len(test),
        "HeLa_all_sites": len(hela),
        "cv_folds": 5,
        "selected_model_count": len(bundles),
        "candidate_cv_rows": len(cv),
        "C1_M0_anchor_max_abs_prediction_difference": anchor_max_abs,
        "bootstrap_index_sha256_by_cohort": index_hashes,
        "bootstrap_scope": "fixed_fitted_models",
        "source_sha256": {name: sha256(path) for name, path in INPUTS.items()},
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": {str(p.relative_to(OUT)).replace("\\", "/"): sha256(p)
                          for p in OUT.rglob("*") if p.is_file()},
    }
    with (OUT / "provenance.json").open("x", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"status": provenance["status"], "models": len(bundles),
                      "cv_rows": len(cv), "metric_rows": len(metrics),
                      "contrast_rows": len(contrasts), "anchor_max_abs": anchor_max_abs,
                      "output": str(OUT)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
