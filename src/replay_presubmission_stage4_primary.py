"""Isolated replay of historical stage-08 Ridge and stage-09C frozen transfer."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from replay_presubmission_stage4_models import compare_csv, sha256


ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "results" / "presubmission_stage4" / "primary_replay"
INPUTS = [
    "src/08_run_prediction_modeling.py", "src/09c_run_hela_model_transfer.py",
    "config/08_prediction_modeling.yaml", "config/09c_hela_model_transfer.yaml",
    "data/final/06_hek293t_main_analysis_dataset.csv.gz",
    "data/final/08a_hek293t_predicted_structure_features.csv.gz",
    "data/final/09_hela_main_analysis_dataset.csv.gz",
    "data/final/09a_hela_predicted_structure_features.csv.gz",
]
TABLES = [
    "08_hek293t_development_cv.csv", "08_hek293t_group_assignments.csv",
    "08_hek293t_test_metrics.csv", "08_hek293t_test_predictions.csv",
    "08_hek293t_feature_manifest.csv", "08_hek293t_split_audit.csv",
    "08_hek293t_primary_model_comparison.csv",
    "09c_hela_transfer_metrics.csv", "09c_hela_transfer_model_comparison.csv",
    "09c_hela_transfer_predictions.csv", "09c_hela_unseen_categorical_levels.csv",
    "09c_hela_transfer_strict_subset_summary.csv",
]


def main() -> None:
    if DEST.exists():
        raise FileExistsError(f"Refusing to overwrite {DEST}")
    started = time.perf_counter()
    sources = {}
    for rel in INPUTS:
        source, target = ROOT / rel, DEST / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        sources[rel] = sha256(source)
        if sha256(target) != sources[rel]:
            raise AssertionError(f"Input copy changed: {rel}")
    runs = []
    for name in ("08_run_prediction_modeling.py", "09c_run_hela_model_transfer.py"):
        result = subprocess.run([sys.executable, str(DEST / "src" / name)], cwd=DEST,
                                capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError(f"{name} exited {result.returncode}: {result.stdout[-3000:]} {result.stderr[-3000:]}")
        runs.append({"script": name, "stdout_tail": result.stdout[-1000:]})
    comparisons = [compare_csv(ROOT / "results" / "tables" / name,
                               DEST / "results" / "tables" / name) for name in TABLES]
    models = sorted((DEST / "results" / "models").glob("08_hek293t_*_ridge.joblib"))
    if len(models) != 5:
        raise AssertionError("Expected five frozen Ridge models")
    status = {"status": "PASS", "level": "B_processed_input_to_historical_Ridge_and_frozen_HeLa_transfer",
              "python": sys.executable, "elapsed_seconds": time.perf_counter() - started,
              "input_sha256": sources, "models": [p.name for p in models],
              "tables": comparisons, "runs": runs}
    (DEST / "replay_validation.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "elapsed_seconds": status["elapsed_seconds"],
                      "models": len(models), "tables": len(comparisons),
                      "max_abs_numeric_difference": max(x["max_abs_numeric_difference"] for x in comparisons)}))


if __name__ == "__main__":
    main()
