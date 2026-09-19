"""Isolated B-level replay of the fixed stage-2 modeling comparison."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "results" / "presubmission_stage4" / "replay"
INPUTS = [
    "src/presubmission_stage2_matched_models.py",
    "src/phase1_common.py",
    "src/08_run_prediction_modeling.py",
    "config/17_nonlinear_models.yaml",
    "data/final/06_hek293t_main_analysis_dataset.csv.gz",
    "data/final/09_hela_main_analysis_dataset.csv.gz",
    "data/final/08a_hek293t_predicted_structure_features.csv.gz",
    "data/final/09a_hela_predicted_structure_features.csv.gz",
    "results/tables/08_hek293t_group_assignments.csv",
    "results/tables/08_hek293t_test_predictions.csv",
    "results/tables/09c_hela_transfer_predictions.csv",
    "metadata/manifests/presubmission_phase0/analysis_spec.json",
    "metadata/manifests/presubmission_phase0/freeze_manifest.json",
    "metadata/manifests/presubmission_phase0/hela_site_cohort_flags.csv",
]
TABLES = ["development_cv_candidates.csv", "selected_models.csv", "evaluation_metrics.csv",
          "paired_contrasts.csv", "HEK_reused_test_predictions.csv", "HeLa_predictions.csv"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_csv(reference: Path, replay: Path) -> dict:
    left, right = pd.read_csv(reference, low_memory=False), pd.read_csv(replay, low_memory=False)
    if left.shape != right.shape or list(left.columns) != list(right.columns):
        raise AssertionError(f"Table schema changed: {reference.name}")
    max_abs = 0.0
    for name in left:
        if pd.api.types.is_numeric_dtype(left[name]):
            a, b = left[name].to_numpy(dtype=float), right[name].to_numpy(dtype=float)
            if not np.array_equal(np.isnan(a), np.isnan(b)):
                raise AssertionError(f"Missing numeric values differ: {reference.name}/{name}")
            delta = np.abs(a - b)
            delta = delta[np.isfinite(delta)]
            if len(delta):
                max_abs = max(max_abs, float(delta.max()))
            if not np.allclose(a, b, atol=1e-8, rtol=1e-8, equal_nan=True):
                raise AssertionError(f"Numeric values differ: {reference.name}/{name}")
        elif not left[name].fillna("<NA>").astype(str).equals(right[name].fillna("<NA>").astype(str)):
            raise AssertionError(f"String values differ: {reference.name}/{name}")
    return {"table": reference.name, "rows": len(left), "max_abs_numeric_difference": max_abs,
            "reference_sha256": sha256(reference), "replay_sha256": sha256(replay)}


def main() -> None:
    if DEST.exists():
        raise FileExistsError(f"Refusing to overwrite {DEST}")
    started = time.perf_counter()
    source_hashes = {rel: sha256(ROOT / rel) for rel in INPUTS}
    for rel in INPUTS:
        target = DEST / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, target)
        if sha256(target) != source_hashes[rel]:
            raise AssertionError(f"Copy hash mismatch: {rel}")
    run = subprocess.run([sys.executable, str(DEST / "src" / "presubmission_stage2_matched_models.py")],
                         cwd=DEST, capture_output=True, text=True, check=False)
    if run.returncode:
        raise RuntimeError(f"Isolated replay failed ({run.returncode}): {run.stdout[-3000:]} {run.stderr[-3000:]}")
    reference = ROOT / "results" / "presubmission_stage2"
    replay = DEST / "results" / "presubmission_stage2"
    comparisons = [compare_csv(reference / name, replay / name) for name in TABLES]
    ref_models = sorted((reference / "models").glob("*.joblib"))
    new_models = sorted((replay / "models").glob("*.joblib"))
    if len(ref_models) != len(new_models) or [p.name for p in ref_models] != [p.name for p in new_models]:
        raise AssertionError("Saved candidate model set differs")
    status = {"status": "PASS", "level": "B_processed_input_to_matched_models_and_metrics",
              "python": sys.executable, "elapsed_seconds": time.perf_counter() - started,
              "inputs_sha256": source_hashes, "tables": comparisons,
              "saved_model_count": len(new_models), "replay_stdout": run.stdout.strip()}
    (DEST / "replay_validation.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "elapsed_seconds": status["elapsed_seconds"],
                      "tables": len(comparisons), "models": len(new_models),
                      "max_abs_numeric_difference": max(x["max_abs_numeric_difference"] for x in comparisons)}))


if __name__ == "__main__":
    main()
