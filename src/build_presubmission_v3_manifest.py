"""Hash the revision's code, inputs, cohorts, models, figures and PDFs."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from importlib import metadata


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "metadata" / "manifests" / "release" / "presubmission_v3_release_manifest.csv"
META = ROOT / "metadata" / "manifests" / "release" / "presubmission_v3_release_metadata.json"
ENV = ROOT / "results" / "presubmission_stage4" / "runtime_environment.json"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args: str) -> str:
    result = subprocess.run(["git", "-c", "safe.directory=E:/my_python/RNAmodStruct", *args],
                            cwd=ROOT, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def main() -> None:
    if any(p.exists() for p in (OUTPUT, META, ENV)):
        raise FileExistsError("Refusing to overwrite existing release manifest or environment record")
    import psutil
    versions = {name: metadata.version(name) for name in
                ("numpy", "pandas", "scipy", "scikit-learn", "statsmodels", "patsy", "matplotlib", "joblib", "PyYAML")}
    runtime = {"recorded_utc": datetime.now(timezone.utc).isoformat(), "python_executable": sys.executable,
               "python_version": sys.version, "platform": platform.platform(),
               "logical_cpu_count": psutil.cpu_count(logical=True),
               "physical_cpu_count": psutil.cpu_count(logical=False),
               "memory_bytes": psutil.virtual_memory().total, "packages": versions,
               "tex": "TeX Live 2026; pdflatex from C:/Users/ROG/texlive/2026/bin/windows"}
    ENV.parent.mkdir(parents=True, exist_ok=True)
    ENV.write_text(json.dumps(runtime, indent=2), encoding="utf-8")

    entries: list[tuple[str, Path]] = []
    def add(category: str, rel: str) -> None:
        path = ROOT / rel
        if not path.is_file():
            raise FileNotFoundError(path)
        entries.append((category, path))
    def add_tree(category: str, rel: str, allowed=None) -> None:
        base = ROOT / rel
        if not base.is_dir():
            raise FileNotFoundError(base)
        for path in sorted(base.rglob("*")):
            if path.is_file() and (allowed is None or path.suffix.lower() in allowed):
                entries.append((category, path))

    for rel in ("environment.yml", ".gitmodules", "metadata/provenance/environment.txt",
                "metadata/manifests/source/source_manifest.csv",
                "metadata/manifests/source/02_reference_manifest.csv",
                "metadata/manifests/download/download_checksums.csv",
                "config/08_prediction_modeling.yaml", "config/09c_hela_model_transfer.yaml",
                "config/09b_hela_association.yaml", "config/17_nonlinear_models.yaml",
                "config/19_invivo_invitro_pairing.yaml"):
        add("environment_or_config", rel)
    for rel in ("src/08_run_prediction_modeling.py", "src/09c_run_hela_model_transfer.py",
                "src/phase1_common.py", "src/freeze_presubmission_phase0.py",
                "src/presubmission_stage1_diagnostics.py", "src/presubmission_stage2_matched_models.py",
                "src/presubmission_stage3_association_context.py",
                "src/prepare_presubmission_manuscript_v3.py", "src/verify_presubmission_stage4_split.py",
                "src/replay_presubmission_stage4_models.py", "src/replay_presubmission_stage4_primary.py",
                "src/validate_presubmission_v3.py", "src/build_presubmission_v3_manifest.py"):
        add("code", rel)
    for rel in ("data/final/06_hek293t_main_analysis_dataset.csv.gz",
                "data/final/09_hela_main_analysis_dataset.csv.gz",
                "data/final/08a_hek293t_predicted_structure_features.csv.gz",
                "data/final/09a_hela_predicted_structure_features.csv.gz",
                "data/final/19_hek293t_invivo_invitro_paired.csv.gz",
                "results/tables/08_hek293t_group_assignments.csv",
                "results/tables/08_hek293t_feature_manifest.csv",
                "results/tables/08_hek293t_test_predictions.csv",
                "results/tables/09c_hela_transfer_predictions.csv"):
        add("processed_input_split_feature_or_frozen_prediction", rel)
    add_tree("phase0_cohort_and_plan", "metadata/manifests/presubmission_phase0")
    add_tree("historical_frozen_model", "results/models", {".joblib"})
    for n in (1, 2, 3):
        add_tree(f"presubmission_stage{n}", f"results/presubmission_stage{n}")
    for rel in ("results/presubmission_stage4/split_reproduction.json",
                "results/presubmission_stage4/replay/replay_validation.json",
                "results/presubmission_stage4/primary_replay/replay_validation.json",
                "results/presubmission_stage4/evidence_validation_v3.json",
                "results/presubmission_stage4/runtime_environment.json"):
        add("acceptance_record", rel)
    for rel in ("docs/manuscript_v3/main_v3.tex", "docs/manuscript_v3/supplement_v3.tex",
                "docs/manuscript_v3/main_v3.pdf", "docs/manuscript_v3/supplement_v3.pdf",
                "docs/manuscript_v3/references.bib", "docs/manuscript_v3/build_assets.py",
                "docs/manuscript_v3/build_v3.ps1", "docs/manuscript_v3/README.md",
                "docs/manuscript_v3/source_data_manifest.json"):
        add("revision_package", rel)
    add_tree("revision_figure", "docs/manuscript_v3/figures", {".pdf", ".svg", ".png", ".tiff"})
    add_tree("revision_table", "docs/manuscript_v3/tables", {".tex"})
    add_tree("revision_frozen_source_copy", "docs/manuscript_v3/source_data", {".csv"})

    unique = {}
    for category, path in entries:
        rel = path.relative_to(ROOT).as_posix()
        if rel in unique:
            raise AssertionError(f"Duplicate release entry: {rel}")
        unique[rel] = {"category": category, "path": rel, "bytes": path.stat().st_size,
                       "sha256": digest(path)}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["category", "path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(unique.values())
    meta = {"created_utc": datetime.now(timezone.utc).isoformat(),
            "git_head": git("rev-parse", "HEAD"),
            "working_tree_clean": not bool(git("status", "--short")),
            "submodule_status": git("submodule", "status"),
            "entry_count": len(unique), "manifest_sha256": digest(OUTPUT),
            "validation_role": "presubmission_revision_with_posthoc_sensitivities",
            "historical_primary_models": "five stage08 joblib bundles preserved under results/models",
            "model_count": 21}
    META.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "entries": len(unique), "sha256": meta["manifest_sha256"]}))


if __name__ == "__main__":
    main()
