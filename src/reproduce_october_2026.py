"""Replay the October 2026 paper's association and prediction tables."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from zipfile import ZIP_STORED, ZipFile

from replay_presubmission_stage4_models import (
    INPUTS as MATCHED_INPUTS,
    TABLES as MATCHED_TABLES,
    compare_csv,
    sha256,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = Path("data/final")
TABLES = Path("results/tables")
MANIFEST = "processed_input_sha256.json"
GEO_TRACKS = (
    "GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz",
    "GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities.txt.gz",
)
INPUTS = (
    "06_hek293t_site_level_dataset.csv.gz",
    "06_hek293t_main_analysis_dataset.csv.gz",
    "09_hela_site_level_dataset.csv.gz",
    "09_hela_main_analysis_dataset.csv.gz",
    "08a_hek293t_predicted_structure_features.csv.gz",
    "09a_hela_predicted_structure_features.csv.gz",
)
STAGES = {
    "07_run_association_analysis": (
        "07_association_analysis",
        ("07_hek293t_primary_association", "07_hek293t_sensitivity_associations"),
    ),
    "09b_run_hela_association": (
        "09b_hela_association",
        ("09b_hela_primary_association", "09b_hela_sensitivity_associations"),
    ),
    "12_run_meta_analysis": ("12_meta_analysis", ("12_meta_analysis",)),
    "08_run_prediction_modeling": (
        "08_prediction_modeling",
        ("08_hek293t_test_metrics", "08_hek293t_split_audit"),
    ),
    "09c_run_hela_model_transfer": (
        "09c_hela_model_transfer",
        ("09c_hela_transfer_metrics", "09c_hela_transfer_model_comparison"),
    ),
    "17_run_nonlinear_models": (
        "17_nonlinear_models",
        ("17_nonlinear_development_cv", "17_nonlinear_external_metrics"),
    ),
    "19_run_invivo_invitro_pairing": (
        "19_invivo_invitro_pairing",
        ("19_paired_condition_associations", "19_invivo_invitro_attrition"),
    ),
}


def pack(archive: Path) -> None:
    files = [DATA / name for name in INPUTS]
    checksums = {file.as_posix(): sha256(ROOT / file) for file in files}
    with ZipFile(archive, "x", compression=ZIP_STORED) as bundle:
        for file in files:
            bundle.write(ROOT / file, file.as_posix())
        bundle.writestr(MANIFEST, json.dumps(checksums, indent=2) + "\n")
    print(f"{sha256(archive)}  {archive}")


def unpack_inputs(archive: Path, destination: Path) -> None:
    expected = {f"{DATA.as_posix()}/{name}" for name in INPUTS}
    with ZipFile(archive) as bundle:
        checksums = json.loads(bundle.read(MANIFEST))
        if set(checksums) != expected or sorted(bundle.namelist()) != sorted(expected | {MANIFEST}):
            raise ValueError("Archive contents do not match the expected processed-input file list.")
        (destination / DATA).mkdir(parents=True)
        for relative, checksum in checksums.items():
            target = destination / relative
            with bundle.open(relative) as source, target.open("xb") as output:
                shutil.copyfileobj(source, output)
            if sha256(target) != checksum:
                raise ValueError(f"SHA-256 mismatch: {relative}")


def replay(archive: Path, destination: Path, scope: str, geo_directory: Path | None) -> None:
    stage_names = {
        "association": ("07_run_association_analysis", "09b_run_hela_association", "12_run_meta_analysis"),
        "prediction": ("08_run_prediction_modeling", "09c_run_hela_model_transfer"),
        "nonlinear": ("17_run_nonlinear_models",),
        "paired": ("19_run_invivo_invitro_pairing",),
    }
    destination.mkdir(parents=True, exist_ok=False)
    unpack_inputs(archive, destination)
    (destination / "src").mkdir()
    (destination / "config").mkdir()
    (destination / "results/reports").mkdir(parents=True)
    shutil.copy2(ROOT / "src/phase1_common.py", destination / "src")
    if scope == "paired":
        raw = destination / "data/raw_processed"
        raw.mkdir()
        for track in GEO_TRACKS:
            shutil.copy2(geo_directory / track, raw)
        shutil.copy2(ROOT / "src/pipeline_common.py", destination / "src")
    comparisons = []
    if scope == "matched":
        for relative in map(Path, MATCHED_INPUTS):
            if relative.parts[:2] != DATA.parts:
                target = destination / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT / relative, target)
        stage = "presubmission_stage2_matched_models"
        print(f"Running {stage}", flush=True)
        subprocess.run([sys.executable, str(destination / "src" / f"{stage}.py")], cwd=destination, check=True)
        for table in MATCHED_TABLES:
            relative = Path("results/presubmission_stage2") / table
            comparisons.append(compare_csv(ROOT / relative, destination / relative))
    for stage in stage_names.get(scope, ()):
        config, tables = STAGES[stage]
        script = destination / "src" / f"{stage}.py"
        shutil.copy2(ROOT / "src" / script.name, script)
        shutil.copy2(ROOT / "config" / f"{config}.yaml", destination / "config")
        print(f"Running {stage}", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=destination, check=True)
        for table in tables:
            relative = TABLES / f"{table}.csv"
            comparisons.append(compare_csv(ROOT / relative, destination / relative))
    summary = {
        "status": "PASS",
        "scope": scope,
        "archive_sha256": sha256(archive),
        "comparison_tolerance": {"relative": 1e-8, "absolute": 1e-8},
        "comparisons": comparisons,
    }
    with (destination / "validation.json").open("x", encoding="utf-8") as output:
        json.dump(summary, output, indent=2)
        output.write("\n")
    print(f"PASS: {len(comparisons)} tables. See {destination / 'validation.json'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    pack_parser = subcommands.add_parser("pack")
    pack_parser.add_argument("archive", type=Path)
    replay_parser = subcommands.add_parser("replay")
    replay_parser.add_argument("archive", type=Path)
    replay_parser.add_argument("destination", type=Path)
    replay_parser.add_argument(
        "--scope", choices=("association", "prediction", "nonlinear", "paired", "matched"), default="association"
    )
    replay_parser.add_argument("--geo-directory", type=Path, help="Directory holding the two GSE74353 tracks")
    args = parser.parse_args()
    if args.command == "pack":
        pack(args.archive)
    else:
        if args.scope == "paired" and args.geo_directory is None:
            replay_parser.error("--scope paired requires --geo-directory")
        replay(args.archive.resolve(), args.destination.resolve(), args.scope, args.geo_directory)


if __name__ == "__main__":
    main()
