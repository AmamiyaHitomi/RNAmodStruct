"""Freeze auditable HeLa cohorts for the main_v2 presubmission analyses.

Creates new files only. Existing outputs are never overwritten.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "metadata" / "manifests" / "presubmission_phase0"
SOURCES = {
    "hek_main": ROOT / "data/final/06_hek293t_main_analysis_dataset.csv.gz",
    "hela_main": ROOT / "data/final/09_hela_main_analysis_dataset.csv.gz",
    "hek_assignments": ROOT / "results/tables/08_hek293t_group_assignments.csv",
    "hela_predictions": ROOT / "results/tables/09c_hela_transfer_predictions.csv",
    "hek_config": ROOT / "config/08_prediction_modeling.yaml",
    "hela_config": ROOT / "config/09c_hela_model_transfer.yaml",
    "nonlinear_config": ROOT / "config/17_nonlinear_models.yaml",
}
COHORTS = {
    "all_qualifying": "all_qualifying_site_ids.csv",
    "development_overlap_excluded": "development_overlap_excluded_site_ids.csv",
    "hek_main_4409_overlap_excluded": "hek_main_4409_overlap_excluded_site_ids.csv",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_sha256(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def id_sha256(ids: list[str]) -> str:
    """Hash sorted site IDs, one UTF-8 ID plus LF per line, including final LF."""
    return hashlib.sha256("".join(f"{item}\n" for item in ids).encode("utf-8")).hexdigest()


def read_main(path: Path) -> pd.DataFrame:
    cols = [
        "site_id", "analysis_gene", "sequence_window_201", "sequence_window_201_full",
        "main_analysis_included", "model_dataset_included",
    ]
    frame = pd.read_csv(path, usecols=cols, low_memory=False)
    if frame["site_id"].isna().any() or not frame["site_id"].is_unique:
        raise ValueError(f"Missing or duplicate site_id in {path}")
    if frame["analysis_gene"].isna().any():
        raise ValueError(f"Missing analysis_gene in {path}")
    return frame


def gene_and_sequence_sets(frame: pd.DataFrame) -> tuple[set[str], set[str]]:
    genes = set(frame["analysis_gene"].astype(str))
    full = frame[frame["sequence_window_201_full"].eq(True)]
    sequences = full["sequence_window_201"].astype(str)
    if not sequences.str.len().eq(201).all():
        raise ValueError("A full sequence is not 201 nt")
    return genes, {sequence_sha256(item) for item in sequences}


def write_ids(path: Path, ids: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["site_id"])
        writer.writerows((item,) for item in ids)


def main() -> None:
    planned = [OUT / "hela_site_cohort_flags.csv", OUT / "freeze_manifest.json"]
    planned.extend(OUT / name for name in COHORTS.values())
    existing = [str(path) for path in planned if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite: " + ", ".join(existing))

    hek = read_main(SOURCES["hek_main"])
    hela = read_main(SOURCES["hela_main"])
    assignments = pd.read_csv(SOURCES["hek_assignments"])
    predictions = pd.read_csv(SOURCES["hela_predictions"])
    if len(hek) != 4409 or int(hek["model_dataset_included"].sum()) != 4096:
        raise ValueError("Unexpected HEK main/model cohort size")
    if len(hela) != 25996 or int(hela["model_dataset_included"].sum()) != 24960:
        raise ValueError("Unexpected HeLa main/model cohort size")
    if len(assignments) != 4096 or not assignments["site_id"].is_unique:
        raise ValueError("Unexpected HEK assignment size or duplicate site_id")
    if len(predictions) != 24960 or not predictions["site_id"].is_unique:
        raise ValueError("Unexpected HeLa prediction size or duplicate site_id")
    if set(assignments["site_id"].astype(str)) != set(hek.loc[hek["model_dataset_included"], "site_id"].astype(str)):
        raise ValueError("HEK model population disagrees with group assignments")
    if set(predictions["site_id"].astype(str)) != set(hela.loc[hela["model_dataset_included"], "site_id"].astype(str)):
        raise ValueError("HeLa model population disagrees with frozen predictions")

    dev = assignments[assignments["split"].eq("development")]
    if len(dev) != 3271:
        raise ValueError("Unexpected HEK development size")
    dev_genes = set(dev["analysis_gene"].astype(str))
    dev_sequences = set(dev["sequence_sha256"].astype(str))
    model_genes, model_sequences = gene_and_sequence_sets(hek[hek["model_dataset_included"]])
    main_genes, main_sequences = gene_and_sequence_sets(hek[hek["main_analysis_included"]])

    frame = predictions[["site_id", "analysis_gene", "joint_group", "sequence_sha256", "overlaps_development"]].copy()
    frame["site_id"] = frame["site_id"].astype(str)
    frame.sort_values("site_id", inplace=True)
    frame.reset_index(drop=True, inplace=True)
    frame["dev_gene_overlap"] = frame["analysis_gene"].astype(str).isin(dev_genes)
    frame["dev_sequence_overlap"] = frame["sequence_sha256"].astype(str).isin(dev_sequences)
    frame["hek_model_gene_overlap"] = frame["analysis_gene"].astype(str).isin(model_genes)
    frame["hek_model_sequence_overlap"] = frame["sequence_sha256"].astype(str).isin(model_sequences)
    frame["hek_main_gene_overlap"] = frame["analysis_gene"].astype(str).isin(main_genes)
    frame["hek_main_sequence_overlap"] = frame["sequence_sha256"].astype(str).isin(main_sequences)
    observed_dev = frame["dev_gene_overlap"] | frame["dev_sequence_overlap"]
    if not observed_dev.eq(frame["overlaps_development"].astype(bool)).all():
        raise ValueError("Development-overlap flags disagree with frozen transfer table")
    hela_model = hela.loc[hela["model_dataset_included"], ["site_id", "sequence_window_201"]].copy()
    merged = frame.merge(hela_model, on="site_id", how="left", validate="one_to_one")
    recomputed = merged["sequence_window_201"].astype(str).map(sequence_sha256)
    if not recomputed.eq(frame["sequence_sha256"]).all():
        raise ValueError("HeLa sequence hashes disagree with the source dataset")

    frame["all_qualifying"] = True
    frame["development_overlap_excluded"] = ~observed_dev
    frame["hek_main_4409_overlap_excluded"] = ~(
        frame["hek_main_gene_overlap"] | frame["hek_main_sequence_overlap"]
    )
    if int(frame["development_overlap_excluded"].sum()) != 22256:
        raise ValueError("Existing strict cohort count did not reproduce")
    if not frame.loc[frame["hek_main_4409_overlap_excluded"], "development_overlap_excluded"].all():
        raise ValueError("HEK main exclusion is not nested within development exclusion")

    OUT.mkdir(parents=True, exist_ok=False)
    flags_path = OUT / "hela_site_cohort_flags.csv"
    with flags_path.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, index=False, lineterminator="\n")

    cohort_manifest = {}
    for key, name in COHORTS.items():
        ids = frame.loc[frame[key], "site_id"].tolist()
        path = OUT / name
        write_ids(path, ids)
        cohort_manifest[key] = {
            "site_count": len(ids),
            "joint_group_count": int(frame.loc[frame[key], "joint_group"].nunique()),
            "site_ids_sha256": id_sha256(ids),
            "csv_sha256": file_sha256(path),
            "csv_path": str(path.relative_to(ROOT)).replace("\\", "/"),
        }

    overlap_summary = {}
    for prefix in ["dev", "hek_model", "hek_main"]:
        gene = frame[f"{prefix}_gene_overlap"]
        seq = frame[f"{prefix}_sequence_overlap"]
        overlap_summary[prefix] = {
            "gene": int(gene.sum()),
            "identical_201nt_sequence": int(seq.sum()),
            "both": int((gene & seq).sum()),
            "union_excluded": int((gene | seq).sum()),
        }
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Phase 0 local presubmission freeze; new HEK-main exclusion is post hoc",
        "hash_rule": "Sort site_id lexicographically; hash UTF-8 IDs, each followed by LF, including the final LF",
        "source_sha256": {key: file_sha256(path) for key, path in SOURCES.items()},
        "source_paths": {key: str(path.relative_to(ROOT)).replace("\\", "/") for key, path in SOURCES.items()},
        "flags_csv_sha256": file_sha256(flags_path),
        "cohorts": cohort_manifest,
        "overlap_summary": overlap_summary,
        "bootstrap_plan": {
            "unit": "HeLa joint_group from the frozen all-qualifying transfer table",
            "resamples": 2000,
            "seed": 20260912,
            "pairing": "Within a cohort, use the same sampled group-index vectors for every prespecified contrast",
            "interval": "percentile 95%",
            "CI_scope": "fixed_fitted_models",
        },
    }
    with (OUT / "freeze_manifest.json").open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(json.dumps({"cohorts": cohort_manifest, "overlap_summary": overlap_summary}, indent=2))


if __name__ == "__main__":
    main()
