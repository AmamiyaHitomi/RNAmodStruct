"""Compute deterministic partition-function structure features for stage 08A."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import RNA
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "08a_predicted_structure.yaml"
INPUT = ROOT / "data" / "final" / "06_hek293t_main_analysis_dataset.csv.gz"
CHECKPOINT = ROOT / "data" / "interim" / "08a_unique_sequence_predicted_structure.csv"
OUTPUT = ROOT / "data" / "final" / "08a_hek293t_predicted_structure_features.csv.gz"
STATUS = ROOT / "results" / "08a_hek293t_predicted_structure_run_status.csv"
LOG = ROOT / "results" / "logs" / "08a_compute_predicted_structure.log"


FIELDS = [
    "sequence_sha256", "sequence_length", "mfe_per_nt", "ensemble_free_energy_per_nt",
    "ensemble_diversity", "mean_pairing_state_entropy", "predicted_unpaired_mean_up10",
    "predicted_unpaired_mean_down10", "predicted_unpaired_mean_flank10",
    "predicted_unpaired_mean_far", "predicted_unpaired_probabilities_json",
]


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def fold_sequence(sequence: str, temperature: float, dangles: int) -> dict:
    md = RNA.md()
    md.temperature = temperature
    md.dangles = dangles
    md.noLP = 0
    compound = RNA.fold_compound(sequence, md)
    _, mfe = compound.mfe()
    _, ensemble_energy = compound.pf()
    bpp = compound.bpp()
    paired = np.zeros(len(sequence), dtype=float)
    for left in range(1, len(sequence) + 1):
        for right in range(left + 1, len(sequence) + 1):
            probability = float(bpp[left][right])
            paired[left - 1] += probability
            paired[right - 1] += probability
    unpaired = np.clip(1.0 - paired, 0.0, 1.0)
    epsilon = np.finfo(float).eps
    p = np.clip(paired, epsilon, 1.0 - epsilon)
    entropy = -(p * np.log(p) + (1.0 - p) * np.log(1.0 - p))
    up10 = unpaired[90:100]
    down10 = unpaired[101:111]
    far = np.concatenate([unpaired[50:90], unpaired[111:151]])
    return {
        "sequence_sha256": sequence_hash(sequence),
        "sequence_length": len(sequence),
        "mfe_per_nt": float(mfe / len(sequence)),
        "ensemble_free_energy_per_nt": float(ensemble_energy / len(sequence)),
        "ensemble_diversity": float(compound.mean_bp_distance()),
        "mean_pairing_state_entropy": float(entropy.mean()),
        "predicted_unpaired_mean_up10": float(up10.mean()),
        "predicted_unpaired_mean_down10": float(down10.mean()),
        "predicted_unpaired_mean_flank10": float(np.concatenate([up10, down10]).mean()),
        "predicted_unpaired_mean_far": float(far.mean()),
        "predicted_unpaired_probabilities_json": json.dumps(
            [round(float(value), 10) for value in unpaired], separators=(",", ":")
        ),
    }


def read_checkpoint() -> dict[str, dict]:
    if not CHECKPOINT.exists():
        return {}
    with CHECKPOINT.open(newline="", encoding="utf-8") as handle:
        return {row["sequence_sha256"]: row for row in csv.DictReader(handle)}


def append_checkpoint(rows: list[dict], write_header: bool) -> None:
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    with CHECKPOINT.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def write_status(row: dict) -> None:
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    with STATUS.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    setup_logging()
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    if RNA.__version__ != config["software"]["required_version"]:
        raise RuntimeError(f"ViennaRNA version {RNA.__version__} does not match frozen version")
    source = pd.read_csv(INPUT, compression="gzip", low_memory=False)
    source = source[source["model_dataset_included"].astype(str).eq("True")].copy()
    sequences = sorted(source["sequence_window_201"].astype(str).str.upper().unique())
    if any(len(sequence) != 201 or set(sequence) - set("ACGT") for sequence in sequences):
        raise ValueError("08A requires unambiguous complete 201-nt sequences")
    completed = read_checkpoint()
    logging.info("Unique sequences=%d checkpointed=%d", len(sequences), len(completed))
    buffer = []
    header_needed = not CHECKPOINT.exists() or CHECKPOINT.stat().st_size == 0
    for index, sequence in enumerate(sequences, 1):
        digest = sequence_hash(sequence)
        if digest not in completed:
            buffer.append(fold_sequence(
                sequence, float(config["software"]["temperature_celsius"]),
                int(config["software"]["dangles"]),
            ))
        if len(buffer) >= 25:
            append_checkpoint(buffer, header_needed)
            header_needed = False
            completed.update({row["sequence_sha256"]: row for row in buffer})
            buffer = []
        if index % 250 == 0:
            logging.info("Processed %d/%d sequences", index, len(sequences))
    if buffer:
        append_checkpoint(buffer, header_needed)
    features = pd.read_csv(CHECKPOINT)
    features = features[features["sequence_sha256"].isin({sequence_hash(s) for s in sequences})].copy()
    features = features.drop_duplicates("sequence_sha256", keep="first").sort_values("sequence_sha256")
    if len(features) != len(sequences):
        raise AssertionError("Checkpoint does not contain every unique sequence")
    mapping = source[["site_id", "sequence_window_201"]].copy()
    mapping["sequence_sha256"] = mapping["sequence_window_201"].astype(str).map(sequence_hash)
    output = mapping.drop(columns="sequence_window_201").merge(features, on="sequence_sha256", how="left", validate="many_to_one")
    if output.isna().any().any() or len(output) != len(source):
        raise AssertionError("Predicted features did not map one-to-one onto all model sites")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, index=False, compression={"method": "gzip", "mtime": 0})
    probability_lengths = output["predicted_unpaired_probabilities_json"].map(lambda value: len(json.loads(value)))
    status = {
        "stage": "08A", "status": "PASS", "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(), "sites": len(output),
        "unique_sequences": len(features), "probability_vector_length_min": int(probability_lengths.min()),
        "probability_vector_length_max": int(probability_lengths.max()), "viennarna_version": RNA.__version__,
        "temperature_celsius": config["software"]["temperature_celsius"],
        "input_sha256": sha256_file(INPUT), "config_sha256": sha256_file(CONFIG_PATH),
        "output_sha256": sha256_file(OUTPUT),
    }
    write_status(status)
    logging.info("08A completed: sites=%d unique_sequences=%d", len(output), len(features))


if __name__ == "__main__":
    main()
