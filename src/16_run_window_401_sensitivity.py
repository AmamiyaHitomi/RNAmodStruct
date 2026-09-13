"""Stage 16: rebuild complete 401-nt sequence/experimental/ViennaRNA windows and model them."""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import math
import tarfile
from array import array
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import RNA
import yaml
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

from phase1_common import TabularEncoder, make_joint_groups, model_metrics, sha256_file, write_rows
from pipeline_common import iter_fasta

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "16_window_401_sensitivity.yaml"
TABLE_CV = ROOT / "results" / "tables" / "16_window_401_development_cv.csv"
TABLE_EXTERNAL = ROOT / "results" / "tables" / "16_window_401_external_metrics.csv"
TABLE_ATTRITION = ROOT / "results" / "tables" / "16_window_401_attrition.csv"
STATUS = ROOT / "results" / "status" / "16_window_401_sensitivity_run_status.csv"
REPORT = ROOT / "results" / "reports" / "16_window_401_sensitivity_report.md"
DATASETS = {name: ROOT / "data" / "final" / f"16_{name.lower()}_401nt_dataset.csv.gz" for name in ["HEK293T", "HeLa"]}
CACHE = ROOT / "data" / "interim" / "16_401nt_vienna_cache.csv"
MISSING = {"NULL", "NA", "N/A", "NaN", "nan", ".", "-999"}


def load_fasta_selected(path: Path, wanted: set[str]) -> dict[str, str]:
    aliases = wanted | {value.split(".", 1)[0] for value in wanted}
    found = {}
    for name, sequence in iter_fasta(path):
        bare = name.split(".", 1)[0]
        if name in aliases or bare in aliases:
            found[name] = sequence
            found.setdefault(bare, sequence)
    return found


def parse_icshape_stream(handle, wanted: set[str]) -> dict[str, tuple[array, str]]:
    result = {}
    aliases = wanted | {value.split(".", 1)[0] for value in wanted}
    for raw in handle:
        line = raw.decode("utf-8") if isinstance(raw, bytes) else raw
        fields = line.rstrip("\n").split("\t")
        transcript = fields[0]
        bare = transcript.split(".", 1)[0]
        if transcript not in aliases and bare not in aliases:
            continue
        values = array("f", (math.nan if value in MISSING else float(value) for value in fields[3:]))
        result[transcript] = (values, fields[2])
        result.setdefault(bare, (values, fields[2]))
    return result


def load_icshape(item: dict, wanted: set[str]) -> dict[str, tuple[array, str]]:
    path = ROOT / item["icshape"]
    if item["icshape_format"] == "gzip_text":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return parse_icshape_stream(handle, wanted)
    with tarfile.open(path, "r") as archive:
        member = archive.extractfile(item["icshape_member"])
        if member is None:
            raise FileNotFoundError(item["icshape_member"])
        with gzip.GzipFile(fileobj=member) as handle:
            return parse_icshape_stream(handle, wanted)


def extract_dataset(dataset: str, item: dict, half: int, minimum: float) -> tuple[pd.DataFrame, list[dict]]:
    source = pd.read_csv(ROOT / item["input"], low_memory=False)
    source = source[source["model_dataset_included"].astype(str).eq("True")].copy()
    wanted = set(source["transcript_id"].astype(str))
    sequences = load_fasta_selected(ROOT / item["transcript_fasta"], wanted)
    reactivities = load_icshape(item, wanted)
    records, reasons = [], {"missing_transcript_sequence": 0, "incomplete_401nt_sequence": 0,
                            "missing_icshape_transcript": 0, "insufficient_401nt_reactivity": 0}
    for _, row in source.iterrows():
        transcript = str(row["transcript_id"])
        bare = transcript.split(".", 1)[0]
        sequence = sequences.get(transcript, sequences.get(bare))
        if sequence is None:
            reasons["missing_transcript_sequence"] += 1; continue
        center = int(row["transcript_pos_1based"]) - 1
        if center - half < 0 or center + half >= len(sequence):
            reasons["incomplete_401nt_sequence"] += 1; continue
        window = sequence[center - half:center + half + 1].upper().replace("U", "T")
        if len(window) != 2 * half + 1 or set(window) - set("ACGT"):
            reasons["incomplete_401nt_sequence"] += 1; continue
        payload = reactivities.get(transcript, reactivities.get(bare))
        if payload is None:
            reasons["missing_icshape_transcript"] += 1; continue
        values, _ = payload
        left = np.asarray(values[center - half:center], dtype=float)
        right = np.asarray(values[center + 1:center + half + 1], dtype=float)
        if len(left) != half or len(right) != half or np.isfinite(left).mean() < minimum or np.isfinite(right).mean() < minimum:
            reasons["insufficient_401nt_reactivity"] += 1; continue
        record = row.to_dict()
        record["sequence_window_401"] = window
        record["sequence_sha256_401"] = hashlib.sha256(window.encode("ascii")).hexdigest()
        record["coverage_up200"] = float(np.isfinite(left).mean())
        record["coverage_down200"] = float(np.isfinite(right).mean())
        record["reactivity_mean_flank200"] = float(np.nanmean(np.concatenate([left, right])))
        records.append(record)
    frame = pd.DataFrame(records)
    attrition = [{"dataset": dataset, "reason": "source_model_population", "sites": len(source)},
                 {"dataset": dataset, "reason": "retained_complete_401", "sites": len(frame)}]
    attrition += [{"dataset": dataset, "reason": key, "sites": value} for key, value in reasons.items()]
    return frame, attrition


def fold_401(sequence: str, temperature: float, dangles: int) -> dict:
    md = RNA.md(); md.temperature = temperature; md.dangles = dangles; md.noLP = 0
    compound = RNA.fold_compound(sequence, md)
    _, mfe = compound.mfe(); _, ensemble = compound.pf(); bpp = compound.bpp()
    paired = np.zeros(len(sequence))
    for left in range(1, len(sequence) + 1):
        for right in range(left + 1, len(sequence) + 1):
            p = float(bpp[left][right]); paired[left - 1] += p; paired[right - 1] += p
    p = np.clip(paired, np.finfo(float).eps, 1 - np.finfo(float).eps)
    entropy = -(p * np.log(p) + (1 - p) * np.log(1 - p))
    unpaired = 1 - paired
    return {"sequence_sha256_401": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "mfe_per_nt_401": float(mfe / len(sequence)), "ensemble_free_energy_per_nt_401": float(ensemble / len(sequence)),
            "ensemble_diversity_401": float(compound.mean_bp_distance()), "mean_pairing_state_entropy_401": float(entropy.mean()),
            "predicted_unpaired_mean_flank200": float(np.concatenate([unpaired[:200], unpaired[201:]]).mean())}


def attach_vienna(frames: list[pd.DataFrame], config: dict) -> list[pd.DataFrame]:
    completed = pd.read_csv(CACHE) if CACHE.exists() else pd.DataFrame()
    known = set(completed["sequence_sha256_401"]) if len(completed) else set()
    sequences = {hashlib.sha256(seq.encode("ascii")).hexdigest(): seq for frame in frames for seq in frame["sequence_window_401"].unique()}
    pending = [(digest, seq) for digest, seq in sorted(sequences.items()) if digest not in known]
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    buffer = []
    for _, sequence in pending:
        buffer.append(fold_401(sequence, float(config["predicted_structure"]["temperature_celsius"]),
                               int(config["predicted_structure"]["dangles"])))
        if len(buffer) >= 50:
            pd.DataFrame(buffer).to_csv(CACHE, mode="a", index=False,
                                        header=not CACHE.exists() or CACHE.stat().st_size == 0)
            buffer = []
    if buffer:
        pd.DataFrame(buffer).to_csv(CACHE, mode="a", index=False,
                                    header=not CACHE.exists() or CACHE.stat().st_size == 0)
    features = pd.read_csv(CACHE)
    features = features.drop_duplicates("sequence_sha256_401", keep="first")
    if not set(sequences).issubset(set(features["sequence_sha256_401"])):
        raise AssertionError("ViennaRNA checkpoint is incomplete")
    return [frame.merge(features, on="sequence_sha256_401", validate="many_to_one") for frame in frames]


class Encoder401(TabularEncoder):
    predicted_names = ["mfe_per_nt_401", "ensemble_free_energy_per_nt_401", "ensemble_diversity_401",
                       "mean_pairing_state_entropy_401", "predicted_unpaired_mean_flank200"]

    def __init__(self):
        super().__init__("sequence_window_401")
        self.vienna_mean = np.zeros(len(self.predicted_names)); self.vienna_scale = np.ones(len(self.predicted_names))
        self.r401_mean = 0.0; self.r401_scale = 1.0

    def fit(self, df):
        super().fit(df)
        raw = df[self.predicted_names].to_numpy(float); self.vienna_mean = raw.mean(0)
        self.vienna_scale = np.where(raw.std(0) > 0, raw.std(0), 1.0)
        self.r401_mean, self.r401_scale = (float(df["reactivity_mean_flank200"].mean()), float(df["reactivity_mean_flank200"].std(ddof=0)) or 1.0)
        return self

    def transform(self, df, model):
        common, common_names = self.common(df); arrays, names = [common], common_names
        if model != "M0":
            sequences = df["sequence_window_401"].astype(str).tolist()
            onehot = np.zeros((len(df), 401 * 4), dtype=np.float32); lookup = {x: i for i, x in enumerate("ACGT")}
            for row, sequence in enumerate(sequences):
                for pos, base in enumerate(sequence): onehot[row, pos * 4 + lookup[base]] = 1
            kmer = (self.kmer_matrix(df) - self.kmer_mean) / self.kmer_scale
            arrays += [onehot, kmer]; names += [f"X401:pos{p-200:+04d}_{b}" for p in range(401) for b in "ACGT"] + [f"X401:3mer_{x}" for x in self.kmers]
        if model in {"M2", "M4"}:
            arrays.append(((df[self.predicted_names].to_numpy(float) - self.vienna_mean) / self.vienna_scale).astype(np.float32)); names += [f"P401:{x}" for x in self.predicted_names]
        if model in {"M3", "M4"}:
            arrays.append(((df["reactivity_mean_flank200"].to_numpy(float) - self.r401_mean) / self.r401_scale)[:, None].astype(np.float32)); names.append("R401:reactivity_mean_flank200")
        return np.hstack(arrays), names


def main() -> None:
    started = datetime.now(timezone.utc); config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    if RNA.__version__ != config["predicted_structure"]["required_version"]:
        raise RuntimeError(f"ViennaRNA version mismatch: {RNA.__version__}")
    frames, attrition = [], []
    for dataset, item in config["datasets"].items():
        frame, rows = extract_dataset(dataset, item, int(config["half_window"]), float(config["minimum_valid_fraction_each_side"]))
        frames.append(frame); attrition += rows
    frames = attach_vienna(frames, config)
    hek, hela = frames
    for dataset, frame in zip(config["datasets"], frames):
        frame.to_csv(DATASETS[dataset], index=False, compression={"method": "gzip", "mtime": 0})
    hek = hek.reset_index(drop=True); hela = hela.reset_index(drop=True)
    hek["joint_group"] = make_joint_groups(hek, "sequence_window_401")
    folds = int(config["learner"]["development_folds"]); alphas = [float(x) for x in config["learner"]["alpha_candidates"]]
    cv_rows, selected = [], {}
    for model in config["models"]:
        for alpha in alphas:
            scores = []
            for fold, (train, valid) in enumerate(GroupKFold(folds).split(hek, groups=hek["joint_group"]), 1):
                encoder = Encoder401().fit(hek.iloc[train]); x_train, _ = encoder.transform(hek.iloc[train], model); x_valid, _ = encoder.transform(hek.iloc[valid], model)
                prediction = np.clip(Ridge(alpha=alpha, solver="lsqr").fit(x_train, hek.iloc[train]["combined_ratio"]).predict(x_valid), 0, 1)
                metrics = model_metrics(hek.iloc[valid]["combined_ratio"].to_numpy(), prediction); scores.append(metrics)
                cv_rows.append({"model": model, "alpha": alpha, "fold": fold, **metrics})
            cv_rows.append({"model": model, "alpha": alpha, "fold": "mean", **{k: float(np.mean([x[k] for x in scores])) for k in scores[0]}})
        selected[model] = min((r for r in cv_rows if r["model"] == model and r["fold"] == "mean"), key=lambda r: r["mae"])["alpha"]
    external = []
    for model in config["models"]:
        encoder = Encoder401().fit(hek); x_hek, _ = encoder.transform(hek, model); x_hela, _ = encoder.transform(hela, model)
        prediction = np.clip(Ridge(alpha=selected[model], solver="lsqr").fit(x_hek, hek["combined_ratio"]).predict(x_hela), 0, 1)
        external.append({"dataset": "HeLa", "model": model, "selected_alpha_from_HEK293T": selected[model], **model_metrics(hela["combined_ratio"].to_numpy(), prediction)})
    write_rows(TABLE_CV, cv_rows); write_rows(TABLE_EXTERNAL, external); write_rows(TABLE_ATTRITION, attrition)
    write_rows(STATUS, [{"stage": 16, "status": "PASS", "started_utc": started.isoformat(), "finished_utc": datetime.now(timezone.utc).isoformat(),
                         "hek293t_sites": len(hek), "hela_sites": len(hela), "models": len(external), "viennarna_version": RNA.__version__,
                         "config_sha256": sha256_file(CONFIG), "cv_sha256": sha256_file(TABLE_CV), "external_sha256": sha256_file(TABLE_EXTERNAL)}])
    m1 = next(x for x in external if x["model"] == "M1"); m4 = next(x for x in external if x["model"] == "M4")
    REPORT.write_text(f"""# Stage 16: 401 nt window sensitivity

## Material Passport

- Origin: phase-1 exploratory extension
- Verification Status: ANALYZED

The full 401 nt window (unpopulated) retains HEK293T {len(hek):,} and HeLa {len(hela):,} sites. Sequence one-hot, 3-mer, ViennaRNA features reconstructed simultaneously with experimental structure windows; grouping based on genes and identical 401 nt sequence connected components. Hyperparameters are only selected in HEK293T grouped cross-validation, and HeLa only performs frozen external review.

HeLa's M1 MAE={m1['mae']:.5f}, M4 MAE={m4['mae']:.5f}, ΔMAE(M1−M4)={m1['mae']-m4['mae']:.6f}. This result is a 401 nt sensitivity analysis only and does not replace the 201 nt master model.""", encoding="utf-8")


if __name__ == "__main__":
    main()
