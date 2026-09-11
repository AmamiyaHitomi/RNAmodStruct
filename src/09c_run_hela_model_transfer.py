"""Transfer the frozen HEK293T stage-08 Ridge models directly onto HeLa (stage 09C)."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import logging
import math
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "09c_hela_model_transfer.yaml"
INPUT = ROOT / "data" / "final" / "09_hela_main_analysis_dataset.csv.gz"
PREDICTED_INPUT = ROOT / "data" / "final" / "09a_hela_predicted_structure_features.csv.gz"
MODELS_DIR = ROOT / "results" / "models"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
LOG = RESULTS / "logs" / "09c_run_hela_model_transfer.log"

ALPHABET = "ACGT"
KMERS = ["".join(parts) for parts in itertools.product("ACGT", repeat=3)]


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows(rows)


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(str(sequence).encode("ascii")).hexdigest()


class UnionFind:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def find(self, value):
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            value, self.parent[value] = self.parent[value], root
        return root

    def union(self, left, right):
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def make_joint_groups(df: pd.DataFrame) -> pd.Series:
    genes = sorted(df["analysis_gene"].astype(str).unique())
    uf = UnionFind(genes)
    for _, values in df.groupby("sequence_window_201")["analysis_gene"]:
        unique = sorted(values.astype(str).unique())
        for gene in unique[1:]:
            uf.union(unique[0], gene)
    roots = {gene: uf.find(gene) for gene in genes}
    ordered = {root: f"G{index:04d}" for index, root in enumerate(sorted(set(roots.values())), 1)}
    return df["analysis_gene"].astype(str).map(lambda gene: ordered[roots[gene]])


class FrozenEncoder:
    """Reconstruct the frozen stage-08 FeatureEncoder from its serialized state."""

    numeric = [
        "gc_fraction_21", "transcript_position_fraction", "log1p_icshape_abundance_rpkm",
        "log1p_combined_agcov", "coverage_flank10", "distance_to_stop_codon_tx",
        "distance_to_nearest_splice_edge_tx",
    ]
    missing = ["distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx"]
    categorical = ["drach_subtype", "transcript_region"]
    predicted_scalars = [
        "mfe_per_nt", "ensemble_free_energy_per_nt", "ensemble_diversity",
        "mean_pairing_state_entropy", "predicted_unpaired_mean_up10",
        "predicted_unpaired_mean_down10", "predicted_unpaired_mean_flank10",
        "predicted_unpaired_mean_far",
    ]

    def __init__(self, state: dict):
        self.medians = {name: float(state["numeric_medians"][name]) for name in self.numeric}
        self.means = {name: float(state["numeric_means"][name]) for name in self.numeric}
        self.scales = {name: float(state["numeric_scales"][name]) for name in self.numeric}
        self.categories = {name: list(state["categorical_levels"][name]) for name in self.categorical}
        self.kmer_means = np.asarray(state["kmer_means"], dtype=float)
        self.kmer_scales = np.asarray(state["kmer_scales"], dtype=float)
        self.structure_mean = float(state["structure_mean"])
        self.structure_scale = float(state["structure_scale"])
        self.predicted_means = np.asarray(state["predicted_means"], dtype=float)
        self.predicted_scales = np.asarray(state["predicted_scales"], dtype=float)

    @staticmethod
    def _numeric_frame(df: pd.DataFrame) -> pd.DataFrame:
        return pd.DataFrame({
            "gc_fraction_21": pd.to_numeric(df["gc_fraction_21"], errors="coerce"),
            "transcript_position_fraction": pd.to_numeric(df["transcript_position_fraction"], errors="coerce"),
            "log1p_icshape_abundance_rpkm": np.log1p(pd.to_numeric(df["icshape_abundance_rpkm"], errors="coerce")),
            "log1p_combined_agcov": np.log1p(pd.to_numeric(df["combined_agcov"], errors="coerce")),
            "coverage_flank10": pd.to_numeric(df["coverage_flank10"], errors="coerce"),
            "distance_to_stop_codon_tx": pd.to_numeric(df["distance_to_stop_codon_tx"], errors="coerce"),
            "distance_to_nearest_splice_edge_tx": pd.to_numeric(df["distance_to_nearest_splice_edge_tx"], errors="coerce"),
        }, index=df.index)

    def common(self, df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
        numeric = self._numeric_frame(df)
        arrays, names = [], []
        for name in self.numeric:
            value = numeric[name].fillna(self.medians[name]).to_numpy(dtype=float)
            arrays.append(((value - self.means[name]) / self.scales[name])[:, None])
            names.append(f"C:z_{name}")
        for name in self.missing:
            arrays.append(numeric[name].isna().astype(float).to_numpy()[:, None])
            names.append(f"C:missing_{name}")
        for name in self.categorical:
            values = df[name].fillna("missing").astype(str).to_numpy()
            for category in self.categories[name]:
                arrays.append((values == category).astype(float)[:, None])
                names.append(f"C:{name}={category}")
        return np.hstack(arrays).astype(np.float32), names

    def _kmer_matrix(self, df: pd.DataFrame) -> np.ndarray:
        sequences = df["sequence_window_201"].astype(str).str.upper().tolist()
        if any(len(seq) != 201 for seq in sequences):
            raise ValueError("All model sequences must contain exactly 201 nt")
        kmer = np.zeros((len(sequences), 64), dtype=np.float32)
        klookup = {value: index for index, value in enumerate(KMERS)}
        for row, seq in enumerate(sequences):
            for pos in range(199):
                value = seq[pos:pos + 3]
                if value in klookup:
                    kmer[row, klookup[value]] += 1.0 / 199.0
        return kmer

    def sequence(self, df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
        sequences = df["sequence_window_201"].astype(str).str.upper().tolist()
        if any(len(seq) != 201 for seq in sequences):
            raise ValueError("All model sequences must contain exactly 201 nt")
        onehot = np.zeros((len(sequences), 201 * 4), dtype=np.float32)
        lookup = {base: index for index, base in enumerate(ALPHABET)}
        for row, seq in enumerate(sequences):
            for pos, base in enumerate(seq):
                if base in lookup:
                    onehot[row, pos * 4 + lookup[base]] = 1.0
        kmer = (self._kmer_matrix(df) - self.kmer_means) / self.kmer_scales
        names = [f"X:pos{pos:+04d}_{base}" for pos in range(-100, 101) for base in ALPHABET]
        names += [f"X:3mer_{value}" for value in KMERS]
        return np.hstack([onehot, kmer]), names

    def structure(self, df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
        value = pd.to_numeric(df["reactivity_mean_flank10"], errors="raise").to_numpy(dtype=float)
        value = (value - self.structure_mean) / self.structure_scale
        return value[:, None].astype(np.float32), ["R:z_reactivity_mean_flank10"]

    def _predicted_matrix(self, df: pd.DataFrame) -> np.ndarray:
        probabilities = np.asarray([
            json.loads(value) for value in df["predicted_unpaired_probabilities_json"]
        ], dtype=np.float32)
        if probabilities.shape != (len(df), 201):
            raise ValueError("Predicted unpaired-probability matrix must have 201 positions")
        scalars = df[self.predicted_scalars].apply(pd.to_numeric, errors="raise").to_numpy(dtype=np.float32)
        return np.hstack([probabilities, scalars])

    def predicted(self, df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
        value = (self._predicted_matrix(df) - self.predicted_means) / self.predicted_scales
        names = [f"P:z_unpaired_pos{pos:+04d}" for pos in range(-100, 101)]
        names += [f"P:z_{name}" for name in self.predicted_scalars]
        return value.astype(np.float32), names

    def transform(self, df: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str]]:
        common, common_names = self.common(df)
        if model == "M0":
            return common, common_names
        sequence, sequence_names = self.sequence(df)
        arrays, names = [common, sequence], common_names + sequence_names
        if model in {"M2", "M4"}:
            predicted, predicted_names = self.predicted(df)
            arrays.append(predicted)
            names += predicted_names
        if model in {"M3", "M4"}:
            structure, structure_names = self.structure(df)
            arrays.append(structure)
            names += structure_names
        return np.hstack(arrays), names


def metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    rho = spearmanr(y, prediction).statistic
    return {
        "mae": float(mean_absolute_error(y, prediction)),
        "rmse": float(math.sqrt(mean_squared_error(y, prediction))),
        "spearman": float(rho),
        "r2": float(r2_score(y, prediction)),
    }


def paired_bootstrap(df: pd.DataFrame, baseline: str, augmented: str, replicates: int, seed: int) -> tuple[np.ndarray, int]:
    grouped = []
    for _, part in df.groupby("joint_group"):
        grouped.append((
            np.abs(part["combined_ratio"] - part[f"prediction_{baseline}"]).to_numpy(),
            np.abs(part["combined_ratio"] - part[f"prediction_{augmented}"]).to_numpy(),
        ))
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(replicates):
        sampled = rng.integers(0, len(grouped), size=len(grouped))
        base = np.concatenate([grouped[index][0] for index in sampled])
        augmented_values = np.concatenate([grouped[index][1] for index in sampled])
        estimates.append(float(base.mean() - augmented_values.mean()))
    return np.asarray(estimates), len(grouped)


def load_development_overlap(path: Path) -> tuple[set[str], set[str]]:
    df = pd.read_csv(path)
    development = df[df["split"].astype(str).eq("development")]
    genes = set(development["analysis_gene"].astype(str))
    seq_hashes = set(development["sequence_sha256"].astype(str))
    return genes, seq_hashes


def main() -> None:
    setup_logging()
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    prefix = config["frozen_model_prefix"]
    model_names = config["frozen_model_names"]

    df = pd.read_csv(INPUT, compression="gzip", low_memory=False)
    df = df[df["model_dataset_included"].astype(str).eq("True")].copy().reset_index(drop=True)
    predicted = pd.read_csv(PREDICTED_INPUT, compression="gzip", low_memory=False)
    df = df.merge(predicted.drop(columns=["sequence_sha256", "sequence_length"]), on="site_id", how="left", validate="one_to_one")
    if df["predicted_unpaired_probabilities_json"].isna().any():
        raise ValueError("09A predicted structure features are incomplete")
    if not df["sequence_window_201_full"].astype(str).eq("True").all():
        raise ValueError("Model population contains incomplete 201-nt sequences")
    df["joint_group"] = make_joint_groups(df)
    df["sequence_sha256"] = df["sequence_window_201"].astype(str).map(sequence_sha256)

    dev_genes, dev_seq_hashes = load_development_overlap(ROOT / config["development_assignment_input"])
    df["overlaps_development"] = df["analysis_gene"].astype(str).isin(dev_genes) | df["sequence_sha256"].isin(dev_seq_hashes)
    strict = df[~df["overlaps_development"]].copy().reset_index(drop=True)
    logging.info("HeLa model population=%d strict_subset=%d (excluded_overlap=%d)", len(df), len(strict), int(df["overlaps_development"].sum()))

    # Load frozen models and one shared frozen encoder state.
    frozen = {}
    encoder_state = None
    for model_name in model_names:
        path = MODELS_DIR / f"{prefix}_{model_name.lower()}_ridge.joblib"
        bundle = joblib.load(path)
        frozen[model_name] = bundle
        if encoder_state is None:
            encoder_state = bundle["encoder"]
        elif bundle["encoder"] != encoder_state:
            raise AssertionError(f"Frozen encoder state differs across models ({model_name})")
    encoder = FrozenEncoder(encoder_state)

    # Audit unseen categorical levels before transfer.
    unseen_rows = []
    for name in encoder.categorical:
        levels = set(encoder.categories[name])
        seen = df[name].fillna("missing").astype(str).isin(levels)
        unseen = (~seen).sum()
        unseen_rows.append({"categorical": name, "unseen_rows": int(unseen), "unseen_fraction": float(unseen / len(df)) if len(df) else ""})

    predictions = df[["site_id", "analysis_gene", "joint_group", "combined_ratio", "sequence_sha256", "overlaps_development"]].copy()
    for model_name in model_names:
        bundle = frozen[model_name]
        estimator = bundle["model"]
        feature_names = bundle["feature_names"]
        clip = bundle["prediction_clip"]
        x, my_names = encoder.transform(df, model_name)
        if my_names != feature_names:
            raise AssertionError(f"Reconstructed feature order differs from frozen model {model_name}")
        prediction = np.clip(estimator.predict(x), clip[0], clip[1])
        predictions[f"prediction_{model_name}"] = prediction
        logging.info("transferred %s: features=%d", model_name, len(feature_names))

    df = predictions
    strict = df[~df["overlaps_development"]].copy().reset_index(drop=True)

    def metric_rows_for(subset: pd.DataFrame) -> list[dict]:
        rows = []
        y = subset["combined_ratio"].to_numpy(dtype=float)
        for model_name in model_names:
            row = metrics(y, subset[f"prediction_{model_name}"].to_numpy(dtype=float))
            rows.append({"model": model_name, "sites": len(subset), "genes": subset["analysis_gene"].nunique(), "groups": subset["joint_group"].nunique(), **row})
        return rows

    def comparison_rows(subset: pd.DataFrame) -> tuple[list[dict], dict]:
        rows = []
        deltas = {}
        comparisons = [("M1", "M3", "primary_experimental_structure_vs_sequence")]
        comparisons += [("M2", "M4", "secondary_complete_matrix_experimental_structure")]
        comparisons += [("M1", "M2", "secondary_predicted_structure_vs_sequence")]
        comparisons += [("M1", "M4", "secondary_joint_structure_vs_sequence")]
        for baseline, augmented, role in comparisons:
            delta_values, bootstrap_groups = paired_bootstrap(
                subset, baseline, augmented,
                int(config["paired_group_bootstrap_replicates"]),
                int(config["bootstrap_seed"]),
            )
            deltas[f"{baseline}_{augmented}"] = delta_values
            base_mae = float(mean_absolute_error(subset["combined_ratio"], subset[f"prediction_{baseline}"]))
            aug_mae = float(mean_absolute_error(subset["combined_ratio"], subset[f"prediction_{augmented}"]))
            difference = base_mae - aug_mae
            rows.append({
                "comparison": f"{augmented}_vs_{baseline}", "role": role,
                "baseline_model": baseline, "augmented_model": augmented,
                "delta_mae_baseline_minus_augmented": difference,
                "delta_mae_percentage_points": difference * 100,
                "bootstrap_ci95_low": float(np.quantile(delta_values, 0.025)),
                "bootstrap_ci95_high": float(np.quantile(delta_values, 0.975)),
                "bootstrap_ci95_low_percentage_points": float(np.quantile(delta_values, 0.025) * 100),
                "bootstrap_ci95_high_percentage_points": float(np.quantile(delta_values, 0.975) * 100),
                "bootstrap_replicates": len(delta_values), "bootstrap_groups": bootstrap_groups,
            })
        return rows, deltas

    all_metrics = metric_rows_for(df)
    all_comparisons, _ = comparison_rows(df)
    strict_metrics = metric_rows_for(strict)
    strict_comparisons, _ = comparison_rows(strict)

    TABLES.mkdir(parents=True, exist_ok=True)
    write_rows(TABLES / "09c_hela_transfer_metrics.csv", [{"subset": "all_qualifying", **row} for row in all_metrics] + [{"subset": "strict", **row} for row in strict_metrics])
    write_rows(TABLES / "09c_hela_transfer_model_comparison.csv", [{"subset": "all_qualifying", **row} for row in all_comparisons] + [{"subset": "strict", **row} for row in strict_comparisons])
    write_rows(TABLES / "09c_hela_unseen_categorical_levels.csv", unseen_rows)
    predictions.to_csv(TABLES / "09c_hela_transfer_predictions.csv", index=False)
    strict_summary = [{
        "subset": "all_qualifying", "sites": len(df), "genes": df["analysis_gene"].nunique(),
        "groups": df["joint_group"].nunique(), "overlaps_development": int(df["overlaps_development"].sum()),
    }, {
        "subset": "strict", "sites": len(strict), "genes": strict["analysis_gene"].nunique(),
        "groups": strict["joint_group"].nunique(), "overlaps_development": 0,
    }]
    write_rows(TABLES / "09c_hela_transfer_strict_subset_summary.csv", strict_summary)

    by_all = {row["model"]: row for row in all_metrics}
    primary_all = next(row for row in all_comparisons if row["comparison"] == "M3_vs_M1")
    primary_strict = next(row for row in strict_comparisons if row["comparison"] == "M3_vs_M1")
    report = f"""# Stage 09C HeLa frozen-model transfer report

## Scope

The five frozen HEK293T Ridge models (M0--M4) and their single frozen preprocessing encoder were applied
directly to the HeLa model population ({len(df):,} sites with complete 201-nt sequence and valid
experimental flank reactivity, {df['analysis_gene'].nunique():,} genes). No HeLa label, feature, or
parameter was used to select or refit anything.

## Population and strict subset

| Subset | Sites | Genes | Joint groups |
|---|---:|---:|---:|
| All qualifying | {len(df):,} | {df['analysis_gene'].nunique():,} | {df['joint_group'].nunique():,} |
| Strict (no HEK293T-development gene/sequence overlap) | {len(strict):,} | {strict['analysis_gene'].nunique():,} | {strict['joint_group'].nunique():,} |

A HeLa site is excluded from the strict subset when its analysis gene or its 201-nt sequence also appears
in the HEK293T development set used to train the frozen models ({int(df['overlaps_development'].sum()):,}
sites excluded).

## Transfer performance (MAE in proportion units; ×100 = percentage points)

| Model | Inputs | All MAE | Strict MAE | All Spearman | Strict Spearman |
|---|---|---:|---:|---:|---:|
| M0 | common background | {by_all['M0']['mae']:.5f} | {next(r for r in strict_metrics if r['model']=='M0')['mae']:.5f} | {by_all['M0']['spearman']:.4f} | {next(r for r in strict_metrics if r['model']=='M0')['spearman']:.4f} |
| M1 | background + sequence | {by_all['M1']['mae']:.5f} | {next(r for r in strict_metrics if r['model']=='M1')['mae']:.5f} | {by_all['M1']['spearman']:.4f} | {next(r for r in strict_metrics if r['model']=='M1')['spearman']:.4f} |
| M2 | background + sequence + predicted structure | {by_all['M2']['mae']:.5f} | {next(r for r in strict_metrics if r['model']=='M2')['mae']:.5f} | {by_all['M2']['spearman']:.4f} | {next(r for r in strict_metrics if r['model']=='M2')['spearman']:.4f} |
| M3 | background + sequence + experimental R_flank10 | {by_all['M3']['mae']:.5f} | {next(r for r in strict_metrics if r['model']=='M3')['mae']:.5f} | {by_all['M3']['spearman']:.4f} | {next(r for r in strict_metrics if r['model']=='M3')['spearman']:.4f} |
| M4 | background + sequence + predicted + experimental | {by_all['M4']['mae']:.5f} | {next(r for r in strict_metrics if r['model']=='M4')['mae']:.5f} | {by_all['M4']['spearman']:.4f} | {next(r for r in strict_metrics if r['model']=='M4')['spearman']:.4f} |

The primary external comparison is ΔMAE = MAE(M1) − MAE(M3): all-qualifying
{primary_all['delta_mae_percentage_points']:.3f} percentage points (95% CI
{primary_all['bootstrap_ci95_low_percentage_points']:.3f} to {primary_all['bootstrap_ci95_high_percentage_points']:.3f});
strict subset {primary_strict['delta_mae_percentage_points']:.3f} percentage points (95% CI
{primary_strict['bootstrap_ci95_low_percentage_points']:.3f} to {primary_strict['bootstrap_ci95_high_percentage_points']:.3f}).
Positive values mean lower error after adding experimental structure.

## Transfer caveats

The frozen encoder imputes HeLa's missing icSHAPE abundance field (literal `*`) with the HEK293T
development median, and one-hot categorical levels discovered only in HEK293T development are applied
as-is; HeLa rows carrying unseen `drach_subtype` motifs contribute all-zero columns for that block. Those
unseen-level counts are recorded in `09c_hela_unseen_categorical_levels.csv`. These results are direct
transfer evidence, not a HeLa-tuned model.
"""
    (RESULTS / "09c_hela_model_transfer_report.md").write_text(report, encoding="utf-8")

    status = [{
        "stage": "09C", "status": "PASS", "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "all_sites": len(df), "strict_sites": len(strict),
        "all_genes": df["analysis_gene"].nunique(), "strict_genes": strict["analysis_gene"].nunique(),
        "models_transferred": len(model_names),
        "development_overlap_excluded": int(df["overlaps_development"].sum()),
        "primary_delta_mae_all_percentage_points": primary_all["delta_mae_percentage_points"],
        "primary_delta_mae_strict_percentage_points": primary_strict["delta_mae_percentage_points"],
    }]
    write_rows(RESULTS / "09c_hela_model_transfer_run_status.csv", status)
    logging.info("09C completed: all=%d strict=%d primary_delta_all=%.4f strict=%.4f", len(df), len(strict), primary_all["delta_mae_percentage_points"], primary_strict["delta_mae_percentage_points"])


if __name__ == "__main__":
    main()
