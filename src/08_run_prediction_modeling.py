"""Fit the frozen stage-08 HEK293T Ridge ablation and evaluate the held-out groups once."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import logging
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import joblib
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "results" / "logs" / "mplconfig"))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "08_prediction_modeling.yaml"
INPUT = ROOT / "data" / "final" / "06_hek293t_main_analysis_dataset.csv.gz"
PREDICTED_INPUT = ROOT / "data" / "final" / "08a_hek293t_predicted_structure_features.csv.gz"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
MODELS = RESULTS / "models"
FIGURES = RESULTS / "figures"
LOG = RESULTS / "logs" / "08_run_prediction_modeling.log"


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows(rows)


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


def assign_split(groups: pd.Series, test_fraction: float, seed: int) -> pd.Series:
    counts = groups.value_counts().sort_index()
    rng = np.random.default_rng(seed)
    order = counts.index.to_numpy()[rng.permutation(len(counts))]
    cumulative = counts.loc[order].cumsum().to_numpy()
    target = len(groups) * test_fraction
    cut = int(np.argmin(np.abs(cumulative - target))) + 1
    test_groups = set(order[:cut])
    return groups.map(lambda value: "test" if value in test_groups else "development")


class FeatureEncoder:
    numeric = [
        "gc_fraction_21", "transcript_position_fraction", "log1p_icshape_abundance_rpkm",
        "log1p_combined_agcov", "coverage_flank10", "distance_to_stop_codon_tx",
        "distance_to_nearest_splice_edge_tx",
    ]
    missing = ["distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx"]
    categorical = ["drach_subtype", "transcript_region"]
    alphabet = "ACGT"
    kmers = ["".join(parts) for parts in itertools.product("ACGT", repeat=3)]
    predicted_scalars = [
        "mfe_per_nt", "ensemble_free_energy_per_nt", "ensemble_diversity",
        "mean_pairing_state_entropy", "predicted_unpaired_mean_up10",
        "predicted_unpaired_mean_down10", "predicted_unpaired_mean_flank10",
        "predicted_unpaired_mean_far",
    ]

    def __init__(self):
        self.medians: dict[str, float] = {}
        self.means: dict[str, float] = {}
        self.scales: dict[str, float] = {}
        self.categories: dict[str, list[str]] = {}
        self.kmer_means = np.zeros(64, dtype=float)
        self.kmer_scales = np.ones(64, dtype=float)
        self.structure_mean = 0.0
        self.structure_scale = 1.0
        self.predicted_means = np.zeros(201 + len(self.predicted_scalars), dtype=float)
        self.predicted_scales = np.ones(201 + len(self.predicted_scalars), dtype=float)

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

    def fit(self, df: pd.DataFrame) -> "FeatureEncoder":
        numeric = self._numeric_frame(df)
        for name in self.numeric:
            median = float(numeric[name].median())
            filled = numeric[name].fillna(median)
            self.medians[name] = median
            self.means[name] = float(filled.mean())
            scale = float(filled.std(ddof=0))
            self.scales[name] = scale if scale > 0 else 1.0
        self.categories = {
            name: sorted(df[name].fillna("missing").astype(str).unique().tolist())
            for name in self.categorical
        }
        raw_kmers = self._kmer_matrix(df)
        self.kmer_means = raw_kmers.mean(axis=0, dtype=np.float64)
        scales = raw_kmers.std(axis=0, dtype=np.float64)
        self.kmer_scales = np.where(scales > 0, scales, 1.0)
        structure = pd.to_numeric(df["reactivity_mean_flank10"], errors="raise").to_numpy(dtype=float)
        self.structure_mean = float(structure.mean())
        scale = float(structure.std(ddof=0))
        self.structure_scale = scale if scale > 0 else 1.0
        if "predicted_unpaired_probabilities_json" in df:
            predicted = self._predicted_matrix(df)
            self.predicted_means = predicted.mean(axis=0, dtype=np.float64)
            scales = predicted.std(axis=0, dtype=np.float64)
            self.predicted_scales = np.where(scales > 0, scales, 1.0)
        return self

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
        klookup = {value: index for index, value in enumerate(self.kmers)}
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
        lookup = {base: index for index, base in enumerate(self.alphabet)}
        for row, seq in enumerate(sequences):
            for pos, base in enumerate(seq):
                if base in lookup:
                    onehot[row, pos * 4 + lookup[base]] = 1.0
        kmer = (self._kmer_matrix(df) - self.kmer_means) / self.kmer_scales
        names = [f"X:pos{pos:+04d}_{base}" for pos in range(-100, 101) for base in self.alphabet]
        names += [f"X:3mer_{value}" for value in self.kmers]
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

    def to_dict(self) -> dict:
        return {
            "numeric_medians": self.medians, "numeric_means": self.means,
            "numeric_scales": self.scales, "categorical_levels": self.categories,
            "kmer_means": self.kmer_means.tolist(), "kmer_scales": self.kmer_scales.tolist(),
            "structure_mean": self.structure_mean, "structure_scale": self.structure_scale,
            "predicted_means": self.predicted_means.tolist(),
            "predicted_scales": self.predicted_scales.tolist(),
            "sequence_encoding": "201nt_position_one_hot_ACGT_plus_normalized_3mer_counts",
        }


def metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    rho = spearmanr(y, prediction).statistic
    return {
        "mae": float(mean_absolute_error(y, prediction)),
        "rmse": float(math.sqrt(mean_squared_error(y, prediction))),
        "spearman": float(rho),
        "r2": float(r2_score(y, prediction)),
    }


def cv_select(df: pd.DataFrame, groups: np.ndarray, model_name: str, alphas: list[float]) -> tuple[float, list[dict]]:
    y = df["combined_ratio"].to_numpy(dtype=float)
    rows = []
    splitter = GroupKFold(n_splits=5)
    for alpha in alphas:
        fold_values = []
        for fold, (train_index, valid_index) in enumerate(splitter.split(df, y, groups), 1):
            train, valid = df.iloc[train_index], df.iloc[valid_index]
            encoder = FeatureEncoder().fit(train)
            x_train, _ = encoder.transform(train, model_name)
            x_valid, _ = encoder.transform(valid, model_name)
            estimator = Ridge(alpha=alpha).fit(x_train, y[train_index])
            prediction = np.clip(estimator.predict(x_valid), 0.0, 1.0)
            score = metrics(y[valid_index], prediction)
            fold_values.append(score)
            rows.append({"model": model_name, "alpha": alpha, "fold": fold, **score})
        rows.append({
            "model": model_name, "alpha": alpha, "fold": "mean",
            **{name: float(np.mean([value[name] for value in fold_values])) for name in fold_values[0]},
        })
    means = [row for row in rows if row["fold"] == "mean"]
    selected = min(means, key=lambda row: (row["mae"], row["alpha"]))
    if selected["alpha"] in {min(alphas), max(alphas)}:
        raise ValueError(
            f"{model_name} selected boundary alpha={selected['alpha']:g}; "
            "expand alpha_candidates before freezing the model"
        )
    return float(selected["alpha"]), rows


def paired_bootstrap(
    df: pd.DataFrame, baseline: str, augmented: str, replicates: int, seed: int
) -> tuple[np.ndarray, int]:
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


def save_figure(metric_rows: list[dict], delta: np.ndarray) -> None:
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.2))
    models = ["M0", "M1", "M2", "M3", "M4"]
    maes = [next(row["mae"] for row in metric_rows if row["model"] == model) * 100 for model in models]
    positions = np.arange(len(models))
    axes[0].bar(positions, maes, color=["#999999", "#4C78A8", "#59A14F", "#E45756", "#B279A2"])
    axes[0].set_xticks(positions, models)
    for position, value in zip(positions, maes):
        axes[0].text(position, value + 0.05, f"{value:.2f}", ha="center", va="bottom", fontsize=7)
    axes[0].set_ylabel("Held-out MAE (percentage points)")
    axes[0].set_title("a  Frozen test performance", loc="left", fontweight="bold")
    axes[1].hist(delta * 100, bins=35, color="#72B7B2", edgecolor="white")
    axes[1].axvline(0, color="black", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel("Paired ΔMAE: M2 − M4 (percentage points)")
    axes[1].set_ylabel("Bootstrap replicates")
    axes[1].set_title("b  Group-bootstrap uncertainty", loc="left", fontweight="bold")
    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(FIGURES / f"08_hek293t_prediction_ablation.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    setup_logging()
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    df = pd.read_csv(INPUT, compression="gzip", low_memory=False)
    df = df[df["model_dataset_included"].astype(str).eq("True")].copy().reset_index(drop=True)
    predicted = pd.read_csv(PREDICTED_INPUT, compression="gzip", low_memory=False)
    df = df.merge(predicted.drop(columns=["sequence_sha256", "sequence_length"]), on="site_id", how="left", validate="one_to_one")
    if df["predicted_unpaired_probabilities_json"].isna().any():
        raise ValueError("08A predicted structure features are incomplete")
    if not df["sequence_window_201_full"].astype(str).eq("True").all():
        raise ValueError("Model population contains incomplete 201-nt sequences")
    df["joint_group"] = make_joint_groups(df)
    df["split"] = assign_split(df["joint_group"], config["split"]["test_fraction"], config["split"]["random_seed"])
    development = df[df["split"].eq("development")].copy().reset_index(drop=True)
    test = df[df["split"].eq("test")].copy().reset_index(drop=True)
    logging.info("Population=%d development=%d test=%d", len(df), len(development), len(test))

    alphas = [float(value) for value in config["learner"]["alpha_candidates"]]
    selected, cv_rows = {}, []
    model_names = ["M0", "M1", "M2", "M3", "M4"]
    for model_name in model_names:
        selected[model_name], rows = cv_select(
            development, development["joint_group"].to_numpy(), model_name, alphas
        )
        cv_rows += rows
        logging.info("%s selected alpha=%s", model_name, selected[model_name])
    write_rows(TABLES / "08_hek293t_development_cv.csv", cv_rows)

    encoder = FeatureEncoder().fit(development)
    metric_rows, manifest_rows = [], []
    MODELS.mkdir(parents=True, exist_ok=True)
    predictions = test[["site_id", "analysis_gene", "joint_group", "combined_ratio"]].copy()
    for model_name in model_names:
        x_development, feature_names = encoder.transform(development, model_name)
        x_test, test_feature_names = encoder.transform(test, model_name)
        if feature_names != test_feature_names:
            raise AssertionError("Feature order changed between development and test")
        estimator = Ridge(alpha=selected[model_name]).fit(
            x_development, development["combined_ratio"].to_numpy(dtype=float)
        )
        prediction = np.clip(estimator.predict(x_test), 0.0, 1.0)
        predictions[f"prediction_{model_name}"] = prediction
        row = metrics(test["combined_ratio"].to_numpy(dtype=float), prediction)
        metric_rows.append({
            "model": model_name, "status": "FIT", "alpha": selected[model_name],
            "features": len(feature_names), "test_sites": len(test),
            "test_genes": test["analysis_gene"].nunique(), "test_groups": test["joint_group"].nunique(), **row,
        })
        joblib.dump({
            "model": estimator, "encoder": encoder.to_dict(), "feature_names": feature_names,
            "prediction_clip": [0.0, 1.0], "training_input_sha256": sha256(INPUT),
            "config_sha256": sha256(CONFIG_PATH), "development_site_ids_sha256": hashlib.sha256(
                "\n".join(sorted(development["site_id"].astype(str))).encode("utf-8")
            ).hexdigest(),
        }, MODELS / f"08_hek293t_{model_name.lower()}_ridge.joblib")
        manifest_rows += [{"model": model_name, "feature_index": index, "feature_name": name}
                          for index, name in enumerate(feature_names)]
    write_rows(TABLES / "08_hek293t_test_metrics.csv", metric_rows)
    write_rows(TABLES / "08_hek293t_feature_manifest.csv", manifest_rows)
    predictions.to_csv(TABLES / "08_hek293t_test_predictions.csv", index=False)
    assignments = df[["site_id", "analysis_gene", "joint_group", "split", "sequence_window_201"]].copy()
    assignments["sequence_sha256"] = assignments["sequence_window_201"].map(
        lambda value: hashlib.sha256(str(value).encode("ascii")).hexdigest()
    )
    assignments.drop(columns="sequence_window_201").to_csv(
        TABLES / "08_hek293t_group_assignments.csv", index=False
    )

    split_rows = []
    for split_name, part in df.groupby("split"):
        split_rows.append({
            "split": split_name, "sites": len(part), "genes": part["analysis_gene"].nunique(),
            "joint_groups": part["joint_group"].nunique(), "unique_sequences": part["sequence_window_201"].nunique(),
            "outcome_mean": part["combined_ratio"].mean(), "outcome_sd": part["combined_ratio"].std(ddof=1),
        })
    leakage = {
        "gene_overlap": len(set(development["analysis_gene"]) & set(test["analysis_gene"])),
        "sequence_overlap": len(set(development["sequence_window_201"]) & set(test["sequence_window_201"])),
        "joint_group_overlap": len(set(development["joint_group"]) & set(test["joint_group"])),
    }
    split_rows.append({"split": "overlap_audit", **leakage})
    write_rows(TABLES / "08_hek293t_split_audit.csv", split_rows)

    by_model = {row["model"]: row for row in metric_rows}
    comparison = []
    deltas = {}
    for baseline, augmented, role in [
        ("M2", "M4", "primary_complete_model_matrix"),
        ("M1", "M3", "secondary_experimental_structure_vs_sequence"),
        ("M1", "M2", "secondary_predicted_structure_vs_sequence"),
        ("M1", "M4", "secondary_joint_structure_vs_sequence"),
    ]:
        delta_values, bootstrap_groups = paired_bootstrap(
            predictions, baseline, augmented,
            config["test_evaluation"]["paired_group_bootstrap_replicates"],
            config["test_evaluation"]["bootstrap_seed"],
        )
        deltas[f"{baseline}_{augmented}"] = delta_values
        difference = by_model[baseline]["mae"] - by_model[augmented]["mae"]
        comparison.append({
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
    write_rows(TABLES / "08_hek293t_primary_model_comparison.csv", comparison)
    delta = deltas["M2_M4"]
    save_figure(metric_rows, delta)

    report = f"""# Stage 08 HEK293T prediction-model report

## Scope and freeze

The model population contains {len(df):,} sites with complete 201-nt sequence windows and valid
experimental flank reactivity. The deterministic joint grouping joins every site from the same gene and
also joins genes connected by an identical 201-nt sequence. It produced {df['joint_group'].nunique():,}
independent groups. A single frozen split assigned {len(development):,} sites to development and
{len(test):,} sites to the held-out test set. Gene, identical-sequence, and joint-group overlap were all zero.

All imputation, scaling, category discovery, and Ridge alpha selection were fitted within development
folds. The candidate alphas were {alphas}; selected values were M0={selected['M0']:g}, M1={selected['M1']:g},
M2={selected['M2']:g}, M3={selected['M3']:g}, and M4={selected['M4']:g}. Predictions were clipped to [0, 1]
for every model.

### Protocol deviation

An initial implementation was evaluated on this same test split before the continuous 3-mer and
experimental-reactivity features were correctly standardized within each training fold. That first run
gave M1 MAE 0.2068117, M3 MAE 0.2068004, and ΔMAE 0.00113 percentage points. The preprocessing defect
was then corrected without changing the population, split, candidate alphas, feature definitions, or
model family, and the results below were generated. Because test outcomes had already been viewed, the
    corrected result is transparently classified as **test-reused after a preprocessing correction**, not as
    a fully untouched confirmatory test. A later acceptance audit found that the original alpha grid ended
    at the development-selected boundary (100); the grid was widened using development CV only, with a
    hard failure now preventing boundary optima from being frozen. No alternative split was searched.
    Independent confirmation comes from the external HeLa stage.

## Held-out result

| Model | Inputs | MAE | RMSE | Spearman | R2 |
|---|---|---:|---:|---:|---:|
| M0 | common background | {metric_rows[0]['mae']:.5f} | {metric_rows[0]['rmse']:.5f} | {metric_rows[0]['spearman']:.4f} | {metric_rows[0]['r2']:.4f} |
| M1 | background + sequence | {metric_rows[1]['mae']:.5f} | {metric_rows[1]['rmse']:.5f} | {metric_rows[1]['spearman']:.4f} | {metric_rows[1]['r2']:.4f} |
| M2 | background + sequence + predicted structure | {metric_rows[2]['mae']:.5f} | {metric_rows[2]['rmse']:.5f} | {metric_rows[2]['spearman']:.4f} | {metric_rows[2]['r2']:.4f} |
| M3 | background + sequence + experimental R_flank10 | {metric_rows[3]['mae']:.5f} | {metric_rows[3]['rmse']:.5f} | {metric_rows[3]['spearman']:.4f} | {metric_rows[3]['r2']:.4f} |
| M4 | background + sequence + predicted + experimental structure | {metric_rows[4]['mae']:.5f} | {metric_rows[4]['rmse']:.5f} | {metric_rows[4]['spearman']:.4f} | {metric_rows[4]['r2']:.4f} |

The primary complete-matrix comparison gives ΔMAE = MAE(M2) − MAE(M4) =
{comparison[0]['delta_mae_percentage_points']:.3f} percentage points (paired {bootstrap_groups}-group
bootstrap 95% CI {comparison[0]['bootstrap_ci95_low_percentage_points']:.3f} to
{comparison[0]['bootstrap_ci95_high_percentage_points']:.3f}; {len(delta):,} replicates). Positive values
mean lower error after adding experimental structure.

Predicted structure did not improve this Ridge sequence baseline: M2 MAE exceeded M1 MAE by
{(by_model['M2']['mae'] - by_model['M1']['mae']) * 100:.3f} percentage points. This comparison is
descriptive within the reused test set and does not justify searching for a favorable alternative split.

## Interpretation boundary

The full Ridge matrix M0--M4 is complete. The M4-vs-M2 comparison estimates the incremental predictive
value of experimental R_flank10 beyond the frozen sequence and partition-function-derived structure
baseline. Predictive gain from experimental structure is not evidence of a causal mechanism, and the held-out interval covers
test-group sampling uncertainty for fixed fitted models rather than batch or training uncertainty.
"""
    (RESULTS / "08_hek293t_prediction_model_report.md").write_text(report, encoding="utf-8")

    status = [{
        "stage": 8, "status": "PASS" if not any(leakage.values()) else "FAIL",
        "started_utc": started.isoformat(), "finished_utc": datetime.now(timezone.utc).isoformat(),
        "population_sites": len(df), "development_sites": len(development), "test_sites": len(test),
        "development_genes": development["analysis_gene"].nunique(), "test_genes": test["analysis_gene"].nunique(),
        "models_fit": 5, "models_pending": 0, "bootstrap_successful": len(delta),
        "alpha_grid_min": min(alphas), "alpha_grid_max": max(alphas),
        "selected_alphas": json.dumps(selected, sort_keys=True),
        "alpha_boundary_check": "PASS",
        "test_status": "REUSED_AFTER_PREPROCESSING_CORRECTION",
        "input_sha256": sha256(INPUT), "config_sha256": sha256(CONFIG_PATH),
    }]
    write_rows(RESULTS / "08_hek293t_prediction_modeling_run_status.csv", status)
    logging.info("Stage 08 completed with status=%s", status[0]["status"])


if __name__ == "__main__":
    main()
