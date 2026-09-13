"""Shared, deterministic helpers for the phase-1 exploratory extensions (stages 12--17)."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


COMMON_NUMERIC = [
    "gc_fraction_21", "transcript_position_fraction", "log1p_icshape_abundance_rpkm",
    "log1p_combined_agcov", "coverage_flank10", "distance_to_stop_codon_tx",
    "distance_to_nearest_splice_edge_tx",
]
COMMON_MISSING = ["distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx"]
COMMON_CATEGORICAL = ["drach_subtype", "transcript_region"]
PREDICTED_SCALARS = [
    "mfe_per_nt", "ensemble_free_energy_per_nt", "ensemble_diversity",
    "mean_pairing_state_entropy", "predicted_unpaired_mean_up10",
    "predicted_unpaired_mean_down10", "predicted_unpaired_mean_flank10",
    "predicted_unpaired_mean_far",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def z_or_one(values: np.ndarray) -> tuple[float, float]:
    if not np.isfinite(values).any():
        return 0.0, 1.0
    mean = float(np.nanmean(values))
    scale = float(np.nanstd(values, ddof=0))
    return mean, scale if scale > 0 else 1.0


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


def make_joint_groups(df: pd.DataFrame, sequence_column: str) -> pd.Series:
    genes = sorted(df["analysis_gene"].fillna("missing_gene").astype(str).unique())
    uf = UnionFind(genes)
    for _, values in df.groupby(sequence_column, sort=False)["analysis_gene"]:
        unique = sorted(values.fillna("missing_gene").astype(str).unique())
        for gene in unique[1:]:
            uf.union(unique[0], gene)
    roots = {gene: uf.find(gene) for gene in genes}
    labels = {root: f"G{index:05d}" for index, root in enumerate(sorted(set(roots.values())), 1)}
    return df["analysis_gene"].fillna("missing_gene").astype(str).map(lambda g: labels[roots[g]])


def add_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    if "combined_agcov" not in result:
        result["combined_agcov"] = pd.to_numeric(result["agcov_rep1"]) + pd.to_numeric(result["agcov_rep2"])
    if "combined_acov" not in result:
        result["combined_acov"] = pd.to_numeric(result["acov_rep1"]) + pd.to_numeric(result["acov_rep2"])
    if "combined_ratio" not in result:
        result["combined_ratio"] = result["combined_acov"] / result["combined_agcov"]
    if "site_id" not in result:
        result["site_id"] = (
            result["hg38_chr"].astype(str) + ":" + result["hg38_pos_1based"].astype(str)
            + ":" + result["hg38_strand"].astype(str)
        )
    if "analysis_gene" not in result:
        gene = result["genes"].mask(result["genes"].isna() | result["genes"].astype(str).eq("ELSE"))
        result["analysis_gene"] = gene.fillna(result["transcript_gene_name"]).fillna(result["transcript_gene_id"])
    return result


class TabularEncoder:
    """Fold-fitted encoder for common, 3-mer, predicted and experimental feature blocks."""

    kmers = ["".join(x) for x in itertools.product("ACGT", repeat=3)]

    def __init__(self, sequence_column: str = "sequence_window_201"):
        self.sequence_column = sequence_column
        self.numeric_stats: dict[str, tuple[float, float, float]] = {}
        self.categories: dict[str, list[str]] = {}
        self.kmer_mean = np.zeros(64)
        self.kmer_scale = np.ones(64)
        self.predicted_mean = np.zeros(len(PREDICTED_SCALARS))
        self.predicted_scale = np.ones(len(PREDICTED_SCALARS))
        self.structure_mean = 0.0
        self.structure_scale = 1.0

    @staticmethod
    def numeric_frame(df: pd.DataFrame) -> pd.DataFrame:
        return pd.DataFrame({
            "gc_fraction_21": pd.to_numeric(df["gc_fraction_21"], errors="coerce"),
            "transcript_position_fraction": pd.to_numeric(df["transcript_position_fraction"], errors="coerce"),
            "log1p_icshape_abundance_rpkm": np.log1p(pd.to_numeric(df["icshape_abundance_rpkm"], errors="coerce")),
            "log1p_combined_agcov": np.log1p(pd.to_numeric(df["combined_agcov"], errors="coerce")),
            "coverage_flank10": pd.to_numeric(df["coverage_flank10"], errors="coerce"),
            "distance_to_stop_codon_tx": pd.to_numeric(df["distance_to_stop_codon_tx"], errors="coerce"),
            "distance_to_nearest_splice_edge_tx": pd.to_numeric(df["distance_to_nearest_splice_edge_tx"], errors="coerce"),
        }, index=df.index)

    def kmer_matrix(self, df: pd.DataFrame) -> np.ndarray:
        sequences = df[self.sequence_column].astype(str).str.upper().tolist()
        lookup = {kmer: index for index, kmer in enumerate(self.kmers)}
        output = np.zeros((len(sequences), len(self.kmers)), dtype=np.float32)
        for row, sequence in enumerate(sequences):
            denominator = max(len(sequence) - 2, 1)
            for pos in range(len(sequence) - 2):
                index = lookup.get(sequence[pos:pos + 3])
                if index is not None:
                    output[row, index] += 1.0 / denominator
        return output

    def fit(self, df: pd.DataFrame) -> "TabularEncoder":
        numeric = self.numeric_frame(df)
        for name in COMMON_NUMERIC:
            median = float(numeric[name].median())
            if not math.isfinite(median):
                median = 0.0
            values = numeric[name].fillna(median).to_numpy(float)
            mean, scale = z_or_one(values)
            self.numeric_stats[name] = (median, mean, scale)
        self.categories = {name: sorted(df[name].fillna("missing").astype(str).unique()) for name in COMMON_CATEGORICAL}
        kmer = self.kmer_matrix(df)
        self.kmer_mean = kmer.mean(axis=0)
        self.kmer_scale = np.where(kmer.std(axis=0) > 0, kmer.std(axis=0), 1.0)
        if set(PREDICTED_SCALARS).issubset(df.columns):
            predicted = df[PREDICTED_SCALARS].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            self.predicted_mean = np.nanmean(predicted, axis=0)
            self.predicted_mean = np.where(np.isfinite(self.predicted_mean), self.predicted_mean, 0.0)
            predicted = np.where(np.isnan(predicted), self.predicted_mean, predicted)
            self.predicted_scale = np.where(predicted.std(axis=0) > 0, predicted.std(axis=0), 1.0)
        structure = pd.to_numeric(df["reactivity_mean_flank10"], errors="coerce").to_numpy(float)
        self.structure_mean, self.structure_scale = z_or_one(structure)
        return self

    def common(self, df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
        numeric = self.numeric_frame(df)
        arrays, names = [], []
        for name in COMMON_NUMERIC:
            median, mean, scale = self.numeric_stats[name]
            arrays.append(((numeric[name].fillna(median).to_numpy(float) - mean) / scale)[:, None])
            names.append(f"C:{name}")
        for name in COMMON_MISSING:
            arrays.append(numeric[name].isna().to_numpy(float)[:, None])
            names.append(f"C:missing_{name}")
        for name in COMMON_CATEGORICAL:
            values = df[name].fillna("missing").astype(str).to_numpy()
            for category in self.categories[name]:
                arrays.append((values == category).astype(float)[:, None])
                names.append(f"C:{name}={category}")
        return np.hstack(arrays).astype(np.float32), names

    def transform(self, df: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str]]:
        arrays, names = list(), list()
        common, common_names = self.common(df)
        arrays.append(common); names += common_names
        if model != "M0":
            kmer = (self.kmer_matrix(df) - self.kmer_mean) / self.kmer_scale
            arrays.append(kmer); names += [f"X:3mer_{x}" for x in self.kmers]
        if model in {"M2", "M4"}:
            raw = df[PREDICTED_SCALARS].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            raw = np.where(np.isnan(raw), self.predicted_mean, raw)
            arrays.append(((raw - self.predicted_mean) / self.predicted_scale).astype(np.float32))
            names += [f"P:{x}" for x in PREDICTED_SCALARS]
        if model in {"M3", "M4"}:
            raw = pd.to_numeric(df["reactivity_mean_flank10"], errors="raise").to_numpy(float)
            arrays.append(((raw - self.structure_mean) / self.structure_scale)[:, None].astype(np.float32))
            names.append("R:reactivity_mean_flank10")
        return np.hstack(arrays), names


def model_metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y, prediction)),
        "rmse": float(math.sqrt(mean_squared_error(y, prediction))),
        "spearman": float(spearmanr(y, prediction).statistic),
        "r2": float(r2_score(y, prediction)),
    }


def association_design(df: pd.DataFrame, predictor: str, adjust: bool, include_r: bool = False) -> pd.DataFrame:
    output = pd.DataFrame({"Intercept": 1.0, predictor: pd.to_numeric(df[predictor])}, index=df.index)
    if include_r and predictor != "reactivity_mean_flank10":
        output["reactivity_mean_flank10"] = pd.to_numeric(df["reactivity_mean_flank10"])
    if not adjust:
        return output
    numeric = TabularEncoder.numeric_frame(df)
    for name in COMMON_NUMERIC:
        values = numeric[name]
        median = float(values.median())
        if not math.isfinite(median):
            median = 0.0
        filled = values.fillna(median)
        mean, scale = z_or_one(filled.to_numpy(float))
        output[f"z_{name}"] = (filled - mean) / scale
        if name in COMMON_MISSING:
            output[f"missing_{name}"] = values.isna().astype(float)
    categorical = pd.get_dummies(df[COMMON_CATEGORICAL].fillna("missing").astype(str), drop_first=True, dtype=float)
    return pd.concat([output, categorical], axis=1).loc[:, lambda x: x.nunique(dropna=False) > 1].assign(Intercept=1.0)


def clustered_association(df: pd.DataFrame, predictor: str, adjust: bool, include_r: bool = False) -> dict:
    design = association_design(df, predictor, adjust, include_r)
    if "Intercept" not in design:
        design.insert(0, "Intercept", 1.0)
    y = pd.to_numeric(df["combined_ratio"]).to_numpy(float)
    ordinary = sm.OLS(y, design).fit()
    fit = sm.OLS(y, design).fit(cov_type="cluster", cov_kwds={"groups": df["analysis_gene"], "use_correction": True})
    beta, se = float(fit.params[predictor]), float(fit.bse[predictor])
    sd = float(pd.to_numeric(df[predictor]).std(ddof=1))
    return {
        "sites": len(df), "genes": df["analysis_gene"].nunique(), "beta_raw": beta,
        "se_cluster": se, "ci95_low_raw": beta - 1.96 * se, "ci95_high_raw": beta + 1.96 * se,
        "p_value": float(fit.pvalues[predictor]), "predictor_sd": sd, "beta_per_sd": beta * sd,
        "se_per_sd": se * sd, "ci95_low_per_sd": (beta - 1.96 * se) * sd,
        "ci95_high_per_sd": (beta + 1.96 * se) * sd, "r_squared": float(ordinary.rsquared),
    }


def bh_adjust(p_values: list[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output = np.empty_like(ranked)
    output[order] = np.clip(ranked, 0, 1)
    return output


def normal_pvalue(z: float) -> float:
    return float(2 * norm.sf(abs(z)))
