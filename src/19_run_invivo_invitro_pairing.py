"""Stage 19: paired HEK293T in-vivo versus in-vitro icSHAPE analysis."""

from __future__ import annotations

import csv
import gzip
import math
from array import array
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml
from scipy.stats import pearsonr, spearmanr

from pipeline_common import MISSING


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "19_invivo_invitro_pairing.yaml"
PAIRED = ROOT / "data" / "final" / "19_hek293t_invivo_invitro_paired.csv.gz"
ATTRITION = ROOT / "results" / "tables" / "19_invivo_invitro_attrition.csv"
SUMMARY = ROOT / "results" / "tables" / "19_structure_condition_summary.csv"
ASSOCIATIONS = ROOT / "results" / "tables" / "19_paired_condition_associations.csv"
STATUS = ROOT / "results" / "status" / "19_invivo_invitro_run_status.csv"
REPORT = ROOT / "results" / "reports" / "19_invivo_invitro_report.md"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)


def load_selected_icshape(path: Path, wanted: set[str]) -> dict[str, tuple[int, float, array]]:
    result: dict[str, tuple[int, float, array]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] not in wanted:
                continue
            values = array("f", (math.nan if value in MISSING else float(value) for value in fields[3:]))
            if len(values) != int(fields[1]):
                raise ValueError(f"Length mismatch for {fields[0]} in {path.name}")
            result[fields[0]] = (int(fields[1]), float(fields[2]), values)
    return result


def flank_metrics(values: array, center1: int, flank: int) -> dict:
    center0 = center1 - 1
    left = values[center0 - flank:center0] if center0 >= flank else array("f")
    right = values[center0 + 1:center0 + flank + 1]
    left_valid = [float(x) for x in left if not math.isnan(x)]
    right_valid = [float(x) for x in right if not math.isnan(x)]
    left_cov = len(left_valid) / flank if len(left) == flank else 0.0
    right_cov = len(right_valid) / flank if len(right) == flank else 0.0
    combined = left_valid + right_valid
    return {
        "coverage_up10": left_cov,
        "coverage_down10": right_cov,
        "coverage_flank10": len(combined) / (2 * flank),
        "reactivity_mean_flank10": float(np.mean(combined)) if combined else math.nan,
    }


class PairedDesign:
    continuous_names = [
        "gc_fraction_21", "transcript_position_fraction", "log1p_mean_icshape_abundance_rpkm",
        "log1p_combined_agcov", "mean_structure_coverage_flank10",
        "distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx",
    ]

    def __init__(self, df: pd.DataFrame):
        raw = self._continuous(df)
        self.stats = {}
        for name, series in raw.items():
            median = float(series.median()) if series.notna().any() else 0.0
            filled = series.fillna(median)
            sd = float(filled.std(ddof=1))
            self.stats[name] = (median, float(filled.mean()), sd if sd > 0 else 1.0)
        categories = pd.get_dummies(df[["drach_subtype", "transcript_region"]].astype(str), prefix=["drach", "region"], drop_first=True, dtype=float)
        self.categories = list(categories.columns)

    @staticmethod
    def _continuous(df: pd.DataFrame) -> dict[str, pd.Series]:
        return {
            "gc_fraction_21": pd.to_numeric(df["gc_fraction_21"], errors="coerce"),
            "transcript_position_fraction": pd.to_numeric(df["transcript_position_fraction"], errors="coerce"),
            "log1p_mean_icshape_abundance_rpkm": np.log1p(pd.to_numeric(df["mean_icshape_abundance_rpkm"], errors="coerce")),
            "log1p_combined_agcov": np.log1p(pd.to_numeric(df["combined_agcov"], errors="coerce")),
            "mean_structure_coverage_flank10": pd.to_numeric(df["mean_structure_coverage_flank10"], errors="coerce"),
            "distance_to_stop_codon_tx": pd.to_numeric(df["distance_to_stop_codon_tx"], errors="coerce"),
            "distance_to_nearest_splice_edge_tx": pd.to_numeric(df["distance_to_nearest_splice_edge_tx"], errors="coerce"),
        }

    def build(self, df: pd.DataFrame, predictor: str) -> pd.DataFrame:
        design = pd.DataFrame(index=df.index)
        design["Intercept"] = 1.0
        design[predictor] = pd.to_numeric(df[predictor], errors="raise")
        for name, series in self._continuous(df).items():
            median, mean, sd = self.stats[name]
            design[f"z_{name}"] = (series.fillna(median) - mean) / sd
            if name.startswith("distance_to_"):
                design[f"missing_{name}"] = series.isna().astype(float)
        categories = pd.get_dummies(df[["drach_subtype", "transcript_region"]].astype(str), prefix=["drach", "region"], dtype=float)
        for name in self.categories:
            design[name] = categories[name] if name in categories else 0.0
        constant = [name for name in design if name != "Intercept" and design[name].nunique(dropna=False) <= 1]
        return design.drop(columns=constant).astype(float)


def fit_association(df: pd.DataFrame, builder: PairedDesign, predictor: str, label: str) -> tuple[dict, pd.DataFrame, pd.Series]:
    design = builder.build(df, predictor)
    y = pd.to_numeric(df["combined_ratio"]).astype(float)
    ordinary = sm.OLS(y, design).fit()
    clustered = sm.OLS(y, design).fit(cov_type="cluster", cov_kwds={"groups": df["analysis_gene"], "use_correction": True})
    beta = float(clustered.params[predictor])
    se = float(clustered.bse[predictor])
    low, high = clustered.conf_int().loc[predictor]
    sd = float(pd.to_numeric(df[predictor]).std(ddof=1))
    return ({
        "row_type": "condition_effect", "analysis": label, "outcome": "combined_ratio", "predictor": predictor,
        "sites": len(df), "genes": df["analysis_gene"].nunique(), "beta_raw": beta, "se_cluster": se,
        "ci95_low_raw": float(low), "ci95_high_raw": float(high), "p_value": float(clustered.pvalues[predictor]),
        "predictor_sd": sd, "beta_per_sd": beta * sd, "ci95_low_per_sd": float(low) * sd,
        "ci95_high_per_sd": float(high) * sd, "r_squared": float(ordinary.rsquared),
    }, design, y)


def paired_bootstrap_difference(df: pd.DataFrame, vivo_design: pd.DataFrame, vitro_design: pd.DataFrame, y: pd.Series, replicates: int, seed: int) -> np.ndarray:
    groups, inverse = np.unique(df["analysis_gene"].astype(str), return_inverse=True)
    rng = np.random.default_rng(seed)
    vivo_index = list(vivo_design).index("reactivity_mean_flank10_invivo")
    vitro_index = list(vitro_design).index("reactivity_mean_flank10_invitro")
    vivo_sd = float(df["reactivity_mean_flank10_invivo"].std(ddof=1))
    vitro_sd = float(df["reactivity_mean_flank10_invitro"].std(ddof=1))
    xv, xt, outcome = vivo_design.to_numpy(), vitro_design.to_numpy(), y.to_numpy()
    estimates = []
    for _ in range(replicates):
        sampled = rng.integers(0, len(groups), len(groups))
        weights = np.bincount(sampled, minlength=len(groups))[inverse].astype(float)
        keep = weights > 0
        root = np.sqrt(weights[keep])
        try:
            beta_vivo = np.linalg.lstsq(xv[keep] * root[:, None], outcome[keep] * root, rcond=None)[0][vivo_index] * vivo_sd
            beta_vitro = np.linalg.lstsq(xt[keep] * root[:, None], outcome[keep] * root, rcond=None)[0][vitro_index] * vitro_sd
            estimates.append(beta_vivo - beta_vitro)
        except np.linalg.LinAlgError:
            continue
    return np.asarray(estimates)


def main() -> None:
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    base = pd.read_csv(ROOT / config["base_dataset"], low_memory=False)
    eligible = base[base["detected_both_replicates"].astype(str).eq("True")].copy()
    wanted = set(eligible["transcript_id"].astype(str))
    vivo = load_selected_icshape(ROOT / config["invivo"], wanted)
    vitro = load_selected_icshape(ROOT / config["invitro"], wanted)
    flank = int(config["flank"])
    minimum = float(config["minimum_valid_fraction_each_side"])
    records, counts = [], {"transcript_missing_either": 0, "invalid_invivo": 0, "invalid_invitro": 0}
    for row in eligible.itertuples(index=False):
        tid = str(row.transcript_id)
        if tid not in vivo or tid not in vitro:
            counts["transcript_missing_either"] += 1
            continue
        vm = flank_metrics(vivo[tid][2], int(row.transcript_pos_1based), flank)
        tm = flank_metrics(vitro[tid][2], int(row.transcript_pos_1based), flank)
        valid_vivo = vm["coverage_up10"] >= minimum and vm["coverage_down10"] >= minimum
        valid_vitro = tm["coverage_up10"] >= minimum and tm["coverage_down10"] >= minimum
        if not valid_vivo:
            counts["invalid_invivo"] += 1
        if not valid_vitro:
            counts["invalid_invitro"] += 1
        if not (valid_vivo and valid_vitro):
            continue
        record = row._asdict()
        record.update({f"{key}_invivo": value for key, value in vm.items()})
        record.update({f"{key}_invitro": value for key, value in tm.items()})
        record["reactivity_delta_invivo_minus_invitro"] = vm["reactivity_mean_flank10"] - tm["reactivity_mean_flank10"]
        record["mean_structure_coverage_flank10"] = (vm["coverage_flank10"] + tm["coverage_flank10"]) / 2
        record["mean_icshape_abundance_rpkm"] = (vivo[tid][1] + vitro[tid][1]) / 2
        records.append(record)
    paired = pd.DataFrame(records)
    if paired.empty or paired["site_id"].duplicated().any():
        raise RuntimeError("Paired population is empty or not site-unique")
    PAIRED.parent.mkdir(parents=True, exist_ok=True)
    paired.to_csv(PAIRED, index=False, compression={"method": "gzip", "mtime": 0})
    attrition = [
        {"step": 1, "stage": "HEK293T site-level mapped population", "sites": len(base), "rule": "Frozen stage-06 site table"},
        {"step": 2, "stage": "Detected in both GLORI replicates", "sites": len(eligible), "rule": "Same m6A outcome for both structure conditions"},
        {"step": 3, "stage": "Paired valid in-vivo and in-vitro flanks", "sites": len(paired), "rule": f"Same transcript and >= {minimum:.0%} valid values on each side in both conditions"},
    ]
    write_csv(ATTRITION, attrition)
    x, z = paired["reactivity_mean_flank10_invivo"], paired["reactivity_mean_flank10_invitro"]
    summary = [
        {"metric": "paired_sites", "value": len(paired), "detail": "site-unique"},
        {"metric": "paired_genes", "value": paired["analysis_gene"].nunique(), "detail": "analysis_gene"},
        {"metric": "mean_invivo", "value": x.mean(), "detail": "flank10"},
        {"metric": "mean_invitro", "value": z.mean(), "detail": "flank10"},
        {"metric": "mean_delta_invivo_minus_invitro", "value": (x - z).mean(), "detail": "paired"},
        {"metric": "pearson_conditions", "value": pearsonr(x, z).statistic, "detail": "paired sites"},
        {"metric": "spearman_conditions", "value": spearmanr(x, z).statistic, "detail": "paired sites"},
    ]
    write_csv(SUMMARY, summary)
    builder = PairedDesign(paired)
    vivo_row, vivo_design, y = fit_association(paired, builder, "reactivity_mean_flank10_invivo", "in_vivo")
    vitro_row, vitro_design, _ = fit_association(paired, builder, "reactivity_mean_flank10_invitro", "in_vitro")
    delta_row, _, _ = fit_association(paired, builder, "reactivity_delta_invivo_minus_invitro", "delta_invivo_minus_invitro")
    boot = paired_bootstrap_difference(paired, vivo_design, vitro_design, y, int(config["inference"]["paired_gene_bootstrap_replicates"]), int(config["inference"]["bootstrap_seed"]))
    point = vivo_row["beta_per_sd"] - vitro_row["beta_per_sd"]
    contrast = {"row_type": "paired_effect_contrast", "analysis": "in_vivo_minus_in_vitro", "outcome": "combined_ratio", "predictor": "condition_specific_reactivity_per_sd", "sites": len(paired), "genes": paired["analysis_gene"].nunique(), "beta_per_sd": point, "ci95_low_per_sd": float(np.quantile(boot, .025)), "ci95_high_per_sd": float(np.quantile(boot, .975)), "bootstrap_replicates_completed": len(boot)}
    write_csv(ASSOCIATIONS, [vivo_row, vitro_row, delta_row, contrast])
    status = "PASS" if len(paired) >= 500 and len(boot) >= int(config["inference"]["paired_gene_bootstrap_replicates"]) * .95 else "FAIL"
    write_csv(STATUS, [{"stage": 19, "status": status, "started_utc": started.isoformat(), "finished_utc": datetime.now(timezone.utc).isoformat(), "paired_sites": len(paired), "paired_genes": paired["analysis_gene"].nunique(), "bootstrap_completed": len(boot), **counts}])
    REPORT.write_text(f"""# Stage 19: HEK293T in vivo-in vitro icSHAPE pairwise comparison

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {datetime.now(timezone.utc).date().isoformat()}
- Verification Status: {'ANALYZED' if status == 'PASS' else 'UNVERIFIED'}
- Version Label: phase2_stage19_v1

Pair both structural conditions with the same representative transcript and the same transcript position in the frozen HEK293T GLORI site; ±10 nt on both sides are required for effective coverage of {minimum:.0%} in both conditions. Finally, {len(paired):,} loci and {paired['analysis_gene'].nunique():,} genes were retained.

The associated effect per SD was {vivo_row['beta_per_sd']:.6f} (95% CI {vivo_row['ci95_low_per_sd']:.6f}–{vivo_row['ci95_high_per_sd']:.6f}) in vivo; {vitro_row['beta_per_sd']:.6f} (95% CI {vitro_row['ci95_low_per_sd']:.6f}–{vitro_row['ci95_high_per_sd']:.6f}) in vitro. The paired gene bootstrap difference in vivo minus in vitro was {point:.6f} (95% CI {contrast['ci95_low_per_sd']:.6f}–{contrast['ci95_high_per_sd']:.6f}).

The present results belong to an observational paired analysis initiated post hoc. Conditional differences can be used to distinguish consistency of structural environment but cannot alone demonstrate a causal role for intracellular environment or RNA structure on m6A.""", encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
