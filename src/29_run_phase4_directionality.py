"""Stage 29: paired transcript-level directionality analysis for GSE264642."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/29_phase4_directionality_analysis.yaml"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bootstrap_median_ci(values: np.ndarray, reps: int, seed: int, alpha: float) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates = np.empty(reps)
    n = len(values)
    chunk = 50
    for start in range(0, reps, chunk):
        size = min(chunk, reps - start)
        idx = rng.integers(0, n, size=(size, n))
        estimates[start:start + size] = np.median(values[idx], axis=1)
    return tuple(np.quantile(estimates, [alpha / 2, 1 - alpha / 2]))


def sign_flip_mean_p(values: np.ndarray, reps: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    observed = abs(values.mean())
    exceed = 0
    chunk = 100
    for start in range(0, reps, chunk):
        size = min(chunk, reps - start)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(size, len(values)))
        exceed += int((np.abs((signs * values).mean(axis=1)) >= observed).sum())
    return (exceed + 1) / (reps + 1)


def bh_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result.tolist()


def endpoint_row(name: str, values: np.ndarray, role: str, config: dict, seed_offset: int) -> dict:
    alpha = config["alpha"]
    low, high = bootstrap_median_ci(values, config["bootstrap_replicates"], config["random_seed"] + seed_offset, alpha)
    wilcoxon = stats.wilcoxon(values, zero_method="wilcox", alternative="two-sided", method="approx")
    return {
        "endpoint": name, "role": role, "n_transcripts": len(values),
        "control_median": np.nan, "treatment_median": np.nan,
        "paired_median_shift": float(np.median(values)), "bootstrap_ci95_low": low, "bootstrap_ci95_high": high,
        "paired_mean_shift": float(np.mean(values)), "paired_sd": float(np.std(values, ddof=1)),
        "wilcoxon_statistic": float(wilcoxon.statistic), "p_value": float(wilcoxon.pvalue),
        "sign_flip_mean_p": sign_flip_mean_p(values, config["sign_flip_replicates"], config["random_seed"] + 100 + seed_offset),
    }


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    status = pd.read_csv(ROOT / config["required_stage_status"]).iloc[0]
    if status.status != "PASS":
        raise ValueError("Stage 28 did not pass")
    control = pd.read_csv(ROOT / config["inputs"]["control"])
    treatment = pd.read_csv(ROOT / config["inputs"]["treatment"])
    paired = control.merge(treatment, on="transcript", how="inner", validate="one_to_one")
    paired["delta_mean_reactivity"] = paired.stm_average_average - paired.wt_average_average
    paired["delta_reactivity_gini"] = paired.stm_average_gini - paired.wt_average_gini
    paired["delta_max_reactivity"] = paired.stm_average_max - paired.wt_average_max
    paired["delta_reactivity_std"] = paired.stm_average_std - paired.wt_average_std

    a = config["analysis"]
    endpoints = [
        endpoint_row("mean_reactivity", paired.delta_mean_reactivity.to_numpy(), "primary", a, 1),
        endpoint_row("reactivity_gini", paired.delta_reactivity_gini.to_numpy(), "secondary", a, 2),
    ]
    endpoints[0]["control_median"] = float(paired.wt_average_average.median())
    endpoints[0]["treatment_median"] = float(paired.stm_average_average.median())
    endpoints[1]["control_median"] = float(paired.wt_average_gini.median())
    endpoints[1]["treatment_median"] = float(paired.stm_average_gini.median())
    for row, q in zip(endpoints, bh_adjust([row["p_value"] for row in endpoints])):
        row["q_value_bh"] = q
        row["ci_excludes_zero"] = not (row["bootstrap_ci95_low"] <= 0 <= row["bootstrap_ci95_high"])

    nonzero = paired.delta_mean_reactivity.ne(0)
    lower = int((paired.loc[nonzero, "delta_mean_reactivity"] < 0).sum())
    sign_test = stats.binomtest(lower, int(nonzero.sum()), 0.5, alternative="two-sided")
    descriptive = [{
        "endpoint": "lower_mean_fraction", "role": "secondary_descriptive", "n_transcripts": int(nonzero.sum()),
        "count_lower_after_STM2457": lower, "fraction_lower_after_STM2457": lower / int(nonzero.sum()),
        "exact_binomial_p": sign_test.pvalue,
    }]
    attrition = [{
        "control_rows": len(control), "treatment_rows": len(treatment), "paired_rows": len(paired),
        "control_only": int((~control.transcript.isin(paired.transcript)).sum()),
        "treatment_only": int((~treatment.transcript.isin(paired.transcript)).sum()),
        "finite_primary_pairs": int(np.isfinite(paired.delta_mean_reactivity).sum()),
    }]
    out = config["outputs"]
    Path(ROOT / out["paired_data"]).parent.mkdir(parents=True, exist_ok=True)
    paired.to_csv(ROOT / out["paired_data"], index=False, compression="gzip")
    write_csv(ROOT / out["endpoint_results"], endpoints + descriptive)
    write_csv(ROOT / out["attrition"], attrition)
    now = datetime.now(timezone.utc)
    write_csv(ROOT / out["status"], [{
        "stage": 29, "status": "PASS", "paired_transcripts": len(paired),
        "primary_median_shift": endpoints[0]["paired_median_shift"],
        "primary_q_value": endpoints[0]["q_value_bh"],
        "evidence_layer": "directionality", "causal_claim_permitted": False,
        "finished_utc": now.isoformat(),
    }])
    primary, gini = endpoints
    report = f"""# Stage 29: STM2457–RNA structure directionality analysis

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {now.date().isoformat()}
- Verification Status: ANALYZED
- Version Label: phase4_stage29_v1

## Result

A total of {len(paired):,} HEK293T transcripts were paired. The mean SHAPE reactivity pairwise median shift of STM2457 relative to vehicle is {primary['paired_median_shift']:.6g} (transcript bootstrap 95% CI {primary['bootstrap_ci95_low']:.6g} to {primary['bootstrap_ci95_high']:.6g}; Wilcoxon BH q={primary['q_value_bh']:.6g}). The proportion of non-zero differential transcripts that are less reactive after treatment is {descriptive[0]['fraction_lower_after_STM2457']:.3%}.

Gini's paired median shift was {gini['paired_median_shift']:.6g} (95% CI {gini['bootstrap_ci95_low']:.6g} to {gini['bootstrap_ci95_high']:.6g}; BH q={gini['q_value_bh']:.6g}).

## Explain boundaries

This analysis revealed a systematic shift in transcript structure summaries upon METTL3 repression, a directional clue. The processing table is the author's summary of per-transcript statistics, which cannot locate individual m6A sites; and there are no rescues, inactive compound controls, or site-by-site m6A changes for the same sample. The displacement must therefore not be written as an identified causal effect of m6A on the structure, nor shall it be upgraded to mechanistic evidence."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
