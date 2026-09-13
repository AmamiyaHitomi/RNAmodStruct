"""Stage 12: pool the two frozen per-SD association estimates without rewriting them."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import chi2, norm

from phase1_common import sha256_file, write_rows

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "12_meta_analysis.yaml"
TABLE = ROOT / "results" / "tables" / "12_meta_analysis.csv"
STATUS = ROOT / "results" / "status" / "12_meta_analysis_run_status.csv"
REPORT = ROOT / "results" / "reports" / "12_meta_analysis_report.md"


def pool(effects: np.ndarray, ses: np.ndarray) -> tuple[list[dict], dict]:
    weights = 1.0 / ses**2
    fixed = float(np.sum(weights * effects) / weights.sum())
    q = float(np.sum(weights * (effects - fixed) ** 2))
    df = len(effects) - 1
    c = float(weights.sum() - np.sum(weights**2) / weights.sum())
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    random_weights = 1.0 / (ses**2 + tau2)
    random = float(np.sum(random_weights * effects) / random_weights.sum())
    rows = []
    for model, estimate, se in [
        ("fixed_effect", fixed, float(np.sqrt(1 / weights.sum()))),
        ("random_effect_DL", random, float(np.sqrt(1 / random_weights.sum()))),
    ]:
        rows.append({"row_type": "pooled", "dataset": "pooled", "model": model, "effect_per_sd": estimate,
                     "standard_error": se, "ci95_low": estimate - 1.96 * se, "ci95_high": estimate + 1.96 * se,
                     "z": estimate / se, "p_value": float(2 * norm.sf(abs(estimate / se)))})
    heterogeneity = {"Q": q, "Q_df": df, "Q_p_value": float(chi2.sf(q, df)),
                     "I2_percent": max(0.0, (q - df) / q * 100) if q > 0 else 0.0, "tau2": tau2}
    return rows, heterogeneity


def main() -> None:
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    rows, effects, ses = [], [], []
    for dataset, relative in config["inputs"].items():
        source = pd.read_csv(ROOT / relative).iloc[0]
        effect = float(source["beta_per_sd"])
        se = float(source["se_cluster"] * source["predictor_sd"])
        effects.append(effect); ses.append(se)
        rows.append({"row_type": "dataset", "dataset": dataset, "model": "dataset_specific",
                     "effect_per_sd": effect, "standard_error": se, "ci95_low": effect - 1.96 * se,
                     "ci95_high": effect + 1.96 * se, "z": effect / se,
                     "p_value": float(2 * norm.sf(abs(effect / se)))})
    pooled, heterogeneity = pool(np.asarray(effects), np.asarray(ses))
    for row in pooled:
        row.update(heterogeneity)
    rows += pooled
    write_rows(TABLE, rows)
    write_rows(STATUS, [{"stage": 12, "status": "PASS", "started_utc": started.isoformat(),
                          "finished_utc": datetime.now(timezone.utc).isoformat(), "datasets": 2,
                          "config_sha256": sha256_file(CONFIG), "output_sha256": sha256_file(TABLE)}])
    fixed, random = pooled
    REPORT.write_text(f"""# Stage 12: Summary of dual cell line effects

## Material Passport

- Origin: phase-1 exploratory extension of the frozen HEK293T–HeLa baseline
- Verification Status: ANALYZED
- Generated UTC: {datetime.now(timezone.utc).isoformat()}

The effects per SD are positive for both sets of data. Fixed effects are summarized as {fixed['effect_per_sd']:.5f} (95% CI {fixed['ci95_low']:.5f}–{fixed['ci95_high']:.5f}) and DerSimonian–Laird random effects are summarized as {random['effect_per_sd']:.5f} (95% CI {random['ci95_low']:.5f}–{random['ci95_high']:.5f}). Heterogeneity Q={heterogeneity['Q']:.3f}, I²={heterogeneity['I2_percent']:.1f}%, τ²={heterogeneity['tau2']:.6g}.

With only two cell lines, the heterogeneity parameters are very unstable; per-dataset effects remain in the machine table and the pooled estimate is not used as the final true value. This stage is an exploratory expansion after the main results are known.""", encoding="utf-8")


if __name__ == "__main__":
    main()
