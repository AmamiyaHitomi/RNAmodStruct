"""Measure the incremental fit of reactivity in the archived primary models."""

from __future__ import annotations

import argparse
import csv
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml

from replay_presubmission_stage4_models import sha256


MODELS = (
    ("HEK293T", "07_run_association_analysis", "07_hek293t_primary_association.csv"),
    ("HeLa", "09b_run_hela_association", "09b_hela_primary_association.csv"),
)


def summarize(cell_line: str, module_name: str, frozen_table: str) -> dict:
    stage = importlib.import_module(module_name)
    config = yaml.safe_load(stage.CONFIG_PATH.read_text(encoding="utf-8"))["primary_model"]
    predictor, outcome = config["predictor"], config["outcome"]
    sites = stage.load_csv(stage.INPUT)
    if sites["site_id"].duplicated().any():
        raise ValueError(f"{cell_line}: repeated site IDs in the primary input")

    full, clustered, design, response = stage.fit_model(
        sites, stage.DesignBuilder(sites), outcome, predictor
    )
    reduced = sm.OLS(response, design.drop(columns=predictor), missing="raise").fit()
    r2_delta = full.rsquared - reduced.rsquared
    partial_r2 = (reduced.ssr - full.ssr) / reduced.ssr
    np.testing.assert_allclose(partial_r2, r2_delta / (1 - reduced.rsquared), atol=1e-12, rtol=0)

    reference = stage.TABLES / frozen_table
    archived = pd.read_csv(reference).iloc[0]
    counts = [len(sites), sites["analysis_gene"].nunique()]
    if (predictor, outcome) != (archived["predictor"], archived["outcome"]):
        raise ValueError(f"{cell_line}: predictor or outcome differs from the archived model")
    np.testing.assert_array_equal(counts, archived[["sites", "genes"]].to_numpy(dtype=int))
    predictor_sd = sites[predictor].std(ddof=1)
    low, high = clustered.conf_int().loc[predictor]
    effect = clustered.params[predictor] * predictor_sd
    low, high = low * predictor_sd, high * predictor_sd
    np.testing.assert_allclose(
        [full.rsquared, effect, low, high],
        archived[["r_squared", "beta_per_sd", "ci95_low_per_sd", "ci95_high_per_sd"]].to_numpy(dtype=float),
        atol=1e-12,
        rtol=0,
        err_msg=f"{cell_line}: the reconstructed full model differs from the archived model",
    )

    return {
        "cell_line": cell_line,
        "sites": counts[0],
        "genes": counts[1],
        "outcome": outcome,
        "predictor": predictor,
        "beta_pp_per_sd": effect * 100,
        "ci95_low_pp_per_sd": low * 100,
        "ci95_high_pp_per_sd": high * 100,
        "r2_full": full.rsquared,
        "r2_covariates_only": reduced.rsquared,
        "delta_r2": r2_delta,
        "partial_r2": partial_r2,
        "partial_r2_percent": 100 * partial_r2,
        "rss_full": full.ssr,
        "rss_covariates_only": reduced.ssr,
        "input_sha256": sha256(stage.INPUT),
        "config_sha256": sha256(stage.CONFIG_PATH),
        "archived_primary_table_sha256": sha256(reference),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_csv", type=Path, help="A new CSV path. Existing files are preserved.")
    args = parser.parse_args()
    rows = [summarize(*model) for model in MODELS]
    with args.output_csv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: archived full-model results verified for {len(rows)} cell lines. {args.output_csv}")


if __name__ == "__main__":
    main()
