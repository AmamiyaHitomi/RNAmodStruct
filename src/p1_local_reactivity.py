"""Compare proximal and distal icSHAPE associations in the archived P1 cohorts."""

from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.linalg import qr
from statsmodels.stats.multitest import multipletests

from p1_coverage_and_within_gene import GENE, OUTCOME, PREDICTOR, coefficient_row, design_matrix
from p1_reactivity_partial_r2 import MODELS
from replay_presubmission_stage4_models import sha256


DISTAL = f"{PREDICTOR}_exclude_m2_p2"
REGIONS = (
    ("upstream_distal", -10, -3, 6),
    ("upstream_proximal", -2, -1, 2),
    ("downstream_proximal", 1, 2, 2),
    ("downstream_distal", 3, 10, 6),
)


def cluster_fit(response, design, frame):
    return sm.OLS(response, design, missing="raise").fit(
        cov_type="cluster", cov_kwds={"groups": frame[GENE], "use_correction": True}
    )


def paired_exclusion(stage, frame, sd, replicates, seed):
    matched = stage.exclude_center_neighborhood_predictor(frame, .70)
    design = design_matrix(stage, matched)
    distal_design = design.assign(**{PREDICTOR: matched[DISTAL]})
    y = matched[OUTCOME]
    designs = {"full_window": design, "exclude_minus2_to_plus2": distal_design}
    rows, samples = [], []
    for label, matrix in designs.items():
        fit = cluster_fit(y, matrix, matched)
        rows.append(coefficient_row(fit, matched, sd, window=label))
        # The same clusters, order, and RNG seed give paired resamples.
        draw = stage.cluster_bootstrap_beta(matrix, y, matched[GENE], PREDICTOR, replicates, seed)
        if len(draw) != replicates:
            raise RuntimeError("A bootstrap replicate failed; pairing cannot be guaranteed")
        samples.append(100 * sd * draw)
    delta = samples[1] - samples[0]
    lo, hi = np.quantile(delta, [.025, .975])
    loss = rows[1]["beta_pp_per_reference_sd"] - rows[0]["beta_pp_per_reference_sd"]
    contrast = {
        "sites": len(matched), "genes": matched[GENE].nunique(), "reference_sd": sd,
        "delta_excluded_minus_full_pp": loss, "ci95_low_pp": lo, "ci95_high_pp": hi,
        "bootstrap_replicates": replicates, "bootstrap_seed": seed,
    }
    return rows, contrast


def local_windows(frame):
    windows = np.asarray([json.loads(item) for item in frame.reactivity_window_m50_p50_json], dtype=float)
    regions = {}
    for label, start, stop, minimum in REGIONS:
        values = windows[:, 50 + start : 51 + stop]
        count = np.isfinite(values).sum(axis=1)
        regions[label] = np.divide(np.nansum(values, axis=1), count,
                                   out=np.full(len(count), np.nan), where=count >= minimum)
    return windows, regions


def local_association(frame, fixed_design, raw_reactivity, sd, **labels):
    covered = np.isfinite(raw_reactivity)
    subset = frame.loc[covered]
    mean_reactivity = raw_reactivity[covered]
    covariates = fixed_design.loc[covered]
    # In the complete-flank cohort, the site coverage covariate is constant.
    # Keep a full-rank basis for the same nuisance column space in each cohort.
    _, r, order = qr(covariates.to_numpy(), mode="economic", pivoting=True)
    covariates = covariates.iloc[:, order[:np.linalg.matrix_rank(r)]]
    forward = cluster_fit(subset[OUTCOME], covariates.assign(**{PREDICTOR: mean_reactivity}), subset)
    reverse = cluster_fit(mean_reactivity, covariates.assign(**{OUTCOME: subset[OUTCOME]}), subset)
    rlo, rhi = reverse.conf_int().loc[OUTCOME]
    return {
        **coefficient_row(forward, subset, sd, **labels),
        "pct_of_primary_sites": 100 * len(subset) / len(frame),
        "p_value": forward.pvalues[PREDICTOR],
        "reactivity_change_per_10pp_m6a": reverse.params[OUTCOME] / 10,
        "reactivity_change_ci95_low": rlo / 10,
        "reactivity_change_ci95_high": rhi / 10,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Fresh results directory")
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    results = {"paired_models": [], "paired_contrasts": [], "regions": [], "positions": [],
               "positions_complete_flanks": []}
    hashes = {}
    for i, (cell, module, _) in enumerate(MODELS):
        stage = importlib.import_module(module)
        frame = stage.load_csv(stage.INPUT)
        if frame.site_id.duplicated().any():
            raise ValueError(f"{cell}: duplicate primary sites")
        sd = frame[PREDICTOR].std(ddof=1)
        print(f"{cell}: {args.bootstrap} paired gene bootstraps and local regressions", flush=True)
        paired_rows, contrast = paired_exclusion(stage, frame, sd, args.bootstrap, 2026112104 + i)
        results["paired_models"].extend(dict(cell_line=cell, **row) for row in paired_rows)
        results["paired_contrasts"].append(dict(cell_line=cell, **contrast))
        windows, regions = local_windows(frame)
        reference = design_matrix(stage, frame).drop(columns=PREDICTOR)
        common = np.isfinite(np.column_stack(list(regions.values()))).all(axis=1)
        for name, start, stop, minimum in REGIONS:
            for population, values in (
                ("region_available", regions[name]),
                ("common_to_four_regions", np.where(common, regions[name], np.nan)),
            ):
                results["regions"].append(
                    local_association(frame, reference, values, sd, cell_line=cell, region=name,
                                      first_nt=start, last_nt=stop, min_valid=minimum, population=population)
                )
        complete = np.isfinite(windows[:, 40:50]).all(axis=1) & np.isfinite(windows[:, 51:61]).all(axis=1)
        for position in (*range(-10, 0), *range(1, 11)):
            reactivity = windows[:, 50 + position]
            results["positions"].append(
                local_association(frame, reference, reactivity, sd, cell_line=cell, position_nt=position)
            )
            results["positions_complete_flanks"].append(
                local_association(frame, reference, np.where(complete, reactivity, np.nan), sd,
                                  cell_line=cell, position_nt=position)
            )
        hashes[cell] = {path.name: sha256(path) for path in (stage.INPUT, stage.CONFIG_PATH)}
    for name, rows in results.items():
        table = pd.DataFrame(rows)
        if name.startswith("positions"):
            table["p_holm"] = multipletests(table.p_value, method="holm")[1]
        elif name == "regions":
            for _, subset in table.groupby("population"):
                table.loc[subset.index, "p_holm"] = multipletests(subset.p_value, method="holm")[1]
        table.to_csv(args.output / f"{name}.csv", mode="x", index=False)
    with (args.output / "input_sha256.json").open("x", encoding="utf-8") as handle:
        json.dump(hashes, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
