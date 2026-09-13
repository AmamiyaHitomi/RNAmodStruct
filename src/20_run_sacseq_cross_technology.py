"""Stage 20: cross-m6A-technology replication using m6A-SAC-seq."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

from phase1_common import TabularEncoder, clustered_association, make_joint_groups, model_metrics


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "20_sacseq_cross_technology.yaml"
JOINED = ROOT / "data" / "final" / "20_sacseq_structure_joined.csv.gz"
ATTRITION = ROOT / "results" / "tables" / "20_sacseq_overlap_attrition.csv"
FIELDS = ROOT / "results" / "tables" / "20_sacseq_field_audit.csv"
ASSOCIATIONS = ROOT / "results" / "tables" / "20_cross_technology_associations.csv"
CV = ROOT / "results" / "tables" / "20_sacseq_nested_cv_metrics.csv"
PREDICTIONS = ROOT / "results" / "tables" / "20_sacseq_oof_predictions.csv.gz"
STATUS = ROOT / "results" / "status" / "20_sacseq_cross_technology_run_status.csv"
REPORT = ROOT / "results" / "reports" / "20_sacseq_cross_technology_report.md"


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)


def load_sacseq(path: Path) -> pd.DataFrame:
    names = ["sac_chrom", "sac_start0", "sac_end1", "sac_source_record", "reported_mean_mutation_ratio_pct", "sac_strand", "reported_calibrated_m6a_fraction_pct"]
    frame = pd.read_csv(path, sep="\t", header=None, names=names)
    if not (frame["sac_end1"] - frame["sac_start0"]).eq(1).all():
        raise ValueError(f"Non-single-base BED interval in {path.name}")
    if not set(frame["sac_strand"]).issubset({"+", "-"}):
        raise ValueError(f"Invalid strand in {path.name}")
    if frame.duplicated(["sac_chrom", "sac_start0", "sac_end1", "sac_strand"]).any():
        raise ValueError(f"Duplicate stranded coordinate in {path.name}")
    for column in ["reported_mean_mutation_ratio_pct", "reported_calibrated_m6a_fraction_pct"]:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
        if not frame[column].between(0, 100).all():
            raise ValueError(f"Out-of-range percentage in {column}")
    frame["sac_fraction"] = frame["reported_calibrated_m6a_fraction_pct"] / 100.0
    return frame


def join_dataset(label: str, item: dict) -> tuple[pd.DataFrame, list[dict], list[dict]]:
    structure = pd.read_csv(ROOT / item["structure_dataset"], low_memory=False)
    sac = load_sacseq(ROOT / item["sacseq"])
    joined = structure.merge(sac, left_on=["hg38_chr", "hg38_pos_1based", "hg38_strand"], right_on=["sac_chrom", "sac_end1", "sac_strand"], validate="one_to_one")
    main = joined[joined["main_analysis_included"].astype(str).eq("True")].copy()
    predicted = pd.read_csv(ROOT / item["predicted_structure"], low_memory=False)
    model = main[main["model_dataset_included"].astype(str).eq("True")].merge(predicted, on="site_id", validate="one_to_one")
    model["dataset"] = label
    model["match_class"] = item["match_class"]
    attrition = [
        {"dataset": label, "step": 1, "population": "SAC-seq submitted sites", "sites": len(sac)},
        {"dataset": label, "step": 2, "population": "Exact stranded overlap with existing structure site table", "sites": len(joined)},
        {"dataset": label, "step": 3, "population": "Existing valid structure main population", "sites": len(main)},
        {"dataset": label, "step": 4, "population": "Model-ready with predicted structure", "sites": len(model)},
    ]
    field_rows = []
    for column, meaning in [
        ("reported_mean_mutation_ratio_pct", "submitted BED column 5; mutation-ratio percentage reported by depositor"),
        ("reported_calibrated_m6a_fraction_pct", "submitted BED column 7; calibrated m6A-fraction percentage reported by depositor"),
    ]:
        field_rows.append({"dataset": label, "field": column, "minimum": sac[column].min(), "median": sac[column].median(), "maximum": sac[column].max(), "missing": int(sac[column].isna().sum()), "meaning": meaning})
    return model, attrition, field_rows


def association_rows(frame: pd.DataFrame, label: str, match_class: str) -> list[dict]:
    rows = []
    correlation = frame[["combined_ratio", "sac_fraction"]].corr(method="spearman").iloc[0, 1]
    for outcome, technology in [("combined_ratio", "GLORI"), ("sac_fraction", "m6A-SAC-seq")]:
        analysis = frame.copy()
        analysis["combined_ratio"] = pd.to_numeric(analysis[outcome])
        result = clustered_association(analysis, "reactivity_mean_flank10", adjust=True)
        rows.append({"dataset": label, "match_class": match_class, "outcome_technology": technology, "outcome": outcome, "same_site_glori_sac_spearman": correlation, **result})
    return rows


def nested_group_cv(frame: pd.DataFrame, config: dict) -> tuple[list[dict], pd.DataFrame]:
    frame = frame.reset_index(drop=True).copy()
    frame["joint_group"] = make_joint_groups(frame, "sequence_window_201")
    groups = frame["joint_group"].to_numpy()
    y = frame["sac_fraction"].to_numpy(float)
    outer = GroupKFold(int(config["learner"]["outer_group_folds"]))
    alpha_candidates = [float(value) for value in config["learner"]["alpha_candidates"]]
    metric_rows, prediction_rows = [], []
    for model_name in config["models"]:
        for outer_fold, (train_idx, test_idx) in enumerate(outer.split(frame, groups=groups), 1):
            train = frame.iloc[train_idx]
            train_groups = groups[train_idx]
            inner = GroupKFold(int(config["learner"]["inner_group_folds"]))
            alpha_scores = []
            for alpha in alpha_candidates:
                fold_mae = []
                for inner_train, inner_valid in inner.split(train, groups=train_groups):
                    encoder = TabularEncoder().fit(train.iloc[inner_train])
                    x_train, _ = encoder.transform(train.iloc[inner_train], model_name)
                    x_valid, _ = encoder.transform(train.iloc[inner_valid], model_name)
                    prediction = np.clip(Ridge(alpha=alpha, solver="lsqr").fit(x_train, train.iloc[inner_train]["sac_fraction"]).predict(x_valid), 0, 1)
                    fold_mae.append(model_metrics(train.iloc[inner_valid]["sac_fraction"].to_numpy(), prediction)["mae"])
                alpha_scores.append((float(np.mean(fold_mae)), alpha))
            selected_alpha = min(alpha_scores)[1]
            encoder = TabularEncoder().fit(train)
            x_train, _ = encoder.transform(train, model_name)
            x_test, _ = encoder.transform(frame.iloc[test_idx], model_name)
            prediction = np.clip(Ridge(alpha=selected_alpha, solver="lsqr").fit(x_train, y[train_idx]).predict(x_test), 0, 1)
            metrics = model_metrics(y[test_idx], prediction)
            metric_rows.append({"model": model_name, "row_type": "outer_fold", "outer_fold": outer_fold, "selected_alpha": selected_alpha, "sites": len(test_idx), **metrics})
            for idx, value in zip(test_idx, prediction):
                prediction_rows.append({"site_id": frame.iloc[idx]["site_id"], "joint_group": groups[idx], "outer_fold": outer_fold, "model": model_name, "observed_sac_fraction": y[idx], "predicted_sac_fraction": float(value), "selected_alpha": selected_alpha})
        model_predictions = pd.DataFrame([row for row in prediction_rows if row["model"] == model_name]).sort_values("site_id")
        overall = model_metrics(model_predictions["observed_sac_fraction"].to_numpy(), model_predictions["predicted_sac_fraction"].to_numpy())
        metric_rows.append({"model": model_name, "row_type": "overall_oof", "outer_fold": "all", "selected_alpha": "nested", "sites": len(model_predictions), **overall})
    return metric_rows, pd.DataFrame(prediction_rows)


def main() -> None:
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    frames, attrition, fields, associations = [], [], [], []
    for label, item in config["datasets"].items():
        frame, dataset_attrition, dataset_fields = join_dataset(label, item)
        frames.append(frame)
        attrition += dataset_attrition
        fields += dataset_fields
        associations += association_rows(frame, label, item["match_class"])
    combined = pd.concat(frames, ignore_index=True)
    JOINED.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(JOINED, index=False, compression={"method": "gzip", "mtime": 0})
    write_csv(ATTRITION, attrition)
    write_csv(FIELDS, fields)
    write_csv(ASSOCIATIONS, associations)
    primary = frames[0]
    cv_rows, predictions = nested_group_cv(primary, config)
    write_csv(CV, cv_rows)
    predictions.to_csv(PREDICTIONS, index=False, compression={"method": "gzip", "mtime": 0})
    overall = pd.DataFrame(cv_rows)
    overall = overall[overall["row_type"].eq("overall_oof")].set_index("model")
    structure_gain = float(overall.loc["M1", "mae"] - overall.loc["M3", "mae"])
    full_gain = float(overall.loc["M2", "mae"] - overall.loc["M4", "mae"])
    primary_associations = [row for row in associations if row["dataset"] == "HeLa"]
    status = "PASS" if len(primary) >= 500 and set(overall.index) == set(config["models"]) else "FAIL"
    write_csv(STATUS, [{"stage": 20, "status": status, "started_utc": started.isoformat(), "finished_utc": datetime.now(timezone.utc).isoformat(), "hela_primary_sites": len(primary), "hek293_near_match_sites": len(frames[1]), "models": len(overall), "m1_minus_m3_mae_gain": structure_gain, "m2_minus_m4_mae_gain": full_gain}])
    glori = next(row for row in primary_associations if row["outcome_technology"] == "GLORI")
    sac = next(row for row in primary_associations if row["outcome_technology"] == "m6A-SAC-seq")
    REPORT.write_text(f"""# Stage 20: SAC-seq cross-m6A measurement technology review

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {datetime.now(timezone.utc).date().isoformat()}
- Verification Status: {'ANALYZED' if status == 'PASS' else 'UNVERIFIED'}
- Version Label: phase2_stage20_v1

The HeLa primary review was performed on {len(primary):,} sites with both valid icSHAPE, complete 201 nt sequence, predicted structure, and SAC-seq quantification. It is a cross-study review of matched cell lines, not the same batch of samples. Only the {len(frames[1]):,} sites of HEK293/HEK293T were analyzed for close matching sensitivity.

In the same HeLa locus, the adjusted effect of experimental structure on GLORI was {glori['beta_per_sd']:.6f}/SD (95% CI {glori['ci95_low_per_sd']:.6f}–{glori['ci95_high_per_sd']:.6f}) and on the SAC-seq calibration score was {sac['beta_per_sd']:.6f}/SD (95% CI {sac['ci95_low_per_sd']:.6f}–{sac['ci95_high_per_sd']:.6f}). The Spearman correlation between the two m6A values ​​at these sites is {sac['same_site_glori_sac_spearman']:.3f}.

In nested gene-identical sequence group cross-validation, the MAE improvement by adding the experimental structure is M1−M3={structure_gain:.6f}, and the improvement based on the existing predicted structure is M2−M4={full_gain:.6f}. All encoding, scaling, and alpha selections are restricted to the outer training fold.

Columns 5 and 7 of SAC-seq BED are saved as `reported_mean_mutation_ratio_pct` and `reported_calibrated_m6a_fraction_pct` respectively, emphasizing that they are the results of the submitter's processing; the main outcome is column 7 divided by 100. Method source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9378555/ .""", encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
