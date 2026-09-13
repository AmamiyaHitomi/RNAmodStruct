"""Stage 15: compare four transcript-mapping estimands while keeping each genomic site intact."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from phase1_common import add_derived_columns, clustered_association, sha256_file, write_rows

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "15_isoform_sensitivity.yaml"
TABLE = ROOT / "results" / "tables" / "15_isoform_sensitivity_associations.csv"
ESTIMATES = ROOT / "results" / "tables" / "15_isoform_site_estimates.csv.gz"
STATUS = ROOT / "results" / "status" / "15_isoform_sensitivity_run_status.csv"
REPORT = ROOT / "results" / "reports" / "15_isoform_sensitivity_report.md"


def build_estimands(raw: pd.DataFrame) -> list[pd.DataFrame]:
    data = add_derived_columns(raw)
    data = data[data["detected_both_replicates"].astype(str).eq("True")].copy()
    valid = data[data["main_window_valid"].astype(str).eq("True") & data["reactivity_mean_flank10"].notna()].copy()
    representative = valid[valid["representative_isoform"].astype(str).eq("True")].copy()
    representative["estimand"] = "deterministic_representative_isoform"
    counts = data.groupby("site_id").size()
    single_ids = set(counts[counts == 1].index)
    single = representative[representative["site_id"].isin(single_ids)].copy()
    single["estimand"] = "single_isoform_sites"
    base = data.sort_values(["site_id", "isoform_selection_rank"]).drop_duplicates("site_id").set_index("site_id")
    equal_r = valid.groupby("site_id")["reactivity_mean_flank10"].mean()
    equal = base.loc[equal_r.index].copy()
    equal["reactivity_mean_flank10"] = equal_r
    equal = equal.reset_index(); equal["estimand"] = "isoform_equal_weight"
    weighted_r = valid.groupby("site_id").apply(
        lambda x: float(np.average(x["reactivity_mean_flank10"], weights=np.maximum(x["coverage_flank10"], 1e-12))),
        include_groups=False,
    )
    weighted = base.loc[weighted_r.index].copy()
    weighted["reactivity_mean_flank10"] = weighted_r
    weighted = weighted.reset_index(); weighted["estimand"] = "available_structure_coverage_weight"
    return [representative, single, equal, weighted]


def main() -> None:
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    rows, estimates = [], []
    for dataset, item in config["datasets"].items():
        raw = pd.read_csv(ROOT / item["mappings"], low_memory=False)
        for frame in build_estimands(raw):
            if len(frame) < 2:
                raise ValueError(f"{dataset} {frame['estimand'].iloc[0]} has too few sites")
            frame = frame.copy()
            frame["dataset"] = dataset
            estimates.append(frame[["dataset", "estimand", "site_id", "analysis_gene", "transcript_id",
                                    "reactivity_mean_flank10", "combined_ratio", "coverage_flank10"]])
            result = clustered_association(frame, "reactivity_mean_flank10", adjust=True)
            rows.append({"dataset": dataset, "estimand": frame["estimand"].iloc[0],
                         "isoform_abundance_interpretation_permitted": False, **result})
    write_rows(TABLE, rows)
    pd.concat(estimates).to_csv(ESTIMATES, index=False, compression={"method": "gzip", "mtime": 0})
    write_rows(STATUS, [{"stage": 15, "status": "PASS", "started_utc": started.isoformat(),
                          "finished_utc": datetime.now(timezone.utc).isoformat(), "estimand_rows": len(rows),
                          "config_sha256": sha256_file(CONFIG), "table_sha256": sha256_file(TABLE),
                          "estimates_sha256": sha256_file(ESTIMATES)}])
    lines = ["# Stage 15: Isoform Sensitivity", "", "## Material Passport", "",
             "- Origin: phase-1 exploratory extension", "- Verification Status: ANALYZED", "",
             "All four mapping calibers use genomic sites as analysis/resampling boundaries. Equal weighting and structure coverage weighting are sensitivity weights only; short-read data are not interpreted as true isoform abundances.", ""]
    for row in rows:
        lines.append(f"- {row['dataset']} / {row['estimand']}：n={row['sites']}，β/SD={row['beta_per_sd']:.5f}（95% CI {row['ci95_low_per_sd']:.5f}–{row['ci95_high_per_sd']:.5f}）。")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
