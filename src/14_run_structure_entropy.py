"""Stage 14: test pairing-state entropy as an independent, non-mechanistic hypothesis."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

from phase1_common import PREDICTED_SCALARS, bh_adjust, clustered_association, sha256_file, write_rows

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "14_structure_entropy.yaml"
TABLE = ROOT / "results" / "tables" / "14_structure_entropy_associations.csv"
STATUS = ROOT / "results" / "status" / "14_structure_entropy_run_status.csv"
REPORT = ROOT / "results" / "reports" / "14_structure_entropy_report.md"


def main() -> None:
    started = datetime.now(timezone.utc)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    rows = []
    specifications = [("unadjusted", False, False), ("covariate_adjusted", True, False),
                      ("covariate_adjusted_plus_R_flank10", True, True)]
    for dataset, item in config["datasets"].items():
        source = pd.read_csv(ROOT / item["input"], low_memory=False)
        source = source[source["model_dataset_included"].astype(str).eq("True")].copy()
        predicted = pd.read_csv(ROOT / item["predicted"], low_memory=False)
        data = source.merge(predicted[["site_id", *PREDICTED_SCALARS]], on="site_id", validate="one_to_one")
        for model, adjust, include_r in specifications:
            start = len(rows)
            for predictor in PREDICTED_SCALARS:
                result = clustered_association(data, predictor, adjust, include_r)
                rows.append({"dataset": dataset, "role": item["role"], "model": model,
                             "predictor": predictor, "is_primary_entropy_hypothesis": predictor == config["primary_predictor"],
                             **result})
            adjusted = bh_adjust([row["p_value"] for row in rows[start:]])
            for row, value in zip(rows[start:], adjusted):
                row["p_value_bh"] = float(value)
    write_rows(TABLE, rows)
    write_rows(STATUS, [{"stage": 14, "status": "PASS", "started_utc": started.isoformat(),
                          "finished_utc": datetime.now(timezone.utc).isoformat(), "tests": len(rows),
                          "config_sha256": sha256_file(CONFIG), "output_sha256": sha256_file(TABLE)}])
    primary = [row for row in rows if row["is_primary_entropy_hypothesis"]]
    lines = ["# Stage 14: Structural entropy independence assumption", "", "## Material Passport", "",
             "- Origin: phase-1 exploratory extension", "- Verification Status: ANALYZED", "",
             "`mean_pairing_state_entropy` is BH corrected as the same comparative family as the remaining predicted structural features; results are observational, post hoc extensions, and not subject to mechanistic interpretation.", ""]
    for row in primary:
        lines.append(f"- {row['dataset']} / {row['model']}：β/SD={row['beta_per_sd']:.5f}，95% CI {row['ci95_low_per_sd']:.5f}–{row['ci95_high_per_sd']:.5f}，BH q={row['p_value_bh']:.4g}。")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
