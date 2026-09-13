"""Stage 21: finalize conditional admission decisions for deferred phase-2 sources."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "21_conditional_source_admission.yaml"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    audit = pd.read_csv(ROOT / config["input_audit"], keep_default_na=False)
    audit = audit.set_index("record_id", drop=False)
    rows: list[dict] = []

    for record_id, candidate in config["candidates"].items():
        if record_id not in audit.index:
            raise ValueError(f"Missing stage-18 audit record: {record_id}")
        source = audit.loc[record_id]
        authoritative = bool(source["source_url"])
        cell_identified = bool(source["cell_line"])
        coordinate_resolved = bool(source["reference"]) and "must be resolved" not in source["reference"]
        structure_semantics = bool(source["value_contract"])
        matched_m6a = source["admission"] not in {
            "AUDIT_ONLY_NO_MATCHED_QUANTITATIVE_M6A",
            "NOT_ADMITTED",
        }
        no_pooling = candidate["raw_structure_pooling_allowed"] is False
        admitted = all(
            [authoritative, cell_identified, coordinate_resolved, structure_semantics, matched_m6a, no_pooling]
        )
        failed = []
        if not coordinate_resolved:
            failed.append("coordinate_reference_resolved")
        if not matched_m6a:
            failed.append("matched_quantitative_m6a_available")
        rows.append(
            {
                "record_id": record_id,
                "accession": source["accession"],
                "cell_line": source["cell_line"],
                "structure_technology": candidate["structure_technology"],
                "authoritative_source_recorded": authoritative,
                "cell_condition_identified": cell_identified,
                "coordinate_reference_resolved": coordinate_resolved,
                "quantitative_structure_semantics_identified": structure_semantics,
                "matched_quantitative_m6a_available": matched_m6a,
                "no_cross_technology_raw_value_pooling": no_pooling,
                "decision": "ADMIT_FOR_ASSOCIATION" if admitted else "HOLD_NOT_ADMITTED",
                "failed_gates": ";".join(failed),
                "next_trigger": candidate["required_modification_measurement"],
            }
        )

    out = config["outputs"]
    write_csv(ROOT / out["decisions"], rows)
    expected_holds = all(row["decision"] == "HOLD_NOT_ADMITTED" for row in rows)
    status = "PASS" if expected_holds else "REVIEW_REQUIRED"
    write_csv(
        ROOT / out["status"],
        [{
            "stage": 21,
            "status": status,
            "candidates_reviewed": len(rows),
            "admitted": sum(row["decision"] == "ADMIT_FOR_ASSOCIATION" for row in rows),
            "held": sum(row["decision"] == "HOLD_NOT_ADMITTED" for row in rows),
            "checked_utc": datetime.now(timezone.utc).isoformat(),
        }],
    )
    decisions = {row["record_id"]: row for row in rows}
    report = f"""# Stage 21: Conditional structure source access conclusion

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {datetime.now(timezone.utc).date().isoformat()}
- Verification Status: {'VERIFIED' if status == 'PASS' else 'REVIEW_REQUIRED'}
- Version Label: phase2_stage21_v1

Structural data lacking matching quantitative m6A were not forced into the correlation analysis at this stage. Human RNA MaP is judged to be `{decisions['human_rna_map']['decision']}`; GSE50676 PARS is judged to be `{decisions['pars_gse50676']['decision']}`.

Human RNA MaP has documented U2OS, DMS-TRAM-seq and hg38/GRCh38.106 conditions, but currently lacks quantitative m6A matching conditions. PARS data also require resolution of source-specific reference versions and again lack matching conditions for quantitative m6A. Neither may be merged directly with the icSHAPE original value.

These HOLD determinations are the result of the execution of preset condition thresholds in the second period; future re-admittances will only be made when quantitative m6A data and full coordinate contracts for the same cell conditions are simultaneously in place."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
