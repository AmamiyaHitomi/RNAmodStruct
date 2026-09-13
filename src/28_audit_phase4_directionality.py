"""Stage 28: audit phase-4 perturbation sources and freeze the analysis contract."""

from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/28_phase4_directionality.yaml"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    candidate = config["primary_candidate"]
    files = {arm: ROOT / rel for arm, rel in candidate["files"].items()}
    missing = [str(path.relative_to(ROOT)) for path in files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing phase-4 inputs: {missing}")

    tables = {arm: pd.read_csv(path) for arm, path in files.items()}
    expected = {
        "treatment": {"transcript", "stm_average_max", "stm_average_average", "stm_average_std", "stm_average_gini"},
        "control": {"transcript", "wt_average_max", "wt_average_average", "wt_average_std", "wt_average_gini"},
    }
    audits = []
    for arm, path in files.items():
        table = tables[arm]
        columns_ok = set(table.columns) == expected[arm]
        audits.append({
            "candidate": "GSE264642", "arm": arm, "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size, "sha256": sha256_file(path), "rows": len(table),
            "unique_transcripts": table.transcript.nunique(), "duplicate_transcripts": int(table.transcript.duplicated().sum()),
            "missing_cells": int(table.isna().sum().sum()), "columns_ok": columns_ok,
            "cell_line": "HEK293T", "decision": "ADMIT_DIRECTIONALITY" if columns_ok else "REJECT_SCHEMA",
        })
    legacy = config["legacy_candidate"]
    audits.append({
        "candidate": f"{legacy['modification_accession']}+{legacy['structure_accession']}", "arm": "WT_vs_KO",
        "path": "data/raw_processed", "bytes": "", "sha256": "", "rows": "", "unique_transcripts": "",
        "duplicate_transcripts": "", "missing_cells": "", "columns_ok": True,
        "cell_line": "v6.5_WT_vs_J1_KO", "decision": legacy["decision"],
    })

    common = set(tables["treatment"].transcript) & set(tables["control"].transcript)
    if len(common) != len(tables["treatment"]) or any(not row["columns_ok"] for row in audits):
        raise ValueError("Primary source pairing or schema audit failed")

    dag = [
        {"source": "STM2457", "target": "METTL3_activity", "edge": "inhibits", "status": "design_intended"},
        {"source": "METTL3_activity", "target": "m6A_level", "edge": "increases", "status": "validated_by_meRIP_qPCR"},
        {"source": "m6A_level", "target": "RNA_structure", "edge": "candidate_effect", "status": "target_estimand"},
        {"source": "STM2457", "target": "RNA_structure", "edge": "off_target_path", "status": "uncontrolled"},
        {"source": "RNA_abundance", "target": "observed_transcripts", "edge": "selection_path", "status": "uncontrolled"},
        {"source": "RBP_occupancy", "target": "RNA_structure", "edge": "alternative_mediator", "status": "uncontrolled"},
    ]
    endpoints = [
        {"endpoint": "mean_reactivity", "role": "primary", "contrast": "STM2457_minus_vehicle",
         "estimand": "paired_median_shift", "test": "wilcoxon_signed_rank", "multiplicity_family": "mean_and_gini"},
        {"endpoint": "reactivity_gini", "role": "secondary", "contrast": "STM2457_minus_vehicle",
         "estimand": "paired_median_shift", "test": "wilcoxon_signed_rank", "multiplicity_family": "mean_and_gini"},
        {"endpoint": "lower_mean_fraction", "role": "secondary_descriptive", "contrast": "STM2457_less_than_vehicle",
         "estimand": "paired_fraction", "test": "exact_binomial", "multiplicity_family": "descriptive"},
    ]
    out = config["outputs"]
    write_csv(ROOT / out["source_audit"], audits)
    write_csv(ROOT / out["dag"], dag)
    write_csv(ROOT / out["endpoint_registry"], endpoints)
    now = datetime.now(timezone.utc)
    write_csv(ROOT / out["status"], [{
        "stage": 28, "status": "PASS", "primary_candidate": "GSE264642",
        "paired_transcripts": len(common), "legacy_candidate_decision": legacy["decision"],
        "admitted_evidence_layer": config["admission"]["admitted_evidence_layer"], "finished_utc": now.isoformat(),
    }])
    prereg = f"""# The fourth phase of directional analysis pre-registration and causal diagram contract

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/plan
- Origin Date: {now.date().isoformat()}
- Verification Status: LOCAL_AUDITED
- Version Label: phase4_directionality_prereg_v1

## Research questions and main indicators

Is there a pairwise shift in transcript average SHAPE reactivity relative to the DMSO control after 24 h of 40 μM STM2457 treatment in the same HEK293T background? The primary estimator was the median within-transcript difference of STM2457 minus control, using transcript resampling confidence intervals and the Wilcoxon paired test. Gini displacement was the secondary endpoint, and BH correction was performed for both tests.

## Admission and Exclusion

- GSE264642: Same cell line, drug and vehicle control, 2 biological replicates for each probe condition; m6A reduction was verified by meRIP-qPCR in the same experiment, and the access was directional.
- GSE52662 + GSE60034: The structure WT is v6.5 and KO is J1, which violates the same cell background and is excluded from causal statistics.

## Evidence Boundary

STM2457 is not a randomly assigned multi-batch experiment and lacks site-by-site sample m6A quantification, rescue, and inactive compound controls; processing may also affect structure through RNA abundance, RBP occupancy, or off-target pathways. Therefore, the results are at most directional clues and cannot be written as proof of causal effects or mechanisms. The analysis was a retrospective exploratory review of published results."""
    (ROOT / out["preregistration"]).write_text(prereg, encoding="utf-8")


if __name__ == "__main__":
    main()
