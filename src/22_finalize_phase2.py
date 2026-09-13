"""Stage 22: verify phase-2 exit gates and create a reproducible release manifest."""

from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "22_phase2_completion.yaml"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def status_row(relative: str) -> dict:
    table = pd.read_csv(ROOT / relative)
    if len(table) != 1:
        raise ValueError(f"Expected one status row: {relative}")
    row = table.iloc[0].to_dict()
    if row.get("status") != "PASS":
        raise ValueError(f"Required stage did not pass: {relative}")
    return row


def release_paths() -> list[Path]:
    paths: set[Path] = set()
    for stage in range(18, 23):
        paths.update((ROOT / "config").glob(f"{stage}_*.yaml"))
        paths.update((ROOT / "src").glob(f"{stage}_*.py"))
        paths.update((ROOT / "results" / "reports").glob(f"{stage}_*"))
        paths.update((ROOT / "results" / "status").glob(f"{stage}_*"))
        paths.update((ROOT / "results" / "tables").glob(f"{stage}_*"))
        paths.update((ROOT / "data" / "final").glob(f"{stage}_*"))
    paths.add(ROOT / "metadata" / "manifests" / "source" / "phase2_source_manifest.csv")
    audit = pd.read_csv(ROOT / "metadata" / "manifests" / "source" / "phase2_source_manifest.csv", keep_default_na=False)
    for relative in audit.loc[audit.local_file.ne(""), "local_file"]:
        paths.add(ROOT / relative)
    return sorted(path for path in paths if path.is_file())


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    statuses = {int(status_row(path)["stage"]): status_row(path) for path in config["required_stage_statuses"]}

    stage19 = pd.read_csv(ROOT / "results/tables/19_paired_condition_associations.csv")
    stage20 = pd.read_csv(ROOT / "results/tables/20_cross_technology_associations.csv")
    stage21 = pd.read_csv(ROOT / "results/tables/21_conditional_source_admission.csv")
    audit = pd.read_csv(ROOT / "results/tables/18_public_data_audit.csv", keep_default_na=False)

    evidence = [
        {
            "stage": 18,
            "evidence_axis": "source_and_coordinate_audit",
            "primary_population": "3 downloaded public files",
            "result": f"verified={int(statuses[18]['files_verified'])}; paired_ready={int(statuses[18]['paired_datasets_ready'])}",
            "interpretation_boundary": "Audit establishes provenance and coordinate contracts, not biological replication.",
        },
        {
            "stage": 19,
            "evidence_axis": "structure_condition_contrast",
            "primary_population": f"{int(statuses[19]['paired_sites'])} HEK293T paired sites",
            "result": "in-vivo and in-vitro icSHAPE effects estimated separately with a paired contrast",
            "interpretation_boundary": "Observational condition contrast; not a causal perturbation.",
        },
        {
            "stage": 20,
            "evidence_axis": "cross_m6a_measurement_technology",
            "primary_population": f"{int(statuses[20]['hela_primary_sites'])} HeLa sites",
            "result": "GLORI and SAC-seq outcomes analyzed separately; HEK293/HEK293T retained as near-match sensitivity only",
            "interpretation_boundary": "Matched cell line across studies, not the same biological samples.",
        },
        {
            "stage": 21,
            "evidence_axis": "conditional_source_admission",
            "primary_population": f"{len(stage21)} deferred sources",
            "result": f"held={int(statuses[21]['held'])}; admitted={int(statuses[21]['admitted'])}",
            "interpretation_boundary": "HOLD prevents unmatched structure-only sources from being presented as m6A replication.",
        },
    ]

    downloaded = audit[audit.local_file.ne("")]
    gate_checks = {
        "authoritative_source_and_coordinate_audit": bool(
            len(downloaded) == 3
            and downloaded.source_url.ne("").all()
            and downloaded.integrity.eq("GZIP_FULL_STREAM_OK").all()
        ),
        "invivo_invitro_or_cross_technology_replication": bool(
            statuses[19]["paired_sites"] > 0 or statuses[20]["hela_primary_sites"] > 0
        ),
        "technology_stratified_effects_reported": bool(
            {"in_vivo", "in_vitro"}.issubset(set(stage19.analysis))
            and {"GLORI", "m6A-SAC-seq"}.issubset(set(stage20.outcome_technology))
        ),
        "incompatible_raw_scales_not_pooled": bool(
            stage21.no_cross_technology_raw_value_pooling.all()
            and set(stage20.outcome_technology) == {"GLORI", "m6A-SAC-seq"}
        ),
    }
    declared = config["exit_gates"]
    if set(declared) != set(gate_checks):
        raise ValueError("Configured and implemented exit gates differ")
    gate_rows = [
        {
            "gate": gate,
            "required": bool(declared[gate]),
            "observed": observed,
            "status": "PASS" if observed else "FAIL",
        }
        for gate, observed in gate_checks.items()
    ]

    out = config["outputs"]
    write_csv(ROOT / out["evidence_summary"], evidence)
    write_csv(ROOT / out["exit_gates"], gate_rows)
    complete = all(row["status"] == "PASS" for row in gate_rows)
    status = "PASS" if complete else "FAIL"
    now = datetime.now(timezone.utc)
    write_csv(
        ROOT / out["status"],
        [{
            "stage": 22,
            "status": status,
            "required_stages_passed": len(statuses),
            "exit_gates_passed": sum(row["status"] == "PASS" for row in gate_rows),
            "exit_gates_total": len(gate_rows),
            "finished_utc": now.isoformat(),
        }],
    )

    invivo = stage19[stage19.analysis.eq("in_vivo")].iloc[0]
    invitro = stage19[stage19.analysis.eq("in_vitro")].iloc[0]
    hela_glori = stage20[(stage20.dataset.eq("HeLa")) & (stage20.outcome_technology.eq("GLORI"))].iloc[0]
    hela_sac = stage20[(stage20.dataset.eq("HeLa")) & (stage20.outcome_technology.eq("m6A-SAC-seq"))].iloc[0]
    report = f"""# Stage 22: Second Phase Completion and Exit Threshold Report

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {now.date().isoformat()}
- Verification Status: {'VERIFIED' if complete else 'UNVERIFIED'}
- Version Label: phase2_stage22_v1

The status of stages 18-21 of the second phase is all PASS, and the exit thresholds for {sum(row['status'] == 'PASS' for row in gate_rows)}/{len(gate_rows)} items have passed. Frozen 06–10 Mainline not being written back.

## Hierarchical results

- HEK293T isosite in-vivo / in-vitro control has a total of {int(statuses[19]['paired_sites']):,} sites. The structural effects were {invivo.beta_per_sd:.6f}/SD (95% CI {invivo.ci95_low_per_sd:.6f}–{invivo.ci95_high_per_sd:.6f}) and {invitro.beta_per_sd:.6f}/SD (95% CI {invitro.ci95_low_per_sd:.6f}–{invitro.ci95_high_per_sd:.6f}), respectively.
- HeLa cross-m6A technology review of a total of {int(statuses[20]['hela_primary_sites']):,} sites. The effects of GLORI and SAC-seq on the same experimental structural features are {hela_glori.beta_per_sd:.6f}/SD and {hela_sac.beta_per_sd:.6f}/SD respectively; the two outcomes are modeled separately without merging the original dimensions.
- Human RNA MaP and GSE50676 PARS remain HOLD due to lack of quantitative m6A matching conditions. This result is the implementation of conditional access thresholds and does not mean that the data itself is of unqualified quality.

## Evidence Boundary

Stage 19 is an observational comparison of structural conditions; Stage 20 is a cross-study and cross-measurement technology review, not an experimental replication of the same biological sample. The current results support the cautious statement of “small effect and technique/condition dependent” and do not support causal or mechanistic conclusions."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")

    manifest_rows = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in release_paths()
        if path != ROOT / out["manifest"]
    ]
    write_csv(ROOT / out["manifest"], manifest_rows)
    if not complete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
