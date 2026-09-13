"""Stage 27: verify phase-3 exit gates and freeze its release manifest."""

from __future__ import annotations

import csv
import hashlib
import math
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "27_phase3_completion.yaml"


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
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
    for stage in range(23, 28):
        paths.update((ROOT / "config").glob(f"{stage}_*.yaml"))
        paths.update((ROOT / "src").glob(f"{stage}_*.py"))
        paths.update((ROOT / "results" / "reports").glob(f"{stage}_*"))
        paths.update((ROOT / "results" / "status").glob(f"{stage}_*"))
        paths.update((ROOT / "results" / "tables").glob(f"{stage}_*"))
        paths.update((ROOT / "data" / "interim").glob(f"{stage}_*"))
        paths.update((ROOT / "data" / "final").glob(f"{stage}_*"))
    paths.update((ROOT / "tests").glob("test_phase3_*.py"))
    paths.update((ROOT / "metadata" / "manifests").rglob("phase3_*.csv"))
    download_manifest = pd.read_csv(ROOT / "metadata" / "manifests" / "download" / "phase3_download_manifest.csv")
    for relative in download_manifest.path:
        paths.add(ROOT / relative)
    return sorted(path for path in paths if path.is_file())


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    statuses = {int(status_row(path)["stage"]): status_row(path) for path in config["required_stage_statuses"]}
    source_cards = pd.read_csv(ROOT / "results/tables/23_multimodification_source_cards.csv")
    file_audit = pd.read_csv(ROOT / "results/tables/24_multimodification_file_audit.csv")
    overlap = pd.read_csv(ROOT / "results/tables/24_multimodification_overlap_audit.csv")
    attrition = pd.read_csv(ROOT / "results/tables/25_multimodification_structure_attrition.csv")
    nm_coord = pd.read_csv(ROOT / "results/tables/25_nm_coordinate_sensitivity.csv")
    associations = pd.read_csv(ROOT / "results/tables/26_stratified_primary_associations.csv")
    diagnostics = pd.read_csv(ROOT / "results/tables/26_stratified_model_diagnostics.csv")
    low_n = pd.read_csv(ROOT / "results/tables/26_descriptive_low_n_tracks.csv")
    matching = pd.read_csv(ROOT / "results/tables/26_nm_matching_qc.csv")

    primary = associations[associations.evidence_tier.eq("primary")].copy()
    held = source_cards[source_cards.decision.str.startswith("HOLD") | source_cards.decision.str.startswith("CONDITIONAL")]
    evidence = [
        {
            "stage": 23, "evidence_axis": "multi_modification_source_scope",
            "population": f"{len(source_cards)} candidate modifications",
            "result": f"download-ready=3; held-or-conditional={len(held)}",
            "interpretation_boundary": "Source scouting is eligibility triage, not biological evidence.",
        },
        {
            "stage": 24, "evidence_axis": "download_coordinate_and_overlap_audit",
            "population": f"{len(file_audit)} downloaded files; {len(overlap)} normalized tracks",
            "result": "all hashes verified; m5C/m7G/Nm admitted with modality-specific contracts",
            "interpretation_boundary": "Nm corrected GEO files are treated as hg38; m7G remains interval-level.",
        },
        {
            "stage": 25, "evidence_axis": "structure_feature_extraction",
            "population": f"{int(attrition.main_valid_records.sum())} valid track-records",
            "result": f"inferential-ready tracks={int(attrition.track_admission.eq('INFERENTIAL_READY').sum())}; descriptive-only={int(attrition.track_admission.eq('DESCRIPTIVE_ONLY_LOW_N').sum())}",
            "interpretation_boundary": "HEK293T m7G n=18 is descriptive only; HEK Nm is a near-cell-label sensitivity.",
        },
        {
            "stage": 26, "evidence_axis": "technology_stratified_associations",
            "population": f"{len(primary)} HeLa primary analyses; {int(matching.matched_pairs.sum())} Nm matched pairs",
            "result": f"primary q range={primary.q_value_primary.min():.6g}-{primary.q_value_primary.max():.6g}; all primary CIs include zero={all((primary.ci95_low <= 0) & (primary.ci95_high >= 0))}",
            "interpretation_boundary": "No clear association was detected; this does not establish exact zero effects or causality.",
        },
    ]

    claim_policy = config["claim_policy"]
    claims = [
        {
            "claim_id": "P3_C1", "topic": "HeLa_m5C_structure",
            "permitted_statement": "Within the audited HeLa data, the adjusted m5C estimate was small and its confidence interval included zero.",
            "forbidden_statement": "RNA structure has no relationship with m5C.",
            "causal_claim_permitted": False, "exact_zero_claim_permitted": False,
            "evidence_tier": "primary_observational",
        },
        {
            "claim_id": "P3_C2", "topic": "HeLa_m7G_structure",
            "permitted_statement": "No clear interval-level association between peak structure and m7G enrichment was detected.",
            "forbidden_statement": "Single-nucleotide m7G structure effects were tested or ruled out.",
            "causal_claim_permitted": False, "exact_zero_claim_permitted": False,
            "evidence_tier": "primary_observational_interval",
        },
        {
            "claim_id": "P3_C3", "topic": "HeLa_Nm_structure",
            "permitted_statement": "Nm sites and matched same-transcript same-base controls showed a small difference with an interval crossing zero.",
            "forbidden_statement": "Nm does not affect or depend on RNA structure.",
            "causal_claim_permitted": False, "exact_zero_claim_permitted": False,
            "evidence_tier": "primary_observational_matched",
        },
        {
            "claim_id": "P3_C4", "topic": "HEK_sensitivity_tracks",
            "permitted_statement": "HEK Nm is a near-cell-label sensitivity and HEK293T m7G is descriptive because n=18.",
            "forbidden_statement": "HEK and HEK293T provide independent confirmatory replication for both modifications.",
            "causal_claim_permitted": False, "exact_zero_claim_permitted": False,
            "evidence_tier": "sensitivity_or_descriptive",
        },
        {
            "claim_id": "P3_C5", "topic": "cross_modification_synthesis",
            "permitted_statement": "Across separate technology-specific analyses, all estimated directions were small and negative with uncertainty spanning zero.",
            "forbidden_statement": "A pooled universal multi-modification effect size was estimated.",
            "causal_claim_permitted": False, "exact_zero_claim_permitted": False,
            "evidence_tier": "narrative_synthesis_no_raw_pooling",
        },
    ]

    gate_checks = {
        "source_scope_and_download_integrity": bool(
            statuses[23]["sources_scouted"] == 7 and statuses[24]["files_verified"] == 5
            and len(file_audit) == 5 and file_audit.sha256_match.all()
        ),
        "coordinate_and_reference_contracts_resolved": bool(
            overlap.coordinate_success_fraction.eq(1.0).all()
            and overlap[overlap.modification.eq("m5C")].reference_base_match_fraction.eq(1.0).all()
            and nm_coord.stability_pass.all()
        ),
        "structure_feature_extraction_completed": bool(
            statuses[25]["tracks_built"] == 5 and statuses[25]["distinct_modifications_admitted"] == 3
            and set(attrition.modification) == {"m5C", "m7G", "Nm"}
        ),
        "technology_stratified_inference_completed": bool(
            statuses[26]["primary_analyses_completed"] == 3
            and set(primary.track_id) == {"m5c_hela", "m7g_hela", "nm_hela"}
            and not bool(statuses[26]["raw_scales_pooled"])
        ),
        "multiple_testing_and_uncertainty_reported": bool(
            len(primary) == 3 and primary.q_value_primary.map(math.isfinite).all()
            and primary.ci95_low.map(math.isfinite).all() and primary.ci95_high.map(math.isfinite).all()
            and len(diagnostics) == 2
        ),
        "low_n_near_match_and_claim_boundaries_preserved": bool(
            len(low_n) == 1 and low_n.analysis_role.eq("DESCRIPTIVE_ONLY_LOW_N").all()
            and associations[associations.track_id.eq("nm_hek")].evidence_tier.eq("near_cell_match_sensitivity").all()
            and not any(row["causal_claim_permitted"] for row in claims)
            and not any(row["exact_zero_claim_permitted"] for row in claims)
            and not claim_policy["cross_modification_raw_effect_pooling_permitted"]
        ),
    }
    if set(gate_checks) != set(config["exit_gates"]):
        raise ValueError("Configured and implemented phase-3 gates differ")
    gates = [{
        "gate": gate, "required": bool(config["exit_gates"][gate]), "observed": observed,
        "status": "PASS" if observed else "FAIL",
    } for gate, observed in gate_checks.items()]
    complete = all(row["status"] == "PASS" for row in gates)
    now = datetime.now(timezone.utc)
    out = config["outputs"]
    write_csv(ROOT / out["evidence_summary"], evidence)
    write_csv(ROOT / out["exit_gates"], gates)
    write_csv(ROOT / out["claim_boundaries"], claims)
    write_csv(ROOT / out["status"], [{
        "stage": 27, "status": "PASS" if complete else "FAIL",
        "required_stages_passed": len(statuses), "exit_gates_passed": sum(row["status"] == "PASS" for row in gates),
        "exit_gates_total": len(gates), "primary_findings_fdr_significant": int((primary.q_value_primary < 0.05).sum()),
        "next_phase_gate": config["next_phase"]["gate"] if complete else "BLOCKED",
        "finished_utc": now.isoformat(),
    }])

    rows = {row.track_id: row for row in associations.itertuples(index=False)}
    report = f"""# Stage 27: Third Phase Completion and Exit Threshold Report

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {now.date().isoformat()}
- Verification Status: {'VERIFIED' if complete else 'UNVERIFIED'}
- Version Label: phase3_stage27_v1

Phases 23–26 of the third phase are all PASS, and the exit thresholds for {sum(row['status'] == 'PASS' for row in gates)}/{len(gates)} are passed. The official status of Issue 3 is `{'COMPLETE' if complete else 'INCOMPLETE'}`; the original multi-modifier dimensions were never merged.

## Main conclusions

- **m5C/HeLa**: n={rows['m5c_hela'].n}, each 1 SD structural change corresponds to a proportional change in m5C {rows['m5c_hela'].effect_per_predictor_sd:.6f} (95% CI {rows['m5c_hela'].ci95_low:.6f} to {rows['m5c_hela'].ci95_high:.6f}; FDR q={rows['m5c_hela'].q_value_primary:.6g}).
- **m7G/HeLa**: n={rows['m7g_hela'].n}, interval-level log2 enrichment change {rows['m7g_hela'].effect_per_predictor_sd:.6f}/structure SD (95% CI {rows['m7g_hela'].ci95_low:.6f} to {rows['m7g_hela'].ci95_high:.6f}; q={rows['m7g_hela'].q_value_primary:.6g}). This result cannot be explained by a single-base m7G effect.
- **Nm/HeLa**: {rows['nm_hela'].n} iso-transcripts, iso-base-paired, poorly structured {rows['nm_hela'].effect_per_predictor_sd:.6f} (cluster bootstrap 95% CI {rows['nm_hela'].ci95_low:.6f} to {rows['nm_hela'].ci95_high:.6f}; q={rows['nm_hela'].q_value_primary:.6g}).
- None of the three main HeLa tests passed FDR 0.05, and the confidence intervals all spanned zero. The safest summary is: **Under the current public data, technology stratification model and coverage threshold, no clear structure-modification association has been detected; this does not prove that the true effect is strictly zero. **

## Keep restrictions

- m7G/HEK293T has only 18 qualified intervals, which are only described and not inferred.
- Nm/HEK is an approximate cell label for HEK versus HEK293T and can only be used as a sensitivity result.
- m5C, m7G, Nm are derived from cross-study observational data; detection efficiency, expression levels and local sequence may still cause residual confounding.
- m1A, Ψ, ac4C and A-to-I stillMaintain HOLD/CONDITIONAL status for technical validity, file size, control reanalysis, or cell matching.

## Next stage

`{config['next_phase']['gate']}`: The fourth phase only authorizes entry into perturbation or verification plan planning with controls, which does not mean that there is causal evidence, nor does it automatically authorize large-scale data downloads."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")

    manifest = [{
        "path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256_file(path),
    } for path in release_paths() if path != ROOT / out["manifest"]]
    write_csv(ROOT / out["manifest"], manifest)
    if not complete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
