"""Stage 31: verify phase-4 gates and freeze the release manifest."""

from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/31_phase4_completion.yaml"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def release_paths() -> list[Path]:
    paths: set[Path] = set()
    for stage in range(28, 32):
        paths.update((ROOT / "config").glob(f"{stage}_*.yaml"))
        paths.update((ROOT / "src").glob(f"{stage}_*.py"))
        paths.update((ROOT / "results" / "reports").glob(f"{stage}_*"))
        paths.update((ROOT / "results" / "status").glob(f"{stage}_*"))
        paths.update((ROOT / "results/tables").glob(f"{stage}_*"))
    paths.update((ROOT / "tests").glob("test_phase4_*.py"))
    paths.update((ROOT / "metadata" / "decisions").glob("phase4_*.md"))
    paths.update((ROOT / "data/raw_processed").glob("GSE264642_*.csv.gz"))
    return sorted(path for path in paths if path.is_file())


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    statuses = [pd.read_csv(ROOT / path).iloc[0] for path in config["required_stage_statuses"]]
    if any(row.status != "PASS" for row in statuses):
        raise ValueError("A required phase-4 stage did not pass")
    source = pd.read_csv(ROOT / "results/tables/28_phase4_perturbation_source_audit.csv")
    endpoints = pd.read_csv(ROOT / "results/tables/29_phase4_directionality_endpoints.csv")
    layers = pd.read_csv(ROOT / "results/tables/30_phase4_evidence_gradient.csv")
    claims = pd.read_csv(ROOT / "results/tables/30_phase4_claim_boundaries.csv")
    model_gates = pd.read_csv(ROOT / "results/tables/30_phase4_advanced_model_gate.csv")
    checks = {
        "perturbation_source_audited": bool((source.candidate == "GSE264642").sum() == 2 and source[source.candidate == "GSE264642"].columns_ok.all()),
        "directionality_analysis_completed": bool({"mean_reactivity", "reactivity_gini", "lower_mean_fraction"}.issubset(set(endpoints.endpoint))),
        "evidence_layers_separated": bool(set(layers.layer) == {"Association", "Directionality", "Causal evidence", "Mechanism"}),
        "causal_and_mechanism_overclaim_prevented": bool(not claims.causal.any() and not claims.mechanism.any()),
        "advanced_model_gate_enforced": bool(not model_gates["pass"].all()),
    }
    gates = [{"gate": key, "required": True, "observed": value, "status": "PASS" if value else "FAIL"} for key, value in checks.items()]
    complete = all(checks.values())
    primary = endpoints[endpoints.endpoint.eq("mean_reactivity")].iloc[0]
    evidence = [
        {"evidence_axis": "source_admission", "result": "GSE264642 admitted for directionality; GSE52662+GSE60034 rejected for causal analysis", "evidence_ceiling": "directionality"},
        {"evidence_axis": "paired_structure_shift", "result": f"n={int(primary.n_transcripts)}; median shift={primary.paired_median_shift:.8g}; BH q={primary.q_value_bh:.6g}", "evidence_ceiling": "directionality"},
        {"evidence_axis": "causal_mechanism", "result": "not identified by available design", "evidence_ceiling": "not_permitted"},
        {"evidence_axis": "advanced_modeling", "result": "not started because corpus and external-test gates failed", "evidence_ceiling": "gate_decision_only"},
    ]
    out = config["outputs"]
    write_csv(ROOT / out["exit_gates"], gates)
    write_csv(ROOT / out["evidence_summary"], evidence)
    now = datetime.now(timezone.utc)
    write_csv(ROOT / out["status"], [{
        "stage": 31, "status": "PASS" if complete else "FAIL",
        "phase4_scope_status": "COMPLETE_DIRECTIONALITY_ONLY" if complete else "INCOMPLETE",
        "exit_gates_passed": sum(checks.values()), "exit_gates_total": len(checks),
        "highest_evidence_layer": "Directionality", "causal_evidence_identified": False,
        "mechanism_identified": False, "advanced_models_started": False, "finished_utc": now.isoformat(),
    }])
    report = f"""# Stage 31: Fourth Stage Completion and Exit Threshold Report

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run+validate
- Origin Date: {now.date().isoformat()}
- Verification Status: {'VERIFIED' if complete else 'UNVERIFIED'}
- Version Label: phase4_stage31_v1

Phases 28–30 of the fourth period are all PASS, and the {sum(checks.values())}/{len(checks)} item scope exit threshold is passed. The scope status is `{'COMPLETE_DIRECTIONALITY_ONLY' if complete else 'INCOMPLETE'}`.

## Core conclusion

- GSE264642 provides structural comparison of STM2457 and vehicle in the same HEK293T background, with access to the directional layer; a total of {int(primary.n_transcripts):,} transcripts are paired.
- The mean SHAPE reactive pairwise median shift for STM2457 minus vehicle is {primary.paired_median_shift:.6g} (bootstrap 95% CI {primary.bootstrap_ci95_low:.6g} to {primary.bootstrap_ci95_high:.6g}; BH q={primary.q_value_bh:.6g}).
- GSE52662/GSE60034 did not enter the causal statistics due to inconsistent cell backgrounds of structures WT=v6.5 and KO=J1.
- The current highest level of evidence is Directionality; no defensible Causal evidence or Mechanism has been identified.
- Transformer, basic model and large-scale hyperparameter search are closed according to the preset threshold, and "not started" is not mistakenly written as negative performance.

## Final boundary

can report "systematic pairwise shifts in the aggregate amount of transcript structure observed upon METTL3 inhibition". It is not possible to report that "m6A deletion causes this structural change", nor to locate a site-level mechanism. If the site-by-site m6A quantification, rescue/inactive compound control and independent review data of the same sample are obtained in the future, it should be regarded as a new confirmatory stage rather than writing back the exploratory results of this time."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")
    manifest_path = ROOT / out["manifest"]
    manifest = [{"path": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size, "sha256": sha256_file(p)} for p in release_paths() if p != manifest_path]
    write_csv(manifest_path, manifest)
    if not complete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
