"""Stage 30: enforce the four-layer evidence gradient and model gates."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/30_phase4_evidence_and_model_gate.yaml"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    stage29 = pd.read_csv(ROOT / config["required_stage_status"]).iloc[0]
    if stage29.status != "PASS":
        raise ValueError("Stage 29 did not pass")
    results = pd.read_csv(ROOT / "results/tables/29_phase4_directionality_endpoints.csv")
    primary = results[results.endpoint.eq("mean_reactivity")].iloc[0]
    layers = [
        {"layer": "Association", "status": "COMPLETE", "evidence": "Frozen HEK293T-HeLa observational baseline and phase-3 stratified analyses", "claim_ceiling": "association"},
        {"layer": "Directionality", "status": "COMPLETE", "evidence": f"GSE264642 STM2457 paired shift; n={int(primary.n_transcripts)}", "claim_ceiling": "directionality_clue"},
        {"layer": "Causal evidence", "status": "NOT_IDENTIFIED", "evidence": "No rescue/inactive-compound control and no site-level same-sample m6A change", "claim_ceiling": "not_permitted"},
        {"layer": "Mechanism", "status": "NOT_IDENTIFIED", "evidence": "No site-directed mutation plus orthogonal structure/biochemical validation generated here", "claim_ceiling": "not_permitted"},
    ]
    req = config["advanced_model_requirements"]
    obs = config["observed_corpus"]
    gates = [
        {"gate": "multiple_modifications", "required": req["minimum_modifications"], "observed": obs["perturbation_modifications"], "pass": obs["perturbation_modifications"] >= req["minimum_modifications"]},
        {"gate": "independent_datasets", "required": req["minimum_independent_datasets"], "observed": obs["perturbation_datasets"], "pass": obs["perturbation_datasets"] >= req["minimum_independent_datasets"]},
        {"gate": "conditions", "required": req["minimum_conditions"], "observed": obs["treatment_conditions"], "pass": obs["treatment_conditions"] >= req["minimum_conditions"]},
        {"gate": "independent_external_test", "required": True, "observed": obs["independent_external_test"], "pass": bool(obs["independent_external_test"])},
    ]
    claims = [
        {"topic": "STM2457_structure_shift", "permitted": "STM2457 exposure preceded a paired transcript-level structure shift in this HEK293T dataset.", "forbidden": "m6A loss causally remodels RNA structure transcriptome-wide.", "causal": False, "mechanism": False},
        {"topic": "site_specificity", "permitted": "The source supports transcript-level structural summaries and same-experiment m6A reduction validation.", "forbidden": "The analysis identifies which m6A sites caused each structural change.", "causal": False, "mechanism": False},
        {"topic": "advanced_models", "permitted": "Foundation-model, Transformer, and hyperparameter-search branches remain gated off.", "forbidden": "Advanced modeling was negative after a fair external benchmark.", "causal": False, "mechanism": False},
    ]
    out = config["outputs"]
    write_csv(ROOT / out["evidence_gradient"], layers)
    write_csv(ROOT / out["advanced_model_gate"], gates)
    write_csv(ROOT / out["claim_boundaries"], claims)
    now = datetime.now(timezone.utc)
    write_csv(ROOT / out["status"], [{
        "stage": 30, "status": "PASS", "directionality_layer_complete": True,
        "causal_layer_identified": False, "mechanism_layer_identified": False,
        "advanced_model_gates_passed": sum(bool(row["pass"]) for row in gates),
        "advanced_model_gates_total": len(gates), "advanced_models_started": False,
        "finished_utc": now.isoformat(),
    }])
    failed = [row["gate"] for row in gates if not row["pass"]]
    report = f"""# Stage 30: Four-layer evidence gradient and advanced model threshold

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/validate
- Origin Date: {now.date().isoformat()}
- Verification Status: VERIFIED
- Version Label: phase4_stage30_v1

The fourth issue has advanced the evidence to **Directionality**, but neither Causal evidence nor Mechanism are identified in the current public data. The "unidentified" here is not a negative causal conclusion, but a design capability boundary.

Base model/Transformer and large-scale hyperparameter search are not started. The failure threshold is: {', '.join(failed)}. A single modification, a single perturbed data set, and two treatment conditions are not sufficient to answer the cross-technology/cross-modification generalization question, and there is no truly independent external test set; therefore not doing advanced modeling is the result of the implementation of the preset threshold, not the negative model performance."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
