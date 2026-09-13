"""Stage 23: materialize the phase-3 multi-modification source scout."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "23_multimodification_source_scout.yaml"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def valid_source_url(url: str) -> bool:
    parsed = urlparse(url)
    return parsed.scheme in {"http", "https"} and bool(parsed.netloc)


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    sources = config["sources"]
    expected = set(config["required_modifications"])
    observed = {row["modification"] for row in sources}
    if expected != observed or len(sources) != len(expected):
        raise ValueError("Each required modification must have exactly one evidence card")
    if len({row["record_id"] for row in sources}) != len(sources):
        raise ValueError("record_id values must be unique")
    if not all(valid_source_url(row["source_url"]) and valid_source_url(row["evidence_url"]) for row in sources):
        raise ValueError("Every source requires resolvable URL syntax")

    queue = config["download_queue"]
    source_ids = {row["record_id"] for row in sources}
    if not all(row["record_id"] in source_ids for row in queue):
        raise ValueError("Download queue references an unknown source")
    if len({row["local_file"] for row in queue}) != len(queue):
        raise ValueError("Download destinations must be unique")

    ready_decisions = {
        "DOWNLOAD_READY_QUANTITATIVE",
        "DOWNLOAD_READY_SEMIQUANTITATIVE",
        "DOWNLOAD_READY_BINARY",
    }
    ready_modifications = {row["modification"] for row in sources if row["decision"] in ready_decisions}
    required_ready = int(config["admission_rule"]["minimum_download_ready_modifications"])
    gate = "GO_TO_DOWNLOAD_AUDIT" if len(ready_modifications) >= required_ready else "NO_GO_INSUFFICIENT_SOURCES"
    checked = datetime.now(timezone.utc).isoformat()

    cards = [{**row, "checked_utc": checked, "local_status": "NOT_DOWNLOADED"} for row in sources]
    queue_rows = [{**row, "local_status": "NOT_DOWNLOADED"} for row in queue]
    out = config["outputs"]
    write_csv(ROOT / out["evidence_cards"], cards)
    write_csv(ROOT / out["source_manifest"], cards)
    write_csv(ROOT / out["download_queue"], queue_rows)
    write_csv(
        ROOT / out["status"],
        [{
            "stage": 23,
            "status": "PASS",
            "sources_scouted": len(sources),
            "download_ready_modifications": len(ready_modifications),
            "queued_files": len(queue),
            "phase3_gate": gate,
            "checked_utc": checked,
        }],
    )

    report = f"""# Stage 23: The third phase of multi-modification public data reconnaissance

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: deep-research/source-scout
- Origin Date: {datetime.now(timezone.utc).date().isoformat()}
- Verification Status: SOURCE_METADATA_VERIFIED
- Version Label: phase3_stage23_v1

At this stage, {len(sources)} candidate modifications were checked, forming a download queue of the first batch of {len(queue)} small processed files. The modification to meet download preparation conditions is {', '.join(sorted(ready_modifications))}, and the stage threshold is `{gate}`. This only authorizes access to post-download field and coordinate auditing, which does not mean that cross-modification modeling can already be carried out.

## First batch priority

1. **m5C/GSE140995**: HeLa bsRNA-seq, post-processing workbook containing candidate sites and non-conversion ratios, prioritized as continuous or semi-continuous track audit.
2. **m7G / GSE112276**: Both HeLa and HEK293T have HighFC files, but it is enrichment/fold change semantics and cannot be used as a stoichiometric ratio.
3. **Nm/GSE90164**: Revised compact BED site table for HeLa and HEK, treated by binary site orbitals.

## Suspension source

- m1A: GSE97909 proposes a cap structure cross-reactivity explanation for a broad range of antibody peaks, not suitable as a general mRNA m1A layer.
- Ψ: GSE255287 matches HEK293T, but the three processed matrices total about 12.5 GB, and the site extraction process needs to be designed first.
- ac4C: GSE162043 matches HeLa, but the mRNA ac4C detection remains methodologically controversial and must be reanalyzed using untreated versus NAT10-knockout controls.
- A-to-I: REDIportal v3 is primarily from GTEx/TCGA and does not strictly match existing HeLa/HEK293T culture conditions; engineered HEK293T edited data is more suitable for the Phase 4 perturbation route.

The next mandatory threshold is file-by-file verification after downloading.Reference version, 0/1 starting point, link direction, field semantics, repeat consistency and actual overlap rate with existing structure tables. Only after at least two modifications pass this threshold can the unified feature space be launched."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
