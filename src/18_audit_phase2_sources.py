"""Stage 18: audit and admit public datasets for phase 2."""

from __future__ import annotations

import csv
import gzip
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "18_public_data_audit.yaml"
TABLE = ROOT / "results" / "tables" / "18_public_data_audit.csv"
MANIFEST = ROOT / "metadata" / "manifests" / "source" / "phase2_source_manifest.csv"
STATUS = ROOT / "results" / "status" / "18_public_data_audit_run_status.csv"
REPORT = ROOT / "results" / "reports" / "18_public_data_audit_report.md"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_profile(path: Path) -> tuple[int, int, list[str]]:
    rows = 0
    columns = None
    first_fields: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if columns is None:
                columns = len(fields)
                first_fields = fields[:7]
            rows += 1
    return rows, int(columns or 0), first_fields


def inspect_source(record_id: str, item: dict) -> dict:
    relative = item.get("file")
    base = {
        "record_id": record_id,
        "accession": item["accession"],
        "cell_line": item["cell_line"],
        "technology": item["technology"],
        "reference": item["reference"],
        "coordinate_contract": item["coordinate_contract"],
        "value_contract": item["value_contract"],
        "source_url": item["source_url"],
        "download_url": item.get("download_url", ""),
        "local_file": relative or "",
        "admission": item["admission"],
        "checked_at": datetime.now(timezone.utc).isoformat(),
    }
    if not relative:
        return {**base, "local_status": "NOT_REQUESTED", "bytes": "", "sha256": "", "rows": "", "columns": "", "integrity": "NOT_APPLICABLE", "coordinate_checks": "NOT_RUN"}
    path = ROOT / relative
    if not path.is_file():
        return {**base, "local_status": "MISSING", "bytes": "", "sha256": "", "rows": "", "columns": "", "integrity": "FAIL", "coordinate_checks": "NOT_RUN"}
    try:
        rows, columns, first = gzip_profile(path)
    except Exception as exc:  # pragma: no cover - exercised by integration audit
        return {**base, "local_status": "PRESENT", "bytes": path.stat().st_size, "sha256": sha256_file(path), "rows": "", "columns": "", "integrity": f"FAIL:{type(exc).__name__}", "coordinate_checks": "NOT_RUN"}
    coordinate_checks = "TRANSCRIPT_TABLE_SHAPE_OK"
    if path.name.endswith(".bed.gz"):
        bed = pd.read_csv(path, sep="\t", header=None)
        width_ok = bool((pd.to_numeric(bed[2]) - pd.to_numeric(bed[1])).eq(1).all())
        strand_ok = set(bed[5].astype(str)).issubset({"+", "-"})
        unique_ok = not bed.duplicated([0, 1, 2, 5]).any()
        coordinate_checks = f"width1={width_ok};strand={strand_ok};unique={unique_ok}"
        if not (width_ok and strand_ok and unique_ok and columns == 7):
            return {**base, "local_status": "PRESENT", "bytes": path.stat().st_size, "sha256": sha256_file(path), "rows": rows, "columns": columns, "integrity": "FAIL", "coordinate_checks": coordinate_checks}
    elif columns <= 3 or len(first) <= 3:
        return {**base, "local_status": "PRESENT", "bytes": path.stat().st_size, "sha256": sha256_file(path), "rows": rows, "columns": columns, "integrity": "FAIL", "coordinate_checks": coordinate_checks}
    return {**base, "local_status": "DOWNLOADED_VERIFIED", "bytes": path.stat().st_size, "sha256": sha256_file(path), "rows": rows, "columns": columns, "integrity": "GZIP_FULL_STREAM_OK", "coordinate_checks": coordinate_checks}


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    rows = [inspect_source(name, item) for name, item in config["sources"].items()]
    write_csv(TABLE, rows)
    write_csv(MANIFEST, rows)
    required = [row for row in rows if row["local_file"]]
    passed = all(row["integrity"] == "GZIP_FULL_STREAM_OK" for row in required)
    ready = sum(row["admission"] == "PAIRED_ANALYSIS_READY" and row["integrity"] == "GZIP_FULL_STREAM_OK" for row in rows)
    status = "PASS" if passed and ready >= 1 else "FAIL"
    write_csv(STATUS, [{"stage": 18, "status": status, "files_verified": sum(row["integrity"] == "GZIP_FULL_STREAM_OK" for row in rows), "paired_datasets_ready": ready, "config": str(CONFIG.relative_to(ROOT)).replace("\\", "/")}])
    REPORT.write_text(f"""# Stage 18: The second phase of public data and coordinate audit

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {datetime.now(timezone.utc).date().isoformat()}
- Verification Status: {'VERIFIED' if status == 'PASS' else 'UNVERIFIED'}
- Version Label: phase2_stage18_v1

This stage performs local container, row and column number, and coordinate contract audits on the second phase candidate sources. {sum(row['integrity'] == 'GZIP_FULL_STREAM_OK' for row in rows)} local files have been fully verified; the number of data sets that can directly enter strict pair analysis is {ready}.

- GSE74353 HEK293T in-vitro icSHAPE: admission stage 19, paired with existing in-vivo files in the same GRCh37.74 transcript coordinate system.
- GSE162356 HeLa SAC-seq: Native seven-column BED with coordinate base contract passed; final semantics for numeric columns continue to retain source fields until stage 20 and are not renamed to stoichiometry without confirmation.
- GSE162356 HEK293: Only close-match sensitivity analysis of HEK293T and cannot be marked as strict cell line replication.
- Human RNA MaP with GSE50676: currently missing quantitative m6A matching cell conditions, remains in audit status, and does not splice raw scores.

Machine table: `results/tables/18_public_data_audit.csv`; independent source manifest: `metadata/manifests/source/phase2_source_manifest.csv`.""", encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
