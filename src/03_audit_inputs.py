"""Reproducible streaming audit of input tables, archives, and reference assets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import logging
import math
import tarfile
from collections import Counter
from datetime import datetime
from pathlib import Path

from pipeline_common import MISSING, parse_icshape_line, read_glori, write_csv


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_processed"
META = ROOT / "metadata"
AUDIT = META / "audits" / "input"
LOG = ROOT / "results" / "logs" / "03_audit_inputs.log"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


class ColumnStats:
    def __init__(self) -> None:
        self.total = self.missing = self.numeric = 0
        self.minimum = math.inf
        self.maximum = -math.inf
        self.values = Counter()

    def add(self, value: str) -> None:
        self.total += 1
        if value in MISSING:
            self.missing += 1
            return
        if len(self.values) <= 30 or value in self.values:
            self.values[value] += 1
        try:
            number = float(value)
        except ValueError:
            return
        self.numeric += 1
        self.minimum = min(self.minimum, number)
        self.maximum = max(self.maximum, number)

    def row(self, dataset: str, column: str) -> dict:
        low_cardinality = self.values if sum(self.values.values()) == self.total - self.missing else {}
        return {
            "dataset": dataset,
            "column": column,
            "rows": self.total,
            "missing": self.missing,
            "missing_fraction": round(self.missing / self.total, 8) if self.total else "",
            "numeric": self.numeric,
            "minimum": "" if self.numeric == 0 else self.minimum,
            "maximum": "" if self.numeric == 0 else self.maximum,
            "value_counts_if_low_cardinality": json.dumps(low_cardinality, ensure_ascii=False, sort_keys=True),
        }


def audit_icshape_stream(dataset: str, handle) -> tuple[dict, list[dict]]:
    stats = {name: ColumnStats() for name in ["transcript_id", "length", "abundance", "reactivity"]}
    rows = bases = malformed = length_mismatch = valid = 0
    transcript_versions = Counter()
    for line in io.TextIOWrapper(handle, encoding="utf-8"):
        rows += 1
        try:
            tid, length, abundance, values = parse_icshape_line(line)
        except Exception:
            malformed += 1
            continue
        stats["transcript_id"].add(tid)
        stats["length"].add(str(length))
        stats["abundance"].add(abundance)
        transcript_versions["versioned" if "." in tid else "unversioned"] += 1
        if len(values) != length:
            length_mismatch += 1
        bases += len(values)
        for value in values:
            if math.isnan(value):
                stats["reactivity"].add("NULL")
            else:
                stats["reactivity"].add(str(value))
                valid += 1
    summary = {
        "dataset": dataset,
        "format": "headerless_icshape_vector",
        "rows": rows,
        "columns": "3+declared_length",
        "malformed_rows": malformed,
        "length_mismatch_rows": length_mismatch,
        "declared_bases": bases,
        "valid_reactivities": valid,
        "valid_fraction": round(valid / bases, 8) if bases else "",
        "chromosomes": "",
        "strands": "",
        "transcript_id_style": json.dumps(transcript_versions, sort_keys=True),
        "formula_failures": "",
    }
    return summary, [value.row(dataset, name) for name, value in stats.items()]


def audit_icshape_gzip(dataset: str, path: Path) -> tuple[dict, list[dict]]:
    with path.open("rb") as raw, gzip.GzipFile(fileobj=raw) as handle:
        return audit_icshape_stream(dataset, handle)


def audit_he_la_tar(path: Path) -> tuple[dict, list[dict]]:
    member_name = "GSM4333258_HeLa.out.txt.gz"
    with tarfile.open(path, "r") as archive:
        members = archive.getnames()
        member = archive.getmember(member_name)
        extracted = archive.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"Cannot read {member_name}")
        with gzip.GzipFile(fileobj=extracted) as nested:
            summary, columns = audit_icshape_stream("hela_icshape", nested)
    summary["archive_members"] = len(members)
    summary["selected_member"] = member_name
    return summary, columns


def audit_glori(dataset: str, path: Path) -> tuple[dict, list[dict]]:
    stats = {}
    rows = malformed = ratio_bad = normalized_bad = 0
    chromosomes = Counter()
    strands = Counter()
    headers = []
    for record in read_glori(path):
        rows += 1
        if not headers:
            headers = list(record)
            stats = {name: ColumnStats() for name in headers}
        if None in record or any(value is None for value in record.values()):
            malformed += 1
            continue
        for name, value in record.items():
            stats[name].add(value)
        chromosomes[record["Chr"]] += 1
        strands[record["Strand"]] += 1
        try:
            ratio = float(record["Ratio"])
            expected = round(int(record["Acov"]) / int(record["AGcov"]), 5)
            if abs(ratio - expected) > 5e-6:
                ratio_bad += 1
            if abs(float(record["NormeRatio"]) - ratio * (1 - float(record["NonCR"]))) > 1e-10:
                normalized_bad += 1
        except (ValueError, ZeroDivisionError):
            ratio_bad += 1
            normalized_bad += 1
    summary = {
        "dataset": dataset,
        "format": "gzip_tsv_with_header",
        "rows": rows,
        "columns": len(headers),
        "malformed_rows": malformed,
        "length_mismatch_rows": "",
        "declared_bases": "",
        "valid_reactivities": "",
        "valid_fraction": "",
        "chromosomes": json.dumps(chromosomes, sort_keys=True),
        "strands": json.dumps(strands, sort_keys=True),
        "transcript_id_style": "",
        "formula_failures": json.dumps({"Ratio": ratio_bad, "NormeRatio": normalized_bad}),
    }
    return summary, [value.row(dataset, name) for name, value in stats.items()]


def audit_manifest(manifest_path: Path, base: Path, kind: str) -> tuple[list[dict], list[dict]]:
    file_rows = []
    reference_rows = []
    with manifest_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            relative = row.get("relative_path") or f"data/raw_processed/{row['candidate_file']}"
            path = base / relative
            expected_hash = row["sha256"].lower()
            exists = path.is_file()
            actual_hash = sha256(path) if exists else ""
            expected_bytes = row.get("bytes", "")
            actual_bytes = path.stat().st_size if exists else ""
            size_ok = "" if not expected_bytes else int(expected_bytes) == actual_bytes
            integrity = "NOT_TESTED"
            error = ""
            if exists:
                try:
                    if path.suffix == ".gz":
                        with gzip.open(path, "rb") as stream:
                            for _ in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                                pass
                        integrity = "GZIP_FULL_STREAM_OK"
                    elif path.suffix == ".tar":
                        with tarfile.open(path, "r") as archive:
                            for member in archive.getmembers():
                                if member.isfile():
                                    extracted = archive.extractfile(member)
                                    if extracted:
                                        for _ in iter(lambda: extracted.read(8 * 1024 * 1024), b""):
                                            pass
                        integrity = "TAR_FULL_STREAM_OK"
                    else:
                        with path.open("rb") as stream:
                            for _ in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                                pass
                        integrity = "FULL_STREAM_OK"
                except Exception as exc:
                    integrity = "FAILED"
                    error = repr(exc)
            output = {
                "asset_id": row.get("asset_id") or row.get("record_id"),
                "kind": kind,
                "path": relative,
                "exists": exists,
                "expected_bytes": expected_bytes,
                "actual_bytes": actual_bytes,
                "size_ok": size_ok,
                "expected_sha256": expected_hash,
                "actual_sha256": actual_hash,
                "sha256_ok": exists and actual_hash == expected_hash,
                "container_integrity": integrity,
                "error": error,
            }
            file_rows.append(output)
            if kind == "reference":
                reference_rows.append({**output, "release_or_version": row["release_or_version"], "assembly": row["assembly"], "source_url": row["source_url"]})
            logging.info("checked %s: exists=%s sha=%s integrity=%s", relative, exists, output["sha256_ok"], integrity)
    return file_rows, reference_rows


def main() -> None:
    argparse.ArgumentParser().parse_args()
    setup_logging()
    started = datetime.now().isoformat(timespec="seconds")
    AUDIT.mkdir(parents=True, exist_ok=True)
    logging.info("audit started at %s", started)

    source_files, _ = audit_manifest(META / "manifests" / "source" / "source_manifest.csv", ROOT, "source")
    reference_files, reference_rows = audit_manifest(META / "manifests" / "source" / "02_reference_manifest.csv", ROOT, "reference")

    summaries = []
    column_rows = []
    datasets = [
        ("hek293t_icshape", RAW / "GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz", "icshape"),
        ("hek293t_glori_rep1", RAW / "GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv.gz", "glori"),
        ("hek293t_glori_rep2", RAW / "GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv.gz", "glori"),
        ("hela_glori_rep1", RAW / "GSM6432595_Hela-1_35bp_m2.totalm6A.FDR.csv.gz", "glori"),
        ("hela_glori_rep2", RAW / "GSM6432596_Hela-2_35bp_m2.totalm6A.FDR.csv.gz", "glori"),
    ]
    for name, path, kind in datasets:
        logging.info("parsing %s", name)
        summary, columns = audit_icshape_gzip(name, path) if kind == "icshape" else audit_glori(name, path)
        summaries.append(summary)
        column_rows.extend(columns)
    logging.info("parsing hela_icshape from TAR")
    summary, columns = audit_he_la_tar(RAW / "GSE145805_RAW.tar")
    summaries.append(summary)
    column_rows.extend(columns)

    write_csv(AUDIT / "03_file_integrity_audit.csv", source_files + reference_files)
    write_csv(AUDIT / "03_reference_integrity_audit.csv", reference_rows)
    write_csv(AUDIT / "03_dataset_structure_audit.csv", summaries)
    write_csv(AUDIT / "03_column_statistics_audit.csv", column_rows)
    failures = [r for r in source_files + reference_files if not (r["exists"] and r["sha256_ok"] and r["container_integrity"].endswith("OK"))]
    parse_failures = sum(int(r["malformed_rows"]) for r in summaries)
    formula_failures = sum(sum(json.loads(r["formula_failures"]).values()) for r in summaries if r["formula_failures"])
    status = "PASS" if not failures and parse_failures == 0 and formula_failures == 0 else "FAIL"
    write_csv(AUDIT / "03_audit_run_status.csv", [{"status": status, "started_at": started, "finished_at": datetime.now().isoformat(timespec="seconds"), "file_failures": len(failures), "malformed_rows": parse_failures, "formula_failures": formula_failures}])
    logging.info("audit finished: %s", status)
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
