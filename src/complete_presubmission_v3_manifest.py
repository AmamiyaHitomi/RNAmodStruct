"""Extend the preliminary release index to every analysis script and config."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "metadata" / "manifests" / "release" / "presubmission_v3_release_manifest.csv"
BASE_META = ROOT / "metadata" / "manifests" / "release" / "presubmission_v3_release_metadata.json"
OUTPUT = ROOT / "metadata" / "manifests" / "release" / "presubmission_v3_release_manifest_full.csv"
META = ROOT / "metadata" / "manifests" / "release" / "presubmission_v3_release_metadata_full.json"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    if OUTPUT.exists() or META.exists():
        raise FileExistsError("Refusing to overwrite complete release manifest")
    with BASE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    index = {row["path"]: row for row in rows}
    if len(index) != len(rows):
        raise AssertionError("Preliminary manifest has duplicate paths")
    for row in rows:
        path = ROOT / row["path"]
        if path.stat().st_size != int(row["bytes"]) or digest(path) != row["sha256"]:
            raise AssertionError(f"Preliminary hash no longer valid: {row['path']}")

    additions = []
    for pattern, category in (("src/*.py", "analysis_code_complete"),
                              ("config/*.yaml", "analysis_config_complete"),
                              ("docs/reports/final/main_v2_phase*_zh.md", "stage_report"),
                              ("docs/reports/final/main_v2_*roadmap_zh.md", "revision_plan"),
                              ("docs/reports/final/main_v2_*review_zh.md", "revision_evidence")):
        additions.extend((category, p) for p in sorted(ROOT.glob(pattern)) if p.is_file())
    additions.extend(("manifest_provenance", p) for p in (BASE, BASE_META))
    for category, path in additions:
        rel = path.relative_to(ROOT).as_posix()
        if rel not in index:
            index[rel] = {"category": category, "path": rel, "bytes": path.stat().st_size,
                          "sha256": digest(path)}
    for required in ("src/09b_run_hela_association.py", "src/19_run_invivo_invitro_pairing.py",
                     "src/complete_presubmission_v3_manifest.py", "config/08_prediction_modeling.yaml",
                     "docs/manuscript_v3/main_v3.pdf", "docs/manuscript_v3/supplement_v3.pdf"):
        if required not in index:
            raise AssertionError(f"Missing required release asset: {required}")
    with OUTPUT.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["category", "path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(index.values())
    base_meta = json.loads(BASE_META.read_text(encoding="utf-8"))
    metadata = {**base_meta, "entry_count": len(index), "manifest_sha256": digest(OUTPUT),
                "base_manifest_sha256": digest(BASE),
                "authority": "full_release_manifest_including_every_src_python_and_config_yaml"}
    META.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "entries": len(index),
                      "sha256": metadata["manifest_sha256"]}))


if __name__ == "__main__":
    main()
