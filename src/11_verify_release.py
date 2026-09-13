"""Create or verify the deterministic RNAmodStruct release-file manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata" / "manifests" / "release" / "release_manifest.csv"
SEARCH_ROOTS = [
    "config",
    "src",
    "tests",
    "docs/reports/final",
    "results",
    "data/final",
    "data/interim",
]
TOP_LEVEL = ["README.md", "environment.yml", "pytest.ini", ".gitignore", ".gitmodules"]
EXCLUDED_PARTS = {"__pycache__", ".pytest_cache", "logs", "qa"}
EXCLUDED_FILES = {MANIFEST.resolve()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def release_files() -> list[Path]:
    files = [ROOT / name for name in TOP_LEVEL]
    for name in SEARCH_ROOTS:
        base = ROOT / name
        if base.exists():
            files.extend(path for path in base.rglob("*") if path.is_file())
    return sorted(
        {
            path.resolve()
            for path in files
            if path.exists()
            and path.resolve() not in EXCLUDED_FILES
            and not EXCLUDED_PARTS.intersection(path.relative_to(ROOT).parts)
            and path.suffix.lower() not in {".pyc", ".log"}
        },
        key=lambda path: path.relative_to(ROOT).as_posix(),
    )


def current_rows() -> list[dict[str, str | int]]:
    return [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in release_files()
    ]


def write_manifest() -> None:
    rows = current_rows()
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"WROTE {MANIFEST.relative_to(ROOT).as_posix()} ({len(rows)} files)")


def verify_manifest() -> None:
    if not MANIFEST.exists():
        raise SystemExit("FAIL: metadata/manifests/release/release_manifest.csv is missing")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        expected = {row["path"]: row for row in csv.DictReader(handle)}
    actual = {str(row["path"]): row for row in current_rows()}
    failures: list[str] = []
    for path in sorted(expected.keys() - actual.keys()):
        failures.append(f"missing: {path}")
    for path in sorted(actual.keys() - expected.keys()):
        failures.append(f"unmanifested: {path}")
    for path in sorted(expected.keys() & actual.keys()):
        if expected[path]["sha256"] != actual[path]["sha256"]:
            failures.append(f"sha256 mismatch: {path}")
        if int(expected[path]["bytes"]) != int(actual[path]["bytes"]):
            failures.append(f"size mismatch: {path}")
    if failures:
        raise SystemExit("FAIL:\n" + "\n".join(failures))
    print(f"PASS: {len(actual)} release files match metadata/manifests/release/release_manifest.csv")


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-manifest", action="store_true")
    mode.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write_manifest:
        write_manifest()
    else:
        verify_manifest()


if __name__ == "__main__":
    main()
