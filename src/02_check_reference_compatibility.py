"""Compare icSHAPE transcript identifiers and lengths with Ensembl references."""

from __future__ import annotations

import argparse
import gzip
import re
import tarfile
from collections import defaultdict
from pathlib import Path


TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
TRANSCRIPT_VERSION = re.compile(r'transcript_version "([^"]+)"')


def read_icshape(path: Path, member: str | None) -> dict[str, int]:
    if member:
        with tarfile.open(path) as archive:
            raw = archive.extractfile(member)
            if raw is None:
                raise FileNotFoundError(member)
            with gzip.open(raw, "rt") as handle:
                return {fields[0]: int(fields[1]) for fields in map(lambda x: x.split("\t", 3), handle)}
    with gzip.open(path, "rt") as handle:
        return {fields[0]: int(fields[1]) for fields in map(lambda x: x.split("\t", 3), handle)}


def gtf_lengths(path: Path, targets: dict[str, int]) -> tuple[dict[str, int], dict[str, int]]:
    lengths: dict[str, int] = defaultdict(int)
    base_lengths: dict[str, int] = defaultdict(int)
    target_bases = {key.split(".", 1)[0] for key in targets}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            transcript = TRANSCRIPT_ID.search(fields[8])
            version = TRANSCRIPT_VERSION.search(fields[8])
            if transcript is None:
                continue
            key = transcript.group(1) + (f".{version.group(1)}" if version else "")
            if key in targets:
                lengths[key] += int(fields[4]) - int(fields[3]) + 1
            if transcript.group(1) in target_bases:
                base_lengths[transcript.group(1)] += int(fields[4]) - int(fields[3]) + 1
    return lengths, base_lengths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("icshape", type=Path)
    parser.add_argument("gtf", nargs="+", type=Path)
    parser.add_argument("--member")
    args = parser.parse_args()
    targets = read_icshape(args.icshape, args.member)
    print("gtf\ttargets\tfound\tlength_exact\tlength_mismatch\tmissing\tbase_found\tbase_length_exact")
    for path in args.gtf:
        found, base_found = gtf_lengths(path, targets)
        exact = sum(length == targets[key] for key, length in found.items())
        base_targets = {key.split(".", 1)[0]: length for key, length in targets.items()}
        base_exact = sum(length == base_targets[key] for key, length in base_found.items())
        print(
            f"{path}\t{len(targets)}\t{len(found)}\t{exact}\t{len(found)-exact}\t"
            f"{len(targets)-len(found)}\t{len(base_found)}\t{base_exact}"
        )


if __name__ == "__main__":
    main()
