"""Shared parsers and coordinate transforms for the RNAmodStruct audit pipeline."""

from __future__ import annotations

import bisect
import csv
import gzip
import io
import math
import re
import tarfile
from array import array
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


MISSING = {"", "NULL", "NA", "N/A", "NaN", "nan", ".", "-999"}


def norm_chr(value: str) -> str:
    value = value.strip()
    if value.startswith("chr"):
        return value
    if value == "MT":
        return "chrM"
    return f"chr{value}"


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+) "([^"]*)";', text)}


@dataclass(frozen=True)
class Exon:
    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str
    gene_id: str
    gene_name: str
    exon_number: int
    tx_start: int = 0
    tx_end: int = 0


def load_gtf_exons(
    path: Path,
    allowed_transcripts: set[str] | None = None,
    versioned: bool = False,
) -> dict[str, list[Exon]]:
    grouped: dict[str, list[Exon]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attributes(fields[8])
            tid = attrs.get("transcript_id")
            if versioned and attrs.get("transcript_version"):
                tid = f"{tid}.{attrs['transcript_version']}"
            if not tid or (allowed_transcripts is not None and tid not in allowed_transcripts):
                continue
            grouped[tid].append(
                Exon(
                    chrom=norm_chr(fields[0]),
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    transcript_id=tid,
                    gene_id=attrs.get("gene_id", ""),
                    gene_name=attrs.get("gene_name", ""),
                    exon_number=int(attrs.get("exon_number", "0") or 0),
                )
            )

    result: dict[str, list[Exon]] = {}
    for tid, exons in grouped.items():
        ordered = sorted(exons, key=lambda e: (e.start, e.end), reverse=exons[0].strand == "-")
        offset = 1
        materialized = []
        for exon in ordered:
            length = exon.end - exon.start + 1
            materialized.append(
                Exon(**{**exon.__dict__, "tx_start": offset, "tx_end": offset + length - 1})
            )
            offset += length
        result[tid] = materialized
    return result


def tx_to_genome(exons: list[Exon], tx_pos: int) -> tuple[str, int, str, int] | None:
    for exon in exons:
        if exon.tx_start <= tx_pos <= exon.tx_end:
            delta = tx_pos - exon.tx_start
            genome_pos = exon.start + delta if exon.strand == "+" else exon.end - delta
            return exon.chrom, genome_pos, exon.strand, exon.exon_number
    return None


def genome_to_tx(exon: Exon, genome_pos: int) -> int:
    delta = genome_pos - exon.start if exon.strand == "+" else exon.end - genome_pos
    return exon.tx_start + delta


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        name = None
        chunks: list[str] = []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def load_selected_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    return {name: seq for name, seq in iter_fasta(path) if name in wanted}


def fetch_fasta_bases(path: Path, queries: dict[str, set[int]]) -> dict[tuple[str, int], str]:
    """Fetch a small/medium set of 1-based bases with one streaming FASTA pass."""
    normalized = {norm_chr(chrom): positions for chrom, positions in queries.items()}
    result: dict[tuple[str, int], str] = {}
    for name, seq in iter_fasta(path):
        chrom = norm_chr(name)
        for pos in normalized.get(chrom, ()):
            if 1 <= pos <= len(seq):
                result[(chrom, pos)] = seq[pos - 1]
    return result


def parse_icshape_line(line: str) -> tuple[str, int, str, array]:
    fields = line.rstrip("\n").split("\t")
    tid, length, abundance = fields[:3]
    values = array("f", (math.nan if x in MISSING else float(x) for x in fields[3:]))
    return tid, int(length), abundance, values


def load_icshape(path: Path) -> dict[str, tuple[int, str, array]]:
    result = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            tid, length, abundance, values = parse_icshape_line(line)
            result[tid] = (length, abundance, values)
    return result


def load_icshape_tar(path: Path, member: str) -> dict[str, tuple[int, str, array]]:
    """Load a headerless icSHAPE reactivity table nested inside a gzip TAR member."""
    result = {}
    with tarfile.open(path, "r") as archive:
        extracted = archive.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"Cannot read TAR member {member}")
        with gzip.GzipFile(fileobj=extracted) as nested, io.TextIOWrapper(nested, encoding="utf-8") as handle:
            for line in handle:
                tid, length, abundance, values = parse_icshape_line(line)
                result[tid] = (length, abundance, values)
    return result


@dataclass(frozen=True)
class ChainBlock:
    source_chrom: str
    source_start: int
    source_end: int
    target_chrom: str
    target_start: int
    target_end: int
    target_size: int
    target_strand: str
    score: int
    chain_id: int


class _BlockIndex:
    def __init__(self, blocks: Iterable[ChainBlock], by_source: bool):
        self.by_source = by_source
        self.blocks: dict[str, list[ChainBlock]] = defaultdict(list)
        for block in blocks:
            chrom = block.source_chrom if by_source else block.target_chrom
            self.blocks[chrom].append(block)
        self.starts: dict[str, list[int]] = {}
        self.max_ends: dict[str, list[int]] = {}
        for chrom, values in self.blocks.items():
            values.sort(key=lambda b: (b.source_start if by_source else b.target_start, -b.score))
            starts = [b.source_start if by_source else b.target_start for b in values]
            ends = [b.source_end if by_source else b.target_end for b in values]
            running = 0
            maxima = []
            for end in ends:
                running = max(running, end)
                maxima.append(running)
            self.starts[chrom] = starts
            self.max_ends[chrom] = maxima

    def candidates(self, chrom: str, pos0: int) -> list[ChainBlock]:
        chrom = norm_chr(chrom)
        values = self.blocks.get(chrom, [])
        if not values:
            return []
        starts = self.starts[chrom]
        maxima = self.max_ends[chrom]
        index = bisect.bisect_right(starts, pos0) - 1
        found = []
        while index >= 0 and maxima[index] > pos0:
            block = values[index]
            end = block.source_end if self.by_source else block.target_end
            start = block.source_start if self.by_source else block.target_start
            if start <= pos0 < end:
                found.append(block)
            index -= 1
        return found


class ChainLiftOver:
    """Bidirectional single-base lookup for a UCSC old-to-new chain file."""

    def __init__(self, path: Path):
        blocks = list(self._read_blocks(path))
        self.forward_index = _BlockIndex(blocks, by_source=True)
        self.reverse_index = _BlockIndex(blocks, by_source=False)

    @staticmethod
    def _read_blocks(path: Path) -> Iterator[ChainBlock]:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            header = None
            source_pos = target_pos = 0
            for raw in handle:
                line = raw.strip()
                if not line:
                    header = None
                    continue
                if line.startswith("chain "):
                    parts = line.split()
                    header = {
                        "score": int(parts[1]),
                        "source_chrom": norm_chr(parts[2]),
                        "source_pos": int(parts[5]),
                        "target_chrom": norm_chr(parts[7]),
                        "target_size": int(parts[8]),
                        "target_strand": parts[9],
                        "target_pos": int(parts[10]),
                        "chain_id": int(parts[12]),
                    }
                    source_pos = header["source_pos"]
                    target_pos = header["target_pos"]
                    continue
                if header is None:
                    continue
                parts = [int(x) for x in line.split()]
                size = parts[0]
                if header["target_strand"] == "+":
                    target_start = target_pos
                else:
                    target_start = header["target_size"] - (target_pos + size)
                yield ChainBlock(
                    source_chrom=header["source_chrom"],
                    source_start=source_pos,
                    source_end=source_pos + size,
                    target_chrom=header["target_chrom"],
                    target_start=target_start,
                    target_end=target_start + size,
                    target_size=header["target_size"],
                    target_strand=header["target_strand"],
                    score=header["score"],
                    chain_id=header["chain_id"],
                )
                if len(parts) == 3:
                    source_pos += size + parts[1]
                    target_pos += size + parts[2]

    @staticmethod
    def _unique(mappings: list[tuple[str, int, str, int]]) -> list[tuple[str, int, str, int]]:
        unique = {}
        for item in mappings:
            unique[item[:3]] = item
        return sorted(unique.values(), key=lambda x: (-x[3], x[0], x[1], x[2]))

    def forward(self, chrom: str, pos1: int, strand: str) -> list[tuple[str, int, str, int]]:
        pos0 = pos1 - 1
        mappings = []
        for block in self.forward_index.candidates(chrom, pos0):
            offset = pos0 - block.source_start
            if block.target_strand == "+":
                target0 = block.target_start + offset
                out_strand = strand
            else:
                target0 = block.target_end - 1 - offset
                out_strand = "+" if strand == "-" else "-"
            mappings.append((block.target_chrom, target0 + 1, out_strand, block.score))
        return self._unique(mappings)

    def reverse(self, chrom: str, pos1: int, strand: str) -> list[tuple[str, int, str, int]]:
        pos0 = pos1 - 1
        mappings = []
        for block in self.reverse_index.candidates(chrom, pos0):
            if block.target_strand == "+":
                offset = pos0 - block.target_start
                out_strand = strand
            else:
                offset = block.target_end - 1 - pos0
                out_strand = "+" if strand == "-" else "-"
            mappings.append((block.source_chrom, block.source_start + offset + 1, out_strand, block.score))
        return self._unique(mappings)


class ExonBinIndex:
    def __init__(self, transcripts: dict[str, list[Exon]], bin_size: int = 16384):
        self.bin_size = bin_size
        self.bins: dict[tuple[str, str, int], list[Exon]] = defaultdict(list)
        for exons in transcripts.values():
            for exon in exons:
                for bin_id in range((exon.start - 1) // bin_size, (exon.end - 1) // bin_size + 1):
                    self.bins[(exon.chrom, exon.strand, bin_id)].append(exon)

    def query(self, chrom: str, pos1: int, strand: str) -> list[Exon]:
        key = (norm_chr(chrom), strand, (pos1 - 1) // self.bin_size)
        return [exon for exon in self.bins.get(key, ()) if exon.start <= pos1 <= exon.end]


def read_glori(path: Path) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def coverage(values: array, center1: int, flank: int) -> tuple[float, float, float]:
    left = values[max(0, center1 - flank - 1) : center1 - 1]
    right = values[center1 : min(len(values), center1 + flank)]
    left_cov = sum(not math.isnan(x) for x in left) / flank if len(left) == flank else math.nan
    right_cov = sum(not math.isnan(x) for x in right) / flank if len(right) == flank else math.nan
    combined = (
        sum(not math.isnan(x) for x in left) + sum(not math.isnan(x) for x in right)
    ) / (2 * flank)
    if len(left) != flank or len(right) != flank:
        combined = math.nan
    return left_cov, right_cov, combined


def write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str] | None = None) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    iterator = iter(rows)
    try:
        first = next(iterator)
    except StopIteration:
        if fieldnames:
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.DictWriter(handle, fieldnames=fieldnames).writeheader()
        return 0
    names = fieldnames or list(first)
    count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerow(first)
        count += 1
        for row in iterator:
            writer.writerow(row)
            count += 1
    return count
