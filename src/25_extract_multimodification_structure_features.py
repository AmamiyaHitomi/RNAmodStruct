"""Stage 25: extract technology-stratified icSHAPE features for phase 3."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import math
import statistics
import sys
import tarfile
from array import array
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pipeline_common import ExonBinIndex, coverage, genome_to_tx, load_gtf_exons, parse_icshape_line  # noqa: E402

CONFIG = ROOT / "config" / "25_extract_multimodification_structure_features.yaml"


def setup_logging(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(path, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_structure_lengths_gzip(path: Path) -> dict[str, int]:
    result = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            tid, length, *_ = line.split("\t", 2)
            result[tid] = int(length)
    return result


def read_structure_lengths_tar(path: Path, member: str) -> dict[str, int]:
    result = {}
    with tarfile.open(path, "r") as archive:
        extracted = archive.extractfile(member)
        if extracted is None:
            raise ValueError(f"Missing TAR member: {member}")
        with gzip.GzipFile(fileobj=extracted) as nested:
            for raw in nested:
                fields = raw.split(b"\t", 2)
                result[fields[0].decode()] = int(fields[1])
    return result


def load_selected_gzip(path: Path, wanted: set[str]) -> dict[str, tuple[int, str, array]]:
    result = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            tid = line.split("\t", 1)[0]
            if tid in wanted:
                parsed = parse_icshape_line(line)
                result[tid] = parsed[1:]
    return result


def load_selected_tar(path: Path, member: str, wanted: set[str]) -> dict[str, tuple[int, str, array]]:
    result = {}
    with tarfile.open(path, "r") as archive:
        extracted = archive.extractfile(member)
        if extracted is None:
            raise ValueError(f"Missing TAR member: {member}")
        with gzip.GzipFile(fileobj=extracted) as nested:
            for raw in nested:
                line = raw.decode("utf-8")
                tid = line.split("\t", 1)[0]
                if tid in wanted:
                    parsed = parse_icshape_line(line)
                    result[tid] = parsed[1:]
    return result


def exact_length_transcripts(gtf: Path, lengths: dict[str, int], versioned: bool) -> dict:
    loaded = load_gtf_exons(gtf, set(lengths), versioned=versioned)
    return {
        tid: exons for tid, exons in loaded.items()
        if sum(exon.end - exon.start + 1 for exon in exons) == lengths[tid]
    }


def interval_hits(index: ExonBinIndex, chrom: str, start1: int, end1: int, strand: str) -> dict[str, list]:
    found: dict[str, dict[tuple, object]] = defaultdict(dict)
    for bin_id in range((start1 - 1) // index.bin_size, (end1 - 1) // index.bin_size + 1):
        for exon in index.bins.get((chrom, strand, bin_id), ()):
            if exon.start <= end1 and exon.end >= start1:
                found[exon.transcript_id][(exon.exon_number, exon.start, exon.end)] = exon
    return {tid: list(exons.values()) for tid, exons in found.items()}


def valid_mean(values) -> float | str:
    kept = [float(value) for value in values if not math.isnan(value)]
    return round(sum(kept) / len(kept), 8) if kept else ""


def window_json(values: array, center1: int, flank: int) -> tuple[str, str, int]:
    scores, mask = [], []
    for pos1 in range(center1 - flank, center1 + flank + 1):
        if pos1 < 1 or pos1 > len(values) or math.isnan(values[pos1 - 1]):
            scores.append(None)
            mask.append(0)
        else:
            scores.append(round(float(values[pos1 - 1]), 6))
            mask.append(1)
    return json.dumps(scores, separators=(",", ":")), json.dumps(mask, separators=(",", ":")), sum(mask)


def point_hit_tids(rows: list[dict], index: ExonBinIndex, offsets: list[int]) -> set[str]:
    tids = set()
    for row in rows:
        for offset in offsets:
            for exon in index.query(row["hg38_chr"], int(row["hg38_pos_1based"]) + offset, row["hg38_strand"]):
                tids.add(exon.transcript_id)
    return tids


def peak_hit_tids(rows: list[dict], index: ExonBinIndex) -> set[str]:
    tids = set()
    for row in rows:
        hits = interval_hits(index, row["hg38_chr"], int(row["hg38_start_0based"]) + 1, int(row["hg38_end_0based_exclusive"]), row["hg38_strand"])
        tids.update(hits)
    return tids


def point_features(row: dict, exon, structures: dict, offset: int, flank: int, display_flank: int, min_cov: float) -> dict:
    tid = exon.transcript_id
    length, abundance, values = structures[tid]
    pos1 = int(row["hg38_pos_1based"]) + offset
    tx_pos = genome_to_tx(exon, pos1)
    up_cov, down_cov, flank_cov = coverage(values, tx_pos, flank)
    left = values[max(0, tx_pos - flank - 1):tx_pos - 1]
    right = values[tx_pos:min(len(values), tx_pos + flank)]
    win_json, mask_json, display_valid = window_json(values, tx_pos, display_flank)
    main_valid = (
        not math.isnan(up_cov) and not math.isnan(down_cov)
        and up_cov >= min_cov and down_cov >= min_cov
    )
    center = values[tx_pos - 1]
    splice_edges = [edge for item in structures_transcripts[tid][:-1] for edge in (item.tx_end, item.tx_end + 1)]
    return {
        **row, "coordinate_offset_nt": offset, "mapped_hg38_pos_1based": pos1,
        "transcript_id": tid, "transcript_pos_1based": tx_pos,
        "transcript_gene_id": exon.gene_id, "transcript_gene_name": exon.gene_name,
        "transcript_length": length, "transcript_position_fraction": round(tx_pos / length, 8),
        "icshape_abundance_rpkm": "" if abundance == "*" else abundance,
        "icshape_center": "" if math.isnan(center) else round(float(center), 6),
        "coverage_up10": up_cov, "coverage_down10": down_cov, "coverage_flank10": flank_cov,
        "valid_count_up10": sum(not math.isnan(value) for value in left),
        "valid_count_down10": sum(not math.isnan(value) for value in right),
        "reactivity_mean_up10": valid_mean(left), "reactivity_mean_down10": valid_mean(right),
        "reactivity_mean_flank10": valid_mean(list(left) + list(right)),
        "main_window_valid": main_valid,
        "reactivity_window_m50_p50_json": win_json, "reactivity_mask_m50_p50_json": mask_json,
        "valid_reactivity_count_m50_p50": display_valid,
        "distance_to_nearest_splice_edge_tx": min((abs(tx_pos - edge) for edge in splice_edges), default=""),
    }


def extract_points(rows: list[dict], index: ExonBinIndex, structures: dict, offset: int, flank: int, display_flank: int, min_cov: float) -> list[dict]:
    selected = []
    for row in rows:
        pos1 = int(row["hg38_pos_1based"]) + offset
        exons = index.query(row["hg38_chr"], pos1, row["hg38_strand"])
        candidates = []
        seen = set()
        for exon in exons:
            if exon.transcript_id in seen or exon.transcript_id not in structures:
                continue
            seen.add(exon.transcript_id)
            candidates.append(point_features(row, exon, structures, offset, flank, display_flank, min_cov))
        if not candidates:
            continue
        ranked = sorted(candidates, key=lambda item: (not item["main_window_valid"], -item["valid_reactivity_count_m50_p50"], item["transcript_id"]))
        chosen = ranked[0]
        chosen["isoform_count"] = len(ranked)
        chosen["representative_isoform"] = True
        chosen["main_analysis_included"] = bool(chosen["main_window_valid"])
        chosen["main_exclusion_reason"] = "" if chosen["main_window_valid"] else "insufficient_icshape_flank10_coverage"
        selected.append(chosen)
    return selected


def extract_peaks(rows: list[dict], index: ExonBinIndex, structures: dict, flank: int, min_side: float, min_peak: float) -> list[dict]:
    selected = []
    for row in rows:
        start1, end1 = int(row["hg38_start_0based"]) + 1, int(row["hg38_end_0based_exclusive"])
        candidates = []
        for tid, exons in interval_hits(index, row["hg38_chr"], start1, end1, row["hg38_strand"]).items():
            if tid not in structures:
                continue
            length, abundance, values = structures[tid]
            tx_positions = set()
            for exon in exons:
                lo, hi = max(start1, exon.start), min(end1, exon.end)
                a, b = genome_to_tx(exon, lo), genome_to_tx(exon, hi)
                tx_positions.update(range(min(a, b), max(a, b) + 1))
            ordered = sorted(tx_positions)
            if not ordered:
                continue
            peak_values = [values[pos - 1] for pos in ordered]
            first, last = ordered[0], ordered[-1]
            upstream = values[max(0, first - flank - 1):first - 1]
            downstream = values[last:min(len(values), last + flank)]
            up_cov = sum(not math.isnan(v) for v in upstream) / flank if len(upstream) == flank else math.nan
            down_cov = sum(not math.isnan(v) for v in downstream) / flank if len(downstream) == flank else math.nan
            peak_valid = sum(not math.isnan(v) for v in peak_values)
            peak_fraction = peak_valid / len(peak_values)
            main_valid = (
                not math.isnan(up_cov) and not math.isnan(down_cov)
                and up_cov >= min_side and down_cov >= min_side and peak_fraction >= min_peak
            )
            splice_edges = [edge for item in structures_transcripts[tid][:-1] for edge in (item.tx_end, item.tx_end + 1)]
            candidates.append({
                **row, "transcript_id": tid, "transcript_gene_id": exons[0].gene_id,
                "transcript_gene_name": exons[0].gene_name, "transcript_length": length,
                "icshape_abundance_rpkm": "" if abundance == "*" else abundance,
                "peak_exonic_nt": len(ordered), "peak_tx_start_1based": first, "peak_tx_end_1based": last,
                "peak_midpoint_tx_1based": round((first + last) / 2, 1),
                "peak_valid_reactivity_count": peak_valid, "peak_valid_reactivity_fraction": round(peak_fraction, 8),
                "reactivity_mean_peak": valid_mean(peak_values),
                "coverage_upstream10": up_cov, "coverage_downstream10": down_cov,
                "reactivity_mean_upstream10": valid_mean(upstream), "reactivity_mean_downstream10": valid_mean(downstream),
                "distance_to_nearest_splice_edge_tx": min((abs(round((first + last) / 2) - edge) for edge in splice_edges), default=""),
                "main_window_valid": main_valid,
            })
        if not candidates:
            continue
        ranked = sorted(candidates, key=lambda item: (not item["main_window_valid"], -item["peak_valid_reactivity_count"], -item["peak_exonic_nt"], item["transcript_id"]))
        chosen = ranked[0]
        chosen["isoform_count"] = len(ranked)
        chosen["representative_isoform"] = True
        chosen["main_analysis_included"] = bool(chosen["main_window_valid"])
        chosen["main_exclusion_reason"] = "" if chosen["main_window_valid"] else "insufficient_peak_or_flank_structure_coverage"
        selected.append(chosen)
    return selected


structures_transcripts: dict = {}


def main() -> None:
    global structures_transcripts
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    out = config["outputs"]
    setup_logging(ROOT / out["log"])
    inputs = {key: pd.read_csv(ROOT / value, keep_default_na=False).to_dict("records") for key, value in config["inputs"].items()}
    refs = config["references"]
    gtf = ROOT / refs["gtf_hg38"]
    hek_lengths = read_structure_lengths_gzip(ROOT / refs["hek293t_structure"])
    hela_lengths = read_structure_lengths_tar(ROOT / refs["hela_structure_tar"], refs["hela_structure_member"])
    logging.info("structure transcript vectors: HeLa=%d HEK293T=%d", len(hela_lengths), len(hek_lengths))
    hela_tx = exact_length_transcripts(gtf, hela_lengths, versioned=True)
    hek_tx = exact_length_transcripts(gtf, hek_lengths, versioned=False)
    hela_index, hek_index = ExonBinIndex(hela_tx), ExonBinIndex(hek_tx)
    logging.info("exact-length hg38 transcripts: HeLa=%d HEK293T=%d", len(hela_tx), len(hek_tx))

    nm_offsets = [int(value) for value in config["coordinate_sensitivity"]["nm_offsets_nt"]]
    hela_wanted = (
        point_hit_tids(inputs["m5c_hela"], hela_index, [0])
        | point_hit_tids(inputs["nm_hela"], hela_index, nm_offsets)
        | peak_hit_tids(inputs["m7g_hela"], hela_index)
    )
    hek_wanted = (
        point_hit_tids(inputs["nm_hek"], hek_index, nm_offsets)
        | peak_hit_tids(inputs["m7g_hek293t"], hek_index)
    )
    hela_struct = load_selected_tar(ROOT / refs["hela_structure_tar"], refs["hela_structure_member"], hela_wanted)
    hek_struct = load_selected_gzip(ROOT / refs["hek293t_structure"], hek_wanted)
    if set(hela_struct) != hela_wanted or set(hek_struct) != hek_wanted:
        raise RuntimeError("Failed to load every selected structure vector")
    logging.info("loaded selected structure vectors: HeLa=%d HEK293T=%d", len(hela_struct), len(hek_struct))

    flank = int(config["structure"]["flank_each_side_nt"])
    display = int(config["structure"]["display_flank_each_side_nt"])
    min_side = float(config["structure"]["minimum_valid_fraction_each_side"])
    min_peak = float(config["structure"]["minimum_peak_valid_fraction"])
    structures_transcripts = hela_tx
    m5c = extract_points(inputs["m5c_hela"], hela_index, hela_struct, 0, flank, display, min_side)
    nm_hela_by_offset = {offset: extract_points(inputs["nm_hela"], hela_index, hela_struct, offset, flank, display, min_side) for offset in nm_offsets}
    m7g_hela = extract_peaks(inputs["m7g_hela"], hela_index, hela_struct, flank, min_side, min_peak)
    structures_transcripts = hek_tx
    nm_hek_by_offset = {offset: extract_points(inputs["nm_hek"], hek_index, hek_struct, offset, flank, display, min_side) for offset in nm_offsets}
    m7g_hek = extract_peaks(inputs["m7g_hek293t"], hek_index, hek_struct, flank, min_side, min_peak)
    nm_hela, nm_hek = nm_hela_by_offset[0], nm_hek_by_offset[0]

    datasets = {"m5c_hela": m5c, "m7g_hela": m7g_hela, "m7g_hek293t": m7g_hek, "nm_hela": nm_hela, "nm_hek": nm_hek}
    for key, rows in datasets.items():
        write_gzip_csv(ROOT / out[key], rows)

    attrition = []
    for key, rows in datasets.items():
        source_n = len(inputs[key])
        valid_n = sum(bool(row["main_analysis_included"]) for row in rows)
        attrition.append({
            "track_id": key, "modification": rows[0]["modification"], "cell_line": rows[0]["cell_line"],
            "input_records": source_n, "mapped_representative_records": len(rows),
            "mapping_fraction": round(len(rows) / source_n, 8), "main_valid_records": valid_n,
            "main_valid_fraction_of_input": round(valid_n / source_n, 8),
            "measurement_scale": "quantitative_site_fraction" if key == "m5c_hela" else ("semiquantitative_interval_peak" if key.startswith("m7g") else "binary_site"),
            "track_admission": "INFERENTIAL_READY" if valid_n >= config["admission"]["minimum_main_records_per_modification"] else "DESCRIPTIVE_ONLY_LOW_N",
        })
    sensitivity = []
    stability_by_track = {}
    for key, groups in (("nm_hela", nm_hela_by_offset), ("nm_hek", nm_hek_by_offset)):
        valid_counts = [sum(bool(row["main_analysis_included"]) for row in groups[offset]) for offset in nm_offsets]
        ratio = min(valid_counts) / max(valid_counts) if max(valid_counts) else 0
        stability_by_track[key] = ratio
        for offset, valid in zip(nm_offsets, valid_counts):
            sensitivity.append({
                "track_id": key, "coordinate_offset_nt": offset, "mapped_records": len(groups[offset]),
                "main_valid_records": valid, "valid_count_stability_ratio_across_offsets": round(ratio, 8),
                "stability_pass": ratio >= config["coordinate_sensitivity"]["minimum_valid_count_stability_ratio"],
            })
    write_csv(ROOT / out["attrition"], attrition)
    write_csv(ROOT / out["nm_sensitivity"], sensitivity)

    minimum = int(config["admission"]["minimum_main_records_per_modification"])
    admitted = set()
    for modification in {row["modification"] for row in attrition}:
        if sum(row["main_valid_records"] for row in attrition if row["modification"] == modification) >= minimum:
            admitted.add(modification)
    nm_stable = all(value >= config["coordinate_sensitivity"]["minimum_valid_count_stability_ratio"] for value in stability_by_track.values())
    if not nm_stable:
        admitted.discard("Nm")
    gate = "GO_TO_STRATIFIED_ASSOCIATION_ANALYSIS" if len(admitted) >= config["admission"]["minimum_distinct_modifications_for_stage26"] else "NO_GO"
    status = "PASS" if gate.startswith("GO") else "FAIL"
    now = datetime.now(timezone.utc)
    write_csv(ROOT / out["status"], [{
        "stage": 25, "status": status, "tracks_built": len(datasets),
        "distinct_modifications_admitted": len(admitted), "nm_coordinate_sensitivity_pass": nm_stable,
        "phase3_gate": gate, "raw_scales_pooled": False, "finished_utc": now.isoformat(),
    }])
    manifest = []
    for key in ["m5c_hela", "m7g_hela", "m7g_hek293t", "nm_hela", "nm_hek", "attrition", "nm_sensitivity", "status"]:
        path = ROOT / out[key]
        manifest.append({"artifact": key, "path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    write_csv(ROOT / out["manifest"], manifest)

    by_id = {row["track_id"]: row for row in attrition}
    report = f"""# Stage 25: Multi-modified structural feature extraction

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {now.date().isoformat()}
- Verification Status: {'VERIFIED' if status == 'PASS' else 'UNVERIFIED'}
- Version Label: phase3_stage25_v1

At this stage, representative transcripts were selected and icSHAPE features were extracted for the 5 technology-layered tracks. The stage status is `{status}` and the threshold is `{gate}`; the original modified dimensions are not merged.

##Analyzable sample size

- m5C/HeLa: Input {by_id['m5c_hela']['input_records']:,}, mapped {by_id['m5c_hela']['mapped_representative_records']:,}, ±10 nt structure window qualified {by_id['m5c_hela']['main_valid_records']:,}.
- m7G / HeLa: Input {by_id['m7g_hela']['input_records']:,} intervals, map {by_id['m7g_hela']['mapped_representative_records']:,}, and the structure coverage within the peak and on both sides is qualified {by_id['m7g_hela']['main_valid_records']:,}.
- m7G / HEK293T: Input {by_id['m7g_hek293t']['input_records']:,} intervals, map {by_id['m7g_hek293t']['mapped_representative_records']:,}, and cover qualified {by_id['m7g_hek293t']['main_valid_records']:,}; below the 20 single-track threshold, only descriptive results are retained, and no independent inference is made.
- Nm/HeLa: Input {by_id['nm_hela']['input_records']:,}, mapped {by_id['nm_hela']['mapped_representative_records']:,}, ±10 nt structural window qualified {by_id['nm_hela']['main_valid_records']:,}.
- Nm/HEK: input {by_id['nm_hek']['input_records']:,}, map {by_id['nm_hek']['mapped_representative_records']:,}, window qualified {by_id['nm_hek']['main_valid_records']:,}; only as HEK293T approximate cell label.

## Coordinate sensitivity and boundaries

The −1/0/+1 nt qualifying sample size stability ratio of Nm is HeLa {stability_by_track['nm_hela']:.4f}, HEK {stability_by_track['nm_hek']:.4f}, and the threshold is {config['coordinate_sensitivity']['minimum_valid_count_stability_ratio']:.2f}. m7G is always analyzed by peak interval, `reactivity_mean_peak` does not represent the local structure of a single m7G base. The next phase can only be modeled hierarchically by modification, technology and cell line."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
