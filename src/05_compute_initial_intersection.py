"""Full preliminary HEK293T GLORI×icSHAPE intersection after validated mapping."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from pipeline_common import (
    ChainLiftOver,
    ExonBinIndex,
    coverage,
    fetch_fasta_bases,
    genome_to_tx,
    load_gtf_exons,
    load_icshape,
    load_selected_fasta,
    norm_chr,
    read_glori,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_processed"
REF = ROOT / "data" / "reference"
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"
LOG = RESULTS / "logs" / "05_compute_initial_intersection.log"


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def load_glori_table(path: Path) -> dict[tuple[str, int, str], dict[str, str]]:
    rows = {}
    for row in read_glori(path):
        key = (norm_chr(row["Chr"]), int(row["Sites"]), row["Strand"])
        if key in rows:
            raise ValueError(f"Duplicate GLORI genomic key in {path.name}: {key}")
        rows[key] = row
    return rows


def write_gzip_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_median(values: list[float]) -> str | float:
    return statistics.median(values) if values else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-flank-coverage", type=float, default=0.70)
    args = parser.parse_args()
    setup_logging()
    started = datetime.now().isoformat(timespec="seconds")

    rep1_path = RAW / "GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv.gz"
    rep2_path = RAW / "GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv.gz"
    logging.info("loading GLORI replicates")
    rep1 = load_glori_table(rep1_path)
    rep2 = load_glori_table(rep2_path)
    keys1, keys2 = set(rep1), set(rep2)
    union_keys = sorted(keys1 | keys2)
    common_keys = keys1 & keys2

    icshape_path = RAW / "GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz"
    gtf_path = REF / "ensembl" / "release-74_GRCh37" / "Homo_sapiens.GRCh37.74.gtf.gz"
    cdna_path = REF / "ensembl" / "release-74_GRCh37" / "Homo_sapiens.GRCh37.74.cdna.all.fa.gz"
    genome38_path = REF / "ensembl" / "release-88_GRCh38" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    chain_path = REF / "ucsc" / "hg19ToHg38.over.chain.gz"

    logging.info("loading icSHAPE vectors and compatible GRCh37.74 annotation")
    icshape = load_icshape(icshape_path)
    ids = set(icshape)
    transcripts_all = load_gtf_exons(gtf_path, ids)
    sequences_all = load_selected_fasta(cdna_path, ids)
    transcripts = {}
    sequences = {}
    exclusion = Counter()
    for tid in sorted(ids):
        if tid not in transcripts_all:
            exclusion["no_gtf_transcript"] += 1
            continue
        if tid not in sequences_all:
            exclusion["no_cdna_sequence"] += 1
            continue
        declared_length = icshape[tid][0]
        exon_length = sum(e.end - e.start + 1 for e in transcripts_all[tid])
        if len(sequences_all[tid]) != declared_length:
            exclusion["cdna_length_mismatch"] += 1
            continue
        if exon_length != declared_length:
            exclusion["gtf_exon_length_mismatch"] += 1
            continue
        transcripts[tid] = transcripts_all[tid]
        sequences[tid] = sequences_all[tid]
    logging.info("eligible icSHAPE transcripts=%d excluded=%s", len(transcripts), dict(exclusion))

    exon_index = ExonBinIndex(transcripts)
    lifter = ChainLiftOver(chain_path)
    status_rows: list[dict] = []
    intersection_rows: list[dict] = []
    reasons = Counter()
    matched_site_queries: dict[str, set[int]] = defaultdict(set)
    matched_genes = set()

    logging.info("mapping %d GLORI union sites back to GRCh37", len(union_keys))
    for counter, key in enumerate(union_keys, start=1):
        chrom38, pos38, strand38 = key
        row1 = rep1.get(key)
        row2 = rep2.get(key)
        genes = sorted({r["Gene"] for r in (row1, row2) if r and r["Gene"] != "ELSE"})
        base_status = {
            "hg38_chr": chrom38,
            "hg38_pos_1based": pos38,
            "hg38_strand": strand38,
            "rep1_present": row1 is not None,
            "rep2_present": row2 is not None,
            "rep1_ratio": row1["Ratio"] if row1 else "",
            "rep2_ratio": row2["Ratio"] if row2 else "",
            "genes": ";".join(genes),
        }
        reverse = lifter.reverse(chrom38, pos38, strand38)
        if not reverse:
            status_rows.append({**base_status, "mapping_status": "liftover_unmapped", "transcript_hits": 0, "valid_main_windows": 0})
            reasons["liftover_unmapped"] += 1
            continue
        if len(reverse) != 1:
            status_rows.append({**base_status, "mapping_status": "liftover_ambiguous", "transcript_hits": 0, "valid_main_windows": 0})
            reasons["liftover_ambiguous"] += 1
            continue
        chrom19, pos19, strand19, chain_score = reverse[0]
        exon_hits = exon_index.query(chrom19, pos19, strand19)
        if not exon_hits:
            status_rows.append({**base_status, "mapping_status": "no_icshape_exon", "hg19_chr": chrom19, "hg19_pos_1based": pos19, "hg19_strand": strand19, "transcript_hits": 0, "valid_main_windows": 0})
            reasons["no_icshape_exon"] += 1
            continue
        valid_hits = 0
        valid_windows = 0
        non_a_hits = 0
        seen_transcripts = set()
        for exon in exon_hits:
            tid = exon.transcript_id
            if tid in seen_transcripts:
                continue
            seen_transcripts.add(tid)
            tx_pos = genome_to_tx(exon, pos19)
            seq = sequences[tid]
            if not 1 <= tx_pos <= len(seq) or seq[tx_pos - 1] != "A":
                non_a_hits += 1
                continue
            valid_hits += 1
            values = icshape[tid][2]
            center_score = values[tx_pos - 1]
            up10, down10, flank10 = coverage(values, tx_pos, 10)
            _, _, flank50 = coverage(values, tx_pos, 50)
            window_valid = (
                not math.isnan(up10)
                and not math.isnan(down10)
                and up10 >= args.min_flank_coverage
                and down10 >= args.min_flank_coverage
            )
            valid_windows += int(window_valid)
            intersection_rows.append(
                {
                    **base_status,
                    "hg19_chr": chrom19,
                    "hg19_pos_1based": pos19,
                    "hg19_strand": strand19,
                    "chain_score": chain_score,
                    "transcript_id": tid,
                    "transcript_pos_1based": tx_pos,
                    "transcript_gene_id": exon.gene_id,
                    "transcript_gene_name": exon.gene_name,
                    "transcript_base": "A",
                    "icshape_center": "" if math.isnan(center_score) else round(float(center_score), 6),
                    "coverage_up10": up10,
                    "coverage_down10": down10,
                    "coverage_flank10": flank10,
                    "coverage_flank50": flank50,
                    "main_window_valid": window_valid,
                }
            )
        if valid_hits:
            mapping_status = "matched_transcript_position"
            reasons[mapping_status] += 1
            matched_site_queries[chrom38].add(pos38)
            matched_genes.update(genes)
        else:
            mapping_status = "center_not_transcript_A"
            reasons[mapping_status] += 1
        status_rows.append(
            {
                **base_status,
                "mapping_status": mapping_status,
                "hg19_chr": chrom19,
                "hg19_pos_1based": pos19,
                "hg19_strand": strand19,
                "transcript_hits": valid_hits,
                "non_a_exon_hits": non_a_hits,
                "valid_main_windows": valid_windows,
            }
        )
        if counter % 50000 == 0:
            logging.info("mapped %d/%d sites", counter, len(union_keys))

    logging.info("checking GRCh38 reference bases for %d matched sites", sum(map(len, matched_site_queries.values())))
    bases38 = fetch_fasta_bases(genome38_path, matched_site_queries)
    center_base_failures = 0
    for row in status_rows:
        if row["mapping_status"] != "matched_transcript_position":
            row["hg38_reference_base"] = ""
            row["hg38_center_base_ok"] = "NOT_APPLICABLE"
            continue
        observed = bases38.get((row["hg38_chr"], row["hg38_pos_1based"]), "")
        expected = "A" if row["hg38_strand"] == "+" else "T"
        row["hg38_reference_base"] = observed
        row["hg38_center_base_ok"] = observed == expected
        center_base_failures += int(observed != expected)

    matched_sites = [r for r in status_rows if r["mapping_status"] == "matched_transcript_position"]
    valid_window_sites = [r for r in matched_sites if int(r["valid_main_windows"]) > 0]
    matched_keys = {(r["hg38_chr"], int(r["hg38_pos_1based"]), r["hg38_strand"]) for r in matched_sites}
    valid_window_keys = {(r["hg38_chr"], int(r["hg38_pos_1based"]), r["hg38_strand"]) for r in valid_window_sites}
    mapped_coverages = [float(r["coverage_flank10"]) for r in intersection_rows if r["coverage_flank10"] != "" and not math.isnan(float(r["coverage_flank10"]))]
    valid_mapping_rows = [r for r in intersection_rows if r["main_window_valid"]]
    center_available_rows = [r for r in intersection_rows if r["icshape_center"] != ""]

    summary_rows = [
        {"section": "glori_replicates", "metric": "rep1_unique_sites", "value": len(keys1), "notes": "Exact (Chr, Sites, Strand) keys"},
        {"section": "glori_replicates", "metric": "rep2_unique_sites", "value": len(keys2), "notes": "Exact (Chr, Sites, Strand) keys"},
        {"section": "glori_replicates", "metric": "replicate_common_sites", "value": len(common_keys), "notes": "Exact genomic-site overlap"},
        {"section": "glori_replicates", "metric": "replicate_union_sites", "value": len(union_keys), "notes": "Exact genomic-site union"},
        {"section": "coordinate_mapping", "metric": "matched_union_sites", "value": len(matched_sites), "notes": "At least one exact HEK293T icSHAPE transcript position"},
        {"section": "coordinate_mapping", "metric": "matched_rep1_sites", "value": len(keys1 & matched_keys), "notes": "Replicate-1 sites with at least one icSHAPE transcript position"},
        {"section": "coordinate_mapping", "metric": "matched_rep2_sites", "value": len(keys2 & matched_keys), "notes": "Replicate-2 sites with at least one icSHAPE transcript position"},
        {"section": "coordinate_mapping", "metric": "matched_common_sites", "value": len(common_keys & matched_keys), "notes": "Replicate-common sites with at least one icSHAPE transcript position"},
        {"section": "coordinate_mapping", "metric": "matched_transcript_rows", "value": len(intersection_rows), "notes": "One genomic site can map to multiple isoforms"},
        {"section": "coordinate_mapping", "metric": "independent_glori_genes", "value": len(matched_genes), "notes": "Non-ELSE GLORI gene labels among matched sites"},
        {"section": "structure_coverage", "metric": "sites_with_valid_flank10_window", "value": len(valid_window_sites), "notes": f"Both ±10 sides separately ≥ {args.min_flank_coverage:.0%}; center excluded"},
        {"section": "structure_coverage", "metric": "site_valid_window_fraction", "value": round(len(valid_window_sites) / len(matched_sites), 8) if matched_sites else "", "notes": "At least one qualifying transcript window per genomic site"},
        {"section": "structure_coverage", "metric": "valid_transcript_windows", "value": len(valid_mapping_rows), "notes": "Transcript-level qualifying windows"},
        {"section": "structure_coverage", "metric": "transcript_window_valid_fraction", "value": round(len(valid_mapping_rows) / len(intersection_rows), 8) if intersection_rows else "", "notes": "Transcript-level denominator"},
        {"section": "structure_coverage", "metric": "rep1_sites_with_valid_flank10_window", "value": len(keys1 & valid_window_keys), "notes": "At least one qualifying transcript window"},
        {"section": "structure_coverage", "metric": "rep2_sites_with_valid_flank10_window", "value": len(keys2 & valid_window_keys), "notes": "At least one qualifying transcript window"},
        {"section": "structure_coverage", "metric": "center_reactivity_available_rows", "value": len(center_available_rows), "notes": "Transcript mappings with a non-missing center score"},
        {"section": "structure_coverage", "metric": "center_reactivity_available_fraction", "value": round(len(center_available_rows) / len(intersection_rows), 8) if intersection_rows else "", "notes": "Transcript-level denominator"},
        {"section": "structure_coverage", "metric": "mean_flank10_coverage", "value": round(sum(mapped_coverages) / len(mapped_coverages), 8) if mapped_coverages else "", "notes": "Only full-length ±10 windows"},
        {"section": "structure_coverage", "metric": "median_flank10_coverage", "value": safe_median(mapped_coverages), "notes": "Across mapped transcript rows"},
        {"section": "validation", "metric": "grch38_center_base_failures", "value": center_base_failures, "notes": "Expected A on + and T on -"},
        {"section": "reference", "metric": "eligible_icshape_transcripts", "value": len(transcripts), "notes": "Exact length agreement in icSHAPE, cDNA and GTF"},
    ]
    for reason, count in sorted(reasons.items()):
        summary_rows.append({"section": "mapping_status", "metric": reason, "value": count, "notes": "GLORI union-site count"})
    for reason, count in sorted(exclusion.items()):
        summary_rows.append({"section": "transcript_exclusion", "metric": reason, "value": count, "notes": "icSHAPE transcript count"})

    INTERIM.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    status_fields = list(dict.fromkeys(key for row in status_rows for key in row))
    intersection_fields = list(dict.fromkeys(key for row in intersection_rows for key in row)) if intersection_rows else []
    write_gzip_csv(INTERIM / "05_hek293t_glori_site_mapping_status.csv.gz", status_rows, status_fields)
    if intersection_rows:
        write_gzip_csv(INTERIM / "05_hek293t_glori_icshape_intersection.csv.gz", intersection_rows, intersection_fields)
    write_csv(RESULTS / "05_hek293t_initial_intersection_summary.csv", summary_rows)
    write_csv(
        RESULTS / "05_hek293t_unmatched_reason_counts.csv",
        [
            {"mapping_status": reason, "sites": count, "fraction_of_union": round(count / len(union_keys), 8)}
            for reason, count in sorted(reasons.items())
            if reason != "matched_transcript_position"
        ],
    )
    run_status = "PASS" if matched_sites and center_base_failures == 0 else "FAIL"
    write_csv(
        RESULTS / "05_hek293t_initial_intersection_run_status.csv",
        [{"status": run_status, "started_at": started, "finished_at": datetime.now().isoformat(timespec="seconds"), "glori_union_sites": len(union_keys), "matched_sites": len(matched_sites), "center_base_failures": center_base_failures, "parameters": json.dumps({"min_flank_coverage": args.min_flank_coverage}, sort_keys=True)}],
    )
    logging.info("initial intersection finished: status=%s matched_sites=%d valid_window_sites=%d reasons=%s", run_status, len(matched_sites), len(valid_window_sites), dict(reasons))
    if run_status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
