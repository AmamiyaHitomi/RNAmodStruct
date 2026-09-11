"""Stratified HEK293T transcript→hg19→hg38 coordinate mapping prototype."""

from __future__ import annotations

import argparse
import gzip
import logging
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from pipeline_common import (
    ChainLiftOver,
    coverage,
    fetch_fasta_bases,
    load_gtf_exons,
    load_icshape,
    load_selected_fasta,
    norm_chr,
    read_glori,
    tx_to_genome,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_processed"
REF = ROOT / "data" / "reference"
INTERIM = ROOT / "data" / "interim"
COORD = ROOT / "metadata" / "coordinate"
LOG = ROOT / "results" / "logs" / "04_validate_coordinate_mapping.log"


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def load_glori_keys() -> set[tuple[str, int, str]]:
    paths = [
        RAW / "GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv.gz",
        RAW / "GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv.gz",
    ]
    return {
        (norm_chr(row["Chr"]), int(row["Sites"]), row["Strand"])
        for path in paths
        for row in read_glori(path)
    }


def is_splice_boundary_position(exons, pos: int, distance: int = 2) -> bool:
    junctions = []
    for exon in exons[:-1]:
        junctions.extend([exon.tx_end, exon.tx_end + 1])
    return any(abs(pos - boundary) <= distance for boundary in junctions)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--per-stratum", type=int, default=3)
    parser.add_argument("--max-checks", type=int, default=30)
    args = parser.parse_args()
    setup_logging()
    started = datetime.now().isoformat(timespec="seconds")

    icshape_path = RAW / "GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz"
    gtf_path = REF / "ensembl" / "release-74_GRCh37" / "Homo_sapiens.GRCh37.74.gtf.gz"
    cdna_path = REF / "ensembl" / "release-74_GRCh37" / "Homo_sapiens.GRCh37.74.cdna.all.fa.gz"
    genome38_path = REF / "ensembl" / "release-88_GRCh38" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    chain_path = REF / "ucsc" / "hg19ToHg38.over.chain.gz"

    logging.info("loading HEK293T icSHAPE")
    icshape = load_icshape(icshape_path)
    ids = set(icshape)
    logging.info("loading GRCh37.74 exons and cDNA for %d transcripts", len(ids))
    transcripts = load_gtf_exons(gtf_path, ids)
    sequences = load_selected_fasta(cdna_path, ids)
    lifter = ChainLiftOver(chain_path)
    glori_keys = load_glori_keys()
    logging.info("loaded %d GTF transcripts, %d cDNAs, %d GLORI union sites", len(transcripts), len(sequences), len(glori_keys))

    strata: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    extras: list[dict] = []
    failures: list[dict] = []
    scanned_a = 0
    eligible_transcripts = 0
    for tid in sorted(ids):
        if tid not in transcripts or tid not in sequences:
            continue
        length, abundance, values = icshape[tid]
        seq = sequences[tid]
        exons = transcripts[tid]
        if len(seq) != length or sum(e.end - e.start + 1 for e in exons) != length:
            continue
        eligible_transcripts += 1
        for pos, (base, score) in enumerate(zip(seq, values), start=1):
            if base != "A" or math.isnan(score):
                continue
            scanned_a += 1
            mapped19 = tx_to_genome(exons, pos)
            if mapped19 is None:
                continue
            chrom19, pos19, strand19, exon_number = mapped19
            mappings = lifter.forward(chrom19, pos19, strand19)
            boundary = "splice_boundary" if is_splice_boundary_position(exons, pos) else "exon_internal"
            base_row = {
                "transcript_id": tid,
                "transcript_pos_1based": pos,
                "coordinate_basis": "1-based_transcript_oriented",
                "icshape_reactivity": round(float(score), 6),
                "transcript_base": base,
                "gene_id": exons[0].gene_id,
                "gene_name": exons[0].gene_name,
                "exon_number": exon_number,
                "exon_count": len(exons),
                "case_type": boundary,
                "hg19_chr": chrom19,
                "hg19_pos_1based": pos19,
                "hg19_strand": strand19,
                "liftover_mapping_count": len(mappings),
            }
            if len(mappings) != 1:
                if len(failures) < 6:
                    failures.append({**base_row, "mapping_status": "unmapped" if not mappings else "ambiguous"})
                continue
            chrom38, pos38, strand38, chain_score = mappings[0]
            glori_match = (chrom38, pos38, strand38) in glori_keys
            left10, right10, cov10 = coverage(values, pos, 10)
            _, _, cov50 = coverage(values, pos, 50)
            row = {
                **base_row,
                "mapping_status": "unique",
                "hg38_chr": chrom38,
                "hg38_pos_1based": pos38,
                "hg38_strand": strand38,
                "chain_score": chain_score,
                "glori_match": glori_match,
                "coverage_up10": left10,
                "coverage_down10": right10,
                "coverage_flank10": cov10,
                "coverage_flank50": cov50,
            }
            key = (strand19, boundary, "glori_match" if glori_match else "no_glori_match")
            if len(strata[key]) < args.per_stratum:
                strata[key].append(row)
            elif len(extras) < args.max_checks:
                extras.append(row)

    selected = []
    for key in sorted(strata):
        selected.extend(strata[key])
    seen = {(r["transcript_id"], r["transcript_pos_1based"]) for r in selected}
    for row in extras:
        key = (row["transcript_id"], row["transcript_pos_1based"])
        if key not in seen and len(selected) < args.max_checks:
            selected.append(row)
            seen.add(key)
    for row in failures:
        key = (row["transcript_id"], row["transcript_pos_1based"])
        if key not in seen and len(selected) < args.max_checks:
            selected.append(row)
            seen.add(key)

    queries: dict[str, set[int]] = defaultdict(set)
    for row in selected:
        if row["mapping_status"] == "unique":
            queries[row["hg38_chr"]].add(row["hg38_pos_1based"])
    logging.info("streaming GRCh38 FASTA for %d selected center bases", sum(map(len, queries.values())))
    bases38 = fetch_fasta_bases(genome38_path, queries)

    checks = []
    for index, row in enumerate(selected, start=1):
        row = dict(row)
        row["check_id"] = f"HEK_CHECK_{index:03d}"
        if row["mapping_status"] == "unique":
            base38 = bases38.get((row["hg38_chr"], row["hg38_pos_1based"]), "")
            expected = "A" if row["hg38_strand"] == "+" else "T"
            reverse = lifter.reverse(row["hg38_chr"], row["hg38_pos_1based"], row["hg38_strand"])
            roundtrip = any(
                chrom == row["hg19_chr"] and pos == row["hg19_pos_1based"] and strand == row["hg19_strand"]
                for chrom, pos, strand, _ in reverse
            )
            row.update(
                {
                    "hg38_reference_base": base38,
                    "expected_genomic_base": expected,
                    "check_transcript_base_A": row["transcript_base"] == "A",
                    "check_hg38_center_base": base38 == expected,
                    "check_liftover_roundtrip": roundtrip,
                    "check_pass": row["transcript_base"] == "A" and base38 == expected and roundtrip,
                }
            )
        else:
            row.update(
                {
                    "hg38_chr": "",
                    "hg38_pos_1based": "",
                    "hg38_strand": "",
                    "chain_score": "",
                    "glori_match": False,
                    "hg38_reference_base": "",
                    "expected_genomic_base": "",
                    "check_transcript_base_A": row["transcript_base"] == "A",
                    "check_hg38_center_base": "NOT_APPLICABLE",
                    "check_liftover_roundtrip": "NOT_APPLICABLE",
                    "check_pass": row["transcript_base"] == "A",
                }
            )
        checks.append(row)

    COORD.mkdir(parents=True, exist_ok=True)
    INTERIM.mkdir(parents=True, exist_ok=True)
    write_csv(COORD / "04_hek293t_coordinate_checks.csv", checks)
    write_csv(INTERIM / "04_hek293t_coordinate_mapping_prototype.csv", checks)
    stratum_counts = Counter(
        f"{r['hg19_strand']}|{r['case_type']}|{'match' if r.get('glori_match') else 'no_match'}"
        for r in checks if r["mapping_status"] == "unique"
    )
    unique_checks = [r for r in checks if r["mapping_status"] == "unique"]
    failed_assertions = [r for r in unique_checks if not r["check_pass"]]
    summary = [{
        "status": "PASS" if len(unique_checks) >= 20 and not failed_assertions else "FAIL",
        "started_at": started,
        "finished_at": datetime.now().isoformat(timespec="seconds"),
        "eligible_transcripts": eligible_transcripts,
        "measured_a_positions_scanned": scanned_a,
        "selected_checks": len(checks),
        "unique_mapped_checks": len(unique_checks),
        "glori_matched_checks": sum(bool(r.get("glori_match")) for r in unique_checks),
        "failed_assertions": len(failed_assertions),
        "stratum_counts": dict(sorted(stratum_counts.items())),
    }]
    write_csv(COORD / "04_hek293t_coordinate_check_summary.csv", summary)
    logging.info("coordinate prototype finished: %s", summary[0])
    if summary[0]["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
