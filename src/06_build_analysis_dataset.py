"""Build analysis-ready HEK293T site and isoform tables from the stage-05 join."""

from __future__ import annotations

import csv
import gzip
import json
import logging
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import yaml
from scipy.stats import pearsonr, spearmanr

from pipeline_common import (
    genome_to_tx,
    load_gtf_exons,
    load_icshape,
    load_selected_fasta,
    norm_chr,
    parse_gtf_attributes,
    read_glori,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_processed"
INTERIM = ROOT / "data" / "interim"
FINAL = ROOT / "data" / "final"
RESULTS = ROOT / "results"
CONFIG_PATH = ROOT / "config" / "06_analysis_config.yaml"
LOG = RESULTS / "logs" / "06_build_analysis_dataset.log"


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def read_gzip_csv(path: Path):
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def write_gzip_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        raise ValueError(f"Cannot infer fields for empty output: {path}")
    names = fieldnames or list(dict.fromkeys(key for row in rows for key in row))
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_glori(path: Path) -> dict[tuple[str, int, str], dict[str, str]]:
    return {
        (norm_chr(row["Chr"]), int(row["Sites"]), row["Strand"]): row
        for row in read_glori(path)
    }


def json_window(values, center1: int, flank: int) -> tuple[str, str, int]:
    output = []
    mask = []
    for pos1 in range(center1 - flank, center1 + flank + 1):
        if pos1 < 1 or pos1 > len(values) or math.isnan(values[pos1 - 1]):
            output.append(None)
            mask.append(0)
        else:
            output.append(round(float(values[pos1 - 1]), 6))
            mask.append(1)
    return json.dumps(output, separators=(",", ":")), json.dumps(mask, separators=(",", ":")), sum(mask)


def valid_mean(values) -> str | float:
    kept = [float(x) for x in values if not math.isnan(x)]
    return round(sum(kept) / len(kept), 8) if kept else ""


def gene_matches(glori_genes: str, transcript_gene: str) -> bool:
    return bool(transcript_gene) and transcript_gene in {x for x in glori_genes.split(";") if x}


def choose_analysis_gene(glori_genes: str, transcript_gene: str, gene_id: str) -> tuple[str, str]:
    submitted = [x for x in glori_genes.split(";") if x and x != "ELSE"]
    if submitted:
        if transcript_gene in submitted:
            return transcript_gene, "GLORI_gene_matching_transcript"
        return submitted[0], "GLORI_gene"
    if transcript_gene:
        return transcript_gene, "Ensembl_transcript_gene_name"
    return gene_id, "Ensembl_gene_id"


def load_coding_annotations(gtf_path: Path, transcript_ids: set[str], exons_by_transcript):
    cds_features = defaultdict(list)
    stop_features = defaultdict(list)
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"CDS", "stop_codon"}:
                continue
            tid = parse_gtf_attributes(fields[8]).get("transcript_id")
            if tid not in transcript_ids:
                continue
            target = cds_features if fields[2] == "CDS" else stop_features
            target[tid].append((int(fields[3]), int(fields[4])))

    annotations = {}
    for tid in transcript_ids:
        exons = exons_by_transcript.get(tid, [])
        cds_tx = []
        stop_tx = []
        for feature_start, feature_end in cds_features.get(tid, []):
            for genomic_pos in (feature_start, feature_end):
                exon = next((e for e in exons if e.start <= genomic_pos <= e.end), None)
                if exon:
                    cds_tx.append(genome_to_tx(exon, genomic_pos))
        for feature_start, feature_end in stop_features.get(tid, []):
            for genomic_pos in (feature_start, feature_end):
                exon = next((e for e in exons if e.start <= genomic_pos <= e.end), None)
                if exon:
                    stop_tx.append(genome_to_tx(exon, genomic_pos))
        annotations[tid] = {
            "cds_start_tx": min(cds_tx) if cds_tx else None,
            "cds_end_tx": max(cds_tx) if cds_tx else None,
            "stop_codon_start_tx": min(stop_tx) if stop_tx else None,
        }
    return annotations


def transcript_region(tx_pos: int, annotation: dict) -> str:
    start, end = annotation["cds_start_tx"], annotation["cds_end_tx"]
    if start is None or end is None:
        return "noncoding_or_CDS_unavailable"
    if tx_pos < start:
        return "5UTR"
    if tx_pos <= end:
        return "CDS"
    return "3UTR"


def drach_subtype(motif: str) -> str:
    if len(motif) != 5:
        return "window_unavailable"
    is_drach = motif[0] in "AGT" and motif[1] in "AG" and motif[2] == "A" and motif[3] == "C" and motif[4] in "ACT"
    return motif if is_drach else "non_DRACH"


def correlation_rows(main_rows: list[dict]) -> list[dict]:
    x = [float(row["ratio_rep1"]) for row in main_rows]
    y = [float(row["ratio_rep2"]) for row in main_rows]
    if len(x) < 3:
        return []
    pearson = pearsonr(x, y)
    spearman = spearmanr(x, y)
    differences = [abs(a - b) for a, b in zip(x, y)]
    return [
        {"category": "replicate_qc", "metric": "pearson_ratio", "value": pearson.statistic, "n": len(x), "notes": "Ratio rep1 versus rep2 in main dataset"},
        {"category": "replicate_qc", "metric": "spearman_ratio", "value": spearman.statistic, "n": len(x), "notes": "Ratio rep1 versus rep2 in main dataset"},
        {"category": "replicate_qc", "metric": "median_absolute_ratio_difference", "value": statistics.median(differences), "n": len(x), "notes": "Absolute replicate difference"},
    ]


def main() -> None:
    setup_logging()
    started = datetime.now().isoformat(timespec="seconds")
    with CONFIG_PATH.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    flank10 = int(config["structure"]["main_window_each_side_nt"])
    min_coverage = float(config["structure"]["minimum_valid_fraction_each_side"])
    flank50 = int(config["structure"]["display_window_each_side_nt"])
    model_flank = (int(config["sequence"]["model_window_nt"]) - 1) // 2

    rep1 = load_glori(RAW / "GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv.gz")
    rep2 = load_glori(RAW / "GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv.gz")
    logging.info("loaded GLORI replicates: %d and %d sites", len(rep1), len(rep2))

    intersection_path = INTERIM / "05_hek293t_glori_icshape_intersection.csv.gz"
    base_rows = list(read_gzip_csv(intersection_path))
    transcript_ids = {row["transcript_id"] for row in base_rows}
    icshape = load_icshape(RAW / "GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz")
    cdna_path = ROOT / "data" / "reference" / "ensembl" / "release-74_GRCh37" / "Homo_sapiens.GRCh37.74.cdna.all.fa.gz"
    gtf_path = ROOT / "data" / "reference" / "ensembl" / "release-74_GRCh37" / "Homo_sapiens.GRCh37.74.gtf.gz"
    sequences = load_selected_fasta(cdna_path, transcript_ids)
    exons_by_transcript = load_gtf_exons(gtf_path, transcript_ids)
    coding_annotations = load_coding_annotations(gtf_path, transcript_ids, exons_by_transcript)
    logging.info("enriching %d site/transcript rows across %d transcripts", len(base_rows), len(transcript_ids))

    enriched_by_site: dict[tuple[str, int, str], list[dict]] = defaultdict(list)
    for row in base_rows:
        key = (row["hg38_chr"], int(row["hg38_pos_1based"]), row["hg38_strand"])
        tid = row["transcript_id"]
        tx_pos = int(row["transcript_pos_1based"])
        values = icshape[tid][2]
        sequence = sequences[tid]
        left10 = values[max(0, tx_pos - flank10 - 1) : tx_pos - 1]
        right10 = values[tx_pos : min(len(values), tx_pos + flank10)]
        valid_left = sum(not math.isnan(x) for x in left10)
        valid_right = sum(not math.isnan(x) for x in right10)
        main_window_valid = (
            len(left10) == flank10
            and len(right10) == flank10
            and valid_left / flank10 >= min_coverage
            and valid_right / flank10 >= min_coverage
        )
        reaction_json, mask_json, valid_101 = json_window(values, tx_pos, flank50)
        full_sequence = model_flank < tx_pos <= len(sequence) - model_flank
        sequence_201 = sequence[tx_pos - model_flank - 1 : tx_pos + model_flank] if full_sequence else ""
        motif_5 = sequence[tx_pos - 3 : tx_pos + 2] if 3 <= tx_pos <= len(sequence) - 2 else ""
        local_21 = sequence[tx_pos - 11 : tx_pos + 10] if 11 <= tx_pos <= len(sequence) - 10 else ""
        gc_fraction_21 = (local_21.count("G") + local_21.count("C")) / len(local_21) if local_21 else ""
        annotation = coding_annotations[tid]
        region = transcript_region(tx_pos, annotation)
        stop_distance = tx_pos - annotation["stop_codon_start_tx"] if annotation["stop_codon_start_tx"] is not None else ""
        junction_edges = [edge for exon in exons_by_transcript[tid][:-1] for edge in (exon.tx_end, exon.tx_end + 1)]
        splice_distance = min(abs(tx_pos - edge) for edge in junction_edges) if junction_edges else ""
        match_gene = gene_matches(row["genes"], row["transcript_gene_name"])
        enriched = {
            **row,
            "valid_count_up10": valid_left,
            "valid_count_down10": valid_right,
            "reactivity_mean_up10": valid_mean(left10),
            "reactivity_mean_down10": valid_mean(right10),
            "reactivity_mean_flank10": valid_mean(list(left10) + list(right10)),
            "main_window_valid": main_window_valid,
            "reactivity_window_m50_p50_json": reaction_json,
            "reactivity_mask_m50_p50_json": mask_json,
            "valid_reactivity_count_m50_p50": valid_101,
            "sequence_window_201": sequence_201,
            "sequence_window_201_full": full_sequence,
            "sequence_motif_5nt": motif_5,
            "drach_subtype": drach_subtype(motif_5),
            "sequence_window_21": local_21,
            "gc_fraction_21": gc_fraction_21,
            "transcript_length": len(sequence),
            "transcript_position_fraction": round(tx_pos / len(sequence), 8),
            "transcript_region": region,
            "distance_to_stop_codon_tx": stop_distance,
            "distance_to_nearest_splice_edge_tx": splice_distance,
            "icshape_abundance_rpkm": icshape[tid][1],
            "glori_gene_matches_transcript_gene": match_gene,
        }
        enriched_by_site[key].append(enriched)

    all_isoforms: list[dict] = []
    site_rows: list[dict] = []
    main_rows: list[dict] = []
    excluded_exact = set(config["exclusions"]["exact_sites"])
    for key in sorted(enriched_by_site):
        candidates = enriched_by_site[key]
        r1, r2 = rep1.get(key), rep2.get(key)
        both = r1 is not None and r2 is not None
        replicate_fields = {
            "detected_rep1": r1 is not None,
            "detected_rep2": r2 is not None,
            "detected_both_replicates": both,
            "agcov_rep1": r1["AGcov"] if r1 else "",
            "acov_rep1": r1["Acov"] if r1 else "",
            "ratio_rep1": r1["Ratio"] if r1 else "",
            "norme_ratio_rep1": r1["NormeRatio"] if r1 else "",
            "p_adjust_rep1": r1["P_adjust"] if r1 else "",
            "agcov_rep2": r2["AGcov"] if r2 else "",
            "acov_rep2": r2["Acov"] if r2 else "",
            "ratio_rep2": r2["Ratio"] if r2 else "",
            "norme_ratio_rep2": r2["NormeRatio"] if r2 else "",
            "p_adjust_rep2": r2["P_adjust"] if r2 else "",
        }
        ranked = sorted(
            candidates,
            key=lambda row: (
                not bool(row["main_window_valid"]),
                -int(row["valid_reactivity_count_m50_p50"]),
                not bool(row["glori_gene_matches_transcript_gene"]),
                row["transcript_id"],
            ),
        )
        for rank, candidate in enumerate(ranked, start=1):
            candidate.update(replicate_fields)
            candidate["isoform_selection_rank"] = rank
            candidate["representative_isoform"] = rank == 1
            all_isoforms.append(candidate)
        selected = dict(ranked[0])
        selected["isoform_count"] = len(ranked)
        selected["site_id"] = f"{key[0]}:{key[1]}:{key[2]}"
        selected.update(replicate_fields)
        total_ag = sum(int(r["AGcov"]) for r in (r1, r2) if r)
        total_a = sum(int(r["Acov"]) for r in (r1, r2) if r)
        selected["combined_agcov"] = total_ag
        selected["combined_acov"] = total_a
        selected["combined_ratio"] = round(total_a / total_ag, 8)
        analysis_gene, source = choose_analysis_gene(selected["genes"], selected["transcript_gene_name"], selected["transcript_gene_id"])
        selected["analysis_gene"] = analysis_gene
        selected["analysis_gene_source"] = source
        exact_excluded = selected["site_id"] in excluded_exact
        selected["main_analysis_included"] = both and bool(selected["main_window_valid"]) and not exact_excluded
        if exact_excluded:
            selected["main_exclusion_reason"] = "prespecified_exact_site_exclusion"
        elif not both:
            selected["main_exclusion_reason"] = "not_detected_in_both_replicates"
        elif not bool(selected["main_window_valid"]):
            selected["main_exclusion_reason"] = "insufficient_icshape_flank10_coverage"
        else:
            selected["main_exclusion_reason"] = ""
        selected["model_dataset_included"] = bool(selected["main_analysis_included"]) and bool(selected["sequence_window_201_full"])
        site_rows.append(selected)
        if selected["main_analysis_included"]:
            main_rows.append(selected)

    if len({row["site_id"] for row in site_rows}) != len(site_rows):
        raise RuntimeError("Site-level output contains duplicate site IDs")
    if any(not row["main_window_valid"] or not row["detected_both_replicates"] for row in main_rows):
        raise RuntimeError("Main analysis inclusion invariant failed")

    gene_counts = Counter(row["analysis_gene"] for row in main_rows)
    gene_rows = [
        {"analysis_gene": gene, "main_analysis_sites": count}
        for gene, count in sorted(gene_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    exclusion_counts = Counter(row["main_exclusion_reason"] or "included" for row in site_rows)
    attrition = [
        {"step": 1, "stage": "GLORI union", "sites": len(set(rep1) | set(rep2)), "retained_from_previous": "", "rule": "Exact genomic-site union"},
        {"step": 2, "stage": "Matched to icSHAPE transcript", "sites": len(site_rows), "retained_from_previous": round(len(site_rows) / len(set(rep1) | set(rep2)), 8), "rule": "At least one exact transcript-oriented A mapping"},
        {"step": 3, "stage": "Detected in both GLORI replicates", "sites": sum(row["detected_both_replicates"] for row in site_rows), "retained_from_previous": round(sum(row["detected_both_replicates"] for row in site_rows) / len(site_rows), 8), "rule": "Primary replicate rule"},
        {"step": 4, "stage": "Main structure dataset", "sites": len(main_rows), "retained_from_previous": round(len(main_rows) / sum(row["detected_both_replicates"] for row in site_rows), 8), "rule": "Representative isoform has ≥70% valid values on each ±10 side"},
        {"step": 5, "stage": "Model-ready 201-nt sequence", "sites": sum(row["model_dataset_included"] for row in site_rows), "retained_from_previous": round(sum(row["model_dataset_included"] for row in site_rows) / len(main_rows), 8), "rule": "Full unpadded 201-nt transcript sequence"},
    ]

    qc_rows = [
        {"category": "dataset", "metric": "site_level_rows", "value": len(site_rows), "n": len(site_rows), "notes": "One selected representative transcript per mapped genomic site"},
        {"category": "dataset", "metric": "all_isoform_rows", "value": len(all_isoforms), "n": len(all_isoforms), "notes": "All mapping-compatible transcript isoforms"},
        {"category": "dataset", "metric": "main_analysis_rows", "value": len(main_rows), "n": len(main_rows), "notes": "Both replicates plus valid ±10 structure window"},
        {"category": "dataset", "metric": "model_ready_rows", "value": sum(row["model_dataset_included"] for row in site_rows), "n": len(main_rows), "notes": "Main rows with full 201-nt sequence"},
        {"category": "dataset", "metric": "main_analysis_genes", "value": len(gene_counts), "n": len(main_rows), "notes": "Analysis gene assignment"},
        {"category": "isoforms", "metric": "sites_with_multiple_isoforms", "value": sum(int(row["isoform_count"]) > 1 for row in site_rows), "n": len(site_rows), "notes": "Multiplicity preserved in all-isoform table"},
        {"category": "isoforms", "metric": "maximum_isoform_count", "value": max(int(row["isoform_count"]) for row in site_rows), "n": len(site_rows), "notes": "Per genomic site"},
        {"category": "gene_clustering", "metric": "maximum_sites_per_gene", "value": max(gene_counts.values()), "n": len(gene_counts), "notes": "Main dataset"},
    ]
    for reason, count in sorted(exclusion_counts.items()):
        qc_rows.append({"category": "main_inclusion", "metric": reason, "value": count, "n": len(site_rows), "notes": "Representative site-level rows"})
    qc_rows.extend(correlation_rows(main_rows))

    FINAL.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_gzip_csv(FINAL / "06_hek293t_site_level_dataset.csv.gz", site_rows)
    write_gzip_csv(FINAL / "06_hek293t_main_analysis_dataset.csv.gz", main_rows)
    write_gzip_csv(FINAL / "06_hek293t_all_isoform_mappings.csv.gz", all_isoforms)
    write_csv(RESULTS / "06_hek293t_attrition_summary.csv", attrition)
    write_csv(RESULTS / "06_hek293t_dataset_qc_summary.csv", qc_rows)
    write_csv(RESULTS / "06_hek293t_gene_site_counts.csv", gene_rows)
    status = "PASS" if len(main_rows) >= 500 and len(site_rows) == len(enriched_by_site) else "FAIL"
    write_csv(
        RESULTS / "06_hek293t_build_analysis_dataset_run_status.csv",
        [{"status": status, "started_at": started, "finished_at": datetime.now().isoformat(timespec="seconds"), "site_level_rows": len(site_rows), "main_analysis_rows": len(main_rows), "model_ready_rows": sum(row["model_dataset_included"] for row in site_rows), "analysis_genes": len(gene_counts), "config": str(CONFIG_PATH.relative_to(ROOT)).replace("\\", "/")}],
    )
    logging.info("stage 06 finished: status=%s sites=%d main=%d genes=%d exclusions=%s", status, len(site_rows), len(main_rows), len(gene_counts), dict(exclusion_counts))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
