"""Build analysis-ready HeLa site and isoform tables from GLORI x icSHAPE join (GRCh38 direct)."""

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
    ExonBinIndex,
    coverage,
    fetch_fasta_bases,
    genome_to_tx,
    load_gtf_exons,
    load_icshape_tar,
    load_selected_fasta,
    norm_chr,
    parse_gtf_attributes,
    read_glori,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_processed"
REF = ROOT / "data" / "reference"
INTERIM = ROOT / "data" / "interim"
FINAL = ROOT / "data" / "final"
RESULTS = ROOT / "results"
CONFIG_PATH = ROOT / "config" / "09_hela_analysis_config.yaml"
LOG = RESULTS / "logs" / "09_build_hela_dataset.log"

HELA_TAR = "GSE145805_RAW.tar"
HELA_ICSHAPE_MEMBER = "GSM4333258_HeLa.out.txt.gz"


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
    rows = {}
    for row in read_glori(path):
        key = (norm_chr(row["Chr"]), int(row["Sites"]), row["Strand"])
        if key in rows:
            raise ValueError(f"Duplicate GLORI genomic key in {path.name}: {key}")
        rows[key] = row
    return rows


def safe_median(values: list[float]) -> str | float:
    return statistics.median(values) if values else ""


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


def versioned_tid(attrs: dict[str, str]) -> str:
    tid = attrs.get("transcript_id", "")
    if attrs.get("transcript_version"):
        return f"{tid}.{attrs['transcript_version']}"
    return tid


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
            tid = versioned_tid(parse_gtf_attributes(fields[8]))
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


def abundance_or_empty(abundance: str) -> str:
    return "" if abundance == "*" else abundance


def main() -> None:
    setup_logging()
    started = datetime.now().isoformat(timespec="seconds")
    with CONFIG_PATH.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    flank10 = int(config["structure"]["main_window_each_side_nt"])
    min_coverage = float(config["structure"]["minimum_valid_fraction_each_side"])
    flank50 = int(config["structure"]["display_window_each_side_nt"])
    model_flank = (int(config["sequence"]["model_window_nt"]) - 1) // 2

    rep1 = load_glori(RAW / "GSM6432595_Hela-1_35bp_m2.totalm6A.FDR.csv.gz")
    rep2 = load_glori(RAW / "GSM6432596_Hela-2_35bp_m2.totalm6A.FDR.csv.gz")
    logging.info("loaded GLORI replicates: %d and %d sites", len(rep1), len(rep2))
    keys1, keys2 = set(rep1), set(rep2)
    union_keys = sorted(keys1 | keys2)

    icshape = load_icshape_tar(RAW / HELA_TAR, HELA_ICSHAPE_MEMBER)
    gtf_path = REF / "ensembl" / "release-88_GRCh38" / "Homo_sapiens.GRCh38.88.gtf.gz"
    cdna_path = REF / "ensembl" / "release-88_GRCh38" / "Homo_sapiens.GRCh38.cdna.all.fa.gz"
    genome38_path = REF / "ensembl" / "release-88_GRCh38" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    logging.info("loading HeLa icSHAPE and GRCh38.88 versioned annotation")
    ids = set(icshape)
    transcripts_all = load_gtf_exons(gtf_path, ids, versioned=True)
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
    status_rows: list[dict] = []
    intersection_rows: list[dict] = []
    reasons = Counter()
    matched_site_queries: dict[str, set[int]] = defaultdict(set)
    matched_genes = set()

    logging.info("mapping %d GLORI union sites to GRCh38 transcripts", len(union_keys))
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
        exon_hits = exon_index.query(chrom38, pos38, strand38)
        if not exon_hits:
            status_rows.append({**base_status, "mapping_status": "no_icshape_exon", "transcript_hits": 0, "valid_main_windows": 0})
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
            tx_pos = genome_to_tx(exon, pos38)
            seq = sequences[tid]
            if not 1 <= tx_pos <= len(seq) or seq[tx_pos - 1] != "A":
                non_a_hits += 1
                continue
            valid_hits += 1
            values = icshape[tid][2]
            center_score = values[tx_pos - 1]
            up10, down10, flank10_cov = coverage(values, tx_pos, 10)
            _, _, flank50cov = coverage(values, tx_pos, 50)
            window_valid = (
                not math.isnan(up10)
                and not math.isnan(down10)
                and up10 >= min_coverage
                and down10 >= min_coverage
            )
            valid_windows += int(window_valid)
            intersection_rows.append(
                {
                    **base_status,
                    "transcript_id": tid,
                    "transcript_pos_1based": tx_pos,
                    "transcript_gene_id": exon.gene_id,
                    "transcript_gene_name": exon.gene_name,
                    "transcript_base": "A",
                    "icshape_center": "" if math.isnan(center_score) else round(float(center_score), 6),
                    "coverage_up10": up10,
                    "coverage_down10": down10,
                    "coverage_flank10": flank10_cov,
                    "coverage_flank50": flank50cov,
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
    mapped_coverages = [float(r["coverage_flank10"]) for r in intersection_rows if r["coverage_flank10"] != "" and not math.isnan(float(r["coverage_flank10"]))]
    valid_mapping_rows = [r for r in intersection_rows if r["main_window_valid"]]
    center_available_rows = [r for r in intersection_rows if r["icshape_center"] != ""]

    summary_rows = [
        {"section": "glori_replicates", "metric": "rep1_unique_sites", "value": len(keys1), "notes": "Exact (Chr, Sites, Strand) keys"},
        {"section": "glori_replicates", "metric": "rep2_unique_sites", "value": len(keys2), "notes": "Exact (Chr, Sites, Strand) keys"},
        {"section": "glori_replicates", "metric": "replicate_union_sites", "value": len(union_keys), "notes": "Exact genomic-site union"},
        {"section": "coordinate_mapping", "metric": "matched_union_sites", "value": len(matched_sites), "notes": "At least one exact HeLa icSHAPE transcript position"},
        {"section": "coordinate_mapping", "metric": "matched_transcript_rows", "value": len(intersection_rows), "notes": "One genomic site can map to multiple isoforms"},
        {"section": "coordinate_mapping", "metric": "independent_glori_genes", "value": len(matched_genes), "notes": "Non-ELSE GLORI gene labels among matched sites"},
        {"section": "structure_coverage", "metric": "sites_with_valid_flank10_window", "value": len(valid_window_sites), "notes": f"Both ±10 sides separately ≥ {min_coverage:.0%}; center excluded"},
        {"section": "structure_coverage", "metric": "site_valid_window_fraction", "value": round(len(valid_window_sites) / len(matched_sites), 8) if matched_sites else "", "notes": "At least one qualifying transcript window per genomic site"},
        {"section": "structure_coverage", "metric": "valid_transcript_windows", "value": len(valid_mapping_rows), "notes": "Transcript-level qualifying windows"},
        {"section": "structure_coverage", "metric": "transcript_window_valid_fraction", "value": round(len(valid_mapping_rows) / len(intersection_rows), 8) if intersection_rows else "", "notes": "Transcript-level denominator"},
        {"section": "structure_coverage", "metric": "center_reactivity_available_rows", "value": len(center_available_rows), "notes": "Transcript mappings with a non-missing center score"},
        {"section": "structure_coverage", "metric": "mean_flank10_coverage", "value": round(sum(mapped_coverages) / len(mapped_coverages), 8) if mapped_coverages else "", "notes": "Only full-length ±10 windows"},
        {"section": "structure_coverage", "metric": "median_flank10_coverage", "value": safe_median(mapped_coverages), "notes": "Across mapped transcript rows"},
        {"section": "validation", "metric": "grch38_center_base_failures", "value": center_base_failures, "notes": "Expected A on + and T on -"},
        {"section": "reference", "metric": "eligible_icshape_transcripts", "value": len(transcripts), "notes": "Exact length agreement in icSHAPE, cDNA and GTF"},
    ]
    for reason, count in sorted(reasons.items()):
        summary_rows.append({"section": "mapping_status", "metric": reason, "value": count, "notes": "GLORI union-site count"})
    for reason, count in sorted(exclusion.items()):
        summary_rows.append({"section": "transcript_exclusion", "metric": reason, "value": count, "notes": "icSHAPE transcript count"})

    # ---- Enrichment and analysis-dataset build (stage-06 equivalent for HeLa) ----
    base_rows = intersection_rows
    transcript_ids = {row["transcript_id"] for row in base_rows}
    coding_annotations = load_coding_annotations(gtf_path, transcript_ids, transcripts)
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
        junction_edges = [edge for exon in transcripts[tid][:-1] for edge in (exon.tx_end, exon.tx_end + 1)]
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
            "icshape_abundance_rpkm": abundance_or_empty(icshape[tid][1]),
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
    union_total = len(union_keys)
    matched_total = len(site_rows)
    both_total = sum(row["detected_both_replicates"] for row in site_rows)
    attrition = [
        {"step": 1, "stage": "GLORI union", "sites": union_total, "retained_from_previous": "", "rule": "Exact genomic-site union"},
        {"step": 2, "stage": "Matched to HeLa icSHAPE transcript", "sites": matched_total, "retained_from_previous": round(matched_total / union_total, 8), "rule": "At least one exact transcript-oriented A mapping"},
        {"step": 3, "stage": "Detected in both GLORI replicates", "sites": both_total, "retained_from_previous": round(both_total / matched_total, 8), "rule": "Primary replicate rule"},
        {"step": 4, "stage": "Main structure dataset", "sites": len(main_rows), "retained_from_previous": round(len(main_rows) / both_total, 8), "rule": "Representative isoform has ≥70% valid values on each ±10 side"},
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

    INTERIM.mkdir(parents=True, exist_ok=True)
    FINAL.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    status_fields = list(dict.fromkeys(key for row in status_rows for key in row))
    intersection_fields = list(dict.fromkeys(key for row in intersection_rows for key in row)) if intersection_rows else []
    write_gzip_csv(INTERIM / "09_hela_glori_site_mapping_status.csv.gz", status_rows, status_fields)
    if intersection_rows:
        write_gzip_csv(INTERIM / "09_hela_glori_icshape_intersection.csv.gz", intersection_rows, intersection_fields)
    write_gzip_csv(FINAL / "09_hela_site_level_dataset.csv.gz", site_rows)
    write_gzip_csv(FINAL / "09_hela_main_analysis_dataset.csv.gz", main_rows)
    write_gzip_csv(FINAL / "09_hela_all_isoform_mappings.csv.gz", all_isoforms)
    write_csv(RESULTS / "09_hela_initial_intersection_summary.csv", summary_rows)
    write_csv(
        RESULTS / "09_hela_unmatched_reason_counts.csv",
        [
            {"mapping_status": reason, "sites": count, "fraction_of_union": round(count / len(union_keys), 8)}
            for reason, count in sorted(reasons.items())
            if reason != "matched_transcript_position"
        ],
    )
    write_csv(RESULTS / "09_hela_attrition_summary.csv", attrition)
    write_csv(RESULTS / "09_hela_dataset_qc_summary.csv", qc_rows)
    write_csv(RESULTS / "09_hela_gene_site_counts.csv", gene_rows)
    status = "PASS" if len(main_rows) >= 500 and len(site_rows) == len(enriched_by_site) and center_base_failures == 0 else "FAIL"
    write_csv(
        RESULTS / "09_hela_build_analysis_dataset_run_status.csv",
        [{"status": status, "started_at": started, "finished_at": datetime.now().isoformat(timespec="seconds"), "site_level_rows": len(site_rows), "main_analysis_rows": len(main_rows), "model_ready_rows": sum(row["model_dataset_included"] for row in site_rows), "analysis_genes": len(gene_counts), "matched_sites": len(matched_sites), "center_base_failures": center_base_failures, "config": str(CONFIG_PATH.relative_to(ROOT)).replace("\\", "/")}],
    )
    logging.info("stage 09 finished: status=%s sites=%d main=%d genes=%d center_base_failures=%d exclusions=%s", status, len(site_rows), len(main_rows), len(gene_counts), center_base_failures, dict(exclusion_counts))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
