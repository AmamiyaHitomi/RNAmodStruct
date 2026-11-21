"""Exploratory cross-study SAC-seq/icSHAPE sensitivity in HepG2 and HeLa."""

from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from p1_coverage_and_within_gene import GENE, PREDICTOR, ROOT, coefficient_row
from pipeline_common import ExonBinIndex, coverage, genome_to_tx, load_selected_fasta
from replay_presubmission_stage4_models import sha256


RAW = ROOT / "data/raw_processed"
TAR = RAW / "GSE145805_RAW.tar"
ASSEMBLY = ROOT / "data/reference/ensembl/release-88_GRCh38"
GTF = ASSEMBLY / "Homo_sapiens.GRCh38.88.gtf.gz"
FASTA = ASSEMBLY / "Homo_sapiens.GRCh38.cdna.all.fa.gz"
OUTCOME = "sac_fraction"


def map_cells(bed, member, mapping, features):
    fields = ["chrom", "start", "end", "site_description", "mutation_ratio_pct", "strand", "calibrated_fraction_pct"]
    points = pd.read_csv(bed, sep="\t", header=None, names=fields)
    keys = ["chrom", "end", "strand"]
    if (not points.end.sub(points.start).eq(1).all() or points.duplicated(keys).any()
            or not points.strand.isin(["+", "-"]).all() or not points.calibrated_fraction_pct.between(0, 100).all()):
        raise ValueError("Expected unique stranded, single-base GRCh38 BED rows and calibrated percentages")
    lengths = features.read_structure_lengths_tar(TAR, member)
    transcripts = features.exact_length_transcripts(GTF, lengths, versioned=True)
    index = ExonBinIndex(transcripts)
    hits = {(r.chrom, r.end, r.strand): index.query(r.chrom, r.end, r.strand) for r in points.itertuples()}
    wanted = {e.transcript_id for exons in hits.values() for e in exons}
    sequences = load_selected_fasta(FASTA, wanted)
    annotations = mapping.load_coding_annotations(GTF, wanted, transcripts)
    tracks = features.load_selected_tar(TAR, member, wanted)
    counts = dict(input_bed_sites=len(points), sites_overlapping_length_matched_transcript=0,
                  sites_confirmed_transcript_A=0, sites_sufficient_flanks=0)
    rows = []
    for site in points.itertuples():
        exons = hits[(site.chrom, site.end, site.strand)]
        counts["sites_overlapping_length_matched_transcript"] += bool(exons)
        adenosine_count = 0
        candidates = []
        for exon in exons:
            tid = exon.transcript_id
            seq = sequences.get(tid)
            center = genome_to_tx(exon, site.end)
            if seq is None or len(seq) != lengths[tid] or seq[center-1] != "A":
                continue
            adenosine_count += 1
            _, _, vector = tracks[tid]
            left, right, both = coverage(vector, center, 10)
            if not (left >= .70 and right >= .70):
                continue
            values = np.asarray(vector, dtype=float)
            window = values[max(0, center-51):center+50]
            flanks = np.r_[values[center-11:center-1], values[center:center+10]]
            candidates.append((-int(np.isfinite(window).sum()), tid, center, exon.gene_id, both, np.nanmean(flanks)))
        counts["sites_confirmed_transcript_A"] += bool(adenosine_count)
        if not candidates:
            continue
        counts["sites_sufficient_flanks"] += 1
        _, tid, center, gene, cov, x = min(candidates)
        seq = sequences[tid]
        annotation = annotations[tid]
        stop = annotation["stop_codon_start_tx"]
        edges = [edge for exon in transcripts[tid][:-1] for edge in (exon.tx_end, exon.tx_end + 1)]
        local = seq[center-11:center+10]
        rows.append({
            "chrom": site.chrom, "pos_1based": site.end, "strand": site.strand, "transcript_id": tid,
            "transcript_pos_1based": center, GENE: gene, OUTCOME: site.calibrated_fraction_pct / 100,
            PREDICTOR: x, "coverage_flank10": cov, "gc_fraction_21": (local.count("G") + local.count("C")) / 21,
            "transcript_position_fraction": center / len(seq), "transcript_region": mapping.transcript_region(center, annotation),
            "drach_subtype": mapping.drach_subtype(seq[center-3:center+2]),
            "distance_to_stop_codon_tx": np.nan if stop is None else center-stop,
            "distance_to_nearest_splice_edge_tx": min((abs(center-edge) for edge in edges), default=np.nan),
        })
    return pd.DataFrame(rows), counts


def fit_sites(frame, sd, method):
    covariates = frame[[PREDICTOR]].copy()
    if method != "unadjusted":
        for c in ("gc_fraction_21", "transcript_position_fraction", "coverage_flank10",
                  "distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx"):
            value = frame[c]
            if value.isna().any():
                covariates[f"{c}_missing"] = value.isna().astype(float)
            covariates[c] = value.fillna(value.median())
        covariates = covariates.join(pd.get_dummies(frame[["transcript_region", "drach_subtype"]],
                                                    drop_first=True, dtype=float))
    covariates = sm.add_constant(covariates.loc[:, covariates.nunique().gt(1)], has_constant="add")
    model = sm.OLS(frame[OUTCOME], covariates).fit(cov_type="cluster", cov_kwds={"groups": frame[GENE]})
    return coefficient_row(model, frame, sd, analysis=method)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hepg2_bed", type=Path)
    parser.add_argument("output", type=Path, help="Fresh output directory")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    mapping = importlib.import_module("09_build_hela_dataset")
    features = importlib.import_module("25_extract_multimodification_structure_features")
    inputs = (
        ("HepG2", args.hepg2_bed, "GSM4333262_HepG2.out.txt.gz"),
        ("HeLa", RAW / "GSE162356_HeLa_polyA.bed.gz", "GSM4333258_HeLa.out.txt.gz"),
    )
    counts, estimates = [], []
    hashes = {path.name: sha256(path) for path in (TAR, GTF, FASTA)}
    for cell, bed, member in inputs:
        print(f"{cell}: stranded mapping, transcript A check, and paired-flank coverage", flush=True)
        frame, attrition = map_cells(bed, member, mapping, features)
        counts.append(dict(cell_line=cell, icshape_member=member, **attrition))
        sd = frame[PREDICTOR].std(ddof=1)
        for subset, label in ((frame, "unadjusted"), (frame, "adjusted"),
                              (frame.loc[frame.coverage_flank10.eq(1)], "adjusted_complete_flanks")):
            estimates.append({"cell_line": cell, **fit_sites(subset, sd, label)})
        frame.to_csv(args.output / f"{cell}_mapped_sites.csv", index=False, mode="x")
        hashes[bed.name] = sha256(bed)
    pd.DataFrame(counts).to_csv(args.output / "mapping_attrition.csv", index=False, mode="x")
    pd.DataFrame(estimates).to_csv(args.output / "cross_study_associations.csv", index=False, mode="x")
    (args.output / "input_sha256.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
