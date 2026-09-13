"""Stage 26: technology-stratified association analyses for phase 3."""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import io
import logging
import math
import sys
import tarfile
from array import array
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml
from scipy.stats import spearmanr
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pipeline_common import load_selected_fasta, parse_icshape_line  # noqa: E402

CONFIG = ROOT / "config" / "26_stratified_multimodification_associations.yaml"


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
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    scale = float(values.std(ddof=0))
    return (values - float(values.mean())) / (scale if scale > 0 else 1.0)


def fit_clustered_continuous(df: pd.DataFrame, track_id: str, outcome: str, predictor: str, covariates: list[str], cluster: str, evidence_tier: str) -> tuple[dict, dict]:
    columns = [outcome, predictor, cluster, *covariates]
    model_df = df[columns].copy().replace([np.inf, -np.inf], np.nan).dropna()
    design = pd.DataFrame({"Intercept": 1.0, f"z_{predictor}": zscore(model_df[predictor])}, index=model_df.index)
    for name in covariates:
        design[f"z_{name}"] = zscore(model_df[name])
    y = pd.to_numeric(model_df[outcome]).to_numpy(float)
    ordinary = sm.OLS(y, design).fit()
    fallback_groups = pd.Series(model_df.index.astype(str), index=model_df.index)
    groups = model_df[cluster].astype("string").where(model_df[cluster].notna(), fallback_groups).astype(str)
    clustered = sm.OLS(y, design).fit(cov_type="cluster", cov_kwds={"groups": groups, "use_correction": True})
    term = f"z_{predictor}"
    beta, se = float(clustered.params[term]), float(clustered.bse[term])
    outcome_sd = float(np.std(y, ddof=1))
    rho = spearmanr(pd.to_numeric(model_df[predictor]), y)
    influence = ordinary.get_influence()
    cooks = influence.cooks_distance[0]
    bp = het_breuschpagan(ordinary.resid, ordinary.model.exog)
    result = {
        "track_id": track_id, "modification": str(df.modification.iloc[0]), "cell_line": str(df.cell_line.iloc[0]),
        "evidence_tier": evidence_tier, "method": "adjusted_OLS_cluster_robust",
        "outcome": outcome, "predictor": predictor, "n": len(model_df), "clusters": groups.nunique(),
        "effect_per_predictor_sd": beta, "se_cluster": se, "ci95_low": beta - 1.96 * se,
        "ci95_high": beta + 1.96 * se, "p_value": float(clustered.pvalues[term]),
        "fully_standardized_effect": beta / outcome_sd if outcome_sd > 0 else math.nan,
        "spearman_rho_unadjusted": float(rho.statistic), "spearman_p_unadjusted": float(rho.pvalue),
        "q_value_primary": "",
    }
    diagnostics = {
        "track_id": track_id, "n": len(model_df), "clusters": groups.nunique(),
        "r_squared_ordinary": float(ordinary.rsquared), "adjusted_r_squared_ordinary": float(ordinary.rsquared_adj),
        "breusch_pagan_p_value": float(bp[1]), "max_cooks_distance": float(np.max(cooks)),
        "influential_cooks_gt_4_over_n": int(np.sum(cooks > 4 / len(model_df))),
        "design_rank": int(np.linalg.matrix_rank(design.to_numpy(float))), "design_columns": design.shape[1],
    }
    return result, diagnostics


def load_selected_structure_gzip(path: Path, wanted: set[str]) -> dict[str, tuple[int, str, array]]:
    result = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            tid = line.split("\t", 1)[0]
            if tid in wanted:
                parsed = parse_icshape_line(line)
                result[tid] = parsed[1:]
    return result


def load_selected_structure_tar(path: Path, member: str, wanted: set[str]) -> dict[str, tuple[int, str, array]]:
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


def valid_control_positions(sequence: str, values: array, flank: int, min_fraction: float, excluded: set[int]) -> dict[str, list[int]]:
    valid = {base: [] for base in "ACGT"}
    mask = np.fromiter((not math.isnan(value) for value in values), dtype=np.int32)
    prefix = np.concatenate(([0], np.cumsum(mask)))
    required = math.ceil(flank * min_fraction)
    for pos1 in range(flank + 1, len(sequence) - flank + 1):
        if pos1 in excluded:
            continue
        left = int(prefix[pos1 - 1] - prefix[pos1 - flank - 1])
        right = int(prefix[pos1 + flank] - prefix[pos1])
        base = sequence[pos1 - 1]
        if base in valid and left >= required and right >= required:
            valid[base].append(pos1)
    return valid


def nearest_unused(candidates: list[int], pos1: int, minimum_distance: int, used: set[int]) -> int | None:
    pivot = bisect.bisect_left(candidates, pos1)
    left, right = pivot - 1, pivot
    while left >= 0 or right < len(candidates):
        options = []
        if left >= 0:
            options.append(candidates[left])
        if right < len(candidates):
            options.append(candidates[right])
        options.sort(key=lambda value: (abs(value - pos1), value))
        for value in options:
            if abs(value - pos1) >= minimum_distance and value not in used:
                return value
        left -= 1
        right += 1
    return None


def flank_mean(values: array, pos1: int, flank: int) -> float:
    selected = list(values[pos1 - flank - 1:pos1 - 1]) + list(values[pos1:pos1 + flank])
    kept = [float(value) for value in selected if not math.isnan(value)]
    return float(np.mean(kept))


def build_nm_pairs(df: pd.DataFrame, sequences: dict[str, str], structures: dict[str, tuple[int, str, array]], track_id: str, flank: int, min_fraction: float, minimum_distance: int) -> list[dict]:
    positives = df[df.main_analysis_included.astype(bool)].copy()
    excluded_by_tid = {tid: set(group.transcript_pos_1based.astype(int)) for tid, group in positives.groupby("transcript_id")}
    candidates_by_tid = {}
    for tid in sorted(excluded_by_tid):
        if tid not in sequences or tid not in structures:
            continue
        candidates_by_tid[tid] = valid_control_positions(sequences[tid], structures[tid][2], flank, min_fraction, excluded_by_tid[tid])
    used_by_tid: dict[str, set[int]] = {tid: set() for tid in candidates_by_tid}
    pairs = []
    for source in positives.sort_values(["transcript_id", "transcript_pos_1based", "site_id_hg38"]).itertuples(index=False):
        tid, pos1 = source.transcript_id, int(source.transcript_pos_1based)
        if tid not in candidates_by_tid:
            continue
        base = sequences[tid][pos1 - 1]
        control = nearest_unused(candidates_by_tid[tid].get(base, []), pos1, minimum_distance, used_by_tid[tid])
        if control is None:
            continue
        used_by_tid[tid].add(control)
        control_mean = flank_mean(structures[tid][2], control, flank)
        positive_mean = float(source.reactivity_mean_flank10)
        pairs.append({
            "track_id": track_id, "cell_line": source.cell_line, "transcript_id": tid,
            "transcript_gene_id": source.transcript_gene_id, "site_id_hg38": source.site_id_hg38,
            "modified_tx_pos_1based": pos1, "control_tx_pos_1based": control,
            "matched_base": base, "absolute_distance_nt": abs(control - pos1),
            "modified_reactivity_flank10": positive_mean, "control_reactivity_flank10": control_mean,
            "paired_difference_modified_minus_control": positive_mean - control_mean,
        })
    return pairs


def clustered_paired_result(pairs: list[dict], evidence_tier: str, bootstrap_n: int, permutation_n: int, seed: int) -> tuple[dict, dict]:
    frame = pd.DataFrame(pairs)
    diffs = frame.paired_difference_modified_minus_control.to_numpy(float)
    grouped = frame.groupby("transcript_id", sort=True).paired_difference_modified_minus_control.agg(["sum", "count"])
    sums, counts = grouped["sum"].to_numpy(float), grouped["count"].to_numpy(int)
    rng = np.random.default_rng(seed)
    boot = np.empty(bootstrap_n)
    for index in range(bootstrap_n):
        draw = rng.integers(0, len(grouped), len(grouped))
        boot[index] = sums[draw].sum() / counts[draw].sum()
    extreme = 0
    batch = 500
    observed = abs(float(np.mean(diffs)))
    for start in range(0, permutation_n, batch):
        size = min(batch, permutation_n - start)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(size, len(grouped)))
        permuted = np.abs((signs * sums).sum(axis=1) / counts.sum())
        extreme += int(np.sum(permuted >= observed))
    p_value = (extreme + 1) / (permutation_n + 1)
    mean_diff = float(np.mean(diffs))
    sd_diff = float(np.std(diffs, ddof=1))
    result = {
        "track_id": str(frame.track_id.iloc[0]), "modification": "Nm", "cell_line": str(frame.cell_line.iloc[0]),
        "evidence_tier": evidence_tier, "method": "same_transcript_same_base_cluster_sign_permutation",
        "outcome": "modified_minus_matched_control_reactivity", "predictor": "Nm_site_status",
        "n": len(frame), "clusters": len(grouped), "effect_per_predictor_sd": mean_diff,
        "se_cluster": "", "ci95_low": float(np.quantile(boot, 0.025)), "ci95_high": float(np.quantile(boot, 0.975)),
        "p_value": p_value, "fully_standardized_effect": mean_diff / sd_diff if sd_diff > 0 else math.nan,
        "spearman_rho_unadjusted": "", "spearman_p_unadjusted": "", "q_value_primary": "",
    }
    qc = {
        "track_id": str(frame.track_id.iloc[0]), "eligible_positive_sites": len(frame),
        "matched_pairs": len(frame), "unique_transcripts": len(grouped),
        "median_absolute_match_distance_nt": float(frame.absolute_distance_nt.median()),
        "maximum_absolute_match_distance_nt": int(frame.absolute_distance_nt.max()),
        "mean_modified_reactivity": float(frame.modified_reactivity_flank10.mean()),
        "mean_control_reactivity": float(frame.control_reactivity_flank10.mean()),
        "cluster_bootstrap_replicates": bootstrap_n, "cluster_sign_permutations": permutation_n,
    }
    return result, qc


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    out = config["outputs"]
    setup_logging(ROOT / out["log"])
    frames = {key: pd.read_csv(ROOT / path) for key, path in config["inputs"].items()}
    for key, frame in frames.items():
        frames[key] = frame[frame.main_analysis_included.astype(bool)].copy()
    logging.info("main-analysis rows: %s", {key: len(value) for key, value in frames.items()})

    m5c = frames["m5c_hela"]
    m5c["mean_log1p_coverage"] = np.log1p(m5c[["coverage_B", "coverage_C", "coverage_E"]].mean(axis=1))
    m5c["log1p_transcript_length"] = np.log1p(m5c.transcript_length)
    m7g = frames["m7g_hela"]
    m7g["log2_fold_enrichment"] = np.log2(m7g.fold_enrichment)
    m7g["log2_interval_width"] = np.log2(m7g.interval_width_nt)
    m7g["peak_position_fraction"] = m7g.peak_midpoint_tx_1based / m7g.transcript_length
    m7g["log1p_transcript_length"] = np.log1p(m7g.transcript_length)

    associations, diagnostics = [], []
    for key, frame in (("m5c_hela", m5c), ("m7g_hela", m7g)):
        spec = config["continuous_models"][key]
        result, diagnostic = fit_clustered_continuous(frame, key, spec["outcome"], spec["predictor"], spec["covariates"], spec["cluster"], spec["evidence_tier"])
        associations.append(result)
        diagnostics.append(diagnostic)

    refs = config["references"]
    nm_frames = {"nm_hela": frames["nm_hela"], "nm_hek": frames["nm_hek"]}
    hela_tids, hek_tids = set(nm_frames["nm_hela"].transcript_id), set(nm_frames["nm_hek"].transcript_id)
    hela_sequences = load_selected_fasta(ROOT / refs["hela_cdna"], hela_tids)
    hek_sequences = load_selected_fasta(ROOT / refs["hek_cdna"], hek_tids)
    hela_structures = load_selected_structure_tar(ROOT / refs["hela_structure_tar"], refs["hela_structure_member"], hela_tids)
    hek_structures = load_selected_structure_gzip(ROOT / refs["hek_structure"], hek_tids)
    matching = config["nm_matching"]
    all_pairs, matching_qc = [], []
    nm_inputs = [
        ("nm_hela", nm_frames["nm_hela"], hela_sequences, hela_structures, matching["hela_evidence_tier"], config["random_seed"]),
        ("nm_hek", nm_frames["nm_hek"], hek_sequences, hek_structures, matching["hek_evidence_tier"], config["random_seed"] + 1),
    ]
    for key, frame, sequences, structures, tier, seed in nm_inputs:
        pairs = build_nm_pairs(frame, sequences, structures, key, int(matching["flank_each_side_nt"]), float(matching["minimum_valid_fraction_each_side"]), int(matching["minimum_distance_nt"]))
        if len(pairs) < config["admission"]["minimum_matched_pairs"]:
            raise RuntimeError(f"Insufficient matched controls for {key}: {len(pairs)}")
        result, qc = clustered_paired_result(pairs, tier, int(matching["cluster_bootstrap_replicates"]), int(matching["cluster_sign_permutations"]), int(seed))
        associations.append(result)
        matching_qc.append({**qc, "eligible_positive_sites": len(frame), "matching_fraction": round(len(pairs) / len(frame), 8)})
        all_pairs.extend(pairs)
        logging.info("%s matched pairs=%d/%d", key, len(pairs), len(frame))

    primary_indexes = [index for index, row in enumerate(associations) if row["evidence_tier"] == "primary"]
    adjusted = multipletests([associations[index]["p_value"] for index in primary_indexes], method=config["multiple_testing"]["method"])[1]
    for index, q_value in zip(primary_indexes, adjusted):
        associations[index]["q_value_primary"] = float(q_value)

    replicate_rows = []
    for outcome in ["methy_rate_B", "methy_rate_C", "methy_rate_E"]:
        result, _ = fit_clustered_continuous(m5c, "m5c_hela", outcome, "reactivity_mean_flank10", config["continuous_models"]["m5c_hela"]["covariates"], "transcript_gene_id", "replicate_sensitivity")
        replicate_rows.append(result)
    replicate_q = multipletests([row["p_value"] for row in replicate_rows], method="fdr_bh")[1]
    for row, q_value in zip(replicate_rows, replicate_q):
        row["q_value_within_replicates"] = float(q_value)

    low = frames["m7g_hek293t"]
    descriptive = [{
        "track_id": "m7g_hek293t", "modification": "m7G", "cell_line": "HEK293T", "n": len(low),
        "analysis_role": "DESCRIPTIVE_ONLY_LOW_N", "median_fold_enrichment": float(low.fold_enrichment.median()),
        "median_reactivity_mean_peak": float(low.reactivity_mean_peak.median()),
        "spearman_rho": float(spearmanr(low.fold_enrichment, low.reactivity_mean_peak).statistic),
        "inferential_p_value_reported": False,
    }]

    write_csv(ROOT / out["primary_associations"], associations)
    write_csv(ROOT / out["m5c_replicate_sensitivity"], replicate_rows)
    write_gzip_csv(ROOT / out["nm_matched_pairs"], all_pairs)
    write_csv(ROOT / out["nm_matching_qc"], matching_qc)
    write_csv(ROOT / out["model_diagnostics"], diagnostics)
    write_csv(ROOT / out["descriptive_tracks"], descriptive)

    primary_complete = sum(row["evidence_tier"] == "primary" and math.isfinite(float(row["p_value"])) for row in associations)
    gate = "GO_TO_PHASE3_SYNTHESIS" if primary_complete >= config["admission"]["minimum_primary_analyses_completed"] else "NO_GO"
    status = "PASS" if gate.startswith("GO") else "FAIL"
    now = datetime.now(timezone.utc)
    write_csv(ROOT / out["status"], [{
        "stage": 26, "status": status, "primary_analyses_completed": primary_complete,
        "primary_tests_fdr_family_size": len(primary_indexes), "nm_matched_pairs": len(all_pairs),
        "phase3_gate": gate, "raw_scales_pooled": False, "finished_utc": now.isoformat(),
    }])
    manifest = []
    for key in ["primary_associations", "m5c_replicate_sensitivity", "nm_matched_pairs", "nm_matching_qc", "model_diagnostics", "descriptive_tracks", "status"]:
        path = ROOT / out[key]
        manifest.append({"artifact": key, "path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    write_csv(ROOT / out["manifest"], manifest)

    primary = {row["track_id"]: row for row in associations}
    report = f"""# Stage 26: Multi-modification hierarchical correlation analysis

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: {now.date().isoformat()}
- Verification Status: {'VERIFIED' if status == 'PASS' else 'UNVERIFIED'}
- Version Label: phase3_stage26_v1

## Experiment Result

- **ID**: phase3_stage26
- **Type**: analysis
- **Status**: {'completed' if status == 'PASS' else 'failed_gate'}
- **Command**: `E:\\ancd\\envs\\my_pytorch\\python.exe src\\26_run_stratified_multimodification_associations.py`
- **Working Directory**: `E:\\my_python\\RNAmodStruct`
- **Output Gate**: `{gate}`

The three HeLa master tests form the same BH-FDR family; HEK Nm is approximate cell line sensitivity, HEK293T m7G is descriptive only because n={len(low)}. None of the original modification dimensions were pooled across technologies.

## Hierarchical estimation

- **m5C/HeLa**: ±10 nt mean icSHAPE per 1 SD increase in mean m5C proportion change {primary['m5c_hela']['effect_per_predictor_sd']:.6f} (95% CI {primary['m5c_hela']['ci95_low']:.6f} to {primary['m5c_hela']['ci95_high']:.6f}; q={primary['m5c_hela']['q_value_primary']:.6g}; n={primary['m5c_hela']['n']}).
- **m7G/HeLa**: log2 peak enrichment change {primary['m7g_hela']['effect_per_predictor_sd']:.6f} per 1 SD increase in intra-peak mean icSHAPE (95% CI {primary['m7g_hela']['ci95_low']:.6f} to {primary['m7g_hela']['ci95_high']:.6f}; q={primary['m7g_hela']['q_value_primary']:.6g}; n={primary['m7g_hela']['n']}). This is an interval-level association, not a single-base effect.
- **Nm/HeLa**: The structural difference between the modified site and the same transcript and the same base control is {primary['nm_hela']['effect_per_predictor_sd']:.6f} (cluster bootstrap 95% CI {primary['nm_hela']['ci95_low']:.6f} to {primary['nm_hela']['ci95_high']:.6f}; q={primary['nm_hela']['q_value_primary']:.6g}; paired n={primary['nm_hela']['n']}).
- **Nm/HEK**: Approximate cell line sensitivity difference {primary['nm_hek']['effect_per_predictor_sd']:.6f} (95% CI {primary['nm_hek']['ci95_low']:.6f} to {primary['nm_hek']['ci95_high']:.6f}; not included in the HeLa master test FDR family).

## Explain boundaries

These are observational associations across studies. The m5C and m7G models adjust coverage/interval width and transcript position respectively, and robustly estimate clustering by gene; Nm control matching reduces transcript and base composition confounding, but cannot exclude residual confounding by expression, detection efficiency, and local sequence environment. Statistical significance does not constitute evidence of causation or mechanism."""
    (ROOT / out["report"]).write_text(report, encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
