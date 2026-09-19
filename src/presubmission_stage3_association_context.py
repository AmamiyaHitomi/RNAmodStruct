"""Post hoc association context and selection diagnostics for main_v2."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy import build_design_matrices, dmatrix


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "presubmission_stage3"
REPORT = ROOT / "docs" / "reports" / "final" / "main_v2_phase3_association_context_zh.md"
MAIN = ROOT / "data" / "final" / "09_hela_main_analysis_dataset.csv.gz"
SITE = ROOT / "data" / "final" / "09_hela_site_level_dataset.csv.gz"
PAIR = ROOT / "data" / "final" / "19_hek293t_invivo_invitro_paired.csv.gz"
STAGE19 = ROOT / "results" / "tables" / "19_paired_condition_associations.csv"
STAGE09 = ROOT / "results" / "tables" / "09b_hela_primary_association.csv"
STAGE09_CODE = ROOT / "src" / "09b_run_hela_association.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("x", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows(rows)


def load_design_builder():
    spec = importlib.util.spec_from_file_location("stage09b_frozen", STAGE09_CODE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.DesignBuilder


def numerical_balance(df: pd.DataFrame, column: str, transform=None) -> list[dict]:
    values = pd.to_numeric(df[column], errors="coerce")
    if transform:
        values = transform(values)
    groups = [values[df.coverage_group.eq(name)].dropna() for name in ("eligible", "ineligible")]
    a, b = groups
    pooled_sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    smd = (a.mean() - b.mean()) / pooled_sd if pooled_sd > 0 else np.nan
    return [{"variable": column, "level": "", "group": group, "sites": len(part),
             "missing_sites": int(df.coverage_group.eq(group).sum() - len(part)),
             "mean_or_proportion": float(part.mean()), "median": float(part.median()),
             "sd": float(part.std(ddof=1)), "smd_eligible_minus_ineligible": float(smd),
             "scale": "log1p" if transform else "raw"}
            for group, part in zip(("eligible", "ineligible"), groups)]


def categorical_balance(df: pd.DataFrame, column: str) -> list[dict]:
    result = []
    for level in sorted(df[column].fillna("<missing>").astype(str).unique()):
        a = df[df.coverage_group.eq("eligible")][column].fillna("<missing>").astype(str).eq(level)
        b = df[df.coverage_group.eq("ineligible")][column].fillna("<missing>").astype(str).eq(level)
        pa, pb = a.mean(), b.mean()
        pooled = np.sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
        smd = (pa - pb) / pooled if pooled > 0 else np.nan
        for group, flag, proportion in (("eligible", a, pa), ("ineligible", b, pb)):
            result.append({"variable": column, "level": level, "group": group,
                           "sites": int(flag.sum()), "missing_sites": 0,
                           "mean_or_proportion": float(proportion), "median": "", "sd": "",
                           "smd_eligible_minus_ineligible": float(smd), "scale": "indicator"})
    return result


def main() -> None:
    targets = [OUT, REPORT]
    if any(path.exists() for path in targets):
        raise FileExistsError(f"Stage 3 output already exists: {[str(p) for p in targets if p.exists()]}")
    main_df = pd.read_csv(MAIN, compression="gzip", low_memory=False)
    site_df = pd.read_csv(SITE, compression="gzip", low_memory=False)
    pair_df = pd.read_csv(PAIR, compression="gzip", low_memory=False)
    old19 = pd.read_csv(STAGE19)
    old09 = pd.read_csv(STAGE09).iloc[0]
    if len(main_df) != 25996 or main_df.site_id.duplicated().any():
        raise AssertionError("HeLa primary population changed")
    if len(pair_df) != 3705 or pair_df.site_id.duplicated().any():
        raise AssertionError("Stage-19 paired population changed")

    # Stage-19 estimands are copied only after checking the paired data and coefficient scales.
    estimands = []
    for _, row in old19.iterrows():
        item = row.to_dict()
        if row.row_type == "condition_effect":
            observed_sd = float(pair_df[row.predictor].std(ddof=1))
            if not np.isclose(observed_sd, row.predictor_sd, atol=1e-12):
                raise AssertionError(f"Stage-19 predictor SD mismatch: {row.analysis}")
            if not np.isclose(row.beta_raw * observed_sd, row.beta_per_sd, atol=1e-12):
                raise AssertionError(f"Stage-19 standardized coefficient mismatch: {row.analysis}")
        item["analysis_role"] = "historical_stage19_exploratory"
        item["coefficient_difference_interpretation"] = (
            "difference_between_condition_specific_1SD_changes" if row.row_type == "paired_effect_contrast" else ""
        )
        estimands.append(item)
    observed_delta = pair_df.reactivity_mean_flank10_invivo - pair_df.reactivity_mean_flank10_invitro
    if not np.allclose(observed_delta, pair_df.reactivity_delta_invivo_minus_invitro, atol=1e-12):
        raise AssertionError("Stage-19 delta reactivity formula mismatch")
    contrast = old19.set_index("analysis").loc["in_vivo_minus_in_vitro"]
    vivo = old19.set_index("analysis").loc["in_vivo"]
    vitro = old19.set_index("analysis").loc["in_vitro"]
    if not np.isclose(vivo.beta_per_sd - vitro.beta_per_sd, contrast.beta_per_sd, atol=1e-12):
        raise AssertionError("Stage-19 contrast mismatch")

    # Reuse the frozen HeLa design, changing only the response column.
    Builder = load_design_builder()
    builder = Builder(main_df)
    design = builder.build(main_df, "reactivity_mean_flank10")
    constant = [name for name in design if name != "Intercept" and design[name].nunique(dropna=False) <= 1]
    design = design.drop(columns=constant)
    predictor = "reactivity_mean_flank10"
    outcome_rows = []
    for response in ("combined_ratio", "norme_ratio_rep1", "norme_ratio_rep2"):
        y = pd.to_numeric(main_df[response], errors="raise")
        fit = sm.OLS(y, design).fit(cov_type="cluster", cov_kwds={"groups": main_df.analysis_gene,
                                                                 "use_correction": True})
        low, high = fit.conf_int().loc[predictor]
        sd = float(main_df[predictor].std(ddof=1))
        outcome_rows.append({"outcome": response, "analysis_role": "frozen_primary_reproduction" if response == "combined_ratio" else "posthoc_replicate_specific_sensitivity",
                             "sites": len(main_df), "genes": main_df.analysis_gene.nunique(),
                             "beta_raw": float(fit.params[predictor]), "se_cluster": float(fit.bse[predictor]),
                             "ci95_low_raw": float(low), "ci95_high_raw": float(high),
                             "p_value": float(fit.pvalues[predictor]), "predictor_sd": sd,
                             "beta_per_sd": float(fit.params[predictor] * sd),
                             "ci95_low_per_sd": float(low * sd), "ci95_high_per_sd": float(high * sd),
                             "response_definition": "count_pooled_Acov_over_AGcov" if response == "combined_ratio" else "GLORI_replicate_specific_Ratio_times_one_minus_NonCR"})
    if not np.isclose(outcome_rows[0]["beta_raw"], old09.beta_raw, atol=1e-12):
        raise AssertionError("HeLa primary coefficient not reproduced")

    detected = site_df.loc[site_df.detected_both_replicates.astype(str).eq("True")].copy()
    detected["coverage_group"] = np.where(detected.main_window_valid.astype(str).eq("True"), "eligible", "ineligible")
    counts = detected.coverage_group.value_counts().to_dict()
    if counts != {"eligible": 25996, "ineligible": 39054}:
        raise AssertionError(f"HeLa coverage strata changed: {counts}")
    if set(detected.loc[detected.coverage_group.eq("eligible"), "site_id"]) != set(main_df.site_id):
        raise AssertionError("HeLa eligibility does not equal primary population")
    balance = []
    for column in ("combined_ratio", "combined_agcov", "gc_fraction_21", "transcript_position_fraction",
                   "isoform_count", "distance_to_stop_codon_tx", "distance_to_nearest_splice_edge_tx"):
        balance.extend(numerical_balance(detected, column, np.log1p if column == "combined_agcov" else None))
    for column in ("drach_subtype", "transcript_region"):
        balance.extend(categorical_balance(detected, column))
    overview = [{"population": "detected_both_GLORI_replicates_and_mapped_to_icSHAPE", "sites": len(detected),
                 "eligible_sites": counts["eligible"], "ineligible_sites": counts["ineligible"],
                 "eligible_fraction": counts["eligible"] / len(detected),
                 "eligible_genes": detected.loc[detected.coverage_group.eq("eligible"), "analysis_gene"].nunique(),
                 "ineligible_genes": detected.loc[detected.coverage_group.eq("ineligible"), "analysis_gene"].nunique()}]

    # Refit exactly the existing 3-df natural cubic spline and plot a marginal adjusted curve.
    x = pd.to_numeric(main_df[predictor]).to_numpy()
    spline = dmatrix("cr(x, df=3, constraints='center') - 1", {"x": x}, return_type="dataframe")
    spline_design = design.drop(columns=[predictor]).copy()
    spline_names = [f"spline_reactivity_{i+1}" for i in range(spline.shape[1])]
    for i, name in enumerate(spline_names):
        spline_design[name] = spline.iloc[:, i].to_numpy()
    spline_fit = sm.OLS(main_df.combined_ratio, spline_design).fit(cov_type="cluster", cov_kwds={"groups": main_df.analysis_gene, "use_correction": True})
    grid = np.linspace(np.quantile(x, .01), np.quantile(x, .99), 140)
    basis_grid = np.asarray(build_design_matrices([spline.design_info], {"x": grid})[0])
    mean_row = spline_design.mean().to_numpy(dtype=float)
    columns = list(spline_design.columns)
    design_grid = np.repeat(mean_row[None, :], len(grid), axis=0)
    for i, name in enumerate(spline_names):
        design_grid[:, columns.index(name)] = basis_grid[:, i]
    beta = spline_fit.params.to_numpy()
    covariance = spline_fit.cov_params().to_numpy()
    prediction = design_grid @ beta
    se = np.sqrt(np.einsum("ij,jk,ik->i", design_grid, covariance, design_grid))
    curve = [{"reactivity_flank10": float(xx), "adjusted_ratio": float(yy),
              "ci95_low": float(yy - 1.96 * ss), "ci95_high": float(yy + 1.96 * ss),
              "sites": len(main_df), "genes": main_df.analysis_gene.nunique(),
              "definition": "marginal_standardization_over_observed_covariates_spline_df3"}
             for xx, yy, ss in zip(grid, prediction, se)]
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7.2, 5.5), sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1]}, constrained_layout=True)
    top.plot(grid, prediction, color="#177a8a", linewidth=2)
    top.fill_between(grid, prediction - 1.96 * se, prediction + 1.96 * se,
                     color="#177a8a", alpha=.18, linewidth=0)
    top.set_ylabel("Adjusted pooled GLORI Ratio")
    top.set_title("HeLa: adjusted 3-df spline and predictor distribution")
    bottom.hist(x, bins=45, color="#7896a0", edgecolor="white", linewidth=.3)
    bottom.set_ylabel("Sites")
    bottom.set_xlabel("Mean icSHAPE reactivity, ±10-nt flanks")
    top.set_xlim(grid.min(), grid.max())

    OUT.mkdir(parents=True)
    write_csv(OUT / "stage19_estimands.csv", estimands)
    write_csv(OUT / "hela_outcome_sensitivity.csv", outcome_rows)
    write_csv(OUT / "hela_coverage_balance.csv", balance)
    write_csv(OUT / "hela_coverage_population.csv", overview)
    write_csv(OUT / "hela_adjusted_spline_curve.csv", curve)
    fig.savefig(OUT / "hela_adjusted_spline_and_distribution.png", dpi=220)
    plt.close(fig)
    inputs = [MAIN, SITE, PAIR, STAGE19, STAGE09, STAGE09_CODE]
    provenance = {"analysis_role": "posthoc_presubmission_stage3", "inputs_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in inputs},
                  "script_sha256": sha256(Path(__file__)), "sites": {"hela_primary": len(main_df), "hela_both_replicates": len(detected), "stage19_paired": len(pair_df)},
                  "model": "stage09b_frozen_covariate_design_gene_cluster_CR1", "spline": "stage09b_3df_natural_cubic_centered_knots_reused",
                  "normalization": "replicates_analyzed_separately_no_average_or_pooling"}
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")

    balance_df = pd.DataFrame(balance)
    def smd(col):
        return float(balance_df.loc[(balance_df.variable == col) & (balance_df.group == "eligible"), "smd_eligible_minus_ineligible"].iloc[0])
    report = f"""# main_v2 阶段 3：关联估计量与筛选人群补强

状态：完成。全部新增分析为投稿前事后敏感性或描述性诊断；不改变原主分析定义。

## 1. 体内／体外估计量与图 2C 定义

Stage 19 使用同一 HEK293T 位点和代表转录本上的体内、体外 icSHAPE 反应性；两条件的 ±10 nt 两侧各需至少 70% 有效覆盖。配对人群为 **{len(pair_df):,} 位点、{pair_df.analysis_gene.nunique():,} 基因**。GLORI 响应是两重复计数合并后的 `combined_ratio = (Acov₁ + Acov₂)/(AGcov₁ + AGcov₂)`，线性模型按基因采用 CR1 聚类稳健区间，属于观察性关联。

| 预测变量 | 原始斜率（每 1 反应性单位） | 预测变量 SD | 每 SD 斜率（95% CI；修饰比例百分点） |
|---|---:|---:|---:|
| 体内反应性 | {vivo.beta_raw:.4f} [{vivo.ci95_low_raw:.4f}, {vivo.ci95_high_raw:.4f}] | {vivo.predictor_sd:.4f} | {100*vivo.beta_per_sd:.3f} [{100*vivo.ci95_low_per_sd:.3f}, {100*vivo.ci95_high_per_sd:.3f}] |
| 体外反应性 | {vitro.beta_raw:.4f} [{vitro.ci95_low_raw:.4f}, {vitro.ci95_high_raw:.4f}] | {vitro.predictor_sd:.4f} | {100*vitro.beta_per_sd:.3f} [{100*vitro.ci95_low_per_sd:.3f}, {100*vitro.ci95_high_per_sd:.3f}] |
| Δ反应性：体内 − 体外 | {old19.set_index('analysis').loc['delta_invivo_minus_invitro'].beta_raw:.4f} [{old19.set_index('analysis').loc['delta_invivo_minus_invitro'].ci95_low_raw:.4f}, {old19.set_index('analysis').loc['delta_invivo_minus_invitro'].ci95_high_raw:.4f}] | {observed_delta.std(ddof=1):.4f} | {100*old19.set_index('analysis').loc['delta_invivo_minus_invitro'].beta_per_sd:.3f} [{100*old19.set_index('analysis').loc['delta_invivo_minus_invitro'].ci95_low_per_sd:.3f}, {100*old19.set_index('analysis').loc['delta_invivo_minus_invitro'].ci95_high_per_sd:.3f}] |

图 2C 应明确 `Δ reactivity = mean flank10 icSHAPE(in vivo) − mean flank10 icSHAPE(in vitro)`，在相同位点、相同代表转录本上计算；不是两组不配对均值之差。其样本均值为 {observed_delta.mean():.4f}。

**另一个估计量**是体内与体外各自标准化斜率之差：{100*contrast.beta_per_sd:.3f} 个百分点，基因配对 bootstrap 95% CI [{100*contrast.ci95_low_per_sd:.3f}, {100*contrast.ci95_high_per_sd:.3f}]（{int(contrast.bootstrap_replicates_completed)} 次）。它使用各条件自己的 SD（{vivo.predictor_sd:.4f} 与 {vitro.predictor_sd:.4f}），不能解释为同一物理尺度下的斜率差，也不等于以 Δ反应性为单一预测变量的斜率。该对比和 Δ反应性回归必须分别标注。

## 2. HeLa 标签敏感性

主标签继续使用 count-pooled `combined_ratio`。GLORI 的 `NormeRatio = Ratio × (1 − NonCR)` 是重复特异的归一化值，因此对重复 1 和重复 2 **分别**复用 stage 09B 的固定协变量设计（含 GC、转录本位置、测序覆盖、结构覆盖、距离、motif 与区域）和基因聚类 CR1 推断；不取两重复简单均值。共同人群为 {len(main_df):,} 位点、{main_df.analysis_gene.nunique():,} 基因。

| 标签 | 调整斜率/每 SD 反应性（修饰比例百分点；95% CI） | P |
|---|---:|---:|
"""
    for row in outcome_rows:
        report += f"| `{row['outcome']}` | {100*row['beta_per_sd']:.3f} [{100*row['ci95_low_per_sd']:.3f}, {100*row['ci95_high_per_sd']:.3f}] | {row['p_value']:.3g} |\n"
    report += f"""

主模型斜率与原 09B 表精确核对（容差 1e−12）。重复特异归一化标签仅检验方向和量级对 GLORI 字段选择的敏感性，不能取代 count-pooled 主结果，也不增加独立生物重复数。

## 3. HeLa 结构覆盖筛选

在已匹配 icSHAPE 且两个 GLORI 重复均检出的 **{len(detected):,}** 个位点中，±10 nt 两侧各至少 70% 有效覆盖的 **{counts['eligible']:,}** 个（{100*counts['eligible']/len(detected):.1f}%）进入 HeLa 主关联；覆盖不合格 **{counts['ineligible']:,}** 个。合格集合与原 09B 主分析 site_id 完全一致。平衡表采用合格组减不合格组的标准化均差（SMD），并列示缺失数；不将筛选后关联外推到覆盖不足位点。

| 特征 | 合格组均值 | 不合格组均值 | SMD |
|---|---:|---:|---:|
"""
    for col, label in (("combined_ratio", "count-pooled Ratio"), ("combined_agcov", "log1p(合并 AGcov)"),
                       ("gc_fraction_21", "21 nt GC 比例"), ("transcript_position_fraction", "转录本位置比例"),
                       ("isoform_count", "候选 isoform 数")):
        a, b = balance_df[(balance_df.variable == col) & (balance_df.group == "eligible")].iloc[0], balance_df[(balance_df.variable == col) & (balance_df.group == "ineligible")].iloc[0]
        report += f"| {label} | {a.mean_or_proportion:.4f} | {b.mean_or_proportion:.4f} | {smd(col):+.3f} |\n"
    report += "\n完整数值、距离变量缺失和 motif/区域分类比例见 `results/presubmission_stage3/hela_coverage_balance.csv`。这是筛选可观测协变量的描述性平衡，不足以证明缺失结构数据满足随机缺失。\n"
    report += f"""

## 4. 非线性诊断

用原 09B 相同的 3 自由度自然三次样条基函数和主分析协变量重拟合，绘制全体主分析位点的边际调整曲线及反应性直方图：[曲线与分布](../../../results/presubmission_stage3/hela_adjusted_spline_and_distribution.png)。曲线仅绘于反应性第 1–99 百分位；阴影为基因聚类稳健逐点 95% 区间，并非同时置信带。数据表为 `hela_adjusted_spline_curve.csv`。原诊断的样条 − 线性 AIC = −2.699，线性模型 RESET P = 1.11×10⁻¹⁹；这些属于函数形式诊断，不据曲线选择有利亚组或声称因果作用。

## 5. 解释边界和复核

- 本阶段未执行可选的基因内/基因间分解；基因聚类稳健标准误只处理同基因观测误差相关，不能视作完全控制基因层面的混杂。
- Stage 19 原始、标准化斜率和区间已从保存结果提取，并与配对数据的各预测变量 SD、Δ 定义及标准化系数差核对。NormeRatio 模型与原主模型共用位点、协变量矩阵和聚类方式。
- 原始输入及原 09B 设计代码的 SHA-256 存于 `results/presubmission_stage3/provenance.json`；本阶段全部交付物位于 `results/presubmission_stage3/`。下一阶段修订文稿时，应把图 2C 和筛选人群定义逐一映射到正文与图注。
"""
    REPORT.write_text(report, encoding="utf-8", errors="strict")
    print(json.dumps({"status": "PASS", "output": str(OUT), "report": str(REPORT), "outcomes": outcome_rows,
                      "coverage": overview, "spline_range": [float(grid.min()), float(grid.max())]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
