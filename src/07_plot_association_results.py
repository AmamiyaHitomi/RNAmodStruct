"""Render the stage-07 HEK293T association overview with publication-style QA."""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / "logs" / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from audit_panel_alignment import require_matplotlib_panel_alignment


plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42})
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False


FINAL = ROOT / "data" / "final"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
QA = FIGURES / "qa"
BASE = FIGURES / "07_hek293t_association_overview"

BLUE = "#0F4D92"
BLUE_LIGHT = "#8FB8D8"
TEAL = "#42949E"
NEUTRAL = "#767676"
RED = "#B64342"
QUARTILE_COLORS = ["#B8C7D9", "#6E9AC5", "#3775BA", "#B64342"]


def add_panel_label(ax, label: str) -> None:
    from matplotlib.transforms import ScaledTranslation

    offset = ScaledTranslation(-7 / 72, 3 / 72, ax.figure.dpi_scale_trans)
    ax.text(0, 1, label, transform=ax.transAxes + offset, fontsize=8, fontweight="bold", ha="left", va="bottom")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(FINAL / "06_hek293t_main_analysis_dataset.csv.gz", compression="gzip", low_memory=False)
    trend = read_csv(TABLES / "07_hek293t_binned_reactivity_ratio_trend.csv")
    primary = read_csv(TABLES / "07_hek293t_primary_association.csv")
    sensitivity = read_csv(TABLES / "07_hek293t_sensitivity_associations.csv")
    profile = read_csv(TABLES / "07_hek293t_reactivity_profile_by_ratio_quartile.csv")

    fig = plt.figure(figsize=(7.09, 4.13), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.12, 1.0])
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])
    axes = [ax_a, ax_b, ax_c]

    ax_a.scatter(
        data["reactivity_mean_flank10"],
        data["combined_ratio"],
        s=4,
        alpha=0.12,
        color=BLUE_LIGHT,
        edgecolors="none",
        rasterized=True,
    )
    ax_a.plot(trend["reactivity_mean"], trend["ratio_mean"], color=BLUE, marker="o", markersize=2.5, linewidth=1.3)
    ax_a.set_xlabel("Mean icSHAPE reactivity, ±10 nt")
    ax_a.set_ylabel("Pooled m6A proportion")
    ax_a.set_xlim(left=0)
    ax_a.set_ylim(0, 1.02)
    ax_a.text(0.98, 0.97, f"n = {len(data):,} sites", transform=ax_a.transAxes, ha="right", va="top", color=NEUTRAL, fontsize=6)
    add_panel_label(ax_a, "a")

    display_names = {
        "primary_adjusted_OLS": "Primary",
        "exclude_minus2_to_plus2": "Exclude −2…+2",
        "complete_flank10_coverage": "Complete coverage",
        "GLORI_replicate_1": "GLORI replicate 1",
        "GLORI_replicate_2": "GLORI replicate 2",
        "equal_weight_per_gene": "Equal gene weight",
        "single_isoform_sites": "Single isoform",
        "canonical_DRACH_sites": "Canonical DRACH",
    }
    forest = pd.concat([primary.assign(analysis="primary_adjusted_OLS"), sensitivity], ignore_index=True)
    forest["label"] = forest["analysis"].map(display_names)
    forest = forest.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(forest))
    for index, row in forest.iterrows():
        color = BLUE if row["analysis"] == "primary_adjusted_OLS" else NEUTRAL
        ax_b.plot([row["ci95_low_per_sd"], row["ci95_high_per_sd"]], [index, index], color=color, linewidth=1.2)
        ax_b.plot(row["beta_per_sd"], index, marker="o", color=color, markersize=3.5)
    ax_b.axvline(0, color="#B8B8B8", linestyle="--", linewidth=0.8)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(forest["label"], fontsize=6)
    ax_b.set_xlabel("Δ m6A proportion per 1 SD reactivity")
    ax_b.set_ylim(-0.7, len(forest) - 0.3)
    add_panel_label(ax_b, "b")

    for color, quartile in zip(QUARTILE_COLORS, ["Q1", "Q2", "Q3", "Q4"]):
        part = profile[profile["ratio_quartile"] == quartile].sort_values("relative_position")
        x = part["relative_position"].to_numpy()
        mean = part["mean_reactivity"].to_numpy(dtype=float)
        low = part["ci95_low_pointwise"].to_numpy(dtype=float)
        high = part["ci95_high_pointwise"].to_numpy(dtype=float)
        ax_c.fill_between(x, low, high, color=color, alpha=0.13, linewidth=0)
        ax_c.plot(x, mean, color=color, linewidth=1.1, label=quartile)
    ax_c.axvline(0, color="#B8B8B8", linestyle="--", linewidth=0.8)
    ax_c.set_xlabel("Position relative to m6A site (nt)")
    ax_c.set_ylabel("Mean icSHAPE reactivity")
    ax_c.set_xlim(-50, 50)
    ax_c.set_ylim(bottom=0)
    for x_label, color, quartile in zip([0.68, 0.77, 0.86, 0.95], QUARTILE_COLORS, ["Q1", "Q2", "Q3", "Q4"]):
        ax_c.text(
            x_label,
            1.02,
            quartile,
            transform=ax_c.transAxes,
            color=color,
            fontsize=6,
            ha="center",
            va="bottom",
            clip_on=False,
        )
    add_panel_label(ax_c, "c")

    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=["a", "b", "c"],
        json_out=QA / "07_hek293t_association_overview_alignment.json",
        overlay_svg=QA / "07_hek293t_association_overview_alignment_overlay.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )
    fig.savefig(BASE.with_suffix(".svg"))
    fig.savefig(BASE.with_suffix(".pdf"))
    fig.savefig(BASE.with_suffix(".tiff"), dpi=600)
    fig.savefig(BASE.with_suffix(".png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
