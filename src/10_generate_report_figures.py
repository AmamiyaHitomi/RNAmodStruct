"""Generate the two remaining main figures (1: data flow, 5: HeLa replication/transfer)."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "results" / "logs" / "mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "results" / "figures"
TABLES = ROOT / "results" / "tables"

STEPS = ["GLORI\nunion", "Matched\nicSHAPE", "Both GLORI\nreplicates", "Structure\neligible", "Full 201-nt\nwindow"]
HEK_COLOR = "#4C78A8"
HELA_COLOR = "#E45756"


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def style_ax(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def figure1() -> None:
    hek = read_rows(ROOT / "results" / "tables" / "06_hek293t_attrition_summary.csv")
    hela = read_rows(ROOT / "results" / "tables" / "09_hela_attrition_summary.csv")
    hek_counts = [int(row["sites"]) for row in hek]
    hela_counts = [int(row["sites"]) for row in hela]

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), sharey=False)
    plt.rcParams.update({"font.size": 8})
    for ax, counts, title in [(axes[0], hek_counts, "HEK293T (discovery)"), (axes[1], hela_counts, "HeLa (replication)")]:
        ax.bar(np.arange(len(STEPS)), counts, color=HEK_COLOR if title.startswith("HEK") else HELA_COLOR)
        ax.set_xticks(np.arange(len(STEPS)))
        ax.set_xticklabels(STEPS, fontsize=6.2, linespacing=1.05)
        ax.set_yscale("log")
        ax.set_ylim(3_000, max(counts) * 1.45)
        ax.set_ylabel("Sites (log scale)")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=8, pad=12)
        style_ax(ax)
        for x, count in enumerate(counts):
            ax.text(x, count * 1.1, f"{count:,}", ha="center", va="bottom", fontsize=6.5)
    fig.suptitle("Figure 1  Data acquisition and filtering", x=0.0, ha="left", fontweight="bold", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93), w_pad=2.5)
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(FIGURES / f"10_data_overview_figure1.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure5() -> None:
    hek_primary = read_rows(TABLES / "07_hek293t_primary_association.csv")[0]
    hela_primary = read_rows(TABLES / "09b_hela_primary_association.csv")[0]
    metrics = read_rows(TABLES / "09c_hela_transfer_metrics.csv")
    comparisons = read_rows(TABLES / "09c_hela_transfer_model_comparison.csv")

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.0), gridspec_kw={"wspace": 0.62})
    plt.rcParams.update({"font.size": 8})

    # Panel a: association replication
    ax = axes[0]
    labels = ["HEK293T\n(discovery)", "HeLa\n(replication)"]
    betas = [float(hek_primary["beta_per_sd"]) * 100, float(hela_primary["beta_per_sd"]) * 100]
    lows = [float(hek_primary["bootstrap_ci95_low_per_sd"]) * 100, float(hela_primary["bootstrap_ci95_low_per_sd"]) * 100]
    highs = [float(hek_primary["bootstrap_ci95_high_per_sd"]) * 100, float(hela_primary["bootstrap_ci95_high_per_sd"]) * 100]
    errs = [[b - l for b, l in zip(betas, lows)], [h - b for h, b in zip(highs, betas)]]
    xs = np.arange(len(labels))
    for x, beta, low, high, color in zip(xs, betas, lows, highs, [HEK_COLOR, HELA_COLOR]):
        ax.errorbar([x], [beta], yerr=[[beta - low], [high - beta]], fmt="o", color=color, capsize=4, linewidth=1)
    ax.axhline(0, color="black", linewidth=0.7, linestyle="--")
    ax.set_xticks(xs, labels, fontsize=6.5)
    ax.set_ylabel("β per SD reactivity\n(percentage points)")
    ax.set_title("a  Association replication", loc="left", fontweight="bold")
    for index, (x, beta) in enumerate(zip(xs, betas)):
        offset = (6, 7) if index == 0 else (0, 7)
        alignment = "left" if index == 0 else "center"
        ax.annotate(f"{beta:.2f}", (x, beta), xytext=offset, textcoords="offset points", ha=alignment, fontsize=7)
    style_ax(ax)

    # Panel b: transfer MAE
    ax = axes[1]
    models = ["M0", "M1", "M2", "M3", "M4"]
    by = {(r["subset"], r["model"]): r for r in metrics}
    all_mae = [float(by[("all_qualifying", m)]["mae"]) * 100 for m in models]
    strict_mae = [float(by[("strict", m)]["mae"]) * 100 for m in models]
    x = np.arange(len(models))
    width = 0.38
    ax.bar(x - width / 2, all_mae, width, label="All qualifying", color="#B279A2")
    ax.bar(x + width / 2, strict_mae, width, label="Strict subset", color="#72B7B2")
    ax.set_xticks(x, models)
    ax.set_ylabel("MAE (percentage points)")
    ax.set_title("b  HeLa transfer MAE", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=5.5, loc="upper right", ncol=1)
    style_ax(ax)

    # Panel c: transfer increment
    ax = axes[2]
    rows = [r for r in comparisons if r["comparison"] == "M3_vs_M1"]
    labels = ["All\nqualifying", "Strict\nsubset"]
    deltas = [float(r["delta_mae_percentage_points"]) for r in rows]
    lows = [float(r["bootstrap_ci95_low_percentage_points"]) for r in rows]
    highs = [float(r["bootstrap_ci95_high_percentage_points"]) for r in rows]
    errs = [[d - l for d, l in zip(deltas, lows)], [h - d for h, d in zip(highs, deltas)]]
    xs = np.arange(len(labels))
    ax.errorbar(xs, deltas, yerr=errs, fmt="o", color="#F28E2B", capsize=4, linewidth=1)
    ax.axhline(0, color="black", linewidth=0.7, linestyle="--")
    ax.set_xticks(xs, labels, fontsize=6.5)
    ax.set_ylabel("ΔMAE = MAE(M1) − MAE(M3)\n(percentage points)")
    ax.set_title("c  Experimental-structure increment", loc="left", fontweight="bold")
    for x, delta in zip(xs, deltas):
        ax.text(x, delta + 0.002, f"{delta:.3f}", ha="center", fontsize=7)
    style_ax(ax)

    fig.suptitle("Figure 5  HeLa replication and frozen-model transfer", x=0.0, ha="left", fontweight="bold", fontsize=9)
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.20, top=0.78, wspace=0.62)
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(FIGURES / f"10_hela_replication_transfer_figure5.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure1()
    figure5()
    print("generated figure1 and figure5")


if __name__ == "__main__":
    main()
