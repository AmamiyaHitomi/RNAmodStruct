"""Render the October paper's revised attrition and spline figures from P1 CSVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

from audit_panel_alignment import require_matplotlib_panel_alignment


INK = "#233b49"
TEAL = "#009b8c"
TEAL_DARK = "#007269"
SLATE = "#536876"
GRID = "#dce3e8"
CARD = "#f4f6f8"
ORIGINAL_CURVE = "#177a8a"
ORIGINAL_BARS = "#7896a0"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": .7,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
})


def export(fig, basename):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig, json_out=basename.with_suffix(".alignment.json"), strict=True, tolerance_pt=1.5
    )
    for suffix in (".pdf", ".svg", ".png", ".tiff"):
        with basename.with_suffix(suffix).open("xb") as output:
            options = {"pil_kwargs": {"compression": "tiff_lzw"}} if suffix == ".tiff" else {}
            fig.savefig(output, format=suffix[1:], dpi=600, **options)
    plt.close(fig)


def figure_one(data, primary_genes, dest):
    labels = (
        "Sites in either\nGLORI replicate",
        "Matched to structural\ntranscripts",
        "Supported in both\nGLORI replicates",
        "Sufficient structure\non both flanks",
    )
    fig = plt.figure(figsize=(7.15, 4.15))
    grid = fig.add_gridspec(1, 8, left=.045, right=.968, bottom=.08, top=.91, wspace=.65)
    step_ax = fig.add_subplot(grid[0, :2])
    step_ax.set_axis_off()
    centers = [.78, .565, .35, .135]
    step_ax.text(0, .95, "Selection stage", fontsize=8, fontweight="bold", va="center")
    for y, label in zip(centers, labels):
        step_ax.text(0, y, label, va="center", linespacing=1.24)
    for i, (cell, color, last) in enumerate((("HEK293T", INK, INK), ("HeLa", TEAL, TEAL_DARK))):
        group = data.loc[data.cell_line.eq(cell)].sort_values("step")
        ax = fig.add_subplot(grid[0, 2+i*3:5+i*3])
        ax.set_axis_off()
        ax.set_title(cell, fontsize=11, color=last, fontweight="bold", pad=10)
        ax.text(.28, .95, "Of initial", ha="center", va="center", fontsize=7, color=SLATE)
        ax.text(.72, .95, "Of previous", ha="center", va="center", fontsize=7, color=SLATE)
        for j, (y, row) in enumerate(zip(centers, group.itertuples())):
            is_final = j == len(centers)-1
            bg, txt = (last, "white") if is_final else (CARD, INK)
            ax.add_patch(FancyBboxPatch((.035, y-.076), .93, .15, transform=ax.transAxes,
                                       boxstyle="round,pad=0.006,rounding_size=0.026",
                                       facecolor=bg, edgecolor=color, linewidth=.8))
            ax.text(.5, y+.026, f"{row.sites:,}", va="center", ha="center",
                    color=txt, fontsize=11, fontweight="bold")
            ax.text(.28, y-.034, f"{row.fraction_from_initial:.2%}", va="center", ha="center",
                    color=txt, fontsize=7.4)
            previous = "—" if np.isnan(row.fraction_from_previous) else f"{row.fraction_from_previous:.2%}"
            ax.text(.72, y-.034, previous, va="center", ha="center", color=txt, fontsize=7.4)
            if not is_final:
                ax.annotate("", xy=(.5, centers[j+1]+.081), xytext=(.5, y-.083),
                            arrowprops={"arrowstyle": "->", "color": SLATE, "lw": .65})
        ax.text(.5, .008, f"{primary_genes[cell]:,} genes", va="center", ha="center",
                color=SLATE, fontsize=7)
    export(fig, dest / "figure1_attrition_oct2026")


def figure_three(curves, hist, attrition, dest):
    fig = plt.figure(figsize=(7.15, 4.45))
    grid = fig.add_gridspec(2, 2, height_ratios=(3.15, 1), left=.13, right=.967,
                           top=.86, bottom=.16, wspace=.23, hspace=.14)
    low = np.floor(curves.ci95_low.min() * 50) / 50
    high = np.ceil(curves.ci95_high.max() * 50) / 50
    peaks = hist.assign(percent=100 * hist.sites / hist.groupby("cell_line").sites.transform("sum"))
    hist_top = np.ceil(peaks.percent.max() / 2) * 2
    for i, (cell, line_color, bar_color) in enumerate(
        (("HEK293T", INK, "#869aa7"), ("HeLa", ORIGINAL_CURVE, ORIGINAL_BARS))
    ):
        curve = curves.loc[curves.cell_line.eq(cell)]
        distribution = peaks.loc[peaks.cell_line.eq(cell)]
        count = int(attrition.loc[attrition.cell_line.eq(cell), "sites"].iloc[-1])
        ax = fig.add_subplot(grid[0, i])
        density = fig.add_subplot(grid[1, i], sharex=ax)
        ax.set_title(f"{cell}  |  {count:,} sites", fontsize=9.2, fontweight="bold",
                     loc="left", color=line_color, pad=8)
        ax.plot(curve.reactivity_flank10, curve.adjusted_ratio, color=line_color, linewidth=1.8)
        ax.fill_between(curve.reactivity_flank10, curve.ci95_low, curve.ci95_high,
                        color=line_color, linewidth=0, alpha=.18)
        ax.set_ylim(low, high)
        ax.set_xlim(0, 1)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5))
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(1, decimals=0))
        ax.grid(axis="y", linewidth=.5, color=GRID)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        density.bar(distribution.left, distribution.percent, width=distribution.right-distribution.left,
                    align="edge", color=bar_color, edgecolor="white", linewidth=.34)
        density.set_ylim(0, hist_top)
        density.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3, integer=True))
        density.yaxis.set_major_formatter(mticker.PercentFormatter(100, decimals=0))
        density.xaxis.set_major_locator(mticker.FixedLocator([0, .2, .4, .6, .8, 1]))
        density.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
        if i == 0:
            ax.set_ylabel("Adjusted m⁶A level")
            density.set_ylabel("Sites")
        else:
            for axis in (ax, density):
                axis.tick_params(axis="y", labelleft=False)
    fig.supxlabel("Mean icSHAPE reactivity (±10-nt flanks)", y=.048, fontsize=8)
    export(fig, dest / "figure3_adjusted_splines_oct2026")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", type=Path)
    parser.add_argument("output", type=Path, help="New figure-output directory")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    attrition = pd.read_csv(args.tables / "attrition.csv")
    ipw = pd.read_csv(args.tables / "ipw_estimates.csv")
    primary_genes = ipw.loc[ipw.method.eq("unweighted")].set_index("cell_line").genes
    figure_one(attrition, primary_genes, args.output)
    figure_three(pd.read_csv(args.tables / "adjusted_spline.csv"),
                 pd.read_csv(args.tables / "reactivity_histogram.csv"), attrition, args.output)


if __name__ == "__main__":
    main()
