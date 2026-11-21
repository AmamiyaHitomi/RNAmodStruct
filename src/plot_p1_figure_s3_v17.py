"""Draw the supplementary positional association and coverage figure."""

from __future__ import annotations

import argparse
from pathlib import Path

from plot_p1_coverage_and_splines import GRID, INK, ORIGINAL_CURVE, SLATE, export
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["TeX Gyre Termes", "Times New Roman", "DejaVu Serif"],
    "text.color": "#30303d",
    "axes.labelcolor": "#30303d",
    "xtick.color": "#30303d",
    "ytick.color": "#30303d",
})

ESTIMATE = "reactivity_change_per_10pp_m6a"
CI_LOW = "reactivity_change_ci95_low"
CI_HIGH = "reactivity_change_ci95_high"


def reverse_holm(data):
    # Cluster-robust OLS uses normal-reference 95% intervals.
    # Confirm that convention against the reported forward-model p-values.
    quantile = norm.ppf(.975)
    forward_se = (data.ci95_high_pp - data.ci95_low_pp) / (2 * quantile)
    np.testing.assert_allclose(
        2 * norm.sf(abs(data.beta_pp_per_reference_sd / forward_se)),
        data.p_value, rtol=1e-8, atol=1e-13,
    )
    reverse_se = (data[CI_HIGH] - data[CI_LOW]) / (2 * quantile)
    np.testing.assert_allclose(
        (data[CI_HIGH] + data[CI_LOW]) / 2, data[ESTIMATE], atol=1e-14,
    )
    if not np.all(np.isfinite(reverse_se) & (reverse_se > 0)):
        raise ValueError("Invalid reverse-model confidence intervals")
    p = 2 * norm.sf(abs(data[ESTIMATE] / reverse_se))
    return multipletests(p, method="holm")[1]


def mark_center(axis):
    for left, right in ((-2.5, -.5), (.5, 2.5)):
        axis.axvspan(left, right, facecolor=GRID, linewidth=0, alpha=.55, zorder=0)
    axis.axvline(0, color=SLATE, linestyle=(0, (2, 2)), linewidth=.7, zorder=0)


def plot(data, output):
    fig = plt.figure(figsize=(7.15, 4.15))
    grid = fig.add_gridspec(2, 2, height_ratios=(3.4, 1.2), left=.135, right=.963,
                           top=.83, bottom=.19, wspace=.24, hspace=.19)
    step = .005
    limits = (step * np.floor(data[CI_LOW].min() / step) - .001,
              step * np.ceil(data[CI_HIGH].max() / step) + .001)
    coverage_min = max(0, 5 * np.floor(data.pct_of_primary_sites.min() / 5) - 1)
    for i, (cell, color) in enumerate((("HEK293T", "#685580"), ("HeLa", "#176d73"))):
        part = data.query("cell_line == @cell").sort_values("position_nt")
        axis = fig.add_subplot(grid[0, i])
        access = fig.add_subplot(grid[1, i], sharex=axis)
        axis.set_title(cell, loc="left", color=color, fontweight="bold", fontsize=9.5, pad=8)
        for panel in (axis, access):
            mark_center(panel)
            panel.set_axisbelow(True)
            panel.grid(axis="y", linewidth=.5, color=GRID)
        axis.axhline(0, color=SLATE, linewidth=.75, zorder=1)
        axis.vlines(part.position_nt, part[CI_LOW], part[CI_HIGH],
                    color=color, linewidth=.95, alpha=.72, zorder=2)
        selected = part.reverse_p_holm.lt(.05)
        for mask, face in ((selected, color), (~selected, "white")):
            axis.scatter(part.loc[mask, "position_nt"], part.loc[mask, ESTIMATE],
                         facecolor=face, edgecolor=color, s=20, linewidth=.95, zorder=3)
        axis.set_ylim(*limits)
        axis.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5))
        axis.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
        axis.tick_params(axis="x", labelbottom=False, bottom=False)
        for same_side in (part.position_nt.lt(0), part.position_nt.gt(0)):
            counts = part.loc[same_side]
            access.plot(counts.position_nt, counts.pct_of_primary_sites,
                        color=color, marker=".", markersize=3, linewidth=1)
        access.set_ylim(coverage_min, 100.5)
        access.yaxis.set_major_locator(mticker.FixedLocator(np.linspace(coverage_min + 2, 100, 3)))
        access.yaxis.set_major_formatter(mticker.PercentFormatter(100, decimals=0))
        access.xaxis.set_major_locator(mticker.FixedLocator([-10, -5, 0, 5, 10]))
        access.set_xlim(-10.6, 10.6)
        if i == 0:
            axis.set_ylabel("icSHAPE reactivity difference\n(per 10 pp of m6A)", linespacing=1.15)
            access.set_ylabel("Sites with\nvalid values", linespacing=1.15)
        else:
            for panel in (axis, access):
                panel.tick_params(axis="y", labelleft=False)
    key = [
        Line2D([], [], marker="o", color=INK, linestyle="none", markersize=4, label="Holm-adjusted P < 0.05"),
        Line2D([], [], marker="o", color=INK, markerfacecolor="white", linestyle="none",
               markersize=4, label="Holm-adjusted P ≥ 0.05"),
    ]
    fig.legend(handles=key, loc="upper center", bbox_to_anchor=(.57, .99), ncol=2,
               frameon=False, fontsize=7.1, handletextpad=.3, columnspacing=1.35)
    fig.supxlabel("Position relative to the modified A (nt)", fontsize=8, y=.071)
    export(fig, output / "figureS_local_reactivity_oct2026")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("positions_csv", type=Path)
    parser.add_argument("output", type=Path, help="Fresh figure-output directory")
    args = parser.parse_args()
    positions = pd.read_csv(args.positions_csv)
    if len(positions) != 40 or positions.duplicated(["cell_line", "position_nt"]).any():
        raise ValueError("Expected 20 unique flanking positions per cell line")
    args.output.mkdir(parents=True, exist_ok=False)
    positions["reverse_p_holm"] = reverse_holm(positions)
    print(positions.loc[positions.reverse_p_holm.lt(.05), ["cell_line", "position_nt", ESTIMATE]].to_string(index=False))
    plot(positions, args.output)


if __name__ == "__main__":
    main()
