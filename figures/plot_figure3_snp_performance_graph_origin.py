#!/usr/bin/env python3
"""Plot SNP performance stratified by graph type and sample origin."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path("/data/malong/benchmark")
INPUT = ROOT / "paper/result_statistics.tsv"
OUT = ROOT / "paper/png/Figure3_SNP_performance_graph_origin.png"
OUT_SVG = ROOT / "paper/png/Figure3_SNP_performance_graph_origin.svg"

DEPTHS = ["5x", "10x", "15x", "20x"]
METRICS = [("f1", "F1"), ("precision", "Precision"), ("recall", "Recall")]
TOOLS = ["vg", "graphtyper", "pangenie", "bayestyper"]
TOOL_LABELS = {
    "vg": "VG",
    "graphtyper": "GraphTyper2",
    "pangenie": "PanGenie",
    "bayestyper": "BayesTyper",
}
TOOL_COLORS = {
    "vg": "#0066CC",
    "graphtyper": "#009E5A",
    "pangenie": "#F28E00",
    "bayestyper": "#D62728",
}
GROUPS = [
    ("MC", "Graph-in samples"),
    ("MC", "Graph-out samples"),
    ("PGGB", "Graph-in samples"),
    ("PGGB", "Graph-out samples"),
]
SAMPLE_GROUPS = {
    "Graph-in samples": {"W0169", "W1735", "W3033"},
    "Graph-out samples": {"9311", "IR64", "zhonghua11"},
}


def read_values():
    grouped = defaultdict(list)
    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["variant_class"] != "snp" or row["status"] != "complete":
                continue
            if row["graph"] not in {"MC", "PGGB"}:
                continue
            origin = next(
                origin for origin, samples in SAMPLE_GROUPS.items()
                if row["sample"] in samples
            )
            for metric, _ in METRICS:
                value = row[metric]
                if value not in ("", "NA", "."):
                    grouped[(row["graph"], origin, row["tool"], row["depth"], metric)].append(float(value))
    return grouped


def main():
    grouped = read_values()
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 20,
        "axes.titlesize": 23,
        "axes.labelsize": 21,
        "xtick.labelsize": 19,
        "ytick.labelsize": 19,
        "legend.fontsize": 20,
        "svg.fonttype": "none",
    })

    # The middle narrow column creates a clean visual separation between MC and PGGB.
    fig = plt.figure(figsize=(20, 13.5), dpi=200)
    grid = fig.add_gridspec(
        3,
        5,
        width_ratios=[1, 1, 0.18, 1, 1],
        left=0.07,
        right=0.985,
        bottom=0.13,
        top=0.91,
        wspace=0.17,
        hspace=0.26,
    )
    axes = {}
    for row_idx, (metric, _) in enumerate(METRICS):
        for group_idx, (graph, origin) in enumerate(GROUPS):
            actual_col = group_idx if group_idx < 2 else group_idx + 1
            axes[(row_idx, group_idx)] = fig.add_subplot(grid[row_idx, actual_col])

    for row_idx, (metric, metric_label) in enumerate(METRICS):
        for group_idx, (graph, origin) in enumerate(GROUPS):
            ax = axes[(row_idx, group_idx)]
            ax.set_facecolor("white")
            ax.set_ylim(0, 1.02)
            ax.set_xlim(-0.15, len(DEPTHS) - 1 + 0.15)
            ax.set_xticks(range(len(DEPTHS)))
            ax.set_xticklabels(["5×", "10×", "15×", "20×"])
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
            ax.tick_params(axis="both", labelsize=19, length=6, width=1.5, colors="#111111")
            ax.grid(False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("#111111")
            ax.spines["bottom"].set_color("#111111")
            ax.spines["left"].set_linewidth(1.5)
            ax.spines["bottom"].set_linewidth(1.5)

            if row_idx == 0:
                ax.set_title(f"{graph} graph\n{origin}", fontsize=23, fontweight="bold", pad=16)
            if group_idx in (0, 2):
                ax.set_ylabel(metric_label, fontsize=21, fontweight="bold", labelpad=14)
            else:
                ax.tick_params(axis="y", labelleft=False)
            if row_idx < len(METRICS) - 1:
                ax.tick_params(axis="x", labelbottom=False)

            for tool in TOOLS:
                means = []
                present = False
                for depth in DEPTHS:
                    values = grouped.get((graph, origin, tool, depth, metric), [])
                    if values:
                        means.append(float(np.mean(values)))
                        present = True
                    else:
                        means.append(np.nan)
                if not present:
                    continue
                ax.plot(
                    range(len(DEPTHS)),
                    means,
                    marker="o",
                    color=TOOL_COLORS[tool],
                    linewidth=4.4,
                    markersize=11,
                    alpha=0.96,
                    label=TOOL_LABELS[tool],
                )

    legend_handles = [
        Line2D(
            [0], [0], color=TOOL_COLORS[tool], marker="o", linewidth=4.4,
            markersize=11, label=TOOL_LABELS[tool],
        )
        for tool in TOOLS
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=4,
        frameon=False,
        columnspacing=1.8,
        handletextpad=0.6,
    )
    fig.savefig(OUT, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.12)
    fig.savefig(OUT_SVG, format="svg", facecolor="white", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    print(OUT)
    print(OUT_SVG)


if __name__ == "__main__":
    main()
