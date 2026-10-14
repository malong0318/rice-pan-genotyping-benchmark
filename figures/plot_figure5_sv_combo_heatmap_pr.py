#!/usr/bin/env python3
"""Combined SV figure: F1 heatmap and precision-recall panels."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path("/data/malong/benchmark")
INPUT = ROOT / "paper/result_statistics.tsv"
OUT = ROOT / "paper/png/Figure5_SV_performance_heatmap_PR.png"
OUT_SVG = ROOT / "paper/png/Figure5_SV_performance_heatmap_PR.svg"
OUT_SVG_LABELED = ROOT / "paper/png/Figure5_SV_performance_heatmap_PR_labeled.svg"
OUT_PNG_LABELED = ROOT / "paper/png/Figure5_SV_performance_heatmap_PR_labeled.png"

GRAPHS = ["MC", "PGGB", "minigraph"]
TOOLS = ["vg", "graphtyper", "pangenie", "bayestyper"]
TOOL_LABELS = {
    "vg": "VG",
    "graphtyper": "GraphTyper2",
    "pangenie": "PanGenie",
    "bayestyper": "BayesTyper",
}
COLORS = {
    "vg": "#277DA1",
    "graphtyper": "#2E8B57",
    "pangenie": "#E58F17",
    "bayestyper": "#C94C4C",
}
DEPTHS = ["5x", "10x", "15x", "20x"]
DEPTH_MARKERS = {"5x": "o", "10x": "s", "15x": "^", "20x": "D"}
IN_SAMPLES = {"W0169", "W1735", "W3033"}
ORIGINS = ["Graph-in", "Graph-out"]


def load_means():
    values = defaultdict(lambda: defaultdict(list))
    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["variant_class"] != "sv" or row["status"] != "complete":
                continue
            graph, tool, sample, depth = row["graph"], row["tool"], row["sample"], row["depth"]
            if graph not in GRAPHS or tool not in TOOLS or depth not in DEPTHS:
                continue
            origin = "Graph-in" if sample in IN_SAMPLES else "Graph-out"
            key = (graph, tool, origin, depth)
            for metric in ("f1", "precision", "recall"):
                raw = row[metric]
                if raw not in ("", "NA", "."):
                    values[key][metric].append(float(raw))
    return {
        key: {metric: float(np.mean(nums)) for metric, nums in metrics.items()}
        for key, metrics in values.items()
    }


def draw_heatmap(fig, means, blank=False, subplot_spec=121):
    ax = fig.add_subplot(subplot_spec)
    columns = [(origin, depth) for origin in ORIGINS for depth in DEPTHS]
    rows = [(graph, tool) for graph in GRAPHS for tool in TOOLS if not (graph == "PGGB" and tool == "vg")]
    display_rows = []
    row_positions = []
    group_centers = {}
    for graph in GRAPHS:
        graph_rows = [(g, t) for g, t in rows if g == graph]
        start = len(display_rows)
        for row in graph_rows:
            row_positions.append(len(display_rows))
            display_rows.append(row)
        group_centers[graph] = (start + len(display_rows) - 1) / 2
        if graph != GRAPHS[-1]:
            display_rows.append(None)

    matrix = np.full((len(display_rows), len(columns)), np.nan)
    for ri, row in enumerate(display_rows):
        if row is None:
            continue
        graph, tool = row
        for ci, (origin, depth) in enumerate(columns):
            result = means.get((graph, tool, origin, depth))
            if result and "f1" in result:
                matrix[ri, ci] = result["f1"]

    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("white")
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(columns)))
    ax.set_yticks(range(len(rows)))
    if blank:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.tick_params(axis="both", length=0)
    else:
        ax.set_xticklabels([d.replace("x", "×") for _, d in columns], fontsize=17)
        ax.set_yticks(row_positions)
        ax.set_yticklabels([TOOL_LABELS[tool] for _, tool in rows], fontsize=18, fontweight="bold")
        ax.tick_params(axis="both", length=0, pad=8)
        ax.set_xlabel("Sequencing depth", fontsize=20, labelpad=12)
        for tick, (_, tool) in zip(ax.get_yticklabels(), rows):
            tick.set_color(COLORS[tool])

    # Group headings make the repeated depth columns easier to scan.
    if not blank:
        for x, label in [(1.5, "Graph-in"), (5.5, "Graph-out")]:
            ax.text(x, 1.015, label, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=17, fontweight="bold", clip_on=False)
        for label, center in group_centers.items():
            ax.text(-0.30, center, label, transform=ax.get_yaxis_transform(), ha="center", va="center", rotation=90, fontsize=17, fontweight="bold", clip_on=False)
    ax.axvline(3.5, color="white", linewidth=2.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.025)
    if blank:
        cbar.set_ticks([])
        cbar.outline.set_visible(False)
    else:
        cbar.set_label("")
        cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])
        cbar.ax.tick_params(labelsize=17, length=7, width=1.5, colors="#222222")
        cbar.outline.set_visible(False)
    return ax


def draw_pr_panels(fig, means, blank=False, bounds=(0.57, 0.985, 0.28, 0.84)):
    left, right, bottom, top = bounds
    sub = fig.add_gridspec(1, 3, left=left, right=right, bottom=bottom, top=top, wspace=0.22)
    axes = []
    for i, graph in enumerate(GRAPHS):
        ax = fig.add_subplot(sub[0, i])
        axes.append(ax)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("auto")
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        if blank:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.tick_params(axis="both", length=0)
        else:
            ax.tick_params(axis="both", labelsize=17, length=8, width=1.8, colors="#222222")
            ax.set_title(graph, fontsize=17, fontweight="bold", pad=13)
            if i == 0:
                ax.set_ylabel("Recall", fontsize=20, labelpad=12)
            else:
                ax.tick_params(axis="y", labelleft=False)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#222222")
        ax.spines["bottom"].set_color("#222222")
        ax.spines["left"].set_linewidth(1.8)
        ax.spines["bottom"].set_linewidth(1.8)
        # Thin F1 iso-lines are a visual guide to the precision-recall trade-off.
        p = np.linspace(0.03, 1, 300)
        for f1 in (0.25, 0.5, 0.75):
            r = (f1 * p) / (2 * p - f1)
            valid = (2 * p > f1) & (r >= 0) & (r <= 1)
            ax.plot(p[valid], r[valid], color="#A8A8A8", linewidth=1.0, linestyle=(0, (2, 3)), zorder=0)

        for tool in TOOLS:
            for origin in ORIGINS:
                for depth in DEPTHS:
                    result = means.get((graph, tool, origin, depth))
                    if not result or "precision" not in result or "recall" not in result:
                        continue
                    ax.scatter(
                        result["precision"], result["recall"],
                        marker=DEPTH_MARKERS[depth],
                        s=132,
                        facecolors=COLORS[tool] if origin == "Graph-in" else "white",
                        edgecolors=COLORS[tool],
                        linewidths=2.4,
                        zorder=3,
                    )
    return axes


def main():
    means = load_means()
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(21, 10.5), facecolor="white")
    fig.text(0.018, 0.965, "A", fontsize=30, fontweight="bold", ha="left", va="top")
    fig.text(0.545, 0.965, "B", fontsize=30, fontweight="bold", ha="left", va="top")
    draw_heatmap(fig, means)
    draw_pr_panels(fig, means)

    # Keep the PR-panel annotations below panel B, with separate rows so that
    # the labels remain legible in the large print figure.
    fig.text(0.78, 0.215, "Precision", fontsize=20,
             ha="center", va="center")

    tool_handles = [
        Patch(facecolor=COLORS[tool], edgecolor="none", label=TOOL_LABELS[tool])
        for tool in TOOLS
    ]
    origin_handles = [
        Line2D([0], [0], marker="o", linestyle="None",
               markerfacecolor="#555555", markeredgecolor="#555555",
               markersize=14, label="Graph-in"),
        Line2D([0], [0], marker="o", linestyle="None",
               markerfacecolor="white", markeredgecolor="#555555",
               markeredgewidth=2.0, markersize=14, label="Graph-out"),
    ]
    depth_handles = [
        Line2D([0], [0], marker=DEPTH_MARKERS[depth], linestyle="None",
               color="#444444", markerfacecolor="#777777",
               markeredgecolor="#444444", markersize=14,
               label=depth.replace("x", "×"))
        for depth in DEPTHS
    ]
    fig.legend(handles=tool_handles, loc="lower center",
               bbox_to_anchor=(0.78, 0.145), ncol=4, frameon=False,
               fontsize=19, columnspacing=1.6, handletextpad=0.55,
               handlelength=1.5)
    fig.legend(handles=origin_handles, loc="lower center",
               bbox_to_anchor=(0.61, 0.065), ncol=2, frameon=False,
               fontsize=19, columnspacing=1.5, handletextpad=0.55,
               handlelength=1.4)
    fig.legend(handles=depth_handles, loc="lower center",
               bbox_to_anchor=(0.88, 0.065), ncol=4, frameon=False,
               fontsize=19, columnspacing=1.0, handletextpad=0.45,
               handlelength=1.3)

    fig.savefig(OUT_PNG_LABELED, dpi=3000, facecolor="white", bbox_inches="tight", pad_inches=0.18)
    fig.savefig(OUT_SVG_LABELED, format="svg", facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    blank_fig = plt.figure(figsize=(19, 9.2), facecolor="white")
    draw_heatmap(blank_fig, means, blank=True)
    draw_pr_panels(blank_fig, means, blank=True)
    blank_fig.subplots_adjust(left=0.055, right=0.985, bottom=0.045, top=0.955, wspace=0.13)
    blank_fig.savefig(OUT_SVG, format="svg", facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(blank_fig)
    print(OUT)
    print(OUT_SVG)
    print(OUT_SVG_LABELED)
    print(OUT_PNG_LABELED)


if __name__ == "__main__":
    main()
