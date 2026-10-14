#!/usr/bin/env python3
"""Create the combined four-panel Figure 2."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


ROOT = Path("/data/malong/benchmark")
GRAPH_STATS = ROOT / "paper/graph_statistics.tsv"
VARIANT_STATS = ROOT / "paper/graph_candidate_variant_classification.tsv"
OUT = ROOT / "paper/png/Figure2_graph_structure_and_candidate_variants.png"

GRAPH_ORDER = ["MC", "PGGB", "minigraph"]
GRAPH_COLORS = {"MC": "#2E7D5B", "PGGB": "#3568A8", "minigraph": "#D17A22"}
CATEGORY_ORDER = ["snp", "small_indel", "sv", "other"]
CATEGORY_LABELS = {"snp": "SNP", "small_indel": "small indel", "sv": "SV", "other": "other"}
CATEGORY_COLORS = {"snp": "#3B78B4", "small_indel": "#E39A3B", "sv": "#B84A4A", "other": "#A6A6A6"}


def load_table(path: Path):
    with path.open(newline="") as handle:
        return {row["graph"]: row for row in csv.DictReader(handle, delimiter="\t")}


def comma(value: int) -> str:
    return f"{value:,}"


def millions(value: int) -> str:
    return f"{value / 1e6:.1f}M"


def style_axis(ax):
    ax.set_facecolor("white")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    graph = load_table(GRAPH_STATS)
    variants = load_table(VARIANT_STATS)

    file_gb = [int(graph[g]["graph_bytes"]) / 1e9 for g in GRAPH_ORDER]
    nodes = [int(graph[g]["gfa_nodes"]) for g in GRAPH_ORDER]
    edges = [int(graph[g]["gfa_edges"]) for g in GRAPH_ORDER]
    totals = {
        g: sum(int(variants[g][category]) for category in CATEGORY_ORDER)
        for g in GRAPH_ORDER
    }
    fractions = {
        g: [int(variants[g][category]) / totals[g] * 100 if totals[g] else 0 for category in CATEGORY_ORDER]
        for g in GRAPH_ORDER
    }

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.labelsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=150)
    fig.patch.set_facecolor("white")
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # A: uncompressed GFA file size.
    x = list(range(len(GRAPH_ORDER)))
    ax_a.bar(x, file_gb, width=0.62, color=[GRAPH_COLORS[g] for g in GRAPH_ORDER], edgecolor="black", linewidth=0.8)
    for idx, value in enumerate(file_gb):
        ax_a.text(idx, value + 0.65, f"{value:.2f}", ha="center", va="bottom", fontsize=12)
    ax_a.set_xticks(x, GRAPH_ORDER)
    ax_a.set_ylim(0, max(file_gb) * 1.15)
    ax_a.set_ylabel("GB")
    style_axis(ax_a)

    # B: nodes and edges, shown on a log scale.
    width = 0.34
    left = [v - width / 2 for v in x]
    right = [v + width / 2 for v in x]
    bars_nodes = ax_b.bar(left, nodes, width=width, color="#4B78BD", edgecolor="black", linewidth=0.7, label="Nodes")
    bars_edges = ax_b.bar(right, edges, width=width, color="#ED7D31", edgecolor="black", linewidth=0.7, label="Edges")
    for bars, values in [(bars_nodes, nodes), (bars_edges, edges)]:
        for bar, value in zip(bars, values):
            ax_b.text(bar.get_x() + bar.get_width() / 2, value * 1.18, millions(value), ha="center", va="bottom", fontsize=11)
    ax_b.set_xticks(x, GRAPH_ORDER)
    ax_b.set_yscale("log")
    ax_b.set_ylim(6e5, 3.5e8)
    ax_b.set_ylabel("Count")
    ax_b.legend(frameon=False, loc="upper right", fontsize=12)
    style_axis(ax_b)

    # C: total ALT alleles after multiallelic decomposition.
    y = list(range(len(GRAPH_ORDER)))
    for idx, graph_name in enumerate(GRAPH_ORDER):
        value = totals[graph_name]
        ax_c.hlines(idx, 1e5, value, color=GRAPH_COLORS[graph_name], linewidth=3.0, alpha=0.9)
        ax_c.scatter(value, idx, s=100, color=GRAPH_COLORS[graph_name], edgecolor="white", linewidth=0.9, zorder=3)
        ax_c.text(value * 1.16, idx, comma(value), ha="left", va="center", fontsize=11)
    ax_c.set_yticks(y, GRAPH_ORDER)
    ax_c.invert_yaxis()
    ax_c.set_xscale("log")
    ax_c.set_xlim(1e5, max(totals.values()) * 3.0)
    ax_c.set_xticks([1e5, 1e6, 1e7], [r"$10^5$", r"$10^6$", r"$10^7$"])
    ax_c.set_xlabel("Decomposed ALT alleles")
    style_axis(ax_c)

    # D: composition of the same decomposed ALT counts.
    left_values = [0.0] * len(GRAPH_ORDER)
    for category in CATEGORY_ORDER:
        values = [fractions[g][CATEGORY_ORDER.index(category)] for g in GRAPH_ORDER]
        ax_d.barh(y, values, left=left_values, height=0.55, color=CATEGORY_COLORS[category], label=CATEGORY_LABELS[category])
        left_values = [a + b for a, b in zip(left_values, values)]
    ax_d.set_yticks(y, GRAPH_ORDER)
    ax_d.invert_yaxis()
    ax_d.set_xlim(0, 115)
    ax_d.set_xticks([0, 25, 50, 75, 100], ["0", "25", "50", "75", "100%"])
    ax_d.set_xlabel("Proportion")
    for idx, graph_name in enumerate(GRAPH_ORDER):
        ax_d.text(101.5, idx, f"n={comma(totals[graph_name])}", va="center", ha="left", fontsize=11)
    style_axis(ax_d)

    handles = [Patch(facecolor=CATEGORY_COLORS[c], edgecolor="none", label=CATEGORY_LABELS[c]) for c in CATEGORY_ORDER]
    fig.legend(handles=handles, frameon=False, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.015), fontsize=12)
    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.13, top=0.96, wspace=0.28, hspace=0.30)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, facecolor="white")
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
