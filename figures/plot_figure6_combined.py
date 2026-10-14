#!/usr/bin/env python3
"""Draw the aligned four-panel Figure 6 as one vector figure."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpecFromSubplotSpec


PAPER = Path("/data/malong/benchmark/paper")
OUT = PAPER / "png"
sys.path.insert(0, str(PAPER))
from plot_figure6_individual_panels import (  # noqa: E402
    COLORS,
    GRAPH_LABELS,
    GRAPH_MARKERS,
    GRAPHS,
    TOOLS,
    TOOL_LABELS,
    perf_rows,
    read_accuracy,
)


def boxed(ax):
    ax.set_facecolor("white")
    ax.grid(False)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("#222222")
        ax.spines[side].set_linewidth(2.6)
    ax.tick_params(axis="both", colors="#222222", width=2.4, length=10, pad=6)


def bold_ticks(ax, size=17):
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(size)
        label.set_fontweight("bold")


def draw_a(ax, rows, accuracy):
    boxed(ax)
    for row in rows:
        f1 = accuracy.get((row["tool"], row["graph"]))
        if f1 is None:
            continue
        ax.scatter(
            row["runtime_h"], f1,
            s=220,
            marker=GRAPH_MARKERS[row["graph"]],
            color=COLORS[row["tool"]],
            edgecolors="white",
            linewidths=2.4,
            zorder=3,
        )
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g} h"))
    ax.set_ylim(0, max(accuracy.values()) * 1.45)
    ax.set_ylabel("F1", fontsize=25, fontweight="bold", labelpad=12)
    ax.tick_params(axis="both", labelsize=19)
    bold_ticks(ax, 19)
    graph_handles = [
        Line2D(
            [0], [0], marker=GRAPH_MARKERS[g], linestyle="None", color="#555555",
            markersize=12, label=GRAPH_LABELS[g],
        )
        for g in GRAPHS
    ]
    ax.legend(
        handles=graph_handles, loc="upper right", frameon=False,
        fontsize=18, ncol=3, handletextpad=0.5, columnspacing=1.0,
    )


def draw_b(axes, rows):
    ymax = max(r["mean_cores"] for r in rows if r["mean_cores"] is not None) * 1.20
    x = np.arange(len(TOOLS))
    for index, (ax, graph) in enumerate(zip(axes, GRAPHS)):
        boxed(ax)
        present = {
            r["tool"]: r["mean_cores"]
            for r in rows
            if r["graph"] == graph and r["mean_cores"] is not None
        }
        for xi, tool in zip(x, TOOLS):
            value = present.get(tool)
            if value is None:
                continue
            ax.vlines(xi, 0, value, color=COLORS[tool], linewidth=4.0, zorder=2)
            ax.scatter(
                xi, value, s=230, color=COLORS[tool],
                edgecolors="white", linewidths=2.5, zorder=3,
            )
        ax.set_title(GRAPH_LABELS[graph], fontsize=23, fontweight="bold", pad=12)
        ax.set_xticks(x)
        ax.set_xticklabels([])
        ax.tick_params(axis="x", width=2.4, length=10, pad=0)
        ax.tick_params(axis="y", labelsize=17, width=2.4, length=10, pad=5)
        ax.set_ylim(0, ymax)
        if index == 0:
            ax.set_ylabel("CPU cores", fontsize=25, fontweight="bold", labelpad=13)
            bold_ticks(ax, 17)
        else:
            ax.tick_params(axis="y", labelleft=False)


def draw_c(ax, rows):
    boxed(ax)
    offsets = np.linspace(-0.24, 0.24, len(TOOLS))
    for yi, graph in enumerate(GRAPHS):
        present = [
            r for r in rows
            if r["graph"] == graph and r["memory_gb"] is not None
        ]
        present.sort(key=lambda r: TOOLS.index(r["tool"]))
        if len(present) > 1:
            ax.plot(
                [r["memory_gb"] for r in present],
                [yi + offsets[TOOLS.index(r["tool"])] for r in present],
                color="#B7B7B7", linewidth=2.0, zorder=1,
            )
        for row in present:
            ax.scatter(
                row["memory_gb"], yi + offsets[TOOLS.index(row["tool"])],
                s=210, color=COLORS[row["tool"]],
                edgecolors="white", linewidths=2.0, zorder=3,
            )
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g} GB"))
    ax.set_yticks(range(len(GRAPHS)), [GRAPH_LABELS[g] for g in GRAPHS])
    ax.tick_params(axis="both", labelsize=18, width=2.4, length=10, pad=6)
    bold_ticks(ax, 18)
    ax.set_ylim(-0.55, len(GRAPHS) - 0.45)


def draw_d(ax, rows):
    boxed(ax)
    x = np.arange(len(GRAPHS))
    width = 0.18
    offsets = (np.arange(len(TOOLS)) - 1.5) * width
    max_value = max(r["records_n"] for r in rows if r["records_n"] is not None)
    for offset, tool in zip(offsets, TOOLS):
        values = []
        for graph in GRAPHS:
            match = next(
                (r for r in rows
                 if r["tool"] == tool and r["graph"] == graph
                 and r["records_n"] is not None),
                None,
            )
            values.append(np.nan if match is None else match["records_n"])
        values = np.asarray(values, dtype=float)
        valid = ~np.isnan(values)
        ax.bar(
            x[valid] + offset, values[valid], width=width * 0.88,
            color=COLORS[tool], edgecolor="#222222", linewidth=0.8,
        )
    ax.set_xticks(x, [GRAPH_LABELS[g] for g in GRAPHS])
    ax.set_ylabel("VCF records (million)", fontsize=25, fontweight="bold", labelpad=13)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e6:g}"))
    ax.tick_params(axis="both", labelsize=18, width=2.4, length=10, pad=6)
    bold_ticks(ax, 18)
    ax.set_ylim(0, max_value * 1.16)


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
    rows = perf_rows()
    if len(rows) != 11:
        raise RuntimeError(f"Expected 11 complete performance rows, found {len(rows)}")
    accuracy = read_accuracy()

    fig = plt.figure(figsize=(18, 13), facecolor="white")
    gs = fig.add_gridspec(
        2, 2, left=0.115, right=0.985, top=0.955, bottom=0.19,
        wspace=0.36, hspace=0.34,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    bgs = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0, 1], wspace=0.14)
    ax_b = [fig.add_subplot(bgs[0, i]) for i in range(3)]

    draw_a(ax_a, rows, accuracy)
    draw_b(ax_b, rows)
    draw_c(ax_c, rows)
    draw_d(ax_d, rows)

    fig.canvas.draw()
    # Place B and D on one shared absolute x coordinate. Panel B contains
    # internal facets, so its first facet has a different local offset than D.
    def label_position(ax):
        display = ax.transAxes.transform((-0.12, 1.04))
        return fig.transFigure.inverted().transform(display)

    a_x, a_y = label_position(ax_a)
    c_x, c_y = label_position(ax_c)
    _, b_y = label_position(ax_b[0])
    d_x, d_y = label_position(ax_d)
    bd_x = d_x
    for x, y, label in ((a_x, a_y, "A"), (bd_x, b_y, "B"), (c_x, c_y, "C"), (bd_x, d_y, "D")):
        fig.text(
            x, y, label, fontsize=27, fontweight="bold",
            va="bottom", ha="left",
        )

    handles = [Patch(facecolor=COLORS[t], edgecolor="none", label=TOOL_LABELS[t]) for t in TOOLS]
    fig.legend(
        handles=handles, loc="lower center", ncol=4, frameon=False,
        bbox_to_anchor=(0.5, 0.055), fontsize=22,
        handlelength=0.9, handleheight=0.9, columnspacing=2.0,
    )

    svg = OUT / "Figure6_accuracy_and_resources_combined.svg"
    png = OUT / "Figure6_accuracy_and_resources_combined.png"
    fig.savefig(svg, format="svg", facecolor="white", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(png, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print(svg)
    print(png)


if __name__ == "__main__":
    main()
