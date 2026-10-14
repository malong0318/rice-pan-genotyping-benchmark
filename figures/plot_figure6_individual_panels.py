#!/usr/bin/env python3
"""Draw independent panels for Figure 6 performance/resource trade-offs."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path("/data/malong/benchmark")
PERF = ROOT / "paper/performance_statistics.tsv"
RESULTS = ROOT / "paper/result_statistics.tsv"
OUT = ROOT / "paper/png"

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
GRAPHS = ["MC", "PGGB", "minigraph"]
GRAPH_LABELS = {"MC": "MC", "PGGB": "PGGB", "minigraph": "minigraph"}
GRAPH_MARKERS = {"MC": "o", "PGGB": "s", "minigraph": "^"}


def read_performance() -> list[dict[str, str]]:
    with PERF.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_accuracy() -> dict[tuple[str, str], float]:
    values: defaultdict[tuple[str, str], list[float]] = defaultdict(list)
    with RESULTS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] != "9311" or row["depth"] != "5x":
                continue
            if row["status"] != "complete" or row["tool"] not in TOOLS:
                continue
            if row["graph"] not in GRAPHS or row["f1"] in ("", "NA", "."):
                continue
            values[(row["tool"], row["graph"])].append(float(row["f1"]))
    return {key: float(np.mean(vals)) for key, vals in values.items()}


def numeric(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    if value in ("", "NA", "."):
        return None
    return float(value)


def perf_rows() -> list[dict]:
    rows = []
    for row in read_performance():
        if row["status"] != "complete":
            continue
        if row["tool"] not in TOOLS or row["graph"] not in GRAPHS:
            continue
        row = dict(row)
        row["runtime_h"] = float(row["elapsed_seconds"]) / 3600.0
        row["cpu_percent_n"] = numeric(row, "cpu_percent")
        row["mean_cores"] = numeric(row, "average_cpu_cores")
        row["memory_gb"] = numeric(row, "max_rss_gb")
        row["records_n"] = numeric(row, "result_records")
        rows.append(row)
    return rows


def style_axes(ax):
    ax.set_facecolor("white")
    ax.grid(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#222222")
        ax.spines[side].set_linewidth(1.8)
    ax.tick_params(axis="both", labelsize=15, width=1.5, length=6, colors="#222222")


def save(fig, stem: str, dpi: int = 600):
    png = OUT / f"Figure6_{stem}.png"
    svg = OUT / f"Figure6_{stem}.svg"
    fig.savefig(png, dpi=dpi, facecolor="white", bbox_inches="tight", pad_inches=0.12)
    fig.savefig(svg, format="svg", facecolor="white", bbox_inches="tight", pad_inches=0.10)
    plt.close(fig)
    print(png)
    print(svg)


def panel_a(rows, accuracy, blank=False):
    fig, ax = plt.subplots(figsize=(10.5, 7.6), facecolor="white")
    style_axes(ax)
    # Panel A uses a fully boxed plotting area, matching the reference style.
    for side in ("top", "right"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("#222222")
        ax.spines[side].set_linewidth(2.6)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(2.6)
    for row in rows:
        f1 = accuracy.get((row["tool"], row["graph"]))
        if f1 is None:
            continue
        ax.scatter(
            row["runtime_h"], f1,
            s=220,
            marker=GRAPH_MARKERS[row["graph"]],
            c=COLORS[row["tool"]], alpha=1.0,
            edgecolors="white", linewidths=2.4, zorder=3,
        )
    ax.set_xscale("log")
    ax.set_ylim(0, max(accuracy.values()) * 1.45)
    if not blank:
        ax.set_xlabel("Runtime (h)", fontsize=27, fontweight="bold", labelpad=14)
        ax.set_ylabel("F1", fontsize=27, fontweight="bold", labelpad=14)
    ax.tick_params(axis="both", labelsize=22, width=2.4, length=10, pad=6)
    if blank:
        ax.tick_params(axis="both", labelbottom=False, labelleft=False)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(22)
        label.set_fontweight("bold")
    if not blank:
        graph_handles = [Line2D([0], [0], marker=GRAPH_MARKERS[g], linestyle="None", color="#555555", markersize=13, label=GRAPH_LABELS[g]) for g in GRAPHS]
        ax.legend(handles=graph_handles, loc="upper right", frameon=False, fontsize=21, ncol=3, handletextpad=0.6, columnspacing=1.3, borderpad=0.4)
    fig.tight_layout()
    save(fig, f"panel_A_accuracy_runtime{'_blank' if blank else ''}", dpi=3000)


def panel_b(rows, blank=False):
    fig, axes = plt.subplots(
        1, 3,
        figsize=(15.5, 7.6),
        sharey=True,
        facecolor="white",
        gridspec_kw={"wspace": 0.14},
    )
    ymax = max(r["mean_cores"] for r in rows if r["mean_cores"] is not None) * 1.20
    x = np.arange(len(TOOLS))

    for panel_index, (ax, graph) in enumerate(zip(axes, GRAPHS)):
        style_axes(ax)
        # Match Panel A's fully boxed plotting area and heavier visual treatment.
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color("#222222")
            ax.spines[side].set_linewidth(2.6)

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
                xi,
                value,
                s=250,
                color=COLORS[tool],
                edgecolors="white",
                linewidths=2.6,
                zorder=3,
            )

        if not blank:
            ax.set_title(GRAPH_LABELS[graph], fontsize=27, fontweight="bold", pad=16)
        ax.set_xticks(x)
        if not blank:
            ax.set_xticklabels(
                [TOOL_LABELS[t] for t in TOOLS],
                fontsize=17,
                fontweight="bold",
                rotation=35,
                ha="center",
                va="top",
                rotation_mode="anchor",
            )
        ax.tick_params(axis="x", width=2.4, length=10, pad=0)
        ax.tick_params(axis="y", labelsize=22, width=2.4, length=10, pad=7)
        if blank:
            ax.tick_params(axis="both", labelbottom=False, labelleft=False)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight("bold")
        ax.set_ylim(0, ymax)
        if panel_index == 0 and not blank:
            ax.set_ylabel("CPU cores", fontsize=27, fontweight="bold", labelpad=16)
        elif panel_index != 0:
            ax.tick_params(axis="y", labelleft=False)

    fig.subplots_adjust(left=0.11, right=0.99, top=0.88, bottom=0.24)
    save(fig, f"panel_B_mean_cpu_cores{'_blank' if blank else ''}", dpi=3000)


def panel_c(rows, blank=False):
    fig, ax = plt.subplots(figsize=(11.5, 7.6), facecolor="white")
    style_axes(ax)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("#222222")
        ax.spines[side].set_linewidth(2.6)

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
                color="#B7B7B7",
                linewidth=2.0,
                zorder=1,
            )
        for r in present:
            ax.scatter(
                r["memory_gb"],
                yi + offsets[TOOLS.index(r["tool"])],
                s=220,
                color=COLORS[r["tool"]],
                edgecolors="white",
                linewidths=2.0,
                zorder=3,
            )

    ax.set_xscale("log")
    ax.set_yticks(range(len(GRAPHS)))
    if not blank:
        ax.set_yticklabels([GRAPH_LABELS[g] for g in GRAPHS])
        ax.set_xlabel("GB", fontsize=27, fontweight="bold", labelpad=14)
    ax.tick_params(axis="both", labelsize=22, width=2.4, length=10, pad=7)
    if blank:
        ax.tick_params(axis="both", labelbottom=False, labelleft=False)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    ax.set_ylim(-0.55, len(GRAPHS) - 0.45)
    fig.tight_layout()
    save(fig, f"panel_C_peak_memory{'_blank' if blank else ''}", dpi=3000)


def panel_d(rows, blank=False):
    fig, ax = plt.subplots(figsize=(11.5, 7.6), facecolor="white")
    style_axes(ax)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("#222222")
        ax.spines[side].set_linewidth(2.6)

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
            x[valid] + offset,
            values[valid],
            width=width * 0.88,
            color=COLORS[tool],
            label=TOOL_LABELS[tool] if not blank else "_nolegend_",
            edgecolor="#222222",
            linewidth=0.8,
        )

    ax.set_xticks(x)
    if not blank:
        ax.set_xticklabels([GRAPH_LABELS[g] for g in GRAPHS])
        ax.set_ylabel("VCF records (million)", fontsize=27, fontweight="bold", labelpad=16)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value / 1e6:g}"))
    ax.tick_params(axis="both", labelsize=22, width=2.4, length=10, pad=7)
    if blank:
        ax.tick_params(axis="both", labelbottom=False, labelleft=False)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    ax.set_ylim(0, max_value * 1.16)
    fig.tight_layout()
    save(fig, f"panel_D_vcf_record_count{'_blank' if blank else ''}", dpi=3000)


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
    OUT.mkdir(parents=True, exist_ok=True)
    rows = perf_rows()
    if len(rows) != 11:
        raise RuntimeError(f"Expected 11 complete performance rows, found {len(rows)}")
    accuracy = read_accuracy()
    panel_a(rows, accuracy)
    panel_b(rows)
    panel_c(rows)
    panel_d(rows)
    panel_a(rows, accuracy, blank=True)
    panel_b(rows, blank=True)
    panel_c(rows, blank=True)
    panel_d(rows, blank=True)


if __name__ == "__main__":
    main()
