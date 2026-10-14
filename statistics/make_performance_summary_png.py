#!/usr/bin/env python3
"""Create one table-style PNG from the completed performance runs."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path("/data/malong/benchmark")
PERF_ROOT = ROOT / "results_performance"
DEFAULT_OUTPUT = PERF_ROOT / "performance_summary_20261002.png"


def raw_value(text: str, label: str) -> str:
    match = re.search(rf"^\s*{re.escape(label)}:\s*(.+)$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"missing metric: {label}")
    return match.group(1).strip()


def elapsed_seconds(value: str) -> float:
    parts = value.split(":")
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    if len(parts) == 2:
        minutes, seconds = parts
        return int(minutes) * 60 + float(seconds)
    return float(value)


def metric_from_tsv(path: Path, name: str) -> str:
    if not path.exists():
        return "NA"
    for line in path.read_text().splitlines():
        key, _, value = line.partition("\t")
        if key == name:
            return value
    return "NA"


def make_rows():
    graph_order = {"MC": 0, "PGGB": 1, "minigraph": 2}
    rows = []
    for raw in PERF_ROOT.glob("*/[A-Za-z]*/9311/5x/performance.raw.txt"):
        # The glob above is intentionally narrow; use its parent names for labels.
        tool = raw.parents[3].name
        graph = raw.parents[2].name
        text = raw.read_text()
        wall = elapsed_seconds(raw_value(text, "Elapsed (wall clock) time (h:mm:ss or m:ss)"))
        user = float(raw_value(text, "User time (seconds)"))
        system = float(raw_value(text, "System time (seconds)"))
        cpu = float(raw_value(text, "Percent of CPU this job got").rstrip("%"))
        rss_gb = int(raw_value(text, "Maximum resident set size (kbytes)")) / 1024 / 1024
        exit_status = raw_value(text, "Exit status")
        tsv = raw.with_name("performance.tsv")
        records = metric_from_tsv(tsv, "result_records")
        rows.append(
            {
                "tool": tool,
                "graph": graph,
                "threads": metric_from_tsv(tsv, "configured_threads"),
                "wall": wall,
                "cpu": cpu,
                "avg": (user + system) / wall if wall else 0.0,
                "rss": rss_gb,
                "records": records,
                "status": "OK" if exit_status == "0" else f"exit {exit_status}",
            }
        )
    rows.sort(key=lambda row: (row["tool"], graph_order.get(row["graph"], 99)))
    if len(rows) != 11:
        raise RuntimeError(f"expected 11 completed runs, found {len(rows)}")
    return rows


def main() -> None:
    rows = make_rows()
    columns = [
        "Tool",
        "Graph",
        "Threads",
        "Wall time (min)",
        "CPU (%)",
        "Avg CPU cores",
        "Max RSS (GB)",
        "VCF records",
        "Status",
    ]
    data = [
        [
            row["tool"],
            row["graph"],
            row["threads"],
            f'{row["wall"] / 60:.1f}',
            f'{row["cpu"]:.0f}',
            f'{row["avg"]:.2f}',
            f'{row["rss"]:.2f}',
            f'{int(row["records"]):,}' if row["records"] != "NA" else "NA",
            row["status"],
        ]
        for row in rows
    ]

    fig, ax = plt.subplots(figsize=(16, 6.4), dpi=220)
    ax.axis("off")
    table = ax.table(
        cellText=data,
        colLabels=columns,
        cellLoc="center",
        colLoc="center",
        loc="center",
        colWidths=[0.12, 0.11, 0.08, 0.13, 0.09, 0.12, 0.12, 0.14, 0.09],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.65)
    for cell in table.get_celld().values():
        cell.set_edgecolor("#B8C2CC")
        cell.set_linewidth(0.6)
    for col in range(len(columns)):
        cell = table[0, col]
        cell.set_facecolor("#234E70")
        cell.set_text_props(color="white", weight="bold")
    tool_colors = {"vg": "#E8F1F8", "pangenie": "#EDF5E8", "bayestyper": "#FFF1DC", "graphtyper": "#F3E8F7"}
    for row_index, row in enumerate(rows, start=1):
        color = tool_colors.get(row["tool"], "#F7F7F7")
        for col in range(len(columns)):
            table[row_index, col].set_facecolor(color)
    ax.set_title(
        "Genotyping performance benchmark: 9311 5x\n"
        "Wall time, CPU utilization, average CPU cores, maximum resident memory, and output records",
        fontsize=14,
        fontweight="bold",
        pad=18,
    )
    fig.text(
        0.01,
        0.02,
        "Metrics parsed from results_performance/*/*/9311/5x/performance.raw.txt; memory is maximum RSS.",
        fontsize=9,
        color="#555555",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    output = DEFAULT_OUTPUT
    if output.exists():
        output = output.with_name("performance_summary_20261002_rerun.png")
    fig.savefig(output, bbox_inches="tight")
    print(output)


if __name__ == "__main__":
    main()
