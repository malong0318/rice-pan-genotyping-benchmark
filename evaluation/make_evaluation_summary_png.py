#!/usr/bin/env python3
"""Render aligned PNG summary tables for genotyper evaluations."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("/data/malong/benchmark/evaluation")
SAMPLES = ["9311", "IR64", "W0169", "W1735", "W3033", "zhonghua11"]
DEPTHS = ["5x", "10x", "15x", "20x"]
GRAPH_ORDER = {"MC": 0, "PGGB": 1, "minigraph": 2}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sort_key(row: dict[str, str]) -> tuple[int, int, int]:
    sample = SAMPLES.index(row["sample"]) if row["sample"] in SAMPLES else len(SAMPLES)
    depth = DEPTHS.index(row["depth"]) if row["depth"] in DEPTHS else len(DEPTHS)
    graph = GRAPH_ORDER.get(row["graph"], len(GRAPH_ORDER))
    return sample, depth, graph


def int_value(value: object) -> str:
    return str(int(float(value)))


def float_value(value: object) -> str:
    return f"{float(value):.4f}"


def read_rtg_summary(path: Path) -> dict[str, str]:
    with path.open() as handle:
        for line in handle:
            fields = re.split(r"\s+", line.strip())
            if len(fields) >= 8 and fields[0] == "None":
                return {
                    "tp_base": fields[1],
                    "tp_call": fields[2],
                    "fp": fields[3],
                    "fn": fields[4],
                    "precision": fields[5],
                    "recall": fields[6],
                    "f1": fields[7],
                }
    raise ValueError(f"No None-threshold row found in {path}")


def common_fields(row: dict[str, str]) -> tuple[str, str, str, Path]:
    return row["sample"], row["depth"], row["graph"], Path(row["result"])


def build_small_table(kind: str, tool: str, suffix: str) -> tuple[list[str], list[list[str]], Path, str]:
    base = ROOT / kind / tool
    rows = read_tsv(base / "summary.tsv")
    rows.sort(key=sort_key)
    table = []
    for row in rows:
        sample, depth, graph, result = common_fields(row)
        metrics = read_rtg_summary(result / "summary.txt")
        table.append(
            [
                sample,
                depth,
                graph,
                int_value(row["truth_records"]),
                int_value(row["call_records"]),
                int_value(metrics["tp_base"]),
                int_value(metrics["tp_call"]),
                int_value(metrics["fp"]),
                int_value(metrics["fn"]),
                float_value(metrics["precision"]),
                float_value(metrics["recall"]),
                float_value(metrics["f1"]),
            ]
        )
    headers = [
        "Sample", "Depth", "Graph", "Truth n", "Call n", "TP-base",
        "TP-call", "FP", "FN", "Precision", "Recall", "F1",
    ]
    output = base / f"{tool}_{kind}_evaluation_summary{suffix}.png"
    return headers, table, output, f"{tool} {kind.upper()} evaluation summary"


def build_sv_table(tool: str, suffix: str) -> tuple[list[str], list[list[str]], Path, str]:
    base = ROOT / "sv" / tool
    rows = read_tsv(base / "summary.tsv")
    rows.sort(key=sort_key)
    table = []
    for row in rows:
        sample, depth, graph, result = common_fields(row)
        with (result / "summary.json").open() as handle:
            metrics = json.load(handle)
        truth_n = row.get("truth_records_excluding_BND", row.get("truth_records"))
        table.append(
            [
                sample,
                depth,
                graph,
                int_value(truth_n),
                int_value(row["call_records"]),
                int_value(metrics["TP-base"]),
                int_value(metrics["TP-comp"]),
                int_value(metrics["FP"]),
                int_value(metrics["FN"]),
                float_value(metrics["precision"]),
                float_value(metrics["recall"]),
                float_value(metrics["f1"]),
                float_value(metrics["gt_concordance"]),
            ]
        )
    headers = [
        "Sample", "Depth", "Graph", "Truth n", "Call n", "TP-base",
        "TP-comp", "FP", "FN", "Precision", "Recall", "F1", "GT concordance",
    ]
    output = base / f"{tool}_sv_evaluation_summary{suffix}.png"
    return headers, table, output, f"{tool} SV evaluation summary"


def render(headers: list[str], rows: list[list[str]], output: Path, title: str) -> None:
    if not rows:
        raise ValueError(f"No rows available for {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    width = 20 if len(headers) > 12 else 18
    height = max(5.0, 1.2 + len(rows) * 0.34)
    fig, ax = plt.subplots(figsize=(width, height), dpi=200)
    ax.axis("off")
    ax.set_title(title, fontsize=15, fontweight="bold", pad=16)
    widths = [0.095, 0.065, 0.085, 0.085, 0.085, 0.085, 0.085, 0.07, 0.07, 0.095, 0.085, 0.07, 0.12]
    widths = widths[: len(headers)]
    total = sum(widths)
    widths = [value / total for value in widths]
    table = ax.table(
        cellText=rows,
        colLabels=headers,
        colWidths=widths,
        cellLoc="center",
        loc="center",
        bbox=[0.005, 0.01, 0.99, 0.92],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.2 if len(rows) <= 24 else 7.2)
    table.scale(1, 1.15)
    for (row_idx, col_idx), cell in table.get_celld().items():
        cell.set_edgecolor("#B8C0CC")
        cell.set_linewidth(0.45)
        if row_idx == 0:
            cell.set_facecolor("#34495E")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("#F4F7FA" if row_idx % 2 else "white")
            if len(headers) > 2 and rows[row_idx - 1][2] == "minigraph":
                cell.set_facecolor("#FFF5E6" if row_idx % 2 else "#FFF9F0")
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {output} ({len(rows)} rows)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", choices=["vg", "pangenie", "bayestyper", "graphtyper"], required=True)
    parser.add_argument("--kind", choices=["snp", "indel", "sv"], required=True)
    parser.add_argument("--suffix", default="", help="suffix added before .png")
    args = parser.parse_args()
    if args.kind == "sv":
        render(*build_sv_table(args.tool, args.suffix))
    else:
        render(*build_small_table(args.kind, args.tool, args.suffix))


if __name__ == "__main__":
    main()
