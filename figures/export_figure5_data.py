#!/usr/bin/env python3
"""Export sample-level and graph-origin mean SV results for Figure 5."""

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path("/data/malong/benchmark")
INPUT = ROOT / "paper/result_statistics.tsv"
OUTPUT = ROOT / "paper/png/figure_data/figure5_data.tsv"
SAMPLES = ["9311", "IR64", "W0169", "W1735", "W3033", "zhonghua11"]
DEPTHS = ["5x", "10x", "15x", "20x"]
GRAPHS = ["MC", "PGGB", "minigraph"]
TOOLS = ["vg", "graphtyper", "pangenie", "bayestyper"]
IN_SAMPLES = {"W0169", "W1735", "W3033"}
GRAPH_ORDER = {graph: i for i, graph in enumerate(GRAPHS)}
ORIGIN_ORDER = {"Graph-in": 0, "Graph-out": 1}


def main():
    source = {}
    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["variant_class"] != "sv":
                continue
            if row["graph"] not in GRAPHS or row["tool"] not in TOOLS:
                continue
            origin = "Graph-in" if row["sample"] in IN_SAMPLES else "Graph-out"
            key = (row["graph"], row["tool"], row["sample"], row["depth"])
            source[key] = {
                "graph": row["graph"],
                "sample_origin": origin,
                "tool": row["tool"],
                "sample": row["sample"],
                "depth": row["depth"],
                "f1": row["f1"],
                "precision": row["precision"],
                "recall": row["recall"],
                "status": row["status"],
                "reason": row["reason"],
            }

    sample_rows = []
    grouped = defaultdict(list)
    for graph in GRAPHS:
        for tool in TOOLS:
            for sample in SAMPLES:
                for depth in DEPTHS:
                    key = (graph, tool, sample, depth)
                    record = source.get(key)
                    if record is None:
                        continue
                    sample_rows.append({
                        "record_level": "sample",
                        **record,
                        "n_samples": "1" if record["status"] == "complete" else "0",
                        "samples": sample,
                        "note": record["reason"] or "Sample-level result from result_statistics.tsv",
                    })
                    if record["status"] == "complete":
                        grouped[(graph, record["sample_origin"], tool, depth)].append(record)

    mean_rows = []
    for graph in GRAPHS:
        for origin in ("Graph-in", "Graph-out"):
            for tool in TOOLS:
                for depth in DEPTHS:
                    key = (graph, origin, tool, depth)
                    group = grouped.get(key, [])
                    if group:
                        means = {
                            metric: f"{sum(float(row[metric]) for row in group) / len(group):.8f}"
                            for metric in ("f1", "precision", "recall")
                        }
                        status = "complete"
                        note = "Mean across sample-level results; used in Figure 5 panels A and B"
                    else:
                        means = {metric: "NA" for metric in ("f1", "precision", "recall")}
                        status = "NA"
                        note = "No result available; VG-PGGB was not evaluated because graph indexing failed"
                    mean_rows.append({
                        "record_level": "group_mean",
                        "graph": graph,
                        "sample_origin": origin,
                        "tool": tool,
                        "sample": "ALL",
                        "depth": depth,
                        **means,
                        "status": status,
                        "reason": "",
                        "n_samples": str(len(group)),
                        "samples": ",".join(sorted(row["sample"] for row in group)),
                        "note": note,
                    })

    output_rows = sample_rows + mean_rows
    output_rows.sort(key=lambda row: (
        0 if row["record_level"] == "sample" else 1,
        GRAPH_ORDER[row["graph"]],
        TOOLS.index(row["tool"]),
        ORIGIN_ORDER[row["sample_origin"]],
        SAMPLES.index(row["sample"]) if row["sample"] in SAMPLES else -1,
        DEPTHS.index(row["depth"]),
    ))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "record_level", "graph", "sample_origin", "tool", "sample", "depth",
        "f1", "precision", "recall", "n_samples", "samples", "status", "reason", "note",
    ]
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    print(OUTPUT)
    print(f"sample_rows={len(sample_rows)}")
    print(f"group_mean_rows={len(mean_rows)}")
    print(f"complete_group_means={sum(row['status'] == 'complete' for row in mean_rows)}")
    print(f"NA_group_means={sum(row['status'] == 'NA' for row in mean_rows)}")


if __name__ == "__main__":
    main()
