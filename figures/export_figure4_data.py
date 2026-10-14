#!/usr/bin/env python3
"""Export sample-level and graph-origin mean data used for Figure 4."""

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path("/data/malong/benchmark")
INPUT = ROOT / "paper/result_statistics.tsv"
OUTPUT = ROOT / "paper/png/figure_data/figure4_data.tsv"
SAMPLES = ["9311", "IR64", "W0169", "W1735", "W3033", "zhonghua11"]
DEPTHS = ["5x", "10x", "15x", "20x"]
GRAPHS = ["MC", "PGGB"]
TOOLS = ["vg", "graphtyper", "pangenie", "bayestyper"]
IN_SAMPLES = {"W0169", "W1735", "W3033"}
ORIGIN_ORDER = {"Graph-in": 0, "Graph-out": 1}
GRAPH_ORDER = {"MC": 0, "PGGB": 1}


def origin(sample):
    return "Graph-in" if sample in IN_SAMPLES else "Graph-out"


def main():
    rows = []
    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["variant_class"] != "indel" or row["status"] != "complete":
                continue
            if row["graph"] not in GRAPHS or row["tool"] not in TOOLS:
                continue
            rows.append({
                "graph": row["graph"],
                "sample_origin": origin(row["sample"]),
                "tool": row["tool"],
                "sample": row["sample"],
                "depth": row["depth"],
                "f1": row["f1"],
                "precision": row["precision"],
                "recall": row["recall"],
            })

    output_rows = []
    for row in rows:
        output_rows.append({
            "record_level": "sample",
            **row,
            "n_samples": "1",
            "samples": row["sample"],
            "status": "complete",
            "note": "Sample-level result from result_statistics.tsv",
        })

    grouped = defaultdict(list)
    for row in rows:
        key = (row["graph"], row["sample_origin"], row["tool"], row["depth"])
        grouped[key].append(row)

    for (graph, sample_origin, tool, depth), group in grouped.items():
        output_rows.append({
            "record_level": "group_mean",
            "graph": graph,
            "sample_origin": sample_origin,
            "tool": tool,
            "sample": "ALL",
            "depth": depth,
            "f1": f"{sum(float(r['f1']) for r in group) / len(group):.8f}",
            "precision": f"{sum(float(r['precision']) for r in group) / len(group):.8f}",
            "recall": f"{sum(float(r['recall']) for r in group) / len(group):.8f}",
            "n_samples": str(len(group)),
            "samples": ",".join(sorted(r["sample"] for r in group)),
            "status": "complete",
            "note": "Mean across sample-level results; this is the value plotted in Figure 4",
        })

    output_rows.sort(key=lambda r: (
        r["record_level"],
        GRAPH_ORDER[r["graph"]],
        ORIGIN_ORDER[r["sample_origin"]],
        TOOLS.index(r["tool"]),
        SAMPLES.index(r["sample"]) if r["sample"] in SAMPLES else -1,
        DEPTHS.index(r["depth"]),
    ))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "record_level", "graph", "sample_origin", "tool", "sample", "depth",
        "f1", "precision", "recall", "n_samples", "samples", "status", "note",
    ]
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(OUTPUT)
    print(f"sample_rows={sum(r['record_level'] == 'sample' for r in output_rows)}")
    print(f"group_mean_rows={sum(r['record_level'] == 'group_mean' for r in output_rows)}")


if __name__ == "__main__":
    main()
