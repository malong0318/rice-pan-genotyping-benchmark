#!/usr/bin/env python3
"""Export the aggregated data used by the graph-origin SNP figure."""

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path("/data/malong/benchmark")
INPUT = ROOT / "paper/result_statistics.tsv"
OUTPUT = ROOT / "paper/png/figure_data/figure3_data.tsv"

DEPTHS = ["5x", "10x", "15x", "20x"]
TOOLS = ["vg", "graphtyper", "pangenie", "bayestyper"]
SAMPLE_GROUPS = {
    "Graph-in": ["W0169", "W1735", "W3033"],
    "Graph-out": ["9311", "IR64", "zhonghua11"],
}


def main():
    sample_to_group = {
        sample: group
        for group, samples in SAMPLE_GROUPS.items()
        for sample in samples
    }
    grouped = defaultdict(lambda: {"f1": [], "precision": [], "recall": [], "samples": []})

    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["variant_class"] != "snp" or row["status"] != "complete":
                continue
            if row["graph"] not in {"MC", "PGGB"} or row["sample"] not in sample_to_group:
                continue
            group = sample_to_group[row["sample"]]
            key = (row["graph"], group, row["tool"], row["depth"])
            for metric in ("f1", "precision", "recall"):
                if row[metric] not in ("", "NA", "."):
                    grouped[key][metric].append(float(row[metric]))
            grouped[key]["samples"].append(row["sample"])

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "graph", "sample_origin", "tool", "depth", "n_samples", "samples",
        "f1_mean", "precision_mean", "recall_mean",
    ]
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for graph in ("MC", "PGGB"):
            for group in ("Graph-in", "Graph-out"):
                for tool in TOOLS:
                    for depth in DEPTHS:
                        key = (graph, group, tool, depth)
                        values = grouped.get(key)
                        if not values or not values["f1"] or not values["precision"] or not values["recall"]:
                            continue
                        writer.writerow({
                            "graph": graph,
                            "sample_origin": group,
                            "tool": tool,
                            "depth": depth,
                            "n_samples": len(values["f1"]),
                            "samples": ",".join(sorted(values["samples"])),
                            "f1_mean": f"{np.mean(values['f1']):.8f}",
                            "precision_mean": f"{np.mean(values['precision']):.8f}",
                            "recall_mean": f"{np.mean(values['recall']):.8f}",
                        })
    print(OUTPUT)


if __name__ == "__main__":
    main()
