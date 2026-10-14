#!/usr/bin/env python3
"""Export the data used for the final combined Figure 6."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path("/data/malong/benchmark/paper")
PERF = ROOT / "performance_statistics.tsv"
RESULTS = ROOT / "result_statistics.tsv"
OUT = ROOT / "png/figure_data/figure6_data.tsv"
TOOLS = ["vg", "graphtyper", "pangenie", "bayestyper"]
GRAPHS = ["MC", "PGGB", "minigraph"]
CLASSES = ["snp", "indel", "sv"]


def number(value: str) -> float | None:
    if value in ("", "NA", "."):
        return None
    return float(value)


def fmt(value: float | None, digits: int = 6) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def main() -> None:
    performance = {}
    with PERF.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] != "9311" or row["depth"] != "5x":
                continue
            if row["tool"] not in TOOLS or row["graph"] not in GRAPHS:
                continue
            performance[(row["tool"], row["graph"])] = row

    f1_values: defaultdict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    with RESULTS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] != "9311" or row["depth"] != "5x":
                continue
            if row["tool"] not in TOOLS or row["graph"] not in GRAPHS:
                continue
            f1 = number(row["f1"])
            if f1 is not None:
                f1_values[(row["tool"], row["graph"])][row["variant_class"]] = f1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "tool", "graph", "sample", "depth", "status",
        "runtime_h", "f1_mean", "f1_snp", "f1_indel", "f1_sv",
        "mean_cpu_cores", "peak_memory_gb", "vcf_records", "vcf_records_million",
        "configured_threads", "note",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for tool in TOOLS:
            for graph in GRAPHS:
                perf = performance.get((tool, graph), {})
                values = f1_values.get((tool, graph), {})
                mean = sum(values.values()) / len(values) if values else None
                records = number(perf.get("result_records", ""))
                status = perf.get("status", "NA")
                note = ""
                if tool == "vg" and graph == "PGGB":
                    note = "VG-PGGB VG index construction failed; missing values are NA, not zero"
                elif graph == "minigraph":
                    note = "Mean F1 uses available SV evaluation only"
                writer.writerow({
                    "tool": tool,
                    "graph": graph,
                    "sample": "9311",
                    "depth": "5x",
                    "status": status,
                    "runtime_h": fmt(number(perf.get("elapsed_seconds", "")) / 3600 if number(perf.get("elapsed_seconds", "")) is not None else None, 6),
                    "f1_mean": fmt(mean, 6),
                    "f1_snp": fmt(values.get("snp"), 6),
                    "f1_indel": fmt(values.get("indel"), 6),
                    "f1_sv": fmt(values.get("sv"), 6),
                    "mean_cpu_cores": fmt(number(perf.get("average_cpu_cores", "")), 3),
                    "peak_memory_gb": fmt(number(perf.get("max_rss_gb", "")), 3),
                    "vcf_records": fmt(records, 0),
                    "vcf_records_million": fmt(records / 1e6 if records is not None else None, 6),
                    "configured_threads": perf.get("configured_threads", "NA"),
                    "note": note,
                })
    print(OUT)


if __name__ == "__main__":
    main()
