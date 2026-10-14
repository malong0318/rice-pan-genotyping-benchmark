#!/usr/bin/env python3
"""Stream-classify graph-level candidate VCFs for Figure 2C/D."""

import csv
import gzip
import re
from pathlib import Path


DNA = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
ROOT = Path("/data/malong/benchmark")
OUT = ROOT / "paper/graph_candidate_variant_classification.tsv"

SOURCES = {
    "MC": ROOT / "graph/cactus_new/rice36_mc.vcf.gz",
    "PGGB": ROOT / "graph/pggb/rice36_graph/rice36.complete.smooth.final.NIP-T2T.vcf.gz",
    "minigraph": ROOT / "graph/minigraph/rice36.deconstruct/rice36.sv.vcf",
}


def classify(ref: str, alt: str, graph: str) -> str:
    if graph == "minigraph":
        # minigraph's deconstruct VCF uses symbolic graph alleles (for example
        # <X:1>) for its SV bubbles. Keep the minigraph panel intentionally
        # limited to SV versus records that cannot be interpreted as SVs.
        if (alt.startswith("<") and alt.endswith(">")):
            return "sv"
        if DNA.fullmatch(ref) and DNA.fullmatch(alt) and abs(len(ref) - len(alt)) >= 50:
            return "sv"
        return "other"
    if not DNA.fullmatch(ref) or not DNA.fullmatch(alt):
        return "other"
    if len(ref) == 1 and len(alt) == 1 and ref.upper() != alt.upper():
        return "snp"
    if len(ref) != len(alt):
        return "small_indel" if abs(len(ref) - len(alt)) < 50 else "sv"
    return "other"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def count_graph(graph: str, path: Path) -> dict[str, int | str]:
    counts = {"snp": 0, "small_indel": 0, "sv": 0, "other": 0}
    raw_records = 0
    decomposed_alts = 0
    multiallelic_records = 0
    with open_text(path) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n\r").split("\t", 5)
            if len(fields) < 5:
                raise ValueError(f"Malformed VCF record in {path}: {line[:120]}")
            raw_records += 1
            ref, alt_field = fields[3], fields[4]
            alts = alt_field.split(",")
            if len(alts) > 1:
                multiallelic_records += 1
            for alt in alts:
                decomposed_alts += 1
                counts[classify(ref, alt, graph)] += 1
    return {
        "graph": graph,
        "source_vcf": str(path),
        "raw_records": raw_records,
        "decomposed_alt_records": decomposed_alts,
        "multiallelic_raw_records": multiallelic_records,
        **counts,
    }


def main() -> None:
    rows = []
    for graph, path in SOURCES.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        print(f"classifying {graph}: {path}", flush=True)
        row = count_graph(graph, path)
        rows.append(row)
        print(
            f"{graph}: raw={row['raw_records']} split_alts={row['decomposed_alt_records']} "
            f"SNP={row['snp']} small_indel={row['small_indel']} "
            f"SV={row['sv']} other={row['other']}",
            flush=True,
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "graph",
        "source_vcf",
        "raw_records",
        "decomposed_alt_records",
        "multiallelic_raw_records",
        "snp",
        "small_indel",
        "sv",
        "other",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(OUT)


if __name__ == "__main__":
    main()
