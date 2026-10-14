# Benchmark Code

This directory contains the analysis code intended for the benchmark's public code archive. It contains scripts only; reference assemblies, FASTQ/BAM files, graph files, VCFs, evaluation results, software binaries, logs, and private cluster configuration are not included.

## Contents

| Directory | Contents |
|---|---|
| `graph/` | MC Minigraph-Cactus reconstruction script and the historical whole-genome PGGB build script |
| `reads/` | Paired-end read downsampling by target depth |
| `truth/` | Assembly-derived truth-set workflow and chromosome/reference harmonization |
| `vcf_processing/` | Unified variant classification and tool-output evaluation-input preparation |
| `genotyping/` | Full-read calling scripts for VG, PanGenie, GraphTyper2, and BayesTyper |
| `graph_indexes/` | VG index-building scripts for the graph inputs |
| `evaluation/` | RTG vcfeval/Truvari evaluation driver and evaluation-summary plotting helper |
| `performance/` | 11 performance-test launchers and shared resource-recording code |
| `statistics/` | Metadata, accuracy/performance summary generation |
| `figures/` | Scripts used to export figure data and draw Figures 2-6 |

## Inputs and execution

The scripts were developed in the project environment and several retain `/data/malong/benchmark` and `/data/malong/miniconda3` as defaults. Before execution outside that environment, set the paths expected by each script or edit its configuration block. The FASTA/FASTQ/VCF/BAM inputs must be obtained separately from the sources listed in the manuscript and supplementary manifests. Intermediate files and large benchmark outputs are deliberately excluded from this code directory.

The typical workflow is: construct/download the graph inputs and downsample reads; generate and harmonize truth sets; build graph indexes where needed; run the four genotypers; classify calls and prepare evaluation inputs; run RTG vcfeval for SNP/small-indel and Truvari for SV; then generate summary tables and figures. The included VG scripts cover MC and minigraph; VG on PGGB is intentionally absent because its graph index could not be built in this study. The original minigraph graph-construction command was not recovered, so the directory includes its available indexing/deconstruction steps but does not claim to reproduce that graph from source. The MC graph rebuild script takes an existing two-column Cactus seqfile via `SEQFILE`; update its FASTA paths before running.

## MC graph reconstruction

The MC graph used in the completed benchmark is an existing historical output. Its original Cactus version and exact invocation were not preserved, and the 36 source assembly FASTAs are no longer present at their recorded local paths. Therefore `graph/build_mc_cactus_pangenome.sh` is a reconstruction recipe, not a guarantee of byte-identical recovery of the analyzed graph.

Prepare a valid two-column Cactus seqfile whose 36 genome rows include the sample `NIP-T2T`, and update all FASTA paths. The seqfile retained in the project also has an additional `_MINIGRAPH_` row; the reconstruction script permits that row but its role in the historical construction could not be independently confirmed. To reproduce the exact input set, retain it only after verifying the corresponding sequence and intended use. Then run in a fresh output location:

```bash
SEQFILE=/path/to/rice36.seqfile \
OUTDIR=/path/to/new_mc_output \
JOBSTORE=/path/to/new_mc_jobstore \
REFERENCE=NIP-T2T THREADS=128 \
bash graph/build_mc_cactus_pangenome.sh
```

The script requests GFA, VCF, GBZ, and Giraffe outputs, logs the full command, validates the 36 genome rows and input files, and refuses to overwrite a non-empty output/job store. To resume a failed Toil run, use the same paths and `RESTART=1`; it does not delete the job store. Verify Cactus option compatibility with the installed version before starting a large run.

The historical graph's recorded dimensions (about 103.8 million nodes, 143.0 million edges, and 11.84 million candidate VCF records) are useful comparison checks, not strict expected outputs for a reconstruction with a different Cactus version or parameters.

## Software and reproducibility notes

The versions used for downstream analyses are listed in the manuscript's `software_versions.txt` and supplementary materials. The old MC Cactus version is unknown and is intentionally not fabricated here. Performance launchers encode the original benchmark resource design; hardware, scheduler, filesystem, and tool versions affect measured values. Check shell scripts with `bash -n` and Python scripts with `python -m py_compile` before execution, and test on a small subset before submitting full-scale jobs.
