#!/usr/bin/env bash
#SBATCH --job-name=decon_pggb
#SBATCH --cpus-per-task=128
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/malong/benchmark/logs/decon_pggb_%j.out
#SBATCH --error=/data/malong/benchmark/logs/decon_pggb_%j.err

# Convert the merged 12-chromosome PGGB GFA (P-lines for all haplotypes) to a
# VG PackedGraph, then run vg deconstruct once for the complete genome.
# No --mem is specified: memory is derived from cpus-per-task on this cluster.

set -Eeuo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-128}"
BASE=/data/malong/benchmark
VG=/data/malong/miniconda3/envs/benchmark/bin/vg
BGZIP=/data/malong/miniconda3/envs/benchmark/bin/bgzip
BCFTOOLS=/data/malong/miniconda3/envs/benchmark/bin/bcftools

GFA=$BASE/graph/pggb/merged_12chr/rice36.pggb.12chr.merged.gfa
OUTDIR=$BASE/graph/pggb/merged_12chr/vg_deconstruct_full
PG=$OUTDIR/rice36.pggb.12chr.merged.pg
PATHS_TXT=$OUTDIR/imported_paths.txt
CONVERT_LOG=$OUTDIR/rice36.pggb.convert.log
DECON_LOG=$OUTDIR/rice36.pggb.deconstruct.log
VCF=$OUTDIR/rice36.pggb.vgdeconstruct.vcf.gz
DONE=$OUTDIR/.complete

mkdir -p "$OUTDIR" "$BASE/logs"

[[ -x "$VG" ]] || { echo "ERROR: vg not found: $VG" >&2; exit 2; }
[[ -s "$GFA" ]] || { echo "ERROR: input GFA missing/empty: $GFA" >&2; exit 2; }

if [[ ! -s "$PG" ]]; then
    echo "[$(date)] converting merged PGGB GFA to PackedGraph" | tee -a "$CONVERT_LOG"
    tmp_pg="$PG.partial"
    rm -f "$tmp_pg"
    "$VG" convert -g -p -t "$THREADS" "$GFA" \
        > "$tmp_pg" 2>> "$CONVERT_LOG"
    mv "$tmp_pg" "$PG"
else
    echo "[$(date)] PackedGraph already exists; skip conversion" | tee -a "$CONVERT_LOG"
fi

"$VG" paths -L -x "$PG" | sort -u > "$PATHS_TXT"
path_count="$(wc -l < "$PATHS_TXT")"
echo "[$(date)] imported unique paths: $path_count" | tee -a "$CONVERT_LOG"
if (( path_count < 400 )); then
    echo "ERROR: expected roughly 36 samples x 12 chromosomes (>=400 paths), got $path_count" >&2
    exit 3
fi

REF_PATHS=()
for chr in $(seq 1 12); do
    ref="NIP-T2T#1#Chr${chr}"
    grep -Fxq "$ref" "$PATHS_TXT" || { echo "ERROR: missing reference path $ref" >&2; exit 3; }
    REF_PATHS+=( -p "$ref" )
done

if [[ ! -s "$VCF" || ! -s "$VCF.tbi" || ! -f "$DONE" ]]; then
    echo "[$(date)] running vg deconstruct on all 12 chromosomes" | tee -a "$DECON_LOG"
    tmp_vcf="$VCF.partial"
    rm -f "$tmp_vcf" "$tmp_vcf.tbi"
    set -o pipefail
    "$VG" deconstruct "${REF_PATHS[@]}" -t "$THREADS" "$PG" 2> "$DECON_LOG" \
        | "$BGZIP" -@ "$THREADS" -c > "$tmp_vcf"
    mv "$tmp_vcf" "$VCF"
    "$BCFTOOLS" index -f -t "$VCF"
    "$BCFTOOLS" view -h "$VCF" > "$OUTDIR/header.full.txt"
    sed -n '1,30p' "$OUTDIR/header.full.txt" > "$OUTDIR/header.txt"
    nrec="$($BCFTOOLS view -H "$VCF" | wc -l)"
    nsamp="$($BCFTOOLS query -l "$VCF" | wc -l)"
    echo "[$(date)] deconstructed records: $nrec" | tee -a "$DECON_LOG"
    echo "[$(date)] VCF samples: $nsamp" | tee -a "$DECON_LOG"
    (( nrec > 0 )) || { echo "ERROR: deconstruct produced zero records" >&2; exit 4; }
    if (( nsamp < 30 )); then
        echo "ERROR: expected approximately 36 sample columns, got $nsamp" >&2
        exit 5
    fi
    touch "$DONE"
else
    echo "[$(date)] VCF and index already complete; skip deconstruct" | tee -a "$DECON_LOG"
fi

echo "[$(date)] PGGB full VCF: $VCF"
echo "[$(date)] graph: $PG"
