#!/usr/bin/env bash
#SBATCH --job-name=decon_minigraph
#SBATCH --cpus-per-task=128
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/malong/benchmark/logs/decon_minigraph_%j.out
#SBATCH --error=/data/malong/benchmark/logs/decon_minigraph_%j.err

# Convert the minigraph rGFA (SN/SO/SR tags) to a VG PackedGraph while
# retaining every ranked path, then deconstruct the complete genome in one run.
# No --mem is specified: this cluster assigns memory from cpus-per-task.

set -Eeuo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-128}"
BASE=/data/malong/benchmark
VG=/data/malong/miniconda3/envs/benchmark/bin/vg
BGZIP=/data/malong/miniconda3/envs/benchmark/bin/bgzip
BCFTOOLS=/data/malong/miniconda3/envs/benchmark/bin/bcftools

GFA=$BASE/graph/minigraph/rice36.gfa
OUTDIR=$BASE/graph/minigraph/vg_deconstruct_full
PG=$OUTDIR/rice36.minigraph.allpaths.pg
PATHS_TXT=$OUTDIR/imported_paths.txt
VCF=$OUTDIR/rice36.minigraph.vgdeconstruct.vcf.gz
CONVERT_LOG=$OUTDIR/rice36.minigraph.convert.log
DECON_LOG=$OUTDIR/rice36.minigraph.deconstruct.log
DONE=$OUTDIR/.complete

mkdir -p "$OUTDIR" "$BASE/logs"

[[ -x "$VG" ]] || { echo "ERROR: vg not found: $VG" >&2; exit 2; }
[[ -s "$GFA" ]] || { echo "ERROR: input GFA missing/empty: $GFA" >&2; exit 2; }

# The current rice36 rGFA has SR ranks 0..35.  vg convert imports all
# SN-tagged paths whose rank is <= this value; -r 0 would retain only NIP-T2T.
RGFA_MAX_RANK=35

if [[ ! -s "$PG" ]]; then
    echo "[$(date)] converting rGFA to PackedGraph (rank <= $RGFA_MAX_RANK)" | tee -a "$CONVERT_LOG"
    tmp_pg="$PG.partial"
    rm -f "$tmp_pg"
    "$VG" convert -g -r "$RGFA_MAX_RANK" -p -t "$THREADS" "$GFA" \
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
    echo "The graph may have been converted with an insufficient -r value." >&2
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

echo "[$(date)] minigraph full VCF: $VCF"
echo "[$(date)] graph: $PG"
