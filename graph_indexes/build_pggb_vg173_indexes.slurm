#!/usr/bin/env bash
# Build complete-genome PGGB Giraffe indexes with vg 1.73.0.
#
# vg autoindex cannot serialize the default whole-genome distance index for
# this 36-haplotype PGGB graph ("File too large").  This script therefore
# builds GBZ explicitly and uses a bounded, top-level-only distance index.
# If even that distance index is not possible, it falls back to a valid
# no-distance minimizer index; Giraffe can run in that mode.

#SBATCH --job-name=pggb_vg173_idx
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/malong/benchmark/logs/pggb_vg173_idx_%j.out
#SBATCH --error=/data/malong/benchmark/logs/pggb_vg173_idx_%j.err

set -Eeuo pipefail
ulimit -c 0

BASE=/data/malong/benchmark
VG=/data/malong/miniconda3/envs/benchmark/bin/vg
GFA=${PGGB_GFA:-$BASE/graph/pggb/rice36_graph/rice36.complete.fa.a13a25f.37c1a4f.8fffabc.smooth.final.gfa}
OUT=${PGGB_INDEX_DIR:-$BASE/graph/pggb/rice36_graph/vg173_giraffe_index}
PREFIX=$OUT/rice36.pggb
CLEAN=$PREFIX.clean.gfa
PG=$PREFIX.pg
XG=$PREFIX.xg
GBWT=$PREFIX.gbwt
GBZ=$PREFIX.giraffe.gbz
DIST=$PREFIX.dist
MIN=$PREFIX.shortread.withzip.min
ZIP=$PREFIX.shortread.zipcodes
TMP=$OUT/tmp
COMPLETE=$PREFIX.complete
MODE=$PREFIX.minimizer.mode
THREADS=${VG_INDEX_THREADS:-128}
GBWT_JOBS=${VG_GBWT_JOBS:-32}
MIN_THREADS=${VG_MINIMIZER_THREADS:-16}
SNARL_LIMIT=${PGGB_SNARL_LIMIT:-1000}

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
nonempty() { [[ -s "$1" ]] || die "missing or empty file: $1"; }
mark() { printf '%s\n' "$(date -Is)" > "$1"; }

mkdir -p "$OUT" "$BASE/logs"
exec > >(tee -a "$OUT/build.log") 2>&1

cleanup_failed() {
    rc=$?
    if (( rc != 0 )); then
        log "failed (exit $rc); removing partial PGGB index files"
        rm -f "$CLEAN" "$CLEAN.tmp" "$PG" "$PG.tmp" "$XG" "$XG.tmp" \
            "$GBWT" "$GBWT.tmp" "$GBZ" "$GBZ.tmp" "$DIST" "$DIST.tmp" \
            "$MIN" "$MIN.tmp" "$ZIP" "$ZIP.tmp" "$COMPLETE" "$MODE"
        rm -rf "$TMP"
    fi
    exit "$rc"
}
trap cleanup_failed EXIT
trap 'rc=$?; echo "[$(date -Is)] FAILED at line $LINENO (exit $rc)" >&2; exit $rc' ERR

[[ -x "$VG" ]] || die "vg 1.73.0 not found: $VG"
nonempty "$GFA"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "VG_INDEX_THREADS must be positive"
[[ "$GBWT_JOBS" =~ ^[1-9][0-9]*$ ]] || die "VG_GBWT_JOBS must be positive"
[[ "$MIN_THREADS" =~ ^[1-9][0-9]*$ ]] || die "VG_MINIMIZER_THREADS must be positive"
[[ "$SNARL_LIMIT" =~ ^[0-9]+$ ]] || die "PGGB_SNARL_LIMIT must be non-negative"
log "$($VG version | head -n 1)"
log "input PGGB GFA=$GFA"

if [[ -s "$COMPLETE" && -s "$GBZ" && -s "$MIN" && -s "$MODE" ]]; then
    log "final PGGB indexes already complete; nothing to rebuild"
    exit 0
fi

FIRST_PATH=$(grep -m 1 $'^P\t' "$GFA" | cut -f2)
[[ -n "$FIRST_PATH" ]] || die "PGGB GFA has no P-line paths"
log "first PGGB P-line path=$FIRST_PATH"

rm -f "$CLEAN" "$CLEAN.tmp" "$PG" "$PG.tmp" "$XG" "$XG.tmp" \
    "$GBWT" "$GBWT.tmp" "$GBZ" "$GBZ.tmp" "$DIST" "$DIST.tmp" \
    "$MIN" "$MIN.tmp" "$ZIP" "$ZIP.tmp" "$COMPLETE" "$MODE"
rm -rf "$TMP"
mkdir -p "$TMP"

# Keep all P/L/H records and all node sequences.  Remove only S-line metadata
# tags because PGGB may encode the same path both as SN/SO/SR and as P lines.
log "[1/7] making a clean GFA while retaining all P paths"
awk -F '\t' -v OFS='\t' '$1=="S" {print $1,$2,$3; next} {print}' \
    "$GFA" > "$CLEAN.tmp"
nonempty "$CLEAN.tmp"
mv -f "$CLEAN.tmp" "$CLEAN"

log "[2/7] converting clean PGGB GFA to PackedGraph"
"$VG" convert -g -p -t "$THREADS" "$CLEAN" > "$PG.tmp"
nonempty "$PG.tmp"
mv -f "$PG.tmp" "$PG"

log "[3/7] building XG index"
"$VG" index -b "$TMP" -t "$THREADS" -x "$XG.tmp" "$PG"
nonempty "$XG.tmp"
mv -f "$XG.tmp" "$XG"

log "[4/7] indexing every PGGB P path into GBWT"
"$VG" gbwt -G "$CLEAN" -d "$TMP" -o "$GBWT.tmp" \
    --max-node 0 --num-jobs "$GBWT_JOBS" --progress
nonempty "$GBWT.tmp"
mv -f "$GBWT.tmp" "$GBWT"

log "[5/7] combining XG and GBWT into GBZ"
"$VG" gbwt -x "$XG" -g "$GBZ.tmp" --gbz-format "$GBWT"
nonempty "$GBZ.tmp"
mv -f "$GBZ.tmp" "$GBZ"

DIST_MODE=none
log "[6/7] building bounded snarl distance index (snarl-limit=$SNARL_LIMIT)"
if "$VG" index -b "$TMP" -t "$THREADS" -j "$DIST.tmp" \
       --snarl-limit "$SNARL_LIMIT" --no-nested-distance "$GBZ"; then
    if [[ -s "$DIST.tmp" ]]; then
        mv -f "$DIST.tmp" "$DIST"
        DIST_MODE=bounded
        log "bounded distance index completed"
    else
        rm -f "$DIST.tmp"
        log "distance command returned success but produced no file"
    fi
else
    rm -f "$DIST.tmp"
    log "bounded distance index failed; continuing with no-distance minimizer"
fi

MIN_MODE=nodist
if [[ "$DIST_MODE" == bounded ]]; then
    log "[7/7] building distance-aware minimizer and zipcodes"
    if "$VG" minimizer --save-memory -t "$MIN_THREADS" -d "$DIST" \
        -o "$MIN.tmp" -z "$ZIP.tmp" "$GBZ" &&
       [[ -s "$MIN.tmp" && -s "$ZIP.tmp" ]]; then
        mv -f "$MIN.tmp" "$MIN"
        mv -f "$ZIP.tmp" "$ZIP"
        MIN_MODE=bounded
        log "distance-aware minimizer completed"
    else
        rm -f "$MIN.tmp" "$ZIP.tmp"
        log "distance-aware minimizer failed; using no-distance minimizer"
    fi
fi

if [[ "$MIN_MODE" == nodist ]]; then
    log "[7/7] building no-distance minimizer"
    "$VG" minimizer --save-memory --no-dist -t "$MIN_THREADS" \
        -o "$MIN.tmp" "$GBZ"
    nonempty "$MIN.tmp"
    mv -f "$MIN.tmp" "$MIN"
    rm -f "$ZIP" "$ZIP.tmp"
fi

nonempty "$GBZ"
nonempty "$MIN"
"$VG" paths -L -x "$GBZ" | sort -u > "$OUT/path_names.txt.tmp"
nonempty "$OUT/path_names.txt.tmp"
count=$(wc -l < "$OUT/path_names.txt.tmp")
(( count >= 400 )) || die "only $count paths in GBZ; expected >=400"
for chr in $(seq 1 12); do
    grep -Fxq "NIP-T2T#1#Chr${chr}" "$OUT/path_names.txt.tmp" || die "missing reference path Chr${chr}"
done
mv -f "$OUT/path_names.txt.tmp" "$OUT/path_names.txt"

# Retain only Giraffe-consumable final files and provenance.  The original
# PGGB GFA is never modified or removed.
rm -f "$CLEAN" "$PG" "$XG" "$GBWT" "$CLEAN.done" "$PG.done" \
    "$XG.done" "$GBWT.done"
rm -rf "$TMP"
printf '%s\n' "$MIN_MODE" > "$MODE"
mark "$COMPLETE"
printf 'vg\t%s\nGFA\t%s\nGBZ\t%s\nDIST\t%s\nMIN\t%s\nZIP\t%s\nMIN_MODE\t%s\nDIST_MODE\t%s\nP_PATHS\t%s\n' \
    "$($VG version | head -n 1)" "$GFA" "$GBZ" "$DIST" "$MIN" "$ZIP" \
    "$MIN_MODE" "$DIST_MODE" "$count" > "$OUT/INDEX.tsv"
log "completed PGGB vg 1.73.0 indexes; retained $count paths; minimizer mode=$MIN_MODE"
