#!/usr/bin/env bash
# Build the missing vg Giraffe indexes for the externally generated cactus GBZ.
# The GBZ and distance index are inputs; only minimizer and zipcode indexes are
# created here with vg 1.73.0.

#SBATCH --job-name=cactus_vg173_idx
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/malong/benchmark/logs/cactus_vg173_idx_%j.out
#SBATCH --error=/data/malong/benchmark/logs/cactus_vg173_idx_%j.err

set -Eeuo pipefail
ulimit -c 0

BASE=/data/malong/benchmark
VG=/data/malong/miniconda3/envs/benchmark/bin/vg
GBZ=${CACTUS_GBZ:-$BASE/graph/cactus_new/rice36_mc.d2.gbz}
DIST=${CACTUS_DIST:-$BASE/graph/cactus_new/rice36_mc.d2.dist}
OUT=${CACTUS_INDEX_DIR:-$BASE/graph/cactus_new/vg173_giraffe_index}
PREFIX=$OUT/rice36_mc.d2
MIN=$PREFIX.shortread.withzip.min
ZIP=$PREFIX.shortread.zipcodes
THREADS=${VG_MINIMIZER_THREADS:-16}

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
nonempty() { [[ -s "$1" ]] || die "missing or empty file: $1"; }
mark() { printf '%s\n' "$(date -Is)" > "$1"; }

mkdir -p "$OUT" "$BASE/logs"
exec > >(tee -a "$OUT/build.log") 2>&1
trap 'rc=$?; echo "[$(date -Is)] FAILED at line $LINENO (exit $rc)" >&2' ERR

[[ -x "$VG" ]] || die "vg 1.73.0 not found: $VG"
nonempty "$GBZ"
nonempty "$DIST"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "VG_MINIMIZER_THREADS must be positive"
log "$($VG version | head -n 1)"
log "GBZ=$GBZ"
log "DIST=$DIST"
log "output prefix=$PREFIX"

if [[ ! -s "$MIN.done" || ! -s "$MIN" || ! -s "$ZIP" ]]; then
    rm -f "$MIN.tmp" "$ZIP.tmp"
    log "building minimizer and zipcode indexes with $THREADS threads"
    "$VG" minimizer --save-memory -t "$THREADS" -d "$DIST" \
        -o "$MIN.tmp" -z "$ZIP.tmp" "$GBZ"
    nonempty "$MIN.tmp"
    nonempty "$ZIP.tmp"
    mv -f "$MIN.tmp" "$MIN"
    mv -f "$ZIP.tmp" "$ZIP"
    mark "$MIN.done"
else
    log "minimizer and zipcode indexes already complete"
fi

nonempty "$MIN"
nonempty "$ZIP"
printf 'vg\t%s\nGBZ\t%s\nDIST\t%s\nMIN\t%s\nZIP\t%s\n' \
    "$($VG version | head -n 1)" "$GBZ" "$DIST" "$MIN" "$ZIP" > "$OUT/INDEX.tsv"
log "completed cactus vg 1.73.0 indexes"
