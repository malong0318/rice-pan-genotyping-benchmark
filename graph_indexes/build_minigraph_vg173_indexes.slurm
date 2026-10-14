#!/usr/bin/env bash
# Build the complete-genome minigraph Giraffe indexes with vg 1.73.0.
# vg autoindex performs the required node chopping and keeps rGFA paths.

#SBATCH --job-name=minigraph_vg173_idx
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/malong/benchmark/logs/minigraph_vg173_idx_%j.out
#SBATCH --error=/data/malong/benchmark/logs/minigraph_vg173_idx_%j.err

set -Eeuo pipefail
ulimit -c 0

BASE=/data/malong/benchmark
VG=/data/malong/miniconda3/envs/benchmark/bin/vg
GFA=${MINIGRAPH_GFA:-$BASE/graph/minigraph/rice36.gfa}
OUT=${MINIGRAPH_INDEX_DIR:-$BASE/graph/minigraph/vg173_giraffe_index}
PREFIX=$OUT/rice36.minigraph
GBZ=$PREFIX.giraffe.gbz
DIST=$PREFIX.dist
MIN=$PREFIX.shortread.withzip.min
ZIP=$PREFIX.shortread.zipcodes
TMP=$OUT/tmp
COMPLETE=$PREFIX.complete
THREADS=${VG_INDEX_THREADS:-128}
TARGET_MEM=${VG_TARGET_MEM:-450G}

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
nonempty() { [[ -s "$1" ]] || die "missing or empty file: $1"; }
mark() { printf '%s\n' "$(date -Is)" > "$1"; }

mkdir -p "$OUT" "$BASE/logs"
exec > >(tee -a "$OUT/build.log") 2>&1
trap 'rc=$?; echo "[$(date -Is)] FAILED at line $LINENO (exit $rc)" >&2' ERR

[[ -x "$VG" ]] || die "vg 1.73.0 not found: $VG"
nonempty "$GFA"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "VG_INDEX_THREADS must be positive"
log "$($VG version | head -n 1)"
log "input rGFA=$GFA"

if [[ -s "$COMPLETE" && -s "$GBZ" && -s "$DIST" && -s "$MIN" && -s "$ZIP" ]]; then
    log "final Giraffe indexes already complete; nothing to rebuild"
    exit 0
fi

# Remove partial final files so autoindex cannot mistake them for inputs.
rm -f "$GBZ" "$DIST" "$MIN" "$ZIP" "$COMPLETE"
rm -rf "$TMP"
mkdir -p "$TMP"
log "building sr-giraffe indexes with vg autoindex ($THREADS threads, target $TARGET_MEM)"
"$VG" autoindex -w sr-giraffe -g "$GFA" -p "$PREFIX" -t "$THREADS" \
    -M "$TARGET_MEM" -T "$TMP" -n

for f in "$GBZ" "$DIST" "$MIN" "$ZIP"; do nonempty "$f"; done

# Verify that all minigraph haplotypes survived conversion, including all
# twelve NIP-T2T reference paths.
"$VG" paths -L -x "$GBZ" | sort -u > "$OUT/path_names.txt.tmp"
nonempty "$OUT/path_names.txt.tmp"
count=$(wc -l < "$OUT/path_names.txt.tmp")
(( count >= 400 )) || die "only $count paths in GBZ; expected >=400"
for chr in $(seq 1 12); do
    grep -Fxq "NIP-T2T#1#Chr${chr}" "$OUT/path_names.txt.tmp" || die "missing reference path Chr${chr}"
done
mv -f "$OUT/path_names.txt.tmp" "$OUT/path_names.txt"

# autoindex output is sufficient for Giraffe; discard construction scratch.
rm -rf "$TMP"
mark "$COMPLETE"
printf 'vg\t%s\nGFA\t%s\nGBZ\t%s\nDIST\t%s\nMIN\t%s\nZIP\t%s\nPATHS\t%s\n' \
    "$($VG version | head -n 1)" "$GFA" "$GBZ" "$DIST" "$MIN" "$ZIP" "$count" > "$OUT/INDEX.tsv"
log "completed minigraph vg 1.73.0 indexes; retained $count paths"
