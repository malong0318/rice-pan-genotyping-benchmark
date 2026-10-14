#!/usr/bin/env bash
# VG 1.73 Giraffe/pack/call for every sample and depth under Illumina_depth
# against the minigraph-cactus (MC) graph. One array task handles one of the
# 24 sample-depth pairs (6 samples x 4 depths).
#
# Submit:
#   sbatch /data/malong/benchmark/code/03vg/run_all_illumina_mc_vg.slurm
#
# The MC GBZ, distance index, minimizer and zipcode indexes were built before
# this script and are reused. No vg autoindex command is run here.

#SBATCH --job-name=mc_vg_all
# I6230R_384G currently has idle capacity and provides about 384 GB/node.
# This cluster allocates roughly 3 GB per requested CPU.
#SBATCH --partition=I6230R_384G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --time=7-00:00:00
#SBATCH --array=0-23%4
#SBATCH --output=/data/malong/benchmark/results/mc_vg_all_%A_%a.out
#SBATCH --error=/data/malong/benchmark/results/mc_vg_all_%A_%a.err

set -Eeuo pipefail
ulimit -c 0

BASE=/data/malong/benchmark
BIN=/data/malong/miniconda3/envs/benchmark/bin
UTIL=/data/malong/miniconda3/envs/benchmark/bin
READROOT=$BASE/Illumina_depth
GRAPH=${MC_GRAPH:-$BASE/graph/cactus_new/rice36_mc.d2.gbz}
DIST=${MC_DIST:-$BASE/graph/cactus_new/rice36_mc.d2.dist}
INDEXDIR=${MC_INDEX_DIR:-$BASE/graph/cactus_new/vg173_giraffe_index}
MIN=${MC_MIN:-$INDEXDIR/rice36_mc.d2.shortread.withzip.min}
ZIP=${MC_ZIP:-$INDEXDIR/rice36_mc.d2.shortread.zipcodes}
SHARED=$BASE/results/mc_vg_shared
MANIFEST=$SHARED/illumina_depth.manifest.tsv
MANIFEST_LOCK=$SHARED/manifest.lock
THREADS=${SLURM_CPUS_PER_TASK:-128}
MAP_THREADS=${VG_GIRAFFE_THREADS:-32}
PACK_THREADS=${VG_PACK_THREADS:-32}
CALL_THREADS=${VG_CALL_THREADS:-32}
SAMPLE_PL="NIP-T2T#0#Chr"

log(){ printf '[%s] %s\n' "$(date -Is)" "$*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ [[ -x "$1" ]] || die "missing executable: $1"; }
nonempty(){ [[ -s "$1" ]] || die "missing or empty file: $1"; }
positive(){ [[ "$1" =~ ^[1-9][0-9]*$ ]] || die "$2 must be a positive integer"; }

mkdir -p "$SHARED"
need "$BIN/vg"; need "$UTIL/bcftools"; need "$UTIL/bgzip"
[[ -x /usr/bin/flock ]] || die "missing /usr/bin/flock"
for n in MAP_THREADS PACK_THREADS CALL_THREADS; do positive "${!n}" "$n"; done
for f in "$GRAPH" "$DIST" "$MIN" "$ZIP"; do nonempty "$f"; done

# Build a stable mapping from array index to sample/depth/read pair. This is
# shared by all array tasks and regenerated only if absent.
(
  /usr/bin/flock -x 9
  if [[ ! -s "$MANIFEST" ]]; then
    tmp=$MANIFEST.tmp.${SLURM_JOB_ID:-$$}
    : > "$tmp"
    while IFS= read -r r1; do
      sample=$(basename "$(dirname "$(dirname "$r1")")")
      depth=$(basename "$(dirname "$r1")")
      r2=${r1/_R1.fastq.gz/_R2.fastq.gz}
      [[ -s "$r2" ]] || die "missing mate for $r1: $r2"
      printf '%s\t%s\t%s\t%s\n' "$sample" "$depth" "$r1" "$r2" >> "$tmp"
    done < <(find "$READROOT" -mindepth 3 -maxdepth 3 -type f -name '*_R1.fastq.gz' | sort)
    [[ -s "$tmp" ]] || die "no *_R1.fastq.gz files under $READROOT"
    mv "$tmp" "$MANIFEST"
  fi
) 9>"$MANIFEST_LOCK"

mapfile -t ROWS < "$MANIFEST"
(( ${#ROWS[@]} == 24 )) || die "expected 24 sample-depth rows, found ${#ROWS[@]}"
TASK=${SLURM_ARRAY_TASK_ID:-0}
(( TASK < ${#ROWS[@]} )) || die "array task $TASK out of range (rows=${#ROWS[@]})"
IFS=$'\t' read -r SAMPLE DEPTH R1 R2 <<< "${ROWS[$TASK]}"
[[ -n "$SAMPLE" && -n "$DEPTH" && -s "$R1" && -s "$R2" ]] || die "invalid manifest row $TASK"

OUT=$BASE/results/$SAMPLE/$DEPTH/MC/vg
mkdir -p "$OUT"
exec > >(tee -a "$OUT/run.log") 2>&1
trap 'rc=$?; echo "[$(date -Is)] FAILED at line $LINENO (exit $rc)" >&2' ERR

GAM=$OUT/${SAMPLE}.MC.vg.giraffe.gam
PACK=$OUT/${SAMPLE}.MC.vg.pack
VCF=$OUT/${SAMPLE}.MC.vg.vcf
VCFGZ=$OUT/${SAMPLE}.MC.vg.vcf.gz
DONE=$OUT/DONE

log "sample=$SAMPLE depth=$DEPTH task=$TASK"
log "vg=$($BIN/vg version | head -n 1)"
log "graph=$GRAPH"
log "distance=$DIST"
log "minimizer=$MIN"
log "zipcodes=$ZIP"
log "reads=$R1,$R2"
log "output=$OUT"

if [[ ! -s "$GAM" ]]; then
  tmp=$GAM.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp"
  log "mapping reads with vg giraffe (${MAP_THREADS} threads)"
  "$BIN/vg" giraffe -Z "$GRAPH" -d "$DIST" -m "$MIN" -z "$ZIP" \
    -t "$MAP_THREADS" -f "$R1" -f "$R2" -N "$SAMPLE" > "$tmp"
  nonempty "$tmp"
  mv -f "$tmp" "$GAM"
fi
nonempty "$GAM"

if [[ ! -s "$PACK" ]]; then
  tmp=$PACK.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp"
  log "building vg pack (${PACK_THREADS} threads)"
  "$BIN/vg" pack -x "$GRAPH" -g "$GAM" -Q 10 -t "$PACK_THREADS" -o "$tmp"
  nonempty "$tmp"
  mv -f "$tmp" "$PACK"
fi
nonempty "$PACK"

# MC GBZ stores the NIP-T2T reference paths as NIP-T2T#0#Chr1..Chr12.
if [[ ! -s "$VCF" ]]; then
  refs=()
  for c in $(seq 1 12); do refs+=(-p "${SAMPLE_PL}${c}"); done
  tmp=$VCF.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp"
  log "calling graph variants with vg call (${CALL_THREADS} threads)"
  "$BIN/vg" call "$GRAPH" -k "$PACK" -s "$SAMPLE" -t "$CALL_THREADS" \
    "${refs[@]}" > "$tmp"
  nonempty "$tmp"
  grep -q '^#CHROM' "$tmp" || die "vg call output has no VCF header"
  mv -f "$tmp" "$VCF"
fi
nonempty "$VCF"

if [[ ! -s "$VCFGZ" || ! -s "$VCFGZ.tbi" ]]; then
  tmp=$VCFGZ.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp" "$tmp.tbi"
  log "compressing and indexing VCF"
  "$UTIL/bgzip" -@ "$CALL_THREADS" -c "$VCF" > "$tmp"
  nonempty "$tmp"
  "$UTIL/bcftools" index -f -t "$tmp"
  nonempty "$tmp.tbi"
  mv -f "$tmp" "$VCFGZ"
  mv -f "$tmp.tbi" "$VCFGZ.tbi"
fi
nonempty "$VCFGZ"; nonempty "$VCFGZ.tbi"

vcf_sample=$($UTIL/bcftools query -l "$VCFGZ" | paste -sd, -)
[[ "$vcf_sample" == "$SAMPLE" ]] || die "unexpected VCF sample: $vcf_sample (expected $SAMPLE)"
records=$($UTIL/bcftools view -H "$VCFGZ" | wc -l)

{
  printf 'key\tvalue\n'
  printf 'sample\t%s\n' "$SAMPLE"
  printf 'depth\t%s\n' "$DEPTH"
  printf 'task\t%s\n' "$TASK"
  printf 'vg\t%s\n' "$($BIN/vg version | head -n 1)"
  printf 'graph\t%s\n' "$GRAPH"
  printf 'distance\t%s\n' "$DIST"
  printf 'minimizer\t%s\n' "$MIN"
  printf 'zipcodes\t%s\n' "$ZIP"
  printf 'reads_R1\t%s\n' "$R1"
  printf 'reads_R2\t%s\n' "$R2"
  printf 'gam\t%s\n' "$GAM"
  printf 'pack\t%s\n' "$PACK"
  printf 'vcf\t%s\n' "$VCF"
  printf 'vcf_gz\t%s\n' "$VCFGZ"
  printf 'sample_in_vcf\t%s\n' "$vcf_sample"
  printf 'records\t%s\n' "$records"
  printf 'reference_paths\tNIP-T2T#0#Chr1..Chr12\n'
} > "$OUT/INDEX.tsv"

date -Is > "$DONE"
log "completed VG MC genotyping; VCF=$VCFGZ records=$records sample=$vcf_sample"
