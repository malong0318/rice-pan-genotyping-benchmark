#!/usr/bin/env bash
# VG 1.73 Giraffe genotyping on the minigraph rGFA graph for all Illumina
# sample/depth pairs under benchmark/Illumina_depth.
#
# One Slurm array task handles one sample/depth pair. Existing completed runs
# are skipped, so the complete array can safely be submitted for continuation.
# This script uses prebuilt Giraffe indexes and deliberately does not run
# vg autoindex.
#
# Submit all 24 sample/depth pairs:
#   sbatch /data/malong/benchmark/code/03vg/run_all_illumina_minigraph_vg.slurm
#
# To submit only unfinished entries after checking the manifest, use e.g.:
#   sbatch --array=0-23%4 /data/malong/benchmark/code/03vg/run_all_illumina_minigraph_vg.slurm

#SBATCH --job-name=vg_mg_all
# About 3 GB memory is allocated per requested CPU on this cluster.
# 52 CPUs therefore request approximately 156 GB.
#SBATCH --partition=I6230R_384G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --time=7-00:00:00
#SBATCH --array=0-23%4
#SBATCH --output=/data/malong/benchmark/results/vg_minigraph_all_%A_%a.out
#SBATCH --error=/data/malong/benchmark/results/vg_minigraph_all_%A_%a.err

set -Eeuo pipefail
ulimit -c 0

BASE=/data/malong/benchmark
VG=/data/malong/miniconda3/envs/benchmark/bin/vg
BCFTOOLS=/data/malong/miniconda3/envs/benchmark/bin/bcftools
BGZIP=/data/malong/miniconda3/envs/benchmark/bin/bgzip

READROOT=$BASE/Illumina_depth
MANIFEST=$BASE/results/minigraph_vg_shared/illumina_depth.manifest.tsv
MANIFEST_LOCK=$BASE/results/minigraph_vg_shared/manifest.lock

# These are the indexes validated by the successful 9311/5x run. Override
# individual paths at submission time if a relocated copy is required.
GRAPH=${MINIGRAPH_GRAPH:-$BASE/graph/minigraph/vg173_giraffe_index/rice36.minigraph.giraffe.gbz}
DIST=${MINIGRAPH_DIST:-$BASE/graph/minigraph/vg173_giraffe_index/rice36.minigraph.dist}
MIN=${MINIGRAPH_MIN:-$BASE/graph/minigraph/vg173_giraffe_index/rice36.minigraph.shortread.withzip.min}
ZIP=${MINIGRAPH_ZIP:-$BASE/graph/minigraph/vg173_giraffe_index/rice36.minigraph.shortread.zipcodes}

MAP_THREADS=${VG_GIRAFFE_THREADS:-32}
PACK_THREADS=${VG_PACK_THREADS:-32}
CALL_THREADS=${VG_CALL_THREADS:-32}

log(){ printf '[%s] %s\n' "$(date -Is)" "$*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ [[ -x "$1" ]] || die "missing executable: $1"; }
nonempty(){ [[ -s "$1" ]] || die "missing or empty file: $1"; }
positive(){ [[ "$1" =~ ^[1-9][0-9]*$ ]] || die "$2 must be a positive integer"; }

need "$VG"
need "$BCFTOOLS"
need "$BGZIP"
for n in MAP_THREADS PACK_THREADS CALL_THREADS; do positive "${!n}" "$n"; done
for f in "$GRAPH" "$DIST" "$MIN" "$ZIP"; do nonempty "$f"; done

mkdir -p "$(dirname "$MANIFEST")"

# Build a stable, tab-separated manifest once. The array index refers to this
# sorted manifest, so it remains reproducible across re-submissions.
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
    mv -f "$tmp" "$MANIFEST"
  fi
  # Remove a temporary manifest left by a diagnostic/manual invocation.
  find "$(dirname "$MANIFEST")" -maxdepth 1 -type f -name "$(basename "$MANIFEST").tmp.*" -delete
) 9>"$MANIFEST_LOCK"

mapfile -t ROWS < "$MANIFEST"
TASK=${SLURM_ARRAY_TASK_ID:-0}
(( TASK >= 0 && TASK < ${#ROWS[@]} )) || die "array task $TASK out of range (manifest rows=${#ROWS[@]})"
IFS=$'\t' read -r SAMPLE DEPTH R1 R2 <<< "${ROWS[$TASK]}"
[[ -n "$SAMPLE" && -n "$DEPTH" && -s "$R1" && -s "$R2" ]] || die "invalid manifest row $TASK"

OUT=$BASE/results/$SAMPLE/$DEPTH/minigraph/vg
mkdir -p "$OUT"
exec > >(tee -a "$OUT/run.log") 2>&1
trap 'rc=$?; echo "[$(date -Is)] FAILED at line $LINENO (exit $rc)" >&2' ERR

SAMPLE_PREFIX=$SAMPLE
GAM=$OUT/${SAMPLE_PREFIX}.minigraph.vg.giraffe.gam
PACK=$OUT/${SAMPLE_PREFIX}.minigraph.vg.pack
VCF=$OUT/${SAMPLE_PREFIX}.minigraph.vg.vcf
VCFGZ=$OUT/${SAMPLE_PREFIX}.minigraph.vg.vcf.gz

log "sample=$SAMPLE depth=$DEPTH task=$TASK"
log "graph=$GRAPH"
log "distance=$DIST"
log "minimizer=$MIN"
log "zipcodes=$ZIP"
log "reads=$R1,$R2"
log "output=$OUT"

# A prior successful run is complete only when the marker, compressed VCF and
# its tabix index all exist. This also prevents stale partial outputs from
# being mistaken for a finished task.
if [[ -s "$OUT/DONE" && -s "$VCFGZ" && -s "$VCFGZ.tbi" ]]; then
  log "already complete; skipping"
  exit 0
fi

if [[ ! -s "$GAM" ]]; then
  tmp=$GAM.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp"
  log "mapping reads with vg giraffe (${MAP_THREADS} threads)"
  "$VG" giraffe -Z "$GRAPH" -d "$DIST" -m "$MIN" -z "$ZIP" \
    -t "$MAP_THREADS" -f "$R1" -f "$R2" -N "$SAMPLE" > "$tmp"
  nonempty "$tmp"
  mv -f "$tmp" "$GAM"
fi
nonempty "$GAM"

if [[ ! -s "$PACK" ]]; then
  tmp=$PACK.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp"
  log "building read-support pack (${PACK_THREADS} threads)"
  "$VG" pack -x "$GRAPH" -g "$GAM" -Q 10 -t "$PACK_THREADS" -o "$tmp"
  nonempty "$tmp"
  mv -f "$tmp" "$PACK"
fi
nonempty "$PACK"

if [[ ! -s "$VCF" ]]; then
  refs=()
  for c in $(seq 1 12); do refs+=(-p "NIP-T2T#1#Chr${c}"); done
  tmp=$VCF.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp"
  log "calling variants with vg call (${CALL_THREADS} threads)"
  "$VG" call "$GRAPH" -k "$PACK" -s "$SAMPLE" -t "$CALL_THREADS" \
    "${refs[@]}" > "$tmp"
  nonempty "$tmp"
  grep -q '^#CHROM' "$tmp" || die "vg call output has no VCF header"
  mv -f "$tmp" "$VCF"
fi
nonempty "$VCF"

if [[ ! -s "$VCFGZ" || ! -s "$VCFGZ.tbi" ]]; then
  tmp=$VCFGZ.tmp.${SLURM_JOB_ID:-$$}
  rm -f "$tmp" "$tmp.tbi"
  log "compressing/indexing VCF"
  "$BGZIP" -@ "$CALL_THREADS" -c "$VCF" > "$tmp"
  nonempty "$tmp"
  "$BCFTOOLS" index -f -t "$tmp"
  nonempty "$tmp.tbi"
  mv -f "$tmp" "$VCFGZ"
  mv -f "$tmp.tbi" "$VCFGZ.tbi"
fi
nonempty "$VCFGZ"
nonempty "$VCFGZ.tbi"

vcf_sample=$("$BCFTOOLS" query -l "$VCFGZ" | paste -sd, -)
[[ "$vcf_sample" == "$SAMPLE" ]] || die "unexpected VCF sample: $vcf_sample (expected $SAMPLE)"
records=$("$BCFTOOLS" view -H "$VCFGZ" | wc -l)

{
  printf 'key\tvalue\n'
  printf 'sample\t%s\n' "$SAMPLE"
  printf 'depth\t%s\n' "$DEPTH"
  printf 'vg\t%s\n' "$("$VG" version | head -n 1)"
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
  printf 'reference_paths\tNIP-T2T#1#Chr1..Chr12\n'
} > "$OUT/INDEX.tsv"

date -Is > "$OUT/DONE"
log "completed VG minigraph genotyping; VCF=$VCFGZ records=$records sample=$vcf_sample"
