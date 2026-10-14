#!/usr/bin/env bash
#SBATCH --job-name=rice36_pggb_full
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
# AMD_9654 nodes on this cluster expose 128 CPUs per node; requesting 256
# makes Slurm reject the job before the script starts.
#SBATCH --cpus-per-task=128
#SBATCH --requeue
#SBATCH --output=/data/malong/benchmark/logs/rice36_pggb_full_%j.out
#SBATCH --error=/data/malong/benchmark/logs/rice36_pggb_full_%j.err

# Run PGGB once on all 36 rice assemblies (all 12 chromosomes remain in the
# same multi-FASTA).  The input and every expensive PGGB stage are restartable.
# Submit with:
#   sbatch /data/malong/benchmark/code/creat_graph/run_pggb_rice36_full_genome.slurm
#
# Useful overrides:
#   OUTDIR=/path/with/more/space  THREADS=128  POA_THREADS=64  sbatch ...
#   RESUME=0 sbatch ...                 # recompute PGGB stages from scratch
#   FORCE_REBUILD_INPUT=1 sbatch ...    # rebuild the combined FASTA

set -Eeuo pipefail
ulimit -n 65535 2>/dev/null || true

BASE=/data/malong/benchmark
INPUT_DIR="$BASE/graph/pggb/rice36_by_sample"
LIST_ROOT="$BASE/graph/pggb"
SAMPLE_LIST="$INPUT_DIR/sample.list"
FASTA_LIST="$INPUT_DIR/sample_fasta.list"

OUTDIR="${OUTDIR:-$BASE/graph/pggb/rice36_full_genome}"
PGGB_TMP="${PGGB_TMP:-$OUTDIR/.pggb_tmp}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-128}}"
POA_THREADS="${POA_THREADS:-64}"
SEGMENT_LENGTH="${SEGMENT_LENGTH:-10000}"
MAP_PCT_ID="${MAP_PCT_ID:-95}"
MIN_MATCH_LEN="${MIN_MATCH_LEN:-19}"
POA_LENGTH_TARGET="${POA_LENGTH_TARGET:-700,1100}"
RESUME="${RESUME:-1}"
FORCE_REBUILD_INPUT="${FORCE_REBUILD_INPUT:-0}"

ENV=/data/malong/miniconda3/envs/benchmark
BIN="$ENV/bin"
PGGB="$BIN/pggb"
SAMTOOLS="$BIN/samtools"
BCFTOOLS="$BIN/bcftools"
BGZIP="$BIN/bgzip"

mkdir -p "$OUTDIR" "$PGGB_TMP" "$BASE/logs"
exec > >(tee -a "$OUTDIR/run.log") 2>&1

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }
on_exit() {
  local rc=$?
  if (( rc != 0 )); then
    log "FAILED (exit $rc). Existing PGGB intermediates were retained for resume."
  fi
}
trap on_exit EXIT
trap 'log "received termination signal"; exit 143' TERM INT

need_exec() { [[ -x "$1" ]] || die "missing executable: $1"; }
for exe in "$PGGB" "$SAMTOOLS" "$BCFTOOLS" "$BGZIP" \
           "$BIN/wfmash" "$BIN/seqwish" "$BIN/smoothxg" "$BIN/odgi" \
           "$BIN/gfaffix" "$BIN/vg" "$BIN/vcfbub" "$BIN/vcfwave"; do
  need_exec "$exe"
done
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer: $THREADS"
[[ "$POA_THREADS" =~ ^[1-9][0-9]*$ ]] || die "POA_THREADS must be a positive integer: $POA_THREADS"
export PATH="$BIN:$PATH"

[[ -s "$SAMPLE_LIST" ]] || die "missing sample list: $SAMPLE_LIST"
[[ -s "$FASTA_LIST" ]] || die "missing FASTA list: $FASTA_LIST"

INPUT_FA="$OUTDIR/rice36.complete.fa"
INPUT_FAI="$INPUT_FA.fai"
INPUT_DONE="$OUTDIR/.combined_input.complete"
EXPECTED="$OUTDIR/expected_paths.txt"
SOURCE_MANIFEST="$OUTDIR/source_manifest.tsv"

mapfile -t SAMPLES < <(sed 's/\r$//' "$SAMPLE_LIST" | sed '/^[[:space:]]*$/d')
mapfile -t REL_FASTAS < <(sed 's/\r$//' "$FASTA_LIST" | sed '/^[[:space:]]*$/d')
(( ${#SAMPLES[@]} == 36 )) || die "sample.list must contain 36 non-empty sample names (found ${#SAMPLES[@]})"
(( ${#REL_FASTAS[@]} == 36 )) || die "sample_fasta.list must contain 36 non-empty paths (found ${#REL_FASTAS[@]})"

printf '%s\n' "${SAMPLES[@]}" | sort -u > "$OUTDIR/.samples.sorted"
(( $(wc -l < "$OUTDIR/.samples.sorted") == 36 )) || die "sample.list contains duplicate sample names"

build_expected_paths() {
  : > "$EXPECTED.tmp"
  local sample chr
  for sample in "${SAMPLES[@]}"; do
    for chr in $(seq 1 12); do
      printf '%s\n' "$sample#1#Chr$chr" >> "$EXPECTED.tmp"
    done
  done
  sort -u "$EXPECTED.tmp" > "$EXPECTED"
  rm -f "$EXPECTED.tmp"
  (( $(wc -l < "$EXPECTED") == 432 )) || die "internal expected-path construction failed"
}

write_source_manifest() {
  : > "$SOURCE_MANIFEST.tmp"
  local i rel file sample_from_file
  for i in "${!REL_FASTAS[@]}"; do
    rel="${REL_FASTAS[$i]}"
    if [[ "$rel" = /* ]]; then file="$rel"; else file="$LIST_ROOT/$rel"; fi
    [[ -s "$file" ]] || die "missing or empty input FASTA: $file"
    sample_from_file=$(basename "$file")
    sample_from_file=${sample_from_file%.fa}
    [[ "$sample_from_file" == "${SAMPLES[$i]}" ]] || \
      die "sample/list order mismatch: ${SAMPLES[$i]} vs $sample_from_file ($file)"
    stat -c '%n\t%s\t%Y' "$file" >> "$SOURCE_MANIFEST.tmp"
  done
  mv -f "$SOURCE_MANIFEST.tmp" "$SOURCE_MANIFEST"
}

validate_fasta() {
  local file=$1 count duplicates actual sorted_actual
  [[ -s "$file" ]] || die "empty FASTA: $file"
  count=$(awk '/^>/{n++} END{print n+0}' "$file")
  (( count == 432 )) || die "$file has $count sequences; expected 432 (36 samples x 12 chromosomes)"
  duplicates=$(awk '/^>/{h=substr($0,2); n[h]++} END{for(h in n) if(n[h]>1) print h}' "$file" | sort -u)
  [[ -z "$duplicates" ]] || die "duplicate FASTA sequence names found: $duplicates"
  actual="$OUTDIR/.actual_paths.sorted"
  awk '/^>/{print substr($0,2)}' "$file" | sort -u > "$actual"
  sorted_actual=$(diff -u "$EXPECTED" "$actual" || true)
  [[ -z "$sorted_actual" ]] || { printf '%s\n' "$sorted_actual" >&2; die "FASTA names are not exactly the expected 432 PanSN paths"; }
  rm -f "$actual"
}

build_combined_input() {
  build_expected_paths
  write_source_manifest
  if [[ "$FORCE_REBUILD_INPUT" != 1 && -s "$INPUT_FA" && -s "$INPUT_FAI" && -f "$INPUT_DONE" && \
        -s "$SOURCE_MANIFEST" && -f "$OUTDIR/source_manifest.completed.tsv" ]] && \
        cmp -s "$SOURCE_MANIFEST" "$OUTDIR/source_manifest.completed.tsv"; then
    log "combined input and source manifest already complete; validating it"
    validate_fasta "$INPUT_FA"
    return
  fi

  log "building one complete-genome FASTA from all 36 listed assemblies"
  local tmp="$INPUT_FA.partial.$$" rel file
  rm -f "$tmp" "$tmp.fai"
  : > "$tmp"
  for rel in "${REL_FASTAS[@]}"; do
    if [[ "$rel" = /* ]]; then file="$rel"; else file="$LIST_ROOT/$rel"; fi
    cat "$file" >> "$tmp"
    [[ "$(tail -c 1 "$file" | wc -l)" -eq 1 ]] || printf '\n' >> "$tmp"
  done
  validate_fasta "$tmp"
  log "creating FASTA index: $tmp.fai"
  "$SAMTOOLS" faidx "$tmp"
  [[ -s "$tmp.fai" ]] || die "samtools faidx produced no index"
  mv -f "$tmp" "$INPUT_FA"
  mv -f "$tmp.fai" "$INPUT_FAI"
  cp -f "$SOURCE_MANIFEST" "$OUTDIR/source_manifest.completed.tsv"
  touch "$INPUT_DONE"
  log "combined FASTA ready: $(stat -c '%s bytes' "$INPUT_FA"); 432 paths indexed"
}

build_combined_input

if [[ -f "$OUTDIR/.pggb.complete" && "$RESUME" != 0 ]]; then
  log "PGGB completion marker exists; validating final outputs before exiting"
else
  log "running PGGB once over the complete genome (not chromosome-by-chromosome)"
  pggb_args=(
    -i "$INPUT_FA" -o "$OUTDIR" -D "$PGGB_TMP"
    -t "$THREADS" -T "$POA_THREADS"
    -s "$SEGMENT_LENGTH" -p "$MAP_PCT_ID" -g 30 -Y '#'
    -k "$MIN_MATCH_LEN" -G "$POA_LENGTH_TARGET" -n 36
    -V NIP-T2T -v -Z
  )
  if [[ "$RESUME" != 0 ]]; then pggb_args+=( -r ); fi
  log "command: $PGGB ${pggb_args[*]}"
  "$PGGB" "${pggb_args[@]}"
fi

mapfile -t GFA_GZ_FILES < <(find "$OUTDIR" -maxdepth 1 -type f -name '*.smooth.final.gfa.gz' -print | sort)
if (( ${#GFA_GZ_FILES[@]} == 1 )); then
  GFA_FILES=("${GFA_GZ_FILES[0]}")
else
  mapfile -t GFA_FILES < <(find "$OUTDIR" -maxdepth 1 -type f -name '*.smooth.final.gfa' -print | sort)
fi
(( ${#GFA_FILES[@]} == 1 )) || die "expected exactly one final PGGB GFA(.gz), found ${#GFA_FILES[@]}"
[[ -s "${GFA_FILES[0]}" ]] || die "final PGGB graph is empty: ${GFA_FILES[0]}"

mapfile -t VCF_GZ_FILES < <(find "$OUTDIR" -maxdepth 1 -type f -name '*.final.NIP-T2T.vcf.gz' -print | sort)
if (( ${#VCF_GZ_FILES[@]} == 1 )); then
  VCF_FILES=("${VCF_GZ_FILES[0]}")
else
  mapfile -t VCF_FILES < <(find "$OUTDIR" -maxdepth 1 -type f -name '*.final.NIP-T2T.vcf' -print | sort)
fi
(( ${#VCF_FILES[@]} == 1 )) || die "expected exactly one NIP-T2T VCF(.gz), found ${#VCF_FILES[@]}"
VCF="${VCF_FILES[0]}"
[[ -s "$VCF" ]] || die "VCF is empty: $VCF"

# PGGB -Z normally creates a bgzip VCF.  If a plain VCF is returned (for
# example after a version change), create a compressed copy and index that;
# failure to make an optional index is reported but does not discard the VCF.
if [[ "$VCF" != *.gz ]]; then
  VCF_GZ="$VCF.gz"
  if [[ ! -s "$VCF_GZ" ]]; then
    log "compressing plain VCF with bgzip"
    "$BGZIP" -@ "$THREADS" -c "$VCF" > "$VCF_GZ.tmp"
    mv -f "$VCF_GZ.tmp" "$VCF_GZ"
  fi
  VCF="$VCF_GZ"
fi

INDEX_STATUS=not_attempted
if [[ ! -s "$VCF.tbi" && ! -s "$VCF.csi" ]]; then
  log "building optional VCF index"
  if "$BCFTOOLS" index -f -t "$VCF"; then
    INDEX_STATUS="$VCF.tbi"
  else
    log "WARNING: VCF index creation failed; retaining the VCF as requested"
  fi
else
  [[ -s "$VCF.tbi" ]] && INDEX_STATUS="$VCF.tbi"
  [[ -s "$VCF.csi" ]] && INDEX_STATUS="$VCF.csi"
fi

"$BCFTOOLS" view -h "$VCF" | grep -q '^#CHROM' || die "VCF has no #CHROM header: $VCF"
record_count=$("$BCFTOOLS" view -H "$VCF" | wc -l)
(( record_count > 0 )) || die "VCF contains zero variant records: $VCF"

{
  printf 'item\tpath\n'
  printf 'input_fasta\t%s\n' "$INPUT_FA"
  printf 'input_fasta_index\t%s\n' "$INPUT_FAI"
  printf 'final_graph\t%s\n' "${GFA_FILES[0]}"
  printf 'vcf_reference_prefix\tNIP-T2T\n'
  printf 'vcf\t%s\n' "$VCF"
  printf 'vcf_index\t%s\n' "$INDEX_STATUS"
  printf 'vcf_records\t%s\n' "$record_count"
} > "$OUTDIR/RESULTS.tsv"
touch "$OUTDIR/.pggb.complete"
log "completed full-genome PGGB"
log "graph: ${GFA_FILES[0]}"
log "VCF: $VCF"
log "VCF index: $INDEX_STATUS"
