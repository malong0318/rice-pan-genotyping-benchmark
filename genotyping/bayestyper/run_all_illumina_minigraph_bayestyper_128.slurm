#!/usr/bin/env bash
# BayesTyper genotyping of all Illumina sample/depth pairs on the minigraph
# panel.  9311/5x is intentionally skipped because it was already completed
# by run_9311_5x_minigraph_bayestyper_128.slurm.
#
# Submit with:
#   sbatch /public/home/malong/benchmark/code2/04bayestyper/run_all_illumina_minigraph_bayestyper_128.slurm
#
# The shared KMC/Bloom files and normalized candidate VCF are reused.  Each
# array task has its own cluster/genotype work directory, so failed tasks can
# be resubmitted without repeating completed steps.  Per-step wall time,
# CPU time and peak RSS are written to PERFORMANCE.tsv and to the batch
# summary under results2/MINIGRAPH/04bayestyper.

#SBATCH --job-name=bt2_minigraph_all
# BayesTyper and the shared KMC/Bloom assets are x86-64 binaries.
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=5-00:00:00
#SBATCH --array=0-22%2
#SBATCH --output=/public/home/malong/benchmark/logs/bt2_minigraph_all_%A_%a.out
#SBATCH --error=/public/home/malong/benchmark/logs/bt2_minigraph_all_%A_%a.err

set -Eeuo pipefail
ulimit -c 0

ROOT=/public/home/malong/benchmark
CODE=$ROOT/code/04bayestyper
READROOT=$ROOT/Illumina_depth
BASEOUT=$ROOT/results2/MINIGRAPH/04bayestyper
SHARED=$BASEOUT/_shared
MANIFEST=$SHARED/all_except_9311_5x.manifest.tsv
MANIFEST_LOCK=$SHARED/manifest.lock
SUMMARY=$BASEOUT/batch_performance.tsv
SUMMARY_LOCK=$SHARED/summary.lock

BT=$CODE/software/bayesTyper-debug
BCFTOOLS=/public/home/malong/miniconda3/envs/benchmark/bin/bcftools
SAMTOOLS=/public/home/malong/miniconda3/envs/benchmark/bin/samtools
REF=$ROOT/graph/pggb/rice36_by_sample/NIP-T2T.fa
CAND=$ROOT/graph/minigraph/04bayestyper_work/minigraph.bayestyper.sorted.multiallelic.vcf.gz
KMCROOT=$ROOT/04bayestyper_share/kmc
THREADS=${SLURM_CPUS_PER_TASK:-128}

mkdir -p "$ROOT/logs" "$BASEOUT" "$SHARED"

die() { echo "ERROR: $*" >&2; exit 1; }
nonempty() { [[ -s "$1" ]] || die "missing or empty: $1"; }
need_exe() { [[ -x "$1" ]] || die "not executable: $1"; }

need_exe "$BT"
need_exe "$BCFTOOLS"
need_exe "$SAMTOOLS"
nonempty "$REF"
nonempty "$CAND"
nonempty "$CAND.tbi"
if [[ ! -s "$REF.fai" ]]; then
  "$SAMTOOLS" faidx "$REF"
fi
nonempty "$REF.fai"

# Create one deterministic manifest exactly once.  It contains all paired
# FASTQ files except the already completed 9311/5x run.
(
  /usr/bin/flock -x 9
  if [[ ! -s "$MANIFEST" ]]; then
    tmp=$MANIFEST.tmp.${SLURM_JOB_ID:-$$}
    : > "$tmp"
    while IFS= read -r r1; do
      sample=$(basename "$(dirname "$(dirname "$r1")")")
      depth=$(basename "$(dirname "$r1")")
      [[ "$sample/$depth" == "9311/5x" ]] && continue
      r2=${r1%_R1.fastq.gz}_R2.fastq.gz
      [[ -s "$r2" ]] || die "missing mate for $r1: $r2"
      printf '%s\t%s\t%s\t%s\n' "$sample" "$depth" "$r1" "$r2" >> "$tmp"
    done < <(find "$READROOT" -mindepth 3 -maxdepth 3 -type f -name '*_R1.fastq.gz' | sort)
    rows=$(wc -l < "$tmp")
    [[ "$rows" -eq 23 ]] || die "expected 23 rows after excluding 9311/5x, found $rows"
    mv -f "$tmp" "$MANIFEST"
  fi
) 9>"$MANIFEST_LOCK"

mapfile -t ROWS < "$MANIFEST"
TASK=${SLURM_ARRAY_TASK_ID:-0}
(( TASK >= 0 && TASK < ${#ROWS[@]} )) || die "array task $TASK out of range"
IFS=$'\t' read -r SAMPLE DEPTH READ1 READ2 <<< "${ROWS[$TASK]}"
nonempty "$READ1"
nonempty "$READ2"
[[ "$SAMPLE/$DEPTH" != "9311/5x" ]] || die "9311/5x must not be run by this array"

OUT=$BASEOUT/$SAMPLE/$DEPTH
WORK=$OUT/work
LOG=$OUT/log
mkdir -p "$OUT" "$WORK" "$LOG"
exec > >(tee -a "$OUT/run.log") 2>&1
trap 'rc=$?; echo "[$(date -Is)] ERROR line ${LINENO}, status ${rc}" >&2; exit "$rc"' ERR

PREFIX=${SAMPLE}_${DEPTH}.MINIGRAPH.bayestyper
FINAL=$OUT/${PREFIX}.vcf.gz
SNP=$OUT/${PREFIX}.snp.vcf.gz
INDEL=$OUT/${PREFIX}.indel.vcf.gz
SV=$OUT/${PREFIX}.sv50.vcf.gz
PERF=$OUT/PERFORMANCE.tsv

write_perf_header() {
  if [[ ! -s "$PERF" ]]; then
    printf 'sample\tdepth\ttask\tstep\tstatus\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb\tnote\n' > "$PERF"
  fi
}
write_perf_header

record_perf() {
  local step=$1 status=$2 wall=$3 timefile=$4 note=${5:-}
  local user sys rss
  user=$(awk -F: '/User time \(seconds\)/ {gsub(/^ +/,"",$2); print $2}' "$timefile" 2>/dev/null || true)
  sys=$(awk -F: '/System time \(seconds\)/ {gsub(/^ +/,"",$2); print $2}' "$timefile" 2>/dev/null || true)
  rss=$(awk -F: '/Maximum resident set size \(kbytes\)/ {gsub(/^ +/,"",$2); print $2}' "$timefile" 2>/dev/null || true)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$DEPTH" "$TASK" "$step" "$status" "$wall" "${user:-NA}" "${sys:-NA}" "${rss:-NA}" "$note" >> "$PERF"
}

timed() {
  local name=$1
  shift
  local start end rc timefile
  timefile=$LOG/${name}.time
  start=$(date +%s)
  {
    echo "[$(date -Is)] command:"
    printf ' %q' "$@"
    echo
  } >> "$OUT/RUN_INFO.txt"
  set +e
  /usr/bin/time -v -o "$timefile" "$@" > "$LOG/${name}.log" 2>&1
  rc=$?
  set -e
  end=$(date +%s)
  record_perf "$name" "$rc" "$((end-start))" "$timefile"
  echo "[$(date -Is)] $name: exit=$rc; wall_seconds=$((end-start))" | tee -a "$OUT/RUN_INFO.txt"
  (( rc == 0 )) || die "$name failed; see $LOG/${name}.log"
}

if [[ -s "$OUT/COMPLETED" && -s "$FINAL" && -s "$FINAL.tbi" \
      && -s "$SNP" && -s "$SNP.tbi" && -s "$INDEL" && -s "$INDEL.tbi" \
      && -s "$SV" && -s "$SV.tbi" ]]; then
  echo "[$(date -Is)] $SAMPLE/$DEPTH already complete; skipping" | tee -a "$OUT/run.log"
  printf '%s\t%s\t%s\treuse_complete\t0\t0\tNA\tNA\tNA\tfinal outputs already present\n' \
    "$SAMPLE" "$DEPTH" "$TASK" >> "$PERF"
  exit 0
fi

KMC_PREFIX=$KMCROOT/${SAMPLE}_${DEPTH}.kmc
for ext in kmc_pre kmc_suf bloomMeta bloomData; do
  nonempty "$KMC_PREFIX.$ext"
done

{
  echo "开始时间=$(date -Is)"
  echo "样本=$SAMPLE"
  echo "测序深度=$DEPTH"
  echo "图=minigraph"
  echo "任务=$TASK/${#ROWS[@]}"
  echo "线程=$THREADS"
  echo "调度器内存估算=$((THREADS * 3)) GB"
  echo "参考=$REF"
  echo "候选VCF=$CAND"
  echo "共享KMC前缀=$KMC_PREFIX"
  echo "R1=$READ1"
  echo "R2=$READ2"
} > "$OUT/RUN_INFO.txt"
echo "[$(date -Is)] Starting $SAMPLE/$DEPTH minigraph BayesTyper run" | tee "$OUT/STARTED"

# BayesTyper receives a sample table pointing to the shared KMC/Bloom evidence;
# FASTQ files are deliberately not recounted in this job.
SAMPLES=$WORK/samples.tsv
printf '%s\tF\t%s\n' "$SAMPLE" "$KMC_PREFIX" > "$SAMPLES"

CAND_STAMP=$(stat -c '%s:%Y' "$CAND")
CLUSTER_PREFIX=$WORK/minigraph
CLUSTER_DATA=$WORK/minigraph_cluster_data
CLUSTER_STAMP=$WORK/minigraph.cluster.candidate_stamp
CLUSTER_OK=0
if [[ -s "$CLUSTER_DATA/intercluster_regions.txt.gz" && -s "$CLUSTER_STAMP" \
      && "$(cat "$CLUSTER_STAMP")" == "$CAND_STAMP" ]]; then
  CLUSTER_OK=1
fi
if (( CLUSTER_OK == 0 )); then
  # Isolate stale output from a failed or candidate-mismatched cluster.
  stale=0
  [[ -e "$CLUSTER_DATA" ]] && stale=1
  for old in "$WORK"/minigraph_unit_*; do [[ -e "$old" ]] && stale=1; done
  if (( stale )); then
    failed=$WORK/minigraph_cluster_failed_$(date +%Y%m%dT%H%M%S)_$$
    mkdir -p "$failed"
    [[ ! -e "$CLUSTER_DATA" ]] || mv "$CLUSTER_DATA" "$failed/"
    for old in "$WORK"/minigraph_unit_*; do [[ -e "$old" ]] && mv "$old" "$failed/"; done
  fi
  timed cluster "$BT" cluster -v "$CAND" -s "$SAMPLES" -g "$REF" -o "$CLUSTER_PREFIX" -p"$THREADS"
  printf '%s\n' "$CAND_STAMP" > "$CLUSTER_STAMP"
else
  record_perf cluster reused 0 /dev/null completed_cluster_reused
  echo "[$(date -Is)] reusing completed cluster" | tee -a "$OUT/RUN_INFO.txt"
fi
nonempty "$CLUSTER_DATA/intercluster_regions.txt.gz"

shopt -s nullglob
UNITS=("$WORK"/minigraph_unit_*/variant_clusters.bin)
(( ${#UNITS[@]} > 0 )) || die "no BayesTyper units under $WORK"
UNIT_VCFS=()
for unit in "${UNITS[@]}"; do
  udir=$(dirname "$unit")
  uprefix=$udir/bayestyper
  unit_name=$(basename "$udir")
  if [[ ! -s "$uprefix.vcf.gz" ]]; then
    timed "${unit_name}_genotype" "$BT" genotype -v "$unit" -c "$CLUSTER_DATA" \
      -s "$SAMPLES" -g "$REF" -o "$uprefix" -p"$THREADS" -z
  else
    record_perf "${unit_name}_genotype" reused 0 /dev/null completed_unit_reused
  fi
  nonempty "$uprefix.vcf.gz"
  UNIT_VCFS+=("$uprefix.vcf.gz")
done

# BayesTyper may emit ordinary gzip. Convert and index each unit before
# concat; then allow non-contiguous chromosome blocks and globally sort.
CONCAT_VCFS=()
for src in "${UNIT_VCFS[@]}"; do
  bgz=${src%.vcf.gz}.bgzf.vcf.gz
  if [[ ! -s "$bgz" || ! -s "$bgz.tbi" ]]; then
    unit_name=$(basename "$(dirname "$src")")
    timed "${unit_name}_bgzip" "$BCFTOOLS" view -Oz -o "$bgz" "$src"
    timed "${unit_name}_bgzip_index" "$BCFTOOLS" index -t "$bgz"
  else
    record_perf "$(basename "$(dirname "$src")")_bgzip" reused 0 /dev/null completed_bgzf_reused
  fi
  nonempty "$bgz"
  nonempty "$bgz.tbi"
  CONCAT_VCFS+=("$bgz")
done

UNIT_LIST=$WORK/unit_vcfs.list
printf '%s\n' "${CONCAT_VCFS[@]}" | sort -V > "$UNIT_LIST"
if [[ ! -s "$FINAL" || ! -s "$FINAL.tbi" ]]; then
  RAW=$WORK/minigraph.concat.unsorted.vcf.gz
  SORTED=$WORK/minigraph.concat.sorted.vcf.gz
  rm -f "$FINAL" "$FINAL.tbi" "$RAW" "$SORTED" "$SORTED.tbi"
  timed concat "$BCFTOOLS" concat -a -f "$UNIT_LIST" -Oz -o "$RAW"
  nonempty "$RAW"
  SORTTMP=$WORK/minigraph_final_sort_tmp
  mkdir -p "$SORTTMP"
  timed sort_final "$BCFTOOLS" sort -m 300G -T "$SORTTMP" -Oz -o "$SORTED" "$RAW"
  nonempty "$SORTED"
  mv "$SORTED" "$FINAL"
  rm -f "$RAW"
  timed index_final "$BCFTOOLS" index -t "$FINAL"
else
  record_perf final_merge reused 0 /dev/null final_and_index_reused
fi
nonempty "$FINAL"
nonempty "$FINAL.tbi"

if [[ ! -s "$SNP" || ! -s "$SNP.tbi" ]]; then
  timed split_snp "$BCFTOOLS" view -i 'TYPE="snp"' -Oz -o "$SNP" "$FINAL"
  timed index_snp "$BCFTOOLS" index -t "$SNP"
else record_perf split_snp reused 0 /dev/null output_reused; fi
if [[ ! -s "$INDEL" || ! -s "$INDEL.tbi" ]]; then
  timed split_indel "$BCFTOOLS" view -i 'TYPE="indel" && strlen(REF)<50 && strlen(ALT)<50' -Oz -o "$INDEL" "$FINAL"
  timed index_indel "$BCFTOOLS" index -t "$INDEL"
else record_perf split_indel reused 0 /dev/null output_reused; fi
if [[ ! -s "$SV" || ! -s "$SV.tbi" ]]; then
  timed split_sv "$BCFTOOLS" view -i '(strlen(REF)>=50 || strlen(ALT)>=50 || abs(strlen(REF)-strlen(ALT))>=50)' -Oz -o "$SV" "$FINAL"
  timed index_sv "$BCFTOOLS" index -t "$SV"
else record_perf split_sv reused 0 /dev/null output_reused; fi

FINAL_N=$($BCFTOOLS index -n "$FINAL")
SNP_N=$($BCFTOOLS index -n "$SNP")
INDEL_N=$($BCFTOOLS index -n "$INDEL")
SV_N=$($BCFTOOLS index -n "$SV")
{
  echo "完成时间=$(date -Is)"
  echo "sample=$SAMPLE"
  echo "depth=$DEPTH"
  echo "final=$FINAL"
  echo "snp=$SNP"
  echo "indel=$INDEL"
  echo "sv50=$SV"
  echo "final_records=$FINAL_N"
  echo "snp_records=$SNP_N"
  echo "indel_records=$INDEL_N"
  echo "sv50_records=$SV_N"
  echo "threads=$THREADS"
  echo "memory_estimate_GB=$((THREADS * 3))"
} | tee "$OUT/COMPLETED"

# Append a one-line completion summary atomically. Detailed per-step metrics
# remain in the sample's PERFORMANCE.tsv even if this append is interrupted.
(
  /usr/bin/flock -x 9
  if [[ ! -s "$SUMMARY" ]]; then
    printf 'sample\tdepth\ttask\tfinal_records\tsnp_records\tindel_records\tsv50_records\tcompleted_at\toutput\n' > "$SUMMARY"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SAMPLE" "$DEPTH" "$TASK" \
    "$FINAL_N" "$SNP_N" "$INDEL_N" "$SV_N" "$(date -Is)" "$OUT" >> "$SUMMARY"
) 9>"$SUMMARY_LOCK"

# The final VCFs and performance files are the reusable deliverables. Remove
# large sample-specific intermediates only after COMPLETED was written.
rm -rf "$WORK/minigraph_cluster_data" "$WORK"/minigraph_unit_* "$WORK/minigraph_final_sort_tmp"
rm -f "$WORK"/minigraph.concat.*.vcf.gz "$WORK/unit_vcfs.list" "$WORK/samples.tsv"
echo "[$(date -Is)] completed $SAMPLE/$DEPTH; final_records=$FINAL_N" | tee -a "$OUT/RUN_INFO.txt"
