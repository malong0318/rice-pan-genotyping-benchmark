#!/usr/bin/env bash
# Build assembly-derived truth sets for six rice samples.
# Submit: sbatch /data/malong/benchmark/code/call_truth.slurm
# One array task processes one sample and creates SNP/Indel and SV outputs.

#SBATCH --job-name=rice_truth
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60
#SBATCH --time=7-00:00:00
#SBATCH --array=0-5
# Keep Slurm's first-level log directory pre-existing at submission time.
# The detailed per-sample log is written under $TRUTH/$sample/ below.
#SBATCH --output=/data/malong/benchmark/rice_truth_%A_%a.out
#SBATCH --error=/data/malong/benchmark/rice_truth_%A_%a.err

set -Eeuo pipefail

BASE=/data/malong/benchmark
REF="$BASE/graph/pggb/rice36_by_sample/NIP-T2T.fa"
TRUTH="$BASE/truth"

# Use absolute paths so the job does not depend on shell initialization.
TRUTH_ENV=/data/malong/miniconda3/envs/benchmark
SV_ENV=/data/malong/miniconda3/envs/benchmark

THREADS=${SLURM_CPUS_PER_TASK:-48}
# dipcall launches two assembly mappings; keep aggregate memory reasonable.
DIPCALL_THREADS=$(( THREADS > 4 ? THREADS / 2 : THREADS ))
MAKE_JOBS=2

SAMPLES=(9311 IR64 zhonghua11 W0169 W1735 W3033)

declare -A ASSEMBLY=(
  [9311]="$BASE/assembly/9311/9311_pansn.fa"
  [IR64]="$BASE/assembly/IR64/IR64_pansn.fa"
  [zhonghua11]="$BASE/assembly/zhonghua11/zhonghua11_pansn.fa"
  [W0169]="$BASE/graph/pggb/rice36_by_sample/W0169.fa"
  [W1735]="$BASE/graph/pggb/rice36_by_sample/W1735.fa"
  [W3033]="$BASE/graph/pggb/rice36_by_sample/W3033.fa"
)

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]:-}
[[ -n "$sample" ]] || { echo "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-unset}" >&2; exit 2; }
asm=${ASSEMBLY[$sample]}

RUN_DIPCALL="$TRUTH_ENV/bin/run-dipcall"
MINIMAP2="$SV_ENV/bin/minimap2"
SAMTOOLS="$SV_ENV/bin/samtools"
BCFTOOLS="$TRUTH_ENV/bin/bcftools"
BEDTOOLS="$TRUTH_ENV/bin/bedtools"
SVIM="$SV_ENV/bin/svim-asm"

for exe in "$RUN_DIPCALL" "$MINIMAP2" "$SAMTOOLS" "$BCFTOOLS" "$BEDTOOLS" "$SVIM"; do
  [[ -x "$exe" ]] || { echo "Missing executable: $exe" >&2; exit 127; }
done
[[ -s "$REF" ]] || { echo "Missing reference: $REF" >&2; exit 2; }
[[ -s "$asm" ]] || { echo "Missing assembly for $sample: $asm" >&2; exit 2; }

mkdir -p "$TRUTH/logs" "$TRUTH/$sample/dipcall" "$TRUTH/$sample/sv" "$TRUTH/$sample/final"

exec > >(tee "$TRUTH/$sample/${sample}.truth.log") 2>&1
echo "[$(date -Is)] sample=$sample assembly=$asm"
echo "[$(date -Is)] reference=$REF"

# Array tasks start together; serialize creation of the shared FASTA index.
INDEX_LOCK="$TRUTH/reference.index.lock"
exec 9>"$INDEX_LOCK"
flock 9
if [[ ! -s "$REF.fai" ]]; then
  "$SAMTOOLS" faidx "$REF"
fi
flock -u 9
exec 9>&-

echo "=== software versions ==="
"$RUN_DIPCALL" 2>&1 | head -3 || true
"$MINIMAP2" --version 2>&1 | head -1 || true
"$SAMTOOLS" --version 2>&1 | head -1 || true
"$BCFTOOLS" --version 2>&1 | head -1 || true
"$SVIM" --help 2>&1 | head -2 || true

###############################################################################
# 1. Small-variant truth. Each sample is a single consensus assembly. dipcall
# requires two haplotypes, so the same assembly is supplied twice. This gives
# an assembly-derived homozygous/pseudodiploid truth, not a heterozygous truth.
###############################################################################
DIPDIR="$TRUTH/$sample/dipcall"
MAKEFILE="$DIPDIR/$sample.mak"
if [[ ! -s "$DIPDIR/$sample.dip.vcf.gz" || ! -s "$DIPDIR/$sample.dip.bed" ]]; then
  "$RUN_DIPCALL" -a -t "$DIPCALL_THREADS" \
    "$sample" "$REF" "$asm" "$asm" > "$MAKEFILE"
  make -C "$DIPDIR" -f "$MAKEFILE" -j "$MAKE_JOBS"
fi

RAW_VCF="$DIPDIR/$sample.dip.vcf.gz"
RAW_BED="$DIPDIR/$sample.dip.bed"
FINAL="$TRUTH/$sample/final"

"$BCFTOOLS" norm -f "$REF" -m -any "$RAW_VCF" -Ou \
  | "$BCFTOOLS" sort -Oz -o "$FINAL/$sample.norm.vcf.gz"
"$BCFTOOLS" index -f "$FINAL/$sample.norm.vcf.gz"
"$BEDTOOLS" sort -i "$RAW_BED" > "$FINAL/$sample.confident.bed"

"$BCFTOOLS" view -v snps "$FINAL/$sample.norm.vcf.gz" \
  -Oz -o "$FINAL/$sample.snp.vcf.gz"
"$BCFTOOLS" index -f "$FINAL/$sample.snp.vcf.gz"
"$BCFTOOLS" view -v indels "$FINAL/$sample.norm.vcf.gz" \
  -Oz -o "$FINAL/$sample.indel.vcf.gz"
"$BCFTOOLS" index -f "$FINAL/$sample.indel.vcf.gz"
"$BCFTOOLS" stats "$FINAL/$sample.norm.vcf.gz" > "$FINAL/$sample.stats.txt"

###############################################################################
# 2. Assembly-based SV truth. SVIM-asm is run in haploid mode because each
# sample has one consensus assembly. The BAM is retained as audit evidence.
###############################################################################
SVDIR="$TRUTH/$sample/sv"
BAM="$SVDIR/$sample.to_ref.bam"
if [[ ! -s "$BAM" ]]; then
  "$MINIMAP2" -ax asm5 --cs -t "$THREADS" "$REF" "$asm" \
    | "$SAMTOOLS" sort -@ "$THREADS" -o "$BAM" -
  "$SAMTOOLS" index "$BAM"
fi

if [[ ! -s "$SVDIR/calls/variants.vcf" && ! -s "$SVDIR/calls/variants.vcf.gz" ]]; then
  "$SVIM" haploid "$SVDIR/calls" "$BAM" "$REF" --sample "$sample"
fi

if [[ -s "$SVDIR/calls/variants.vcf" ]]; then
  "$BCFTOOLS" sort "$SVDIR/calls/variants.vcf" \
    -Oz -o "$FINAL/$sample.sv.vcf.gz"
elif [[ -s "$SVDIR/calls/variants.vcf.gz" ]]; then
  "$BCFTOOLS" sort "$SVDIR/calls/variants.vcf.gz" \
    -Oz -o "$FINAL/$sample.sv.vcf.gz"
else
  echo "SVIM-asm did not create variants.vcf(.gz)" >&2
  exit 1
fi
"$BCFTOOLS" index -f "$FINAL/$sample.sv.vcf.gz"

cat > "$FINAL/$sample.truth.manifest.tsv" <<EOF
sample	$sample
reference	$REF
assembly	$asm
small_variant_method	dipcall_0.3_single_consensus_supplied_as_two_haplotypes
small_variant_vcf	$FINAL/$sample.norm.vcf.gz
small_variant_callable_bed	$FINAL/$sample.confident.bed
sv_method	SVIM-asm_haploid_assembly_alignment
sv_vcf	$FINAL/$sample.sv.vcf.gz
graph_membership	$(case "$sample" in W0169|W1735|W3033) echo in_graph ;; *) echo out_of_graph ;; esac)
EOF

echo "[$(date -Is)] completed"
ls -lh "$FINAL/$sample.norm.vcf.gz" "$FINAL/$sample.snp.vcf.gz" \
  "$FINAL/$sample.indel.vcf.gz" "$FINAL/$sample.sv.vcf.gz" \
  "$FINAL/$sample.confident.bed"
