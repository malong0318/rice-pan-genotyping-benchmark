#!/usr/bin/env bash
# Run SyRI as an independent assembly-based SV caller for the six rice samples.
# Submit all samples:
#   sbatch /data/malong/benchmark/code/call_syri_sv.slurm
# Submit one sample (0=9311, 1=IR64, 2=zhonghua11, 3=W0169, 4=W1735, 5=W3033):
#   sbatch --array=4 /data/malong/benchmark/code/call_syri_sv.slurm

#SBATCH --job-name=rice_syri
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60
#SBATCH --time=7-00:00:00
#SBATCH --array=0-5
#SBATCH --output=/data/malong/benchmark/rice_syri_%A_%a.out
#SBATCH --error=/data/malong/benchmark/rice_syri_%A_%a.err

set -Eeuo pipefail

BASE=/data/malong/benchmark
REF="$BASE/graph/pggb/rice36_by_sample/NIP-T2T.fa"
TRUTH="$BASE/truth"
SYRI_ENV=/data/malong/miniconda3/envs/benchmark
# SyRI invokes show-snps through /bin/sh internally, so the environment's
# bin directory must be on PATH in addition to using absolute executable paths.
export PATH="$SYRI_ENV/bin:$PATH"

THREADS=${SLURM_CPUS_PER_TASK:-48}
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

NUCMER="$SYRI_ENV/bin/nucmer"
DELTA_FILTER="$SYRI_ENV/bin/delta-filter"
SHOW_COORDS="$SYRI_ENV/bin/show-coords"
SHOW_SNPS="$SYRI_ENV/bin/show-snps"
SYRI="$SYRI_ENV/bin/syri"
BCFTOOLS="$SYRI_ENV/bin/bcftools"
SAMTOOLS="$SYRI_ENV/bin/samtools"

for exe in "$NUCMER" "$DELTA_FILTER" "$SHOW_COORDS" "$SHOW_SNPS" "$SYRI" "$BCFTOOLS" "$SAMTOOLS"; do
  [[ -x "$exe" ]] || { echo "Missing executable: $exe" >&2; exit 127; }
done
[[ -s "$REF" ]] || { echo "Missing reference: $REF" >&2; exit 2; }
[[ -s "$asm" ]] || { echo "Missing assembly for $sample: $asm" >&2; exit 2; }

OUT="$TRUTH/$sample/syri"
mkdir -p "$OUT"
exec > >(tee "$OUT/$sample.syri.log") 2>&1

echo "[$(date -Is)] sample=$sample"
echo "[$(date -Is)] reference=$REF"
echo "[$(date -Is)] assembly=$asm"
echo "[$(date -Is)] SyRI version: $($SYRI --version 2>&1 | head -1 || true)"

# The reference and assemblies use names such as NIP-T2T#1#Chr1 and
# 9311#1#Chr1. SyRI's chromosome matching is more reliable when both files
# use the common Chr1..Chr12 names. The original names are restored in the
# output VCF using chr.rename.tsv.
REF_SYRI="$OUT/reference.Chr.fa"
QRY_SYRI="$OUT/$sample.Chr.fa"
if [[ ! -s "$REF_SYRI" ]]; then
  sed -E 's/^>[^#[:space:]]+#1#/>/' "$REF" > "$REF_SYRI"
fi
if [[ ! -s "$QRY_SYRI" ]]; then
  sed -E 's/^>[^#[:space:]]+#1#/>/' "$asm" > "$QRY_SYRI"
fi

RENAME="$OUT/chr.rename.tsv"
if [[ ! -s "$RENAME" ]]; then
  grep '^>' "$REF" | sed -E 's/^>([^[:space:]]+).*/\1/' \
    | awk -F'#' '{print $3 "\t" $0}' > "$RENAME"
fi

PREFIX="$OUT/$sample"
DELTA="$PREFIX.delta"
FILTERED_DELTA="$PREFIX.filtered.delta"
COORDS="$PREFIX.coords"

# Nucmer --maxmatch plus delta-filter -m retains repeat-aware alignments,
# which is important for identifying rearrangements. The same parameters are
# used for all samples to make the caller comparison reproducible.
if [[ ! -s "$DELTA" ]]; then
  "$NUCMER" --maxmatch -l 100 -c 500 -t "$THREADS" \
    -p "$PREFIX" "$REF_SYRI" "$QRY_SYRI"
fi
if [[ ! -s "$FILTERED_DELTA" ]]; then
  "$DELTA_FILTER" -m -i 90 -l 100 "$DELTA" > "$FILTERED_DELTA"
fi
if [[ ! -s "$COORDS" ]]; then
  "$SHOW_COORDS" -THrd "$FILTERED_DELTA" > "$COORDS"
fi

# SyRI writes several intermediate files and a VCF into --dir. Keep these
# files as audit evidence. --all includes duplication-related calls; --maxsize
# makes very large events symbolic while retaining their coordinates.
if [[ ! -s "$OUT/$sample.syri.vcf" ]]; then
  "$SYRI" -c "$COORDS" -d "$FILTERED_DELTA" \
    -r "$REF_SYRI" -q "$QRY_SYRI" -F T \
    --dir "$OUT/" --prefix "$sample." --samplename "$sample" \
    --all --maxsize 100000 --nc "$THREADS" -k \
    --lf "$sample.syri.log"
fi

RAW_VCF="$OUT/$sample.syri.vcf"
[[ -s "$RAW_VCF" ]] || { echo "SyRI VCF not found: $RAW_VCF" >&2; exit 1; }

# Restore the NIP-T2T reference contig names, sort, and index the complete
# SyRI VCF. The complete file is retained because it contains structural,
# short-variant, and rearrangement annotations.
ALL_VCF="$OUT/$sample.syri.all.vcf.gz"
if [[ ! -s "$ALL_VCF" ]]; then
  "$BCFTOOLS" annotate --rename-chrs "$RENAME" "$RAW_VCF" \
    | "$BCFTOOLS" sort -Oz -o "$ALL_VCF"
fi
"$BCFTOOLS" index -f "$ALL_VCF"

# SyRI encodes structural events with symbolic ALT alleles (<DEL>, <INS>,
# <INV>, <TRANS>, <DUP>, ...). Keep only those records for direct comparison
# with the existing SVIM-asm SV VCF; the complete SyRI VCF remains available.
SV_VCF="$OUT/$sample.syri.sv.vcf.gz"
if [[ ! -s "$SV_VCF" ]]; then
  awk 'BEGIN { OFS="\t" }
       /^#/ { print; next }
       {
         n=split($5, a, ","); keep=0;
         for (i=1; i<=n; i++) if (a[i] ~ /^</) keep=1;
         if (keep) print;
       }' "$RAW_VCF" \
    | "$BCFTOOLS" annotate --rename-chrs "$RENAME" \
    | "$BCFTOOLS" sort -Oz -o "$SV_VCF"
fi
"$BCFTOOLS" index -f "$SV_VCF"

SVIM_VCF="$TRUTH/$sample/final/$sample.sv.vcf.gz"
SUMMARY="$OUT/$sample.syri.summary.tsv"
{
  printf 'sample\tmethod\tvcf\tn_records\n'
  printf '%s\tSyRI_all\t%s\t%s\n' "$sample" "$ALL_VCF" \
    "$($BCFTOOLS view -H "$ALL_VCF" | wc -l)"
  printf '%s\tSyRI_SV_symbolic\t%s\t%s\n' "$sample" "$SV_VCF" \
    "$($BCFTOOLS view -H "$SV_VCF" | wc -l)"
  if [[ -s "$SVIM_VCF" ]]; then
    printf '%s\tSVIM-asm_existing\t%s\t%s\n' "$sample" "$SVIM_VCF" \
      "$($BCFTOOLS view -H "$SVIM_VCF" | wc -l)"
  fi
} > "$SUMMARY"

cat > "$OUT/$sample.syri.manifest.tsv" <<EOF
sample	$sample
reference	$REF
assembly	$asm
alignment	nucmer4_maxmatch_delta_filter_m_i90_l100
caller	SyRI
syri_all_vcf	$ALL_VCF
syri_sv_vcf	$SV_VCF
existing_svim_asm_vcf	$SVIM_VCF
sv_size_policy	structural_symbolic_records; SyRI maxsize=100000; inspect/filter exact size range before scoring
graph_membership	$(case "$sample" in W0169|W1735|W3033) echo in_graph ;; *) echo out_of_graph ;; esac)
EOF

echo "[$(date -Is)] completed"
cat "$SUMMARY"
