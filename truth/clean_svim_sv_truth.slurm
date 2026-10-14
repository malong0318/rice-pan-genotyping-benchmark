#!/usr/bin/env bash
# Filter SVIM-asm VCFs to a strict 50 bp--100 kb SV set.
# The original SVIM-asm VCF is retained unchanged.
#
# Submit all six samples:
#   sbatch /data/malong/benchmark/code/clean_svim_sv_truth.slurm
# Submit one sample (0=9311, 1=IR64, 2=zhonghua11, 3=W0169, 4=W1735, 5=W3033):
#   sbatch --array=4 /data/malong/benchmark/code/clean_svim_sv_truth.slurm

#SBATCH --job-name=clean_svim_sv
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --array=0-5
#SBATCH --output=/data/malong/benchmark/clean_svim_sv_%A_%a.out
#SBATCH --error=/data/malong/benchmark/clean_svim_sv_%A_%a.err

set -Eeuo pipefail

BASE=/data/malong/benchmark
TRUTH="$BASE/truth"
BCFTOOLS=/data/malong/miniconda3/envs/benchmark/bin/bcftools
SAMPLES=(9311 IR64 zhonghua11 W0169 W1735 W3033)
SIZE_MIN=50
SIZE_MAX=100000

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]:-}
[[ -n "$sample" ]] || { echo "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-unset}" >&2; exit 2; }

IN_VCF="$TRUTH/$sample/final/$sample.sv.vcf.gz"
OUT_DIR="$TRUTH/$sample/final/strict"
OUT_VCF="$OUT_DIR/$sample.svim.strict.sv.vcf.gz"
OUT_BED="$OUT_DIR/$sample.svim.strict.sv.bed"
OUT_STATS="$OUT_DIR/$sample.svim.strict.sv.stats.txt"
OUT_SUMMARY="$OUT_DIR/$sample.svim.strict.summary.tsv"
OUT_MANIFEST="$OUT_DIR/$sample.svim.strict.manifest.tsv"
KEEP_IDS="$OUT_DIR/$sample.keep_ids.txt"

[[ -x "$BCFTOOLS" ]] || { echo "Missing executable: $BCFTOOLS" >&2; exit 127; }
[[ -s "$IN_VCF" ]] || { echo "Missing SVIM-asm VCF: $IN_VCF" >&2; exit 2; }
mkdir -p "$OUT_DIR"
exec > >(tee "$OUT_DIR/$sample.clean_svim.log") 2>&1

echo "[$(date -Is)] sample=$sample"
echo "[$(date -Is)] input=$IN_VCF"

# Prefer SVLEN because insertions have END=POS and therefore cannot be sized
# from END-POS. For records without SVLEN (for example inversions), use the
# inclusive reference span END-POS+1. BND records have no defined SV length
# and are excluded from a length-constrained truth set.
if [[ ! -s "$KEEP_IDS" ]]; then
  "$BCFTOOLS" query -f '%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\t%POS\n' "$IN_VCF" \
    | awk -F '\t' -v min="$SIZE_MIN" -v max="$SIZE_MAX" '
      function abs(x) { return x < 0 ? -x : x }
      {
        len=0
        if ($3 != "." && $3 != "0") len=abs($3)
        else if ($4 != "." && $4+0 >= $5+0) len=$4-$5+1
        if ($2 != "BND" && len >= min && len <= max) print $1
      }' > "$KEEP_IDS"
fi

[[ -s "$KEEP_IDS" ]] || { echo "No SV records passed the size filter" >&2; exit 1; }

# IDs are unique in the current SVIM-asm VCFs. Filter by ID so all original
# INFO/FORMAT fields and headers are retained.
if [[ ! -s "$OUT_VCF" ]]; then
  "$BCFTOOLS" view -i "ID=@$KEEP_IDS" -Oz -o "$OUT_VCF" "$IN_VCF"
fi
"$BCFTOOLS" index -f "$OUT_VCF"

# Use a one-base BED interval for insertions (END=POS), otherwise use the
# inclusive VCF interval converted to 0-based, half-open BED coordinates.
if [[ ! -s "$OUT_BED" ]]; then
  "$BCFTOOLS" query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$OUT_VCF" \
    | awk 'BEGIN { OFS="\t" }
      {
        end=$3; if (end=="." || end < $2) end=$2
        print $1, $2-1, end, $4, $5, $6
      }' \
    | sort -k1,1 -k2,2n > "$OUT_BED"
fi

# Validate that every output record has a defined size in the requested range.
invalid=$(
  "$BCFTOOLS" query -f '%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\t%POS\n' "$OUT_VCF" \
    | awk -F '\t' -v min="$SIZE_MIN" -v max="$SIZE_MAX" '
      function abs(x) { return x < 0 ? -x : x }
      {
        len=0
        if ($3 != "." && $3 != "0") len=abs($3)
        else if ($4 != "." && $4+0 >= $5+0) len=$4-$5+1
        if ($2 == "BND" || len < min || len > max) n++
      }
      END { print n+0 }'
)
[[ "$invalid" == 0 ]] || { echo "Invalid records in strict VCF: $invalid" >&2; exit 1; }

"$BCFTOOLS" stats "$OUT_VCF" > "$OUT_STATS"
total=$("$BCFTOOLS" view -H "$OUT_VCF" | wc -l)
printf 'sample\tinput_vcf\tstrict_vcf\tn_records\tmin_bp\tmax_bp\n%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$sample" "$IN_VCF" "$OUT_VCF" "$total" "$SIZE_MIN" "$SIZE_MAX" > "$OUT_SUMMARY"

cat > "$OUT_MANIFEST" <<EOF
sample	$sample
input_vcf	$IN_VCF
strict_vcf	$OUT_VCF
strict_bed	$OUT_BED
caller	SVIM-asm haploid
size_source	abs(INFO/SVLEN), fallback to END-POS+1 when SVLEN is missing
size_policy	${SIZE_MIN}-${SIZE_MAX} bp
excluded	BND records with undefined SV length
coordinate_policy	VCF 1-based inclusive; BED 0-based half-open
graph_membership	$(case "$sample" in W0169|W1735|W3033) echo in_graph ;; *) echo out_of_graph ;; esac)
EOF

rm -f "$KEEP_IDS"
echo "[$(date -Is)] completed"
cat "$OUT_SUMMARY"
