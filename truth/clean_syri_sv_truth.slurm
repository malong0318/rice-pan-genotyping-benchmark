#!/usr/bin/env bash
# Create strict, assembly-based SyRI SV truth VCFs from the completed SyRI runs.
# This is a post-processing job; it does not rerun nucmer or SyRI.
#
# Submit all six samples:
#   sbatch /data/malong/benchmark/code/clean_syri_sv_truth.slurm
# Submit one sample (0=9311, 1=IR64, 2=zhonghua11, 3=W0169, 4=W1735, 5=W3033):
#   sbatch --array=4 /data/malong/benchmark/code/clean_syri_sv_truth.slurm

#SBATCH --job-name=clean_syri_sv
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --array=0-5
#SBATCH --output=/data/malong/benchmark/clean_syri_sv_%A_%a.out
#SBATCH --error=/data/malong/benchmark/clean_syri_sv_%A_%a.err

set -Eeuo pipefail

BASE=/data/malong/benchmark
TRUTH="$BASE/truth"
BCFTOOLS=/data/malong/miniconda3/envs/benchmark/bin/bcftools
SAMPLES=(9311 IR64 zhonghua11 W0169 W1735 W3033)

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]:-}
[[ -n "$sample" ]] || { echo "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-unset}" >&2; exit 2; }

IN_DIR="$TRUTH/$sample/syri"
IN_VCF="$IN_DIR/$sample.syri.sv.vcf.gz"
OUT_DIR="$IN_DIR/strict"
OUT_VCF="$OUT_DIR/$sample.syri.strict.sv.vcf.gz"
OUT_BED="$OUT_DIR/$sample.syri.strict.sv.bed"
OUT_STATS="$OUT_DIR/$sample.syri.strict.sv.stats.txt"
OUT_SUMMARY="$OUT_DIR/$sample.syri.strict.summary.tsv"
OUT_MANIFEST="$OUT_DIR/$sample.syri.strict.manifest.tsv"

[[ -x "$BCFTOOLS" ]] || { echo "Missing executable: $BCFTOOLS" >&2; exit 127; }
[[ -s "$IN_VCF" ]] || { echo "Missing SyRI VCF: $IN_VCF" >&2; exit 2; }
mkdir -p "$OUT_DIR"
exec > >(tee "$OUT_DIR/$sample.clean.log") 2>&1

echo "[$(date -Is)] sample=$sample"
echo "[$(date -Is)] input=$IN_VCF"

# SyRI's VCF contains both structural events and annotation records. For a
# strict SV truth set, retain only primary structural rearrangements:
# DUP, INV, INVDP, INVTR, and TRANS. Exclude SYN/SYNAL, *AL alignment
# annotations, NOTAL/HDR regions, and VarType=ShV short-variant records.
# Restrict all retained events to the same 50 bp--100 kb policy used for SV
# scoring. END is inclusive in the VCF expression and POS is 1-based.
FILTER='INFO/VarType="SR" && (ALT="<DUP>" || ALT="<INV>" || ALT="<INVDP>" || ALT="<INVTR>" || ALT="<TRANS>") && INFO/END>=POS+49 && INFO/END<=POS+99999'

if [[ ! -s "$OUT_VCF" ]]; then
  "$BCFTOOLS" view -i "$FILTER" -Oz -o "$OUT_VCF" "$IN_VCF"
fi
"$BCFTOOLS" index -f "$OUT_VCF"

# Convert the VCF intervals to 0-based, half-open BED intervals for callable
# region and overlap calculations. The VCF END coordinate is inclusive.
if [[ ! -s "$OUT_BED" ]]; then
  "$BCFTOOLS" query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%ALT\n' "$OUT_VCF" \
    | awk 'BEGIN { OFS="\t" } $3 >= $2 && $3 != "." { print $1, $2-1, $3, $4, $5 }' \
    | sort -k1,1 -k2,2n > "$OUT_BED"
fi

# Fail closed if any retained record has an invalid interval or an excluded
# annotation class. This makes the file suitable for downstream scoring.
invalid=$(
  "$BCFTOOLS" query -f '%POS\t%INFO/END\t%INFO/VarType\t%ALT\n' "$OUT_VCF" \
    | awk '$2=="." || $2 < $1 || $3!="SR" || ($4!="<DUP>" && $4!="<INV>" && $4!="<INVDP>" && $4!="<INVTR>" && $4!="<TRANS>") { n++ } END { print n+0 }'
)
[[ "$invalid" == 0 ]] || { echo "Invalid records in strict VCF: $invalid" >&2; exit 1; }

"$BCFTOOLS" stats "$OUT_VCF" > "$OUT_STATS"

total=$("$BCFTOOLS" view -H "$OUT_VCF" | wc -l)
printf 'sample\tinput_vcf\tstrict_vcf\tn_records\n%s\t%s\t%s\t%s\n' \
  "$sample" "$IN_VCF" "$OUT_VCF" "$total" > "$OUT_SUMMARY"

cat > "$OUT_MANIFEST" <<EOF
sample	$sample
input_vcf	$IN_VCF
strict_vcf	$OUT_VCF
strict_bed	$OUT_BED
caller	SyRI 1.8.2
source_filter	INFO/VarType=SR
retained_alt_types	<DUP>,<INV>,<INVDP>,<INVTR>,<TRANS>
excluded_records	SYN/SYNAL,*AL,NOTAL,HDR,VarType=ShV,CPG,CPL,TDM,INS,DEL
size_policy	50-100000 bp; POS-based inclusive END interval
coordinate_policy	VCF 1-based inclusive; BED 0-based half-open
graph_membership	$(case "$sample" in W0169|W1735|W3033) echo in_graph ;; *) echo out_of_graph ;; esac)
EOF

echo "[$(date -Is)] completed"
cat "$OUT_SUMMARY"
