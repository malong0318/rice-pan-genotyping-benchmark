#!/usr/bin/env bash
# Reproducibly generate the final assembly-derived truth sets.
#
# Usage:
#   bash 01generate_truth_sets.sh
#   THREADS=60 bash 01generate_truth_sets.sh
#   ONLY_SAMPLE=9311 bash 01generate_truth_sets.sh

set -Eeuo pipefail

BASE=${BASE:-/data/malong/benchmark}
ENV=${ENV:-/data/malong/miniconda3/envs/benchmark}
ROOT=${ROOT:-$BASE/truth}
REF=${REF:-$BASE/graph/pggb/rice36_by_sample/NIP-T2T.fa}
REF_CHR=${REF_CHR:-$BASE/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa}
THREADS=${THREADS:-60}
DIPCALL_THREADS=${DIPCALL_THREADS:-$(( THREADS > 4 ? THREADS / 2 : THREADS ))}
MAKE_JOBS=${MAKE_JOBS:-2}

RUN_DIPCALL="$ENV/bin/run-dipcall"
MINIMAP2="$ENV/bin/minimap2"
SAMTOOLS="$ENV/bin/samtools"
BCFTOOLS="$ENV/bin/bcftools"
BEDTOOLS="$ENV/bin/bedtools"
SVIM="$ENV/bin/svim-asm"
TABIX="$ENV/bin/tabix"
BGZIP="$ENV/bin/bgzip"
[[ -x "$TABIX" ]] || TABIX=/data/malong/miniconda3/bin/tabix
[[ -x "$BGZIP" ]] || BGZIP=/data/malong/miniconda3/bin/bgzip

SAMPLES=(9311 IR64 zhonghua11 W0169 W1735 W3033)
declare -A ASSEMBLY=(
  [9311]="$BASE/assembly/9311/9311_pansn.fa"
  [IR64]="$BASE/assembly/IR64/IR64_pansn.fa"
  [zhonghua11]="$BASE/assembly/zhonghua11/zhonghua11_pansn.fa"
  [W0169]="$BASE/graph/pggb/rice36_by_sample/W0169.fa"
  [W1735]="$BASE/graph/pggb/rice36_by_sample/W1735.fa"
  [W3033]="$BASE/graph/pggb/rice36_by_sample/W3033.fa"
)

LOG="$ROOT/generate_truth_sets.log"
MAP="$ROOT/truth_to_chr.tsv"
SUMMARY="$ROOT/truth_sets.summary.tsv"
timestamp() { date '+%Y-%m-%dT%H:%M:%S%z'; }
log() { printf '[%s] %s\n' "$(timestamp)" "$*" | tee -a "$LOG"; }
die() { log "ERROR: $*"; exit 1; }
on_error() {
  local rc=$?
  log "FAILED: line=${BASH_LINENO[0]:-unknown}, command=${BASH_COMMAND}, exit=$rc"
  exit "$rc"
}
trap on_error ERR

mkdir -p "$ROOT"
exec > >(tee -a "$LOG") 2>&1
export PATH="$ENV/bin:$PATH"

for exe in "$RUN_DIPCALL" "$MINIMAP2" "$SAMTOOLS" "$BCFTOOLS" "$BEDTOOLS" "$SVIM" "$TABIX" "$BGZIP"; do
  [[ -x "$exe" ]] || die "missing executable: $exe"
done
command -v make >/dev/null 2>&1 || die "make is not available"
[[ -s "$REF" ]] || die "reference is missing or empty: $REF"

if [[ ! -s "$REF.fai" ]]; then
  log "building reference index: $REF.fai"
  "$SAMTOOLS" faidx "$REF"
fi
if [[ ! -s "$REF_CHR" ]]; then
  log "building renamed reference: $REF_CHR"
  tmp="$REF_CHR.tmp.$$"
  awk '/^>/{sub(/^>[^#]+#1#/, ">")}{print}' "$REF" > "$tmp"
  mv "$tmp" "$REF_CHR"
fi
if [[ ! -s "$REF_CHR.fai" ]]; then
  log "building renamed reference index: $REF_CHR.fai"
  "$SAMTOOLS" faidx "$REF_CHR"
fi

if [[ ! -s "$MAP" ]]; then
  awk -F '\t' '{print "NIP-T2T#1#" $1 "\t" $1}' "$REF_CHR.fai" > "$MAP.tmp.$$"
  mv "$MAP.tmp.$$" "$MAP"
fi
for c in Chr{1..12}; do
  grep -F $'NIP-T2T#1#'"$c"$'\t'"$c" "$MAP" >/dev/null || die "mapping missing: $c"
done

log "reference=$REF"
log "evaluation_reference=$REF_CHR"
log "threads=$THREADS dipcall_threads=$DIPCALL_THREADS make_jobs=$MAKE_JOBS"
log "minimap2=$($MINIMAP2 --version | head -1)"
log "$($SAMTOOLS --version | head -1)"
log "$($BCFTOOLS --version | head -1)"
log "$($BEDTOOLS --version 2>&1 | head -1)"
log "$($SVIM --help 2>&1 | head -1)"
log "dipcall pipeline: 0.3 workflow; standalone binary version not recorded"

printf 'sample\ttype\tfinal_records\tfinal_file\tevaluation_file\n' > "$SUMMARY.tmp.$$"
sample_list=("${SAMPLES[@]}")
if [[ -n "${ONLY_SAMPLE:-}" ]]; then sample_list=("$ONLY_SAMPLE"); fi

for sample in "${sample_list[@]}"; do
  asm=${ASSEMBLY[$sample]:-}
  [[ -n "$asm" ]] || die "unknown sample: $sample"
  [[ -s "$asm" ]] || die "assembly is missing or empty: $asm"

  DIPDIR="$ROOT/$sample/dipcall"
  SVDIR="$ROOT/$sample/sv"
  FINAL="$ROOT/$sample/final"
  OUT="$ROOT/$sample/Chr_vcf"
  SORTTMP="$OUT/.sort_tmp"
  mkdir -p "$DIPDIR" "$SVDIR" "$FINAL" "$OUT" "$SORTTMP"
  log "===== START sample=$sample assembly=$asm ====="

  # A. dipcall small-variant truth. The single consensus assembly is supplied
  # as both haplotypes, producing an assembly-derived pseudodiploid truth set.
  MAKEFILE="$DIPDIR/$sample.mak"
  RAW_DIP_VCF="$DIPDIR/$sample.dip.vcf.gz"
  RAW_DIP_BED="$DIPDIR/$sample.dip.bed"
  if [[ ! -s "$RAW_DIP_VCF" || ! -s "$RAW_DIP_BED" ]]; then
    log "[$sample] running dipcall"
    tmp_make="$MAKEFILE.tmp.$$"
    "$RUN_DIPCALL" -a -t "$DIPCALL_THREADS" "$sample" "$REF" "$asm" "$asm" > "$tmp_make"
    mv "$tmp_make" "$MAKEFILE"
    make -C "$DIPDIR" -f "$MAKEFILE" -j "$MAKE_JOBS"
  else
    log "[$sample] reuse dipcall outputs"
  fi
  [[ -s "$RAW_DIP_VCF" && -s "$RAW_DIP_BED" ]] || die "dipcall output incomplete: $sample"

  RAW_NORM="$FINAL/$sample.norm.vcf.gz"
  if [[ ! -s "$RAW_NORM" ]]; then
    log "[$sample] normalize and split dipcall VCF"
    tmp="$RAW_NORM.tmp.$$"
    "$BCFTOOLS" norm -f "$REF" -m -any "$RAW_DIP_VCF" -Ou | "$BCFTOOLS" sort -Oz -o "$tmp"
    mv "$tmp" "$RAW_NORM"
  fi
  "$BCFTOOLS" index -f "$RAW_NORM"
  RAW_SNP="$FINAL/$sample.snp.vcf.gz"
  RAW_INDEL="$FINAL/$sample.indel.vcf.gz"
  [[ -s "$RAW_SNP" ]] || "$BCFTOOLS" view -v snps -Oz -o "$RAW_SNP" "$RAW_NORM"
  [[ -s "$RAW_INDEL" ]] || "$BCFTOOLS" view -v indels -Oz -o "$RAW_INDEL" "$RAW_NORM"
  "$BCFTOOLS" index -f "$RAW_SNP"
  "$BCFTOOLS" index -f "$RAW_INDEL"
  CONF_BED="$FINAL/$sample.confident.bed"
  if [[ ! -s "$CONF_BED" ]]; then
    "$BEDTOOLS" sort -i "$RAW_DIP_BED" > "$CONF_BED"
  fi

  # B. minimap2 plus SVIM-asm haploid SV truth.
  BAM="$SVDIR/$sample.to_ref.bam"
  if [[ ! -s "$BAM" ]]; then
    log "[$sample] align assembly for SVIM-asm"
    "$MINIMAP2" -ax asm5 --cs -t "$THREADS" "$REF" "$asm" | "$SAMTOOLS" sort -@ "$THREADS" -o "$BAM" -
    "$SAMTOOLS" index "$BAM"
  fi
  if [[ ! -s "$SVDIR/calls/variants.vcf" && ! -s "$SVDIR/calls/variants.vcf.gz" ]]; then
    log "[$sample] run SVIM-asm haploid"
    "$SVIM" haploid "$SVDIR/calls" "$BAM" "$REF" --sample "$sample"
  fi
  RAW_SV="$FINAL/$sample.sv.vcf.gz"
  if [[ ! -s "$RAW_SV" ]]; then
    if [[ -s "$SVDIR/calls/variants.vcf" ]]; then
      "$BCFTOOLS" sort "$SVDIR/calls/variants.vcf" -Oz -o "$RAW_SV"
    elif [[ -s "$SVDIR/calls/variants.vcf.gz" ]]; then
      "$BCFTOOLS" sort "$SVDIR/calls/variants.vcf.gz" -Oz -o "$RAW_SV"
    else
      die "SVIM-asm VCF is missing: $sample"
    fi
  fi
  "$BCFTOOLS" index -f "$RAW_SV"

  # C. Common chromosome conversion, header rewrite, normalization, sorting.
  normalize_type() {
    local type=$1 input=$2 stem="$OUT/$sample.$1"
    local renamed="$stem.renamed.vcf.gz" reheadered="$stem.reheadered.vcf.gz"
    local sorted="$stem.sorted.vcf.gz" normed="$stem.norm.vcf.gz" final="$stem.final.vcf.gz"
    if [[ ! -s "$renamed" ]]; then
      if [[ "$type" == sv ]]; then
        local raw="$renamed.raw.tmp.$$" tmp="$renamed.tmp.$$"
        "$BCFTOOLS" annotate --rename-chrs "$MAP" -Oz -o "$raw" "$input"
        "$BGZIP" -cd "$raw" | sed 's/NIP-T2T#1#Chr/Chr/g' | "$BGZIP" -c > "$tmp"
        rm -f "$raw"; mv "$tmp" "$renamed"
      else
        local tmp="$renamed.tmp.$$"
        "$BCFTOOLS" annotate --rename-chrs "$MAP" -Oz -o "$tmp" "$input"
        mv "$tmp" "$renamed"
      fi
    fi
    if [[ ! -s "$reheadered" ]]; then
      local tmp="$reheadered.tmp.$$"
      "$BCFTOOLS" reheader -f "$REF_CHR.fai" -o "$tmp" "$renamed"; mv "$tmp" "$reheadered"
    fi
    if [[ ! -s "$sorted" ]]; then
      local tmp="$sorted.tmp.$$"
      "$BCFTOOLS" sort -T "$SORTTMP" -Oz -o "$tmp" "$reheadered"; mv "$tmp" "$sorted"
    fi
    "$TABIX" -f -p vcf "$sorted"
    if [[ ! -s "$normed" ]]; then
      local tmp="$normed.tmp.$$"
      if [[ "$type" == sv ]]; then
        "$BCFTOOLS" norm -f "$REF_CHR" -c s -m -any -Oz -o "$tmp" "$sorted"
      else
        "$BCFTOOLS" norm -f "$REF_CHR" -c e -m -any -Oz -o "$tmp" "$sorted"
      fi
      mv "$tmp" "$normed"
    fi
    if [[ ! -s "$final" ]]; then
      local tmp="$final.tmp.$$"
      "$BCFTOOLS" sort -T "$SORTTMP" -Oz -o "$tmp" "$normed"; mv "$tmp" "$final"
    fi
    "$TABIX" -f -p vcf "$final"
  }
  normalize_type snp "$RAW_SNP"
  normalize_type indel "$RAW_INDEL"
  normalize_type sv "$RAW_SV"

  # Final true_snp and true_indel used by RTG vcfeval.
  TRUE_SNP="$OUT/$sample.snp.confident.final.vcf.gz"
  TRUE_INDEL="$OUT/$sample.indel.small.confident.final.vcf.gz"
  if [[ ! -s "$TRUE_SNP" ]]; then
    "$BCFTOOLS" view -f . -i 'strlen(REF)=1 && strlen(ALT)=1 && REF!=ALT' -Oz \
      -o "$TRUE_SNP" "$OUT/$sample.snp.final.vcf.gz"
  fi
  if [[ ! -s "$TRUE_INDEL" ]]; then
    "$BCFTOOLS" view -f . -i 'strlen(REF)!=strlen(ALT) && strlen(REF)<50 && strlen(ALT)<50' -Oz \
      -o "$TRUE_INDEL" "$OUT/$sample.indel.final.vcf.gz"
  fi
  "$TABIX" -f -p vcf "$TRUE_SNP"
  "$TABIX" -f -p vcf "$TRUE_INDEL"

  # Final true_sv used by Truvari: PASS records, excluding BND records.
  TRUE_SV="$OUT/$sample.sv.truth.PASS.noBND.vcf.gz"
  if [[ ! -s "$TRUE_SV" ]]; then
    "$BCFTOOLS" view -f PASS -e 'INFO/SVTYPE="BND"' -Oz -o "$TRUE_SV" "$OUT/$sample.sv.final.vcf.gz"
  fi
  "$TABIX" -f -p vcf "$TRUE_SV"

  TRUE_BED="$OUT/$sample.confident.Chr.bed"
  if [[ ! -s "$TRUE_BED" ]]; then
    awk 'BEGIN{OFS="\t"} {$1="Chr" substr($1, index($1, "Chr")+3); print}' "$CONF_BED" > "$TRUE_BED"
  fi

  snp_n=$($BCFTOOLS view -H "$TRUE_SNP" | wc -l)
  indel_n=$($BCFTOOLS view -H "$TRUE_INDEL" | wc -l)
  sv_n=$($BCFTOOLS view -H "$TRUE_SV" | wc -l)
  printf '%s\tsnp\t%s\t%s\t%s\n' "$sample" "$snp_n" "$TRUE_SNP" "$TRUE_SNP" >> "$SUMMARY.tmp.$$"
  printf '%s\tindel\t%s\t%s\t%s\n' "$sample" "$indel_n" "$TRUE_INDEL" "$TRUE_INDEL" >> "$SUMMARY.tmp.$$"
  printf '%s\tsv\t%s\t%s\t%s\n' "$sample" "$sv_n" "$TRUE_SV" "$TRUE_SV" >> "$SUMMARY.tmp.$$"
  printf 'sample\t%s\nreference\t%s\nassembly\t%s\ntrue_snp\t%s\ntrue_indel\t%s\ntrue_sv\t%s\nconfident_bed\t%s\n' \
    "$sample" "$REF_CHR" "$asm" "$TRUE_SNP" "$TRUE_INDEL" "$TRUE_SV" "$TRUE_BED" \
    > "$OUT/$sample.truth.manifest.tsv"
  log "===== DONE sample=$sample SNP=$snp_n indel=$indel_n SV=$sv_n ====="
done

mv "$SUMMARY.tmp.$$" "$SUMMARY"
log "SUCCESS: truth generation completed"
log "summary=$SUMMARY"
