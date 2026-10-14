#!/usr/bin/env bash
set -Eeuo pipefail

TOOL=${1:?usage: $0 vg|pangenie}
case "$TOOL" in
  vg)
    SOURCE_ROOT=/data/malong/benchmark/01classify_vcf/vg
    OUT_ROOT=/data/malong/benchmark/results_vg2
    GRAPHS=(MC minigraph)
    ;;
  pangenie)
    SOURCE_ROOT=/data/malong/benchmark/01classify_vcf/pangenie
    OUT_ROOT=/data/malong/benchmark/results_pangenie2
    GRAPHS=(MC PGGB minigraph)
    ;;
  bayestyper)
    SOURCE_ROOT=/data/malong/benchmark/01classify_vcf/bayestyper
    OUT_ROOT=/data/malong/benchmark/results_bayestyper2
    GRAPHS=(MC PGGB minigraph)
    ;;
  *) printf 'usage: %s vg|pangenie|bayestyper\n' "$0" >&2; exit 2 ;;
esac

ROOT=/data/malong/benchmark
REF=$ROOT/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa
BCF=/data/malong/miniconda3/envs/benchmark/bin/bcftools
TABIX=/data/malong/miniconda3/envs/benchmark/bin/tabix
THREADS=${EVAL_PREP_THREADS:-16}
SORT_MEM=${EVAL_PREP_SORT_MEM:-4G}
ONLY_SAMPLE=${EVAL_PREP_ONLY_SAMPLE:-}
ONLY_DEPTH=${EVAL_PREP_ONLY_DEPTH:-}
ONLY_GRAPH=${EVAL_PREP_ONLY_GRAPH:-}
FORCE=${EVAL_PREP_FORCE:-0}
SAMPLES=(9311 IR64 W0169 W1735 W3033 zhonghua11)
DEPTHS=(5x 10x 15x 20x)

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*"; }
fail() { log "ERROR: $*" >&2; exit 1; }
trap 'rc=$?; log "FAILED line=$LINENO command=$BASH_COMMAND exit=$rc" >&2; exit "$rc"' ERR
need_file() { [[ -s "$1" ]] || fail "missing or empty file: $1"; }
need_bin() { [[ -x "$1" ]] || fail "missing executable: $1"; }
count_records() { "$BCF" index -n "$1"; }

need_bin "$BCF"
need_bin "$TABIX"
need_file "$REF"
need_file "$REF.fai"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || fail "invalid EVAL_PREP_THREADS: $THREADS"

source_for() {
  local sample=$1 depth=$2 graph=$3
  if [[ "$TOOL" == vg ]]; then
    if [[ "$graph" == MC ]]; then
      printf '%s/%s/%s/MC/%s.MC.vg.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample"
    else
      printf '%s/%s/%s/minigraph/%s.minigraph.vg.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample"
    fi
  else
    if [[ "$TOOL" == pangenie ]]; then
      case "$graph" in
        MC) printf '%s/%s/%s/MC/%s.MC.pangenie.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample" ;;
        PGGB) printf '%s/%s/%s/PGGB/%s.PGGB.pangenie.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample" ;;
        minigraph) printf '%s/%s/%s/minigraph/%s.minigraph.pangenie.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample" ;;
      esac
    else
      case "$graph" in
        MC) printf '%s/%s/%s/MC/%s.MC.bayestyper.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample" ;;
        PGGB) printf '%s/%s/%s/PGGB/%s.PGGB.bayestyper.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample" ;;
        minigraph) printf '%s/%s/%s/minigraph/%s.minigraph.bayestyper.PASS.nonref.split.sorted.vcf.gz\n' "$SOURCE_ROOT" "$sample" "$depth" "$sample" ;;
      esac
    fi
  fi
}

write_class() {
  local input=$1 output=$2 expression=$3
  if [[ "$FORCE" == 1 || ! -s "$output" || ! -s "$output.tbi" ]]; then
    local tmp="$output.tmp.$$"
    rm -f "$tmp" "$tmp.tbi"
    "$BCF" view --threads "$THREADS" -i "$expression" -Oz -o "$tmp" "$input"
    "$TABIX" -f -p vcf "$tmp"
    mv -f "$tmp" "$output"
    mv -f "$tmp.tbi" "$output.tbi"
  fi
}

write_complement() {
  local input=$1 output=$2 expression=$3
  if [[ "$FORCE" == 1 || ! -s "$output" || ! -s "$output.tbi" ]]; then
    local tmp="$output.tmp.$$"
    rm -f "$tmp" "$tmp.tbi"
    "$BCF" view --threads "$THREADS" -e "$expression" -Oz -o "$tmp" "$input"
    "$TABIX" -f -p vcf "$tmp"
    mv -f "$tmp" "$output"
    mv -f "$tmp.tbi" "$output.tbi"
  fi
}
prepare_one() {
  local sample=$1 depth=$2 graph=$3 input=$4
  local outdir="$OUT_ROOT/$sample/$depth/$graph/$TOOL"
  local stem="$sample.$graph.$TOOL.PASS"
  local prepared="$outdir/$stem.nonref.split.multi0.completeGT.sorted.vcf.gz"
  local tmpdir="$outdir/.tmp.$$"
  need_file "$input"
  need_file "$input.tbi"
  mkdir -p "$outdir"

  if [[ "$FORCE" == 1 || ! -s "$prepared" || ! -s "$prepared.tbi" ]]; then
    rm -rf "$tmpdir"
    mkdir -p "$tmpdir/sort"
    local tmp="$outdir/.prepared.$$.vcf.gz"
    rm -f "$tmp" "$tmp.tbi"
    log "$TOOL $sample $depth $graph: PASS + complete non-reference GT + normalize + split + sort"
    "$BCF" view --threads "$THREADS" -f PASS -i 'GT="alt" && GT!="mis"' -Ou "$input" |
      "$BCF" norm --threads "$THREADS" -f "$REF" -c e -m -any --multi-overlaps 0 -Ou |
      "$BCF" view --threads "$THREADS" -i 'GT="alt" && GT!="mis"' -Ou |
      "$BCF" sort -T "$tmpdir/sort" -m "$SORT_MEM" -Oz -o "$tmp"
    "$TABIX" -f -p vcf "$tmp"
    mv -f "$tmp" "$prepared"
    mv -f "$tmp.tbi" "$prepared.tbi"
    rm -rf "$tmpdir"
  else
    log "$TOOL $sample $depth $graph: reuse prepared VCF"
  fi

  local dna='REF~"^[ACGTN][ACGTN]*$" && ALT~"^[ACGTN][ACGTN]*$"'
  local snp_expr="$dna && strlen(REF)=1 && strlen(ALT)=1 && REF!=ALT"
  local indel_expr="$dna && strlen(REF)!=strlen(ALT) && (strlen(REF)-strlen(ALT)<50) && (strlen(ALT)-strlen(REF)<50)"
  local sv_expr="$dna && ((strlen(REF)-strlen(ALT)>=50) || (strlen(ALT)-strlen(REF)>=50))"
  local class_union="$snp_expr || $indel_expr || $sv_expr"
  write_class "$prepared" "$outdir/$stem.snp.vcf.gz" "$snp_expr"
  write_class "$prepared" "$outdir/$stem.indel.vcf.gz" "$indel_expr"
  write_class "$prepared" "$outdir/$stem.sv.vcf.gz" "$sv_expr"
  write_complement "$prepared" "$outdir/$stem.other.vcf.gz" "$class_union"

  local all_n snp_n indel_n sv_n other_n
  all_n=$(count_records "$prepared")
  snp_n=$(count_records "$outdir/$stem.snp.vcf.gz")
  indel_n=$(count_records "$outdir/$stem.indel.vcf.gz")
  sv_n=$(count_records "$outdir/$stem.sv.vcf.gz")
  other_n=$(count_records "$outdir/$stem.other.vcf.gz")
  [[ $((snp_n + indel_n + sv_n + other_n)) -eq "$all_n" ]] || fail "classification count mismatch: $prepared"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$depth" "$graph" "$input" "$all_n" "$snp_n" "$indel_n" "$sv_n" "$other_n" >> "$SUMMARY_TMP"
  rm -rf "$tmpdir"
  log "$TOOL $sample $depth $graph: prepared=$all_n SNP=$snp_n small_indel=$indel_n SV=$sv_n other=$other_n"
}

mkdir -p "$OUT_ROOT"
SUMMARY="$OUT_ROOT/summary.tsv"
SUMMARY_TMP="$SUMMARY.tmp.$$"
printf 'sample\tdepth\tgraph\tinput_vcf\tprepared_records\tSNP\tsmall_indel\tSV\tother\n' > "$SUMMARY_TMP"
log "START $TOOL evaluation-input preparation"
log "source=$SOURCE_ROOT output=$OUT_ROOT reference=$REF threads=$THREADS sort_mem=$SORT_MEM"
found=0
for sample in "${SAMPLES[@]}"; do
  [[ -z "$ONLY_SAMPLE" || "$sample" == "$ONLY_SAMPLE" ]] || continue
  for depth in "${DEPTHS[@]}"; do
    [[ -z "$ONLY_DEPTH" || "$depth" == "$ONLY_DEPTH" ]] || continue
    for graph in "${GRAPHS[@]}"; do
      [[ -z "$ONLY_GRAPH" || "$graph" == "$ONLY_GRAPH" ]] || continue
      input=$(source_for "$sample" "$depth" "$graph")
      [[ -s "$input" ]] || fail "missing classified VCF: $input"
      [[ -s "$input.tbi" ]] || fail "missing classified VCF index: $input.tbi"
      prepare_one "$sample" "$depth" "$graph" "$input"
      found=$((found + 1))
    done
  done
done
[[ "$found" -gt 0 ]] || fail "no classified VCF was selected"
mv -f "$SUMMARY_TMP" "$SUMMARY"
log "SUCCESS: processed=$found summary=$SUMMARY"
