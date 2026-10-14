#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/data/malong/benchmark
SOURCE_ROOT=$ROOT/results_graphtyper2
OUT_ROOT=$ROOT/results_graphtyper3
BCF=/data/malong/miniconda3/envs/benchmark/bin/bcftools
TABIX=/data/malong/miniconda3/envs/benchmark/bin/tabix
THREADS=${EVAL_PREP_THREADS:-16}
SORT_MEM=${EVAL_PREP_SORT_MEM:-4G}

SAMPLES=(9311 IR64 W0169 W1735 W3033 zhonghua11)
DEPTHS=(5x 10x 15x 20x)
GRAPHS=(MC PGGB minigraph)

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*"; }
fail() { log "ERROR: $*" >&2; exit 1; }
trap 'rc=$?; log "FAILED line=$LINENO command=$BASH_COMMAND exit=$rc" >&2; exit "$rc"' ERR

[[ -x $BCF && -x $TABIX ]] || fail "bcftools/tabix not found"
[[ $THREADS =~ ^[1-9][0-9]*$ ]] || fail "invalid EVAL_PREP_THREADS=$THREADS"

need_vcf() {
  [[ -s $1 && -s $1.tbi ]] || fail "missing VCF or index: $1"
}

prepare_small() {
  local input=$1 output=$2 work=$3
  local tmp=$output.tmp.$$
  "$BCF" view --threads "$THREADS" -f PASS -i 'GT="alt" && GT!="mis"' -Oz -o "$tmp" "$input"
  "$TABIX" -f -p vcf "$tmp"
  mv -f "$tmp" "$output"
  mv -f "$tmp.tbi" "$output.tbi"
}

prepare_sv() {
  local input=$1 output=$2 work=$3
  local tmp=$output.tmp.$$
  # FORMAT/FT is the sample-level GraphTyper SV genotype filter.
  "$BCF" view --threads "$THREADS" -f PASS \
    -i 'INFO/SVMODEL="AGGREGATED" && GT="alt" && FORMAT/FT="PASS"' -Oz -o "$tmp" "$input"
  "$TABIX" -f -p vcf "$tmp"
  mv -f "$tmp" "$output"
  mv -f "$tmp.tbi" "$output.tbi"
}

mkdir -p "$OUT_ROOT"
printf 'sample\tdepth\tgraph\tSNP\tsmall_indel\tSV\tsv_policy\toutput_dir\n' > "$OUT_ROOT/summary.tsv.tmp"
cat > "$OUT_ROOT/evaluation_input_policy.txt" <<'EOF'
GraphTyper evaluation-input policy (corrected 2026-10-04)
Source: /data/malong/benchmark/results_graphtyper2
SNP/small indel: FILTER=PASS and GT=alt; no extra quality or allele-depth filter.
SV: FILTER=PASS, INFO/SVMODEL=AGGREGATED, GT=alt, and FORMAT/FT=PASS.
FORMAT/FT is the sample-level GraphTyper filter; FAIL1/FAIL2/FAIL3 are excluded.
The previous SV preparation kept failed sample genotypes and was not used here.
EOF

for sample in "${SAMPLES[@]}"; do
  for depth in "${DEPTHS[@]}"; do
    for graph in "${GRAPHS[@]}"; do
      src=$SOURCE_ROOT/$sample/$depth/$graph/graphtyper
      out=$OUT_ROOT/$sample/$depth/$graph/graphtyper
      work=$out/.tmp.$$
      mkdir -p "$out" "$work"
      snp_n=0
      indel_n=0
      if [[ $graph != minigraph ]]; then
        snp=$src/$sample.$graph.graphtyper.PASS.snp.vcf.gz
        indel=$src/$sample.$graph.graphtyper.PASS.indel.vcf.gz
        need_vcf "$snp"
        need_vcf "$indel"
        prepare_small "$snp" "$out/$sample.$graph.graphtyper.PASS.snp.vcf.gz" "$work"
        prepare_small "$indel" "$out/$sample.$graph.graphtyper.PASS.indel.vcf.gz" "$work"
        snp_n=$("$BCF" index -n "$out/$sample.$graph.graphtyper.PASS.snp.vcf.gz")
        indel_n=$("$BCF" index -n "$out/$sample.$graph.graphtyper.PASS.indel.vcf.gz")
      fi
      sv=$src/$sample.$graph.graphtyper.PASS.sv.vcf.gz
      [[ $graph != minigraph ]] || sv=$src/$sample.minigraph.graphtyper.PASS.sv.vcf.gz
      need_vcf "$sv"
      prepare_sv "$sv" "$out/$sample.$graph.graphtyper.PASS.sv.vcf.gz" "$work"
      sv_n=$("$BCF" index -n "$out/$sample.$graph.graphtyper.PASS.sv.vcf.gz")
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$depth" "$graph" "$snp_n" "$indel_n" "$sv_n" 'FILTER=PASS;SVMODEL=AGGREGATED;GT=alt;FORMAT/FT=PASS' "$out" >> "$OUT_ROOT/summary.tsv.tmp"
      rm -rf "$work"
      log "$sample $depth $graph: SNP=$snp_n small_indel=$indel_n SV_FT_PASS=$sv_n"
    done
  done
done
mv -f "$OUT_ROOT/summary.tsv.tmp" "$OUT_ROOT/summary.tsv"
log "SUCCESS: corrected GraphTyper evaluation inputs written to $OUT_ROOT"
