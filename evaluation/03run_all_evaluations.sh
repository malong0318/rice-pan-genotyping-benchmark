#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=${ROOT:-/data/malong/benchmark}
TRUTH_ROOT=$ROOT/truth
REF=$ROOT/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa
RTG=$ROOT/../miniconda3/envs/benchmark/bin/rtg
BCF=$ROOT/../miniconda3/envs/benchmark/bin/bcftools
TABIX=$ROOT/../miniconda3/bin/tabix
PYTHON=$ROOT/../miniconda3/envs/benchmark/bin/python
SDF=$ROOT/evaluation/NIP-T2T_Chr.sdf
OUT_ROOT=${EVAL_OUT_ROOT:-$ROOT/evaluation/paper}
THREADS=${THREADS:-64}
RTG_MEM=${RTG_MEM:-8G}
FORCE=${EVAL_FORCE:-0}
ONLY_TOOL=${EVAL_ONLY_TOOL:-}
ONLY_KIND=${EVAL_ONLY_KIND:-}
ONLY_SAMPLE=${EVAL_ONLY_SAMPLE:-}
ONLY_DEPTH=${EVAL_ONLY_DEPTH:-}
ONLY_GRAPH=${EVAL_ONLY_GRAPH:-}
LOG=$OUT_ROOT/run_all_evaluations.log

SAMPLES=(9311 IR64 W0169 W1735 W3033 zhonghua11)
DEPTHS=(5x 10x 15x 20x)
TOOLS=(vg pangenie bayestyper graphtyper)

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*" | tee -a "$LOG"; }
fail() { log "ERROR: $*" >&2; exit 1; }
trap 'rc=$?; log "FAILED line=$LINENO command=$BASH_COMMAND exit=$rc"; exit "$rc"' ERR
count() { "$BCF" index -n "$1"; }
selected() { [[ -z "$1" || "$1" == "$2" ]]; }

mkdir -p "$OUT_ROOT"
exec > >(tee -a "$LOG") 2>&1
for exe in "$RTG" "$BCF" "$TABIX" "$PYTHON"; do [[ -x "$exe" ]] || fail "missing executable: $exe"; done
[[ -s "$REF" && -s "$REF.fai" ]] || fail "missing reference or index: $REF"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || fail "invalid THREADS=$THREADS"
"$PYTHON" -c 'import truvari' >/dev/null 2>&1 || fail "truvari is not importable"

if [[ ! -d "$SDF" ]]; then
  log "building RTG SDF: $SDF"
  RTG_MEM="$RTG_MEM" "$RTG" format -o "$SDF" "$REF"
fi

cat > "$OUT_ROOT/software_versions.txt" <<EOF
bcftools: $($BCF --version | head -1)
samtools: $(/data/malong/miniconda3/envs/benchmark/bin/samtools --version | head -1)
rtg-tools: 3.13 (conda package)
truvari: $($PYTHON -c 'import truvari; print(getattr(truvari, "__version__", "unknown"))')
python: $($PYTHON --version 2>&1)
reference: $REF
graphTyper_SV_policy: ordinary symbolic-ALT results_graphtyper3; no sequence restoration
EOF

call_root() {
  case "$1" in
    vg) printf '%s\n' "$ROOT/results_vg2" ;;
    pangenie) printf '%s\n' "$ROOT/results_pangenie2" ;;
    bayestyper) printf '%s\n' "$ROOT/results_bayestyper2" ;;
    graphtyper) printf '%s\n' "$ROOT/results_graphtyper3" ;;
    *) fail "unknown tool: $1" ;;
  esac
}

small_graphs() { if [[ "$1" == vg ]]; then printf 'MC\n'; else printf 'MC\nPGGB\n'; fi; }
sv_graphs() { if [[ "$1" == vg ]]; then printf 'MC\nminigraph\n'; else printf 'MC\nPGGB\nminigraph\n'; fi; }

call_file() {
  local tool=$1 sample=$2 depth=$3 graph=$4 kind=$5 root
  root=$(call_root "$tool")
  printf '%s\n' "$root/$sample/$depth/$graph/$tool/$sample.$graph.$tool.PASS.$kind.vcf.gz"
}

prepare_rtg_call() {
  local tool=$1 sample=$2 depth=$3 graph=$4 kind=$5 input=$6
  local out="$OUT_ROOT/$kind/$tool/_work/$sample.$depth.$graph.$tool.PASS.$kind.rtg.vcf.gz"
  mkdir -p "$(dirname "$out")"
  if [[ "$FORCE" == 1 || ! -s "$out" || ! -s "$out.tbi" ]]; then
    local tmp="$out.tmp.$$"
    rm -f "$tmp" "$tmp.tbi"
    log "prepare RTG input: tool=$tool kind=$kind sample=$sample depth=$depth graph=$graph" >&2
    "$BCF" annotate -x '^FORMAT/GT' -Oz -o "$tmp" "$input"
    "$TABIX" -f -p vcf "$tmp"
    mv -f "$tmp" "$out"; mv -f "$tmp.tbi" "$out.tbi"
  fi
  printf '%s\n' "$out"
}

run_small() {
  local tool=$1 kind=$2 outroot summary
  outroot="$OUT_ROOT/$kind/$tool"
  summary="$outroot/summary.tsv"
  mkdir -p "$outroot"
  printf 'sample\tdepth\tgraph\ttruth_records\tcall_records\tresult\n' > "$summary.tmp"
  while read -r graph; do
    [[ -n "$graph" ]] || continue
    selected "$ONLY_GRAPH" "$graph" || continue
    for sample in "${SAMPLES[@]}"; do
      selected "$ONLY_SAMPLE" "$sample" || continue
      truth="$TRUTH_ROOT/$sample/Chr_vcf/$sample.snp.confident.final.vcf.gz"
      [[ "$kind" == indel ]] && truth="$TRUTH_ROOT/$sample/Chr_vcf/$sample.indel.small.confident.final.vcf.gz"
      bed="$TRUTH_ROOT/$sample/Chr_vcf/$sample.confident.Chr.bed"
      [[ -s "$truth" && -s "$truth.tbi" && -s "$bed" ]] || fail "missing truth/BED: $sample $kind"
      for depth in "${DEPTHS[@]}"; do
        selected "$ONLY_DEPTH" "$depth" || continue
        input=$(call_file "$tool" "$sample" "$depth" "$graph" "$kind")
        [[ -s "$input" && -s "$input.tbi" ]] || fail "missing call: $input"
        call=$(prepare_rtg_call "$tool" "$sample" "$depth" "$graph" "$kind" "$input")
        out="$outroot/$graph/$sample/$depth"
        if [[ "$FORCE" == 1 || ! -s "$out/summary.txt" ]]; then
          if [[ "$FORCE" == 1 && -e "$out" ]]; then mv -f "$out" "$out.previous.$(date +%Y%m%dT%H%M%S)"; fi
          runout="$out.running.$$"
          mkdir -p "$(dirname "$out")"
          log "RTG $kind: tool=$tool sample=$sample depth=$depth graph=$graph"
          RTG_MEM="$RTG_MEM" "$RTG" vcfeval -t "$SDF" -b "$truth" -c "$call" \
            -e "$bed" --sample "syndip,$sample" --sample-ploidy 2 --decompose \
            -o "$runout" -T "$THREADS"
          [[ -s "$runout/summary.txt" ]] || fail "RTG output missing: $runout"
          mv -f "$runout" "$out"
        else
          log "reuse RTG result: $out"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$depth" "$graph" \
          "$(count "$truth")" "$(count "$call")" "$out" >> "$summary.tmp"
      done
    done
  done < <(small_graphs "$tool")
  mv -f "$summary.tmp" "$summary"
  log "completed small-variant evaluation: tool=$tool kind=$kind summary=$summary"
}

prepare_sv_truth() {
  local tool=$1 sample=$2 outroot="$OUT_ROOT/sv/$tool"
  local out="$outroot/_work/$sample.truth.PASS.noBND.vcf.gz"
  local input="$TRUTH_ROOT/$sample/Chr_vcf/$sample.sv.final.vcf.gz"
  mkdir -p "$(dirname "$out")"
  [[ -s "$input" && -s "$input.tbi" ]] || fail "missing SV truth: $input"
  if [[ "$FORCE" == 1 || ! -s "$out" || ! -s "$out.tbi" ]]; then
    local tmp="$out.tmp.$$"
    "$BCF" view -f PASS -e 'INFO/SVTYPE="BND"' -Oz -o "$tmp" "$input"
    "$TABIX" -f -p vcf "$tmp"
    mv -f "$tmp" "$out"; mv -f "$tmp.tbi" "$out.tbi"
  fi
  printf '%s\n' "$out"
}

prepare_truvari_call() {
  local tool=$1 sample=$2 depth=$3 graph=$4 input=$5
  local out="$OUT_ROOT/sv/$tool/_work/$sample.$depth.$graph.$tool.PASS.sv.truvari.vcf.gz"
  mkdir -p "$(dirname "$out")"
  if [[ "$FORCE" == 1 || ! -s "$out" || ! -s "$out.tbi" ]]; then
    local tmp="$out.tmp.$$"
    "$BCF" annotate -x '^FORMAT/GT' -Oz -o "$tmp" "$input"
    "$TABIX" -f -p vcf "$tmp"
    mv -f "$tmp" "$out"; mv -f "$tmp.tbi" "$out.tbi"
  fi
  printf '%s\n' "$out"
}

run_sv() {
  local tool=$1 outroot="$OUT_ROOT/sv/$tool" summary
  summary="$outroot/summary.tsv"
  mkdir -p "$outroot"
  printf 'sample\tdepth\tgraph\ttruth_records_excluding_BND\tcall_records\tresult\n' > "$summary.tmp"
  for sample in "${SAMPLES[@]}"; do
    selected "$ONLY_SAMPLE" "$sample" || continue
    truth=$(prepare_sv_truth "$tool" "$sample")
    bed="$TRUTH_ROOT/$sample/Chr_vcf/$sample.confident.Chr.bed"
    [[ -s "$bed" ]] || fail "missing SV confidence BED: $bed"
    while read -r graph; do
      [[ -n "$graph" ]] || continue
      selected "$ONLY_GRAPH" "$graph" || continue
      for depth in "${DEPTHS[@]}"; do
        selected "$ONLY_DEPTH" "$depth" || continue
        input=$(call_file "$tool" "$sample" "$depth" "$graph" sv)
        [[ -s "$input" && -s "$input.tbi" ]] || fail "missing SV call: $input"
        call=$(prepare_truvari_call "$tool" "$sample" "$depth" "$graph" "$input")
        out="$outroot/$graph/$sample/$depth"
        if [[ "$FORCE" == 1 || ! -s "$out/summary.json" ]]; then
          if [[ "$FORCE" == 1 && -e "$out" ]]; then mv -f "$out" "$out.previous.$(date +%Y%m%dT%H%M%S)"; fi
          runout="$out.running.$$"
          mkdir -p "$(dirname "$out")"
          log "Truvari SV: tool=$tool sample=$sample depth=$depth graph=$graph"
          "$PYTHON" -m truvari bench -b "$truth" -c "$call" -f "$REF" -o "$runout" \
            --passonly --bSample "$sample" --cSample "$sample" --includebed "$bed" \
            --sizemin 50 --sizefilt 50 --sizemax -1 --refdist 500 --pctseq 0.7 \
            --pctsize 0.7 --bnddist 100 --no-ref c
          [[ -s "$runout/summary.json" ]] || fail "Truvari output missing: $runout"
          mv -f "$runout" "$out"
        else
          log "reuse Truvari result: $out"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$depth" "$graph" \
          "$(count "$truth")" "$(count "$call")" "$out" >> "$summary.tmp"
      done
    done < <(sv_graphs "$tool")
  done
  mv -f "$summary.tmp" "$summary"
  log "completed SV evaluation: tool=$tool summary=$summary"
}

for tool in "${TOOLS[@]}"; do
  selected "$ONLY_TOOL" "$tool" || continue
  if [[ -z "$ONLY_KIND" || "$ONLY_KIND" == snp ]]; then run_small "$tool" snp; fi
  if [[ -z "$ONLY_KIND" || "$ONLY_KIND" == indel ]]; then run_small "$tool" indel; fi
  if [[ -z "$ONLY_KIND" || "$ONLY_KIND" == sv ]]; then run_sv "$tool"; fi
done

log "SUCCESS: all requested evaluations are complete under $OUT_ROOT"
