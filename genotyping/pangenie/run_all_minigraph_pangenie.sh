#!/usr/bin/env bash
set -Eeuo pipefail
ulimit -c 0

ROOT=/data/malong/benchmark
BIN=/data/malong/miniconda3/envs/benchmark/bin
WORK=$ROOT/graph/minigraph/pangenie_work
READROOT=$ROOT/Illumina_depth/pangenie_fastq
OUTROOT=$ROOT/results2/minigraph/02pangenie
INDEX_DIR=$WORK/index_all
PANEL_DIR=$WORK/vcf/by_chromosome
THREADS=${MINIGRAPH_PANGENIE_THREADS:-128}
MIN_SV_LEN=${MINIGRAPH_PANGENIE_SV_MIN_SIZE:-50}
SKIP_9311_5X=${MINIGRAPH_PANGENIE_SKIP_9311_5X:-1}
SUMMARY=$OUTROOT/run_all_summary.tsv

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
need_file() { [[ -s "$1" ]] || die "missing or empty file: $1"; }
need_bin() { [[ -x "$BIN/$1" ]] || die "missing executable: $BIN/$1"; }

for x in PanGenie bcftools bgzip tabix; do need_bin "$x"; done
[[ -x /usr/bin/time ]] || die 'missing /usr/bin/time'
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "invalid MINIGRAPH_PANGENIE_THREADS: $THREADS"
[[ "$MIN_SV_LEN" =~ ^[1-9][0-9]*$ ]] || die "invalid MINIGRAPH_PANGENIE_SV_MIN_SIZE: $MIN_SV_LEN"
mkdir -p "$OUTROOT"

for c in $(seq 1 12); do
  chr=Chr$c
  prefix=$INDEX_DIR/$chr.pangenie
  need_file "$PANEL_DIR/$chr.vcf"
  need_file "${prefix}_UniqueKmersMap.cereal"
  need_file "${prefix}_path_segments.fasta"
  need_file "${prefix}_${chr}_Graph.cereal"
done

printf 'sample\tdepth\tstatus\tsv_records\telapsed_seconds\tcompleted\n' > "$SUMMARY"

run_one() {
  local sample=$1 depth=$2 reads=$3
  local out=$OUTROOT/$sample/$depth
  local log_dir=$out/logs chr_out=$out/chromosomes perf=$out/performance.tsv
  local all=$out/$sample.$depth.minigraph.pangenie.sv.vcf.gz
  local snp=$out/$sample.$depth.minigraph.pangenie.snp.vcf.gz
  local indel=$out/$sample.$depth.minigraph.pangenie.indel.vcf.gz
  local sv=$out/$sample.$depth.minigraph.pangenie.sv.vcf.gz
  local done=$out/DONE list=$out/chromosomes.list

  need_file "$reads"
  mkdir -p "$out" "$log_dir" "$chr_out"
  if [[ -s "$done" && -s "$sv" && -s "$sv.tbi" ]]; then
    log "reusing completed PanGenie result: $sample $depth"
    printf '%s\t%s\treused\t%s\t%s\t%s\n' "$sample" "$depth" "$("$BIN/bcftools" index -n "$sv")" "$(awk -F= '$1==\"elapsed_seconds\"{print $2}' "$done" 2>/dev/null || printf 0)" "$(awk -F= '$1==\"completed\"{print $2}' "$done" 2>/dev/null || date -Is)" >> "$SUMMARY"
    return
  fi

  printf 'chromosome\tstatus\telapsed_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb\tpage_faults\toutput_records\n' > "$perf"
  : > "$list"
  local start_epoch=$(date +%s)
  for c in $(seq 1 12); do
    local chr=Chr$c
    local prefix=$chr_out/$sample.$depth.minigraph.$chr
    local raw=${prefix}_genotyping.vcf
    local compressed=$raw.gz
    local timefile=$log_dir/$chr.time.txt
    local panlog=$log_dir/$chr.pangenie.log
    if [[ ! -s "$raw" ]]; then
      rm -f "$raw" "$compressed" "$compressed.tbi"
      log "running PanGenie minigraph $sample $depth $chr with $THREADS k-mer threads"
      /usr/bin/time -f 'elapsed_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kb=%M\npage_faults=%F' -o "$timefile" "$BIN/PanGenie" -f "$INDEX_DIR/$chr.pangenie" -i "$reads" -o "$prefix" -s "$sample" -j "$THREADS" -t 1 -u > "$panlog" 2>&1
      need_file "$raw"
    else
      log "reusing PanGenie output $sample $depth $chr"
      if [[ ! -s "$timefile" ]]; then printf 'elapsed_seconds=0\nuser_seconds=0\nsystem_seconds=0\nmax_rss_kb=0\npage_faults=0\n' > "$timefile"; fi
    fi
    if [[ ! -s "$compressed" || ! -s "$compressed.tbi" ]]; then
      local tmp=$compressed.tmp.$$
      "$BIN/bgzip" -@ "$THREADS" -c "$raw" > "$tmp"
      "$BIN/tabix" -f -p vcf "$tmp"
      mv -f "$tmp" "$compressed"
      mv -f "$tmp.tbi" "$compressed.tbi"
    fi
    local n=$("$BIN/bcftools" index -n "$compressed")
    local elapsed=$(awk -F= '$1==\"elapsed_seconds\"{print $2}' "$timefile")
    local user=$(awk -F= '$1==\"user_seconds\"{print $2}' "$timefile")
    local system=$(awk -F= '$1==\"system_seconds\"{print $2}' "$timefile")
    local rss=$(awk -F= '$1==\"max_rss_kb\"{print $2}' "$timefile")
    local faults=$(awk -F= '$1==\"page_faults\"{print $2}' "$timefile")
    printf '%s\tcompleted\t%s\t%s\t%s\t%s\t%s\t%s\n' "$chr" "$elapsed" "$user" "$system" "$rss" "$faults" "$n" >> "$perf"
    printf '%s\n' "$compressed" >> "$list"
    log "$sample $depth $chr done: $n records, elapsed ${elapsed}s, RSS ${rss}KB"
  done

  if [[ ! -s "$all" || ! -s "$all.tbi" ]]; then
    local tmp=$all.tmp.$$
    "$BIN/bcftools" concat -a -Oz -o "$tmp" -f "$list"
    "$BIN/tabix" -f -p vcf "$tmp"
    mv -f "$tmp" "$all"
    mv -f "$tmp.tbi" "$all.tbi"
  fi
  if [[ "$all" != "$sv" ]]; then
    mv -f "$all" "$sv"
    mv -f "$all.tbi" "$sv.tbi"
  fi

  for f in "$sv" "$sv.tbi"; do need_file "$f"; done
  local total_elapsed=$(( $(date +%s) - start_epoch ))
  cat > "$done" <<EOF
sample=$sample
depth=$depth
graph=minigraph
variant_scope=SV_only
panel=$WORK/vcf/rice36_minigraph.pangenie.panel.vcf
reads=$reads
threads=$THREADS
sv_records=$("$BIN/bcftools" index -n "$sv")
elapsed_seconds=$total_elapsed
completed=$(date -Is)
EOF
  printf '%s\t%s\tcompleted\t%s\t%s\t%s\n' "$sample" "$depth" "$("$BIN/bcftools" index -n "$sv")" "$total_elapsed" "$(date -Is)" >> "$SUMMARY"
  log "completed PanGenie: $sample $depth in ${total_elapsed}s"
}

shopt -s nullglob
reads_list=("$READROOT"/*.fastq)
(( ${#reads_list[@]} > 0 )) || die "no uncompressed FASTQ files found in $READROOT"
for reads in "${reads_list[@]}"; do
  base=$(basename "$reads" .fastq)
  sample=${base%_*}
  depth=${base##*_}
  [[ "$sample" != "$base" && "$depth" =~ ^(5|10|15|20)x$ ]] || die "unexpected FASTQ basename: $base"
  if [[ "$sample" == 9311 && "$depth" == 5x && "$SKIP_9311_5X" == 1 ]]; then
    log "skipping already completed dataset: $sample $depth"
    printf '%s\t%s\tskipped_existing\t%s\t%s\t%s\n' "$sample" "$depth" "$("$BIN/bcftools" index -n "$OUTROOT/9311/5x/9311.5x.minigraph.pangenie.sv.vcf.gz" 2>/dev/null || printf 0)" "$(awk -F= '$1==\"elapsed_seconds\"{print $2}' "$OUTROOT/9311/5x/DONE" 2>/dev/null || printf 0)" "$(awk -F= '$1==\"completed\"{print $2}' "$OUTROOT/9311/5x/DONE" 2>/dev/null || date -Is)" >> "$SUMMARY"
    continue
  fi
  run_one "$sample" "$depth" "$reads"
done
log "all minigraph PanGenie sample-depth datasets completed or reused"
cat "$SUMMARY"
