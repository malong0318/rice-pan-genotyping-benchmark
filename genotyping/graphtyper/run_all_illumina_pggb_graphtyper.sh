#!/usr/bin/env bash
# PGGB GraphTyper controller: all 24 Illumina sets, 128 threads/task, 3 concurrent.
set -Eeuo pipefail
ulimit -c 0
BASE=/data/malong/benchmark
READS=$BASE/Illumina_depth
ROOT=$BASE/results2/PGGB/01graphtyper
WORKER=$BASE/code2/01graphtyper/run_9311_5x_pggb_graphtyper.sh
BAMROOT=$BASE/share/graphtyper/bam
LOG=$ROOT/run_all.log
STATUS=$ROOT/run_all.status.tsv
THREADS=${GRAPHTYPER_THREADS:-128}
MAX_CONCURRENT=${GRAPHTYPER_MAX_CONCURRENT:-3}
INTERVAL=${GRAPHTYPER_STATUS_INTERVAL:-60}
FORCE_RERUN=${GRAPHTYPER_FORCE_RERUN:-1}
START_TIME=$(date +%s)

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }
fmt() { local s=$1; printf '%02dh:%02dm:%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60)); }
frags() { { find "$1" -type f -name '*.vcf.gz' -size +0c 2>/dev/null || true; } | wc -l; }
perf() { awk -F= -v k="$2" '$1 == k {print $2; exit}' "$ROOT/$1/performance.raw.txt" 2>/dev/null || true; }

[[ -x "$WORKER" ]] || die "missing worker: $WORKER"
[[ -x /usr/bin/time ]] || die "missing GNU time"
[[ -d "$READS" ]] || die "missing reads: $READS"
(( THREADS > 0 && MAX_CONCURRENT > 0 && INTERVAL > 0 )) || die "invalid numeric setting"
[[ "$FORCE_RERUN" == 0 || "$FORCE_RERUN" == 1 ]] || die "FORCE_RERUN must be 0 or 1"
mkdir -p "$ROOT"
if (( FORCE_RERUN == 1 )); then : > "$LOG"; fi
exec > >(tee -a "$LOG") 2>&1

declare -a TASKS=() RUNNING=()
declare -A SAMPLE=() DEPTH=() PID=() STATE=() START=() ELAPSED=() MSG=()

complete() {
  local t=$1 s=${SAMPLE[$1]} d=${DEPTH[$1]} o=$ROOT/$1
  [[ -s "$o/$s.$d.PGGB.graphtyper.snp_smallindel.vcf.gz" &&
     -s "$o/$s.$d.PGGB.graphtyper.snp_smallindel.vcf.gz.tbi" &&
     -s "$o/$s.$d.PGGB.graphtyper.sv.raw.vcf.gz" &&
     -s "$o/$s.$d.PGGB.graphtyper.sv.raw.vcf.gz.tbi" &&
     -s "$o/$s.$d.PGGB.graphtyper.sv.vcf.gz" &&
     -s "$o/$s.$d.PGGB.graphtyper.sv.vcf.gz.tbi" &&
     -s "$o/performance.tsv" &&
     -e "$o/genotype/.graphtyper_genotype.complete" &&
     -e "$o/genotype_sv/.graphtyper_genotype_sv.complete" ]]
}
perf_summary() {
  local t=$1 r=$ROOT/$1/performance.raw.txt x=$ROOT/$1/performance.tsv.tmp.$$
  [[ -s "$r" ]] || return 0
  { printf 'sample\tdepth\tmetric\tvalue\n'
    while IFS='=' read -r k v; do
      [[ -n "$k" ]] || continue
      printf '%s\t%s\t%s\t%s\n' "${SAMPLE[$t]}" "${DEPTH[$t]}" "$k" "$v"
    done < "$r"
  } > "$x"
  mv -f "$x" "$ROOT/$1/performance.tsv"
}
status() {
  local x=$STATUS.tmp.$$ t o now wall user sys cpu rss fi fo e gt sv
  now=$(date +%s)
  printf 'task\tsample\tdepth\tstate\tpid\telapsed\tgt_fragments\tsv_fragments\twall_seconds\tuser_seconds\tsystem_seconds\tcpu_percent\tmax_rss_kb\tfs_inputs\tfs_outputs\toutput\tmessage\n' > "$x"
  for t in "${TASKS[@]}"; do
    o=$ROOT/$t; gt=$(frags "$o/genotype"); sv=$(frags "$o/genotype_sv")
    e=${ELAPSED[$t]:--}
    [[ ${STATE[$t]} == running ]] && e=$(fmt $((now-START[$t])))
    wall=$(perf "$t" wall_seconds); user=$(perf "$t" user_seconds); sys=$(perf "$t" system_seconds)
    cpu=$(perf "$t" cpu_percent); rss=$(perf "$t" max_rss_kb); fi=$(perf "$t" fs_inputs); fo=$(perf "$t" fs_outputs)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$t" "${SAMPLE[$t]}" "${DEPTH[$t]}" "${STATE[$t]}" "${PID[$t]:--}" "$e" "$gt" "$sv" \
      "${wall:--}" "${user:--}" "${sys:--}" "${cpu:--}" "${rss:--}" "${fi:--}" "${fo:--}" "$o" "${MSG[$t]:-}" >> "$x"
  done
  mv -f "$x" "$STATUS"
}
start_task() {
  local t=$1 s=${SAMPLE[$1]} d=${DEPTH[$1]} o=$ROOT/$1
  mkdir -p "$o"
  log "START $t: PGGB genotype + genotype_sv, threads=$THREADS"
  (
    local rc=0
    GRAPHTYPER_SAMPLE=$s GRAPHTYPER_DEPTH=$d GRAPHTYPER_THREADS=$THREADS \
      /usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu_percent=%P\nmax_rss_kb=%M\nfs_inputs=%I\nfs_outputs=%O\nvoluntary_context_switches=%w\ninvoluntary_context_switches=%c\nexit_status=%x' \
      -o "$o/performance.raw.txt" "$WORKER" > "$o/driver.log" 2>&1 || rc=$?
    (( rc == 0 )) && perf_summary "$t"
    if (( rc == 0 )) && complete "$t"; then touch "$o/.complete"; exit 0; fi
    touch "$o/.failed"; exit "$rc"
  ) &
  PID[$t]=$!; START[$t]=$(date +%s); STATE[$t]=running; MSG[$t]=running; RUNNING+=("$t")
}
poll() {
  local -a n=()
  local t rc now=$(date +%s)
  for t in "${RUNNING[@]}"; do
    if kill -0 "${PID[$t]}" 2>/dev/null; then n+=("$t"); continue; fi
    rc=0; wait "${PID[$t]}" || rc=$?; ELAPSED[$t]=$(fmt $((now-START[$t])))
    if (( rc == 0 )) && complete "$t"; then
      STATE[$t]=completed; MSG[$t]=completed; log "DONE $t elapsed=${ELAPSED[$t]}"
    else
      STATE[$t]=failed; MSG[$t]="worker exit=$rc; inspect $ROOT/$t/driver.log"; log "FAILED $t: ${MSG[$t]}"
    fi
  done
  RUNNING=("${n[@]}")
}
stop_all() {
  local t
  log "termination signal; stopping ${#RUNNING[@]} workers"
  for t in "${RUNNING[@]}"; do kill "${PID[$t]}" 2>/dev/null || true; done
  wait || true; status; exit 143
}
trap stop_all INT TERM

mapfile -t R1S < <(find "$READS" -mindepth 3 -maxdepth 3 -type f -name '*_R1.fastq.gz' -print | sort)
for r1 in "${R1S[@]}"; do
  dir=$(dirname "$r1"); s=$(basename "$(dirname "$dir")"); d=$(basename "$dir"); t=$s/$d
  [[ -s "$dir/${s}_${d}_R2.fastq.gz" ]] || die "missing R2 for $t"
  [[ -s "$BAMROOT/$s/$d/$s.$d.MC.sorted.bam" ]] || die "missing BAM for $t"
  [[ -s "$BAMROOT/$s/$d/$s.$d.MC.sorted.bam.bai" ]] || die "missing BAM index for $t"
  [[ -s "$BAMROOT/$s/$d/bams.list" ]] || die "missing bams.list for $t"
  TASKS+=("$t"); SAMPLE[$t]=$s; DEPTH[$t]=$d
done
(( ${#TASKS[@]} == 24 )) || die "expected 24 tasks, found ${#TASKS[@]}"

for t in "${TASKS[@]}"; do
  if (( FORCE_RERUN == 1 )); then rm -rf "$ROOT/$t"; STATE[$t]=queued; MSG[$t]=forced_rerun
  elif complete "$t"; then STATE[$t]=completed; MSG[$t]=reused
  else STATE[$t]=queued; MSG[$t]=waiting
  fi
done
log "PGGB controller: tasks=${#TASKS[@]} max_concurrent=$MAX_CONCURRENT threads_per_task=$THREADS force_rerun=$FORCE_RERUN"
log "outputs: $ROOT"
status

next=0; last=0
while (( next < ${#TASKS[@]} || ${#RUNNING[@]} > 0 )); do
  poll
  while (( ${#RUNNING[@]} < MAX_CONCURRENT && next < ${#TASKS[@]} )); do
    t=${TASKS[$next]}; next=$((next+1)); [[ ${STATE[$t]} == queued ]] || continue; start_task "$t"
  done
  now=$(date +%s)
  if (( now-last >= INTERVAL )); then
    status; done_n=0; failed_n=0
    for t in "${TASKS[@]}"; do
      [[ ${STATE[$t]} == completed ]] && done_n=$((done_n+1))
      [[ ${STATE[$t]} == failed ]] && failed_n=$((failed_n+1))
    done
    log "STATUS completed=$done_n/${#TASKS[@]} running=${#RUNNING[@]} next=$next failed=$failed_n"; last=$now
  fi
  (( next < ${#TASKS[@]} || ${#RUNNING[@]} > 0 )) && sleep 10
done

status
end=$(date +%s)
printf 'metric\tvalue\ncontroller_wall_seconds\t%s\ntasks_total\t%s\nmax_concurrent\t%s\nthreads_per_task\t%s\n' \
  "$((end-START_TIME))" "${#TASKS[@]}" "$MAX_CONCURRENT" "$THREADS" > "$ROOT/run_all.performance.tsv"
failed_n=0
for t in "${TASKS[@]}"; do [[ ${STATE[$t]} == failed ]] && failed_n=$((failed_n+1)); done
if (( failed_n > 0 )); then log "FINISHED WITH FAILURES: $failed_n"; exit 1; fi
log "COMPLETED all ${#TASKS[@]} PGGB tasks"
