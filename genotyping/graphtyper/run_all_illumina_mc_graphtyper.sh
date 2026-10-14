#!/usr/bin/env bash
# Run MC GraphTyper genotype and genotype_sv for all Illumina data sets.

set -Eeuo pipefail
ulimit -c 0

BASE=/data/malong/benchmark
READS=$BASE/Illumina_depth
ROOT=$BASE/results2/MC/01graphtyper
WORKER=$BASE/code2/01graphtyper/run_9311_5x_mc_graphtyper.sh
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
[[ -x "$WORKER" ]] || die "missing worker: $WORKER"
[[ -x /usr/bin/time ]] || die "missing GNU time"
(( THREADS > 0 && MAX_CONCURRENT > 0 && INTERVAL > 0 )) || die "invalid numeric setting"
[[ "$FORCE_RERUN" == 0 || "$FORCE_RERUN" == 1 ]] || die "FORCE_RERUN must be 0 or 1"

mkdir -p "$ROOT"
if (( FORCE_RERUN == 1 )); then : > "$LOG"; fi
exec > >(tee -a "$LOG") 2>&1

declare -a TASKS=() RUNNING=()
declare -A SAMPLE=() DEPTH=() PID=() STATE=() START=() ELAPSED=() MSG=()

perf_value() {
  local task=$1 key=$2
  awk -F= -v k="$key" '$1 == k {print $2; exit}' "$ROOT/$task/performance.raw.txt" 2>/dev/null || true
}

write_perf_summary() {
  local task=$1 raw=$ROOT/$task/performance.raw.txt tmp=$ROOT/$task/performance.tsv.tmp.$$
  [[ -s "$raw" ]] || return 0
  {
    printf 'sample\tdepth\tmetric\tvalue\n'
    while IFS='=' read -r key value; do
      [[ -n "$key" ]] || continue
      printf '%s\t%s\t%s\t%s\n' "${SAMPLE[$task]}" "${DEPTH[$task]}" "$key" "$value"
    done < "$raw"
  } > "$tmp"
  mv -f "$tmp" "$ROOT/$task/performance.tsv"
}

write_status() {
  local tmp=$STATUS.tmp.$$ task out gt sv wall user system cpu rss input output now elapsed
  now=$(date +%s)
  printf 'task\tsample\tdepth\tstate\tpid\telapsed\tgt_fragments\tsv_fragments\twall_seconds\tuser_seconds\tsystem_seconds\tcpu_percent\tmax_rss_kb\tfs_inputs\tfs_outputs\toutput\tmessage\n' > "$tmp"
  for task in "${TASKS[@]}"; do
    out=$ROOT/$task
    gt=$( { find "$out/genotype" -type f -name '*.vcf.gz' -size +0c 2>/dev/null || true; } | wc -l)
    sv=$( { find "$out/genotype_sv" -type f -name '*.vcf.gz' -size +0c 2>/dev/null || true; } | wc -l)
    elapsed=${ELAPSED[$task]:--}
    if [[ ${STATE[$task]} == running ]]; then elapsed=$(fmt $((now - START[$task]))); fi
    wall=$(perf_value "$task" wall_seconds); user=$(perf_value "$task" user_seconds); system=$(perf_value "$task" system_seconds)
    cpu=$(perf_value "$task" cpu_percent); rss=$(perf_value "$task" max_rss_kb); input=$(perf_value "$task" fs_inputs); output=$(perf_value "$task" fs_outputs)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "${SAMPLE[$task]}" "${DEPTH[$task]}" "${STATE[$task]}" "${PID[$task]:--}" "$elapsed" "$gt" "$sv" "${wall:--}" "${user:--}" "${system:--}" "${cpu:--}" "${rss:--}" "${input:--}" "${output:--}" "$out" "${MSG[$task]:-}" >> "$tmp"
  done
  mv -f "$tmp" "$STATUS"
}

task_complete() {
  local task=$1 s=${SAMPLE[$1]} d=${DEPTH[$1]} out=$ROOT/$1
  [[ -s "$out/$s.$d.MC.graphtyper.snp_smallindel.vcf.gz" && -s "$out/$s.$d.MC.graphtyper.snp_smallindel.vcf.gz.tbi" && -s "$out/$s.$d.MC.graphtyper.sv.vcf.gz" && -s "$out/$s.$d.MC.graphtyper.sv.vcf.gz.tbi" && -s "$out/performance.tsv" && -s "$out/genotype/.graphtyper_genotype.complete" && -s "$out/genotype_sv/.graphtyper_genotype_sv.complete" ]]
}

start_task() {
  local task=$1 s=${SAMPLE[$1]} d=${DEPTH[$1]} out=$ROOT/$1
  mkdir -p "$out"
  log "START $task: MC genotype + genotype_sv, threads=$THREADS"
  (
    rc=0
    GRAPHTYPER_SAMPLE=$s GRAPHTYPER_DEPTH=$d GRAPHTYPER_THREADS=$THREADS \
      /usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu_percent=%P\nmax_rss_kb=%M\nfs_inputs=%I\nfs_outputs=%O\nvoluntary_context_switches=%w\ninvoluntary_context_switches=%c\nexit_status=%x' \
      -o "$out/performance.raw.txt" "$WORKER" > "$out/driver.log" 2>&1 || rc=$?
    if (( rc == 0 )); then write_perf_summary "$task"; fi
    if (( rc == 0 )) && task_complete "$task"; then touch "$out/.complete"; exit 0; fi
    touch "$out/.failed"; exit "$rc"
  ) &
  PID[$task]=$!; START[$task]=$(date +%s); STATE[$task]=running; MSG[$task]=running; RUNNING+=("$task")
}

poll_tasks() {
  local -a next=()
  local task rc now
  now=$(date +%s)
  for task in "${RUNNING[@]}"; do
    if kill -0 "${PID[$task]}" 2>/dev/null; then next+=("$task"); continue; fi
    rc=0; wait "${PID[$task]}" || rc=$?; ELAPSED[$task]=$(fmt $((now - START[$task])))
    if (( rc == 0 )) && task_complete "$task"; then STATE[$task]=completed; MSG[$task]=completed; log "DONE $task elapsed=${ELAPSED[$task]}"; else STATE[$task]=failed; MSG[$task]="worker exit=$rc; inspect $ROOT/$task/driver.log"; log "FAILED $task: ${MSG[$task]}"; fi
  done
  RUNNING=("${next[@]}")
}

on_signal() {
  local task
  log "termination signal; stopping ${#RUNNING[@]} workers"
  for task in "${RUNNING[@]}"; do kill "${PID[$task]}" 2>/dev/null || true; done
  wait || true; write_status; exit 143
}
trap on_signal INT TERM

while IFS= read -r r1; do
  dir=$(dirname "$r1"); s=$(basename "$(dirname "$dir")"); d=$(basename "$dir"); task=$s/$d
  [[ -s "$dir/${s}_${d}_R2.fastq.gz" ]] || die "missing R2 for $task"
  [[ -s "$BASE/share/graphtyper/bam/$s/$d/$s.$d.MC.sorted.bam" ]] || die "missing MC BAM for $task"
  [[ -s "$BASE/share/graphtyper/bam/$s/$d/bams.list" ]] || die "missing BAM list for $task"
  TASKS+=("$task"); SAMPLE[$task]=$s; DEPTH[$task]=$d
done < <(find "$READS" -mindepth 3 -maxdepth 3 -type f -name '*_R1.fastq.gz' -print | sort)
(( ${#TASKS[@]} == 24 )) || die "expected 24 tasks, found ${#TASKS[@]}"

for task in "${TASKS[@]}"; do
  if (( FORCE_RERUN == 1 )); then
    rm -rf "$ROOT/$task"; STATE[$task]=queued; MSG[$task]=forced_rerun
  elif task_complete "$task"; then
    STATE[$task]=completed; MSG[$task]=reused
  else
    STATE[$task]=queued; MSG[$task]=waiting
  fi
done

log "MC controller: tasks=${#TASKS[@]} max_concurrent=$MAX_CONCURRENT threads_per_task=$THREADS force_rerun=$FORCE_RERUN"
write_status

next_index=0
last_status=0
while (( next_index < ${#TASKS[@]} || ${#RUNNING[@]} > 0 )); do
  poll_tasks
  while (( ${#RUNNING[@]} < MAX_CONCURRENT && next_index < ${#TASKS[@]} )); do
    task=${TASKS[$next_index]}; next_index=$((next_index + 1))
    [[ ${STATE[$task]} == queued ]] || continue
    start_task "$task"
  done
  now=$(date +%s)
  if (( now - last_status >= INTERVAL )); then
    write_status
    completed=0; failed=0
    for task in "${TASKS[@]}"; do
      [[ ${STATE[$task]} == completed ]] && completed=$((completed + 1))
      [[ ${STATE[$task]} == failed ]] && failed=$((failed + 1))
    done
    log "STATUS completed=$completed/${#TASKS[@]} running=${#RUNNING[@]} next=$next_index failed=$failed"
    last_status=$now
  fi
  (( next_index < ${#TASKS[@]} || ${#RUNNING[@]} > 0 )) && sleep 10
done
write_status

end_time=$(date +%s)
printf 'metric\tvalue\ncontroller_wall_seconds\t%s\ntasks_total\t%s\nmax_concurrent\t%s\nthreads_per_task\t%s\n' "$((end_time - START_TIME))" "${#TASKS[@]}" "$MAX_CONCURRENT" "$THREADS" > "$ROOT/run_all.performance.tsv"
failed=0
for task in "${TASKS[@]}"; do
  [[ ${STATE[$task]} == failed ]] && failed=$((failed + 1))
done
(( failed == 0 )) || { log "FINISHED WITH FAILURES: $failed"; exit 1; }
log "COMPLETED all ${#TASKS[@]} MC tasks"
