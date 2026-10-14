#!/usr/bin/env bash
# Run minigraph GraphTyper SV genotyping for every Illumina sample/depth pair.
# At most MAX_CONCURRENT workers run at once; each worker uses THREADS threads.

set -Eeuo pipefail
ulimit -c 0

BASE=/data/malong/benchmark
READS=$BASE/Illumina_depth
RESULT_ROOT=$BASE/results2/minigraph/01graphtyper
WORKER=$BASE/code2/01graphtyper/run_9311_5x_minigraph_graphtyper.sh
MASTER_LOG=$RESULT_ROOT/run_all.log
STATUS=$RESULT_ROOT/run_all.status.tsv
THREADS=${GRAPHTYPER_THREADS:-128}
MAX_CONCURRENT=${GRAPHTYPER_MAX_CONCURRENT:-5}
STATUS_INTERVAL=${GRAPHTYPER_STATUS_INTERVAL:-60}
DRY_RUN=${GRAPHTYPER_DRY_RUN:-0}
FORCE_RERUN=${GRAPHTYPER_FORCE_RERUN:-0}
RUN_START=$(date +%s)

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }
fmt_time() {
  local seconds=$1
  printf '%02dh:%02dm:%02ds' $((seconds / 3600)) $(((seconds % 3600) / 60)) $((seconds % 60))
}

mkdir -p "$RESULT_ROOT"
exec > >(tee -a "$MASTER_LOG") 2>&1

[[ -x "$WORKER" ]] || die "missing executable worker script: $WORKER"
(( THREADS > 0 )) || die "GRAPHTYPER_THREADS must be positive"
(( MAX_CONCURRENT > 0 )) || die "GRAPHTYPER_MAX_CONCURRENT must be positive"
(( STATUS_INTERVAL > 0 )) || die "GRAPHTYPER_STATUS_INTERVAL must be positive"
[[ "$FORCE_RERUN" == 0 || "$FORCE_RERUN" == 1 ]] || die "GRAPHTYPER_FORCE_RERUN must be 0 or 1"
[[ -x /usr/bin/time ]] || die "missing GNU time: /usr/bin/time"

declare -a TASKS=() QUEUED=() RUNNING=() DONE=() FAILED=()
declare -A TASK_SAMPLE=() TASK_DEPTH=() TASK_PID=() TASK_START=() TASK_MESSAGE=()
REUSED=0

is_complete() {
  local sample depth out
  sample=$1
  depth=$2
  out=$RESULT_ROOT/$sample/$depth
  [[ -s "$out/$sample.$depth.minigraph.graphtyper.sv.vcf.gz" && \
     -s "$out/$sample.$depth.minigraph.graphtyper.sv.vcf.gz.tbi" && \
     -s "$out/policy.txt" && \
     -s "$out/genotype_sv/.graphtyper_genotype_sv.complete" && \
     -s "$out/performance.tsv" ]]
}

perf_value() {
  local task=$1 key=$2 sample=${TASK_SAMPLE[$task]} depth=${TASK_DEPTH[$task]}
  awk -F= -v key="$key" '$1 == key {print $2; exit}' \
    "$RESULT_ROOT/$sample/$depth/performance.raw.txt" 2>/dev/null || true
}

write_performance_summary() {
  local task=$1 sample=${TASK_SAMPLE[$task]} depth=${TASK_DEPTH[$task]}
  local out=$RESULT_ROOT/$sample/$depth raw=$RESULT_ROOT/$sample/$depth/performance.raw.txt
  local tmp=$out/performance.tsv.tmp.$$
  [[ -s "$raw" ]] || return 0
  {
    printf 'sample\tdepth\tmetric\tvalue\n'
    while IFS='=' read -r key value; do
      [[ -n "$key" ]] || continue
      printf '%s\t%s\t%s\t%s\n' "$sample" "$depth" "$key" "$value"
    done < "$raw"
  } > "$tmp"
  mv -f "$tmp" "$out/performance.tsv"
}

task_state() {
  local task=$1
  case " ${DONE[*]} " in *" $task "*) printf 'completed'; return ;; esac
  case " ${FAILED[*]} " in *" $task "*) printf 'failed'; return ;; esac
  case " ${RUNNING[*]} " in *" $task "*) printf 'running'; return ;; esac
  case " ${QUEUED[*]} " in *" $task "*) printf 'queued'; return ;; esac
  printf 'unknown'
}

write_status() {
  local tmp=$STATUS.tmp.$$ now task sample depth pid elapsed out fragments message
  local wall user system cpu rss read_bytes write_bytes
  now=$(date +%s)
  printf 'task\tsample\tdepth\tstate\tpid\telapsed\tfragments\twall_seconds\tuser_seconds\tsystem_seconds\tcpu_percent\tmax_rss_kb\tfs_inputs\tfs_outputs\toutput\tmessage\n' > "$tmp"
  for task in "${TASKS[@]}"; do
    sample=${TASK_SAMPLE[$task]}
    depth=${TASK_DEPTH[$task]}
    out=$RESULT_ROOT/$sample/$depth
    pid=${TASK_PID[$task]:--}
    elapsed=0
    [[ -n ${TASK_START[$task]:-} ]] && elapsed=$((now - TASK_START[$task]))
    fragments=$( { find "$out/genotype_sv" -type f -name '*.vcf.gz' -size +0c -print 2>/dev/null || true; } | wc -l)
    wall=$(perf_value "$task" wall_seconds)
    user=$(perf_value "$task" user_seconds)
    system=$(perf_value "$task" system_seconds)
    cpu=$(perf_value "$task" cpu_percent)
    rss=$(perf_value "$task" max_rss_kb)
    read_bytes=$(perf_value "$task" fs_inputs)
    write_bytes=$(perf_value "$task" fs_outputs)
    message=${TASK_MESSAGE[$task]:-}
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$task" "$sample" "$depth" "$(task_state "$task")" "$pid" "$(fmt_time "$elapsed")" \
      "$fragments" "${wall:--}" "${user:--}" "${system:--}" "${cpu:--}" "${rss:--}" \
      "${read_bytes:--}" "${write_bytes:--}" "$out" "$message" >> "$tmp"
  done
  mv -f "$tmp" "$STATUS"
}

refresh_workers() {
  local -a next_running=()
  local task pid rc
  for task in "${RUNNING[@]}"; do
    pid=${TASK_PID[$task]}
    if kill -0 "$pid" 2>/dev/null; then
      next_running+=("$task")
      continue
    fi
    rc=0
    wait "$pid" || rc=$?
    if (( rc == 0 )); then
      write_performance_summary "$task"
    fi
    if (( rc == 0 )) && is_complete "${TASK_SAMPLE[$task]}" "${TASK_DEPTH[$task]}"; then
      DONE+=("$task")
      TASK_MESSAGE[$task]=completed
      log "DONE $task; completed=${#DONE[@]} running=$(( ${#RUNNING[@]} - 1 )) queued=${#QUEUED[@]} failed=${#FAILED[@]}"
    else
      FAILED+=("$task")
      TASK_MESSAGE[$task]="worker exit=$rc; inspect $RESULT_ROOT/${TASK_SAMPLE[$task]}/${TASK_DEPTH[$task]}/driver.log"
      log "FAILED $task; ${TASK_MESSAGE[$task]}"
    fi
  done
  RUNNING=("${next_running[@]}")
}

start_workers() {
  local task sample depth out pid
  while (( ${#RUNNING[@]} < MAX_CONCURRENT && ${#QUEUED[@]} > 0 )); do
    task=${QUEUED[0]}
    QUEUED=("${QUEUED[@]:1}")
    sample=${TASK_SAMPLE[$task]}
    depth=${TASK_DEPTH[$task]}
    out=$RESULT_ROOT/$sample/$depth
    mkdir -p "$out"
    log "START $task: GraphTyper genotype_sv, threads=$THREADS, worker_log=$out/driver.log, performance=$out/performance.tsv"
    GRAPHTYPER_SAMPLE=$sample GRAPHTYPER_DEPTH=$depth GRAPHTYPER_THREADS=$THREADS \
      /usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu_percent=%P\nmax_rss_kb=%M\nfs_inputs=%I\nfs_outputs=%O\nvoluntary_context_switches=%w\ninvoluntary_context_switches=%c\nexit_status=%x' \
      -o "$out/performance.raw.txt" "$WORKER" > "$out/driver.log" 2>&1 &
    pid=$!
    TASK_PID[$task]=$pid
    TASK_START[$task]=$(date +%s)
    TASK_MESSAGE[$task]=running
    RUNNING+=("$task")
  done
}

on_signal() {
  local task
  log "received termination signal; stopping ${#RUNNING[@]} running workers"
  for task in "${RUNNING[@]}"; do kill "${TASK_PID[$task]}" 2>/dev/null || true; done
  wait || true
  write_status
  exit 143
}
trap on_signal INT TERM

# R1/R2 pairs are the task manifest. GraphTyper consumes the corresponding
# shared BAM, which has already been produced from that exact read pair.
while IFS= read -r r1; do
  dir=$(dirname "$r1")
  sample=$(basename "$(dirname "$dir")")
  depth=$(basename "$dir")
  r2=$dir/${sample}_${depth}_R2.fastq.gz
  [[ -s "$r2" ]] || die "missing paired read file: $r2"
  [[ -s "$BASE/share/graphtyper/bam/$sample/$depth/$sample.$depth.MC.sorted.bam" ]] || die "missing shared BAM for $sample/$depth"
  [[ -s "$BASE/share/graphtyper/bam/$sample/$depth/bams.list" ]] || die "missing shared BAM list for $sample/$depth"
  task=$sample/$depth
  TASKS+=("$task")
  TASK_SAMPLE[$task]=$sample
  TASK_DEPTH[$task]=$depth
done < <(find "$READS" -mindepth 3 -maxdepth 3 -type f -name '*_R1.fastq.gz' -print | sort)

(( ${#TASKS[@]} == 24 )) || die "expected 24 Illumina sample/depth tasks, found ${#TASKS[@]}"
for task in "${TASKS[@]}"; do
  sample=${TASK_SAMPLE[$task]}
  depth=${TASK_DEPTH[$task]}
  if (( FORCE_RERUN == 1 )); then
    rm -rf "$RESULT_ROOT/$sample/$depth"
    QUEUED+=("$task")
    TASK_MESSAGE[$task]=forced_rerun
  elif is_complete "$sample" "$depth"; then
    DONE+=("$task")
    TASK_MESSAGE[$task]='pre-existing complete result reused'
    REUSED=$((REUSED + 1))
  else
    QUEUED+=("$task")
    TASK_MESSAGE[$task]=waiting
  fi
done

if (( FORCE_RERUN == 1 )); then
  : > "$MASTER_LOG"
fi
log "minigraph GraphTyper controller: tasks=${#TASKS[@]} reused=$REUSED queued=${#QUEUED[@]} max_concurrent=$MAX_CONCURRENT worker_threads=$THREADS total_requested_threads=$((MAX_CONCURRENT * THREADS)) force_rerun=$FORCE_RERUN"
log "status file: $STATUS"
write_status
if [[ "$DRY_RUN" == 1 ]]; then
  log "DRY RUN: validation and status generation completed; no workers started"
  exit 0
fi
last_status=0
while (( ${#QUEUED[@]} > 0 || ${#RUNNING[@]} > 0 )); do
  refresh_workers
  start_workers
  now=$(date +%s)
  if (( now - last_status >= STATUS_INTERVAL )); then
    write_status
    log "STATUS completed=${#DONE[@]}/${#TASKS[@]} running=${#RUNNING[@]} queued=${#QUEUED[@]} failed=${#FAILED[@]}"
    last_status=$now
  fi
  (( ${#QUEUED[@]} > 0 || ${#RUNNING[@]} > 0 )) && sleep 10
done
write_status

RUN_END=$(date +%s)
printf 'metric\tvalue\n' > "$RESULT_ROOT/run_all.performance.tsv"
printf 'controller_wall_seconds\t%s\n' "$((RUN_END - RUN_START))" >> "$RESULT_ROOT/run_all.performance.tsv"
printf 'tasks_total\t%s\n' "${#TASKS[@]}" >> "$RESULT_ROOT/run_all.performance.tsv"
printf 'tasks_reused\t%s\n' "$REUSED" >> "$RESULT_ROOT/run_all.performance.tsv"
printf 'tasks_completed\t%s\n' "${#DONE[@]}" >> "$RESULT_ROOT/run_all.performance.tsv"
printf 'tasks_failed\t%s\n' "${#FAILED[@]}" >> "$RESULT_ROOT/run_all.performance.tsv"
printf 'max_concurrent\t%s\n' "$MAX_CONCURRENT" >> "$RESULT_ROOT/run_all.performance.tsv"
printf 'threads_per_task\t%s\n' "$THREADS" >> "$RESULT_ROOT/run_all.performance.tsv"

if (( ${#FAILED[@]} > 0 )); then
  log "FINISHED WITH FAILURES: completed=${#DONE[@]} failed=${#FAILED[@]}"
  printf '%s\n' "${FAILED[@]}"
  exit 1
fi
log "COMPLETED: all ${#DONE[@]} tasks are complete; reused_pre-existing=$REUSED newly_completed=$((${#DONE[@]} - REUSED))"
