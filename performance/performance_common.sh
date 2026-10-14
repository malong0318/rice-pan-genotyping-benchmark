#!/usr/bin/env bash
set -Eeuo pipefail

perf_init() {
    local out=$1 script=$2
    shift 2
    mkdir -p "$out"
    if [ "$(printenv PERF_TIMED 2>/dev/null || true)" != 1 ]; then
        export PERF_TIMED=1
        exec /usr/bin/time -v -o "$out/performance.raw.txt" "$script" "$@"
    fi
}

perf_log() {
    printf '[%s] %s\n' "$(date -Is)" "$*"
}

perf_die() {
    perf_log "ERROR: $*" >&2
    exit 1
}

perf_need_bin() {
    [ -x "$1" ] || perf_die "missing executable: $1"
}

perf_need_file() {
    [ -s "$1" ] || perf_die "missing or empty file: $1"
}

perf_finish() {
    local out=$1 tool=$2 graph=$3 result=$4 end wall user system rss cpu inputs outputs records avg
    end=$(date +%s)
    wall=$((end - PERF_START))
    user=$(awk -F: '/^User time \(seconds\):/ {gsub(/^ +/,"",$2); print $2}' "$out/performance.raw.txt" 2>/dev/null || true)
    system=$(awk -F: '/^System time \(seconds\):/ {gsub(/^ +/,"",$2); print $2}' "$out/performance.raw.txt" 2>/dev/null || true)
    rss=$(awk -F: '/^Maximum resident set size \(kbytes\):/ {gsub(/^ +/,"",$2); print $2}' "$out/performance.raw.txt" 2>/dev/null || true)
    cpu=$(awk -F: '/^Percent of CPU this job got:/ {gsub(/^ +/,"",$2); print $2}' "$out/performance.raw.txt" 2>/dev/null || true)
    inputs=$(awk -F: '/^File system inputs:/ {gsub(/^ +/,"",$2); print $2}' "$out/performance.raw.txt" 2>/dev/null || true)
    outputs=$(awk -F: '/^File system outputs:/ {gsub(/^ +/,"",$2); print $2}' "$out/performance.raw.txt" 2>/dev/null || true)
    [ -n "$user" ] || user=0
    [ -n "$system" ] || system=0
    [ -n "$rss" ] || rss=0
    [ -n "$cpu" ] || cpu=0
    [ -n "$inputs" ] || inputs=0
    [ -n "$outputs" ] || outputs=0
    avg=$(awk -v u="$user" -v s="$system" -v w="$wall" 'BEGIN {if (w>0) printf "%.3f", (u+s)/w; else print "0"}')
    records=NA
    if [ -s "$result" ]; then
        records=$("$BCFTOOLS" index -n "$result" 2>/dev/null || "$BCFTOOLS" view -H "$result" | wc -l)
    fi
    {
        printf 'metric\tvalue\n'
        printf 'tool\t%s\n' "$tool"
        printf 'graph\t%s\n' "$graph"
        printf 'sample\t%s\n' "$SAMPLE"
        printf 'depth\t%s\n' "$DEPTH"
        printf 'configured_threads\t%s\n' "$THREADS"
        printf 'wall_seconds\t%s\n' "$wall"
        printf 'user_seconds\t%s\n' "$user"
        printf 'system_seconds\t%s\n' "$system"
        printf 'cpu_percent\t%s\n' "$cpu"
        printf 'average_cpu_cores\t%s\n' "$avg"
        printf 'max_rss_kb\t%s\n' "$rss"
        printf 'filesystem_inputs\t%s\n' "$inputs"
        printf 'filesystem_outputs\t%s\n' "$outputs"
        printf 'result_records\t%s\n' "$records"
        printf 'completed\t%s\n' "$(date -Is)"
    } > "$out/performance.tsv"
    perf_log "completed $tool/$graph wall=$wall seconds average_cpu_cores=$avg max_rss_kb=$rss"
}
