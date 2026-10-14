#!/usr/bin/env bash
set -Eeuo pipefail
PERF_MODE=vg_minigraph exec "$(dirname "$0")/run_performance_task.sh" "$@"
