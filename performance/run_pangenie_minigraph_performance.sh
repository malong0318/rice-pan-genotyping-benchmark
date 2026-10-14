#!/usr/bin/env bash
set -Eeuo pipefail
PERF_MODE=pangenie_minigraph exec "$(dirname "$0")/run_performance_task.sh" "$@"
