#!/usr/bin/env bash
set -Eeuo pipefail
PERF_MODE=pangenie_mc exec "$(dirname "$0")/run_performance_task.sh" "$@"
