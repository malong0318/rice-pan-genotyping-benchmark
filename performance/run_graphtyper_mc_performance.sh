#!/usr/bin/env bash
set -Eeuo pipefail
PERF_MODE=graphtyper_mc exec "$(dirname "$0")/run_performance_task.sh" "$@"
