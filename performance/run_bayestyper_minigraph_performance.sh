#!/usr/bin/env bash
set -Eeuo pipefail
PERF_MODE=bayestyper_minigraph exec "$(dirname "$0")/run_performance_task.sh" "$@"
