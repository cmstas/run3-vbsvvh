#!/bin/bash
# Summarize the status of all SLURM tasks matching a given submission date.
#
# Usage:
#   bash slurm/status_all.sh <date>           # e.g. 20260416
#   bash slurm/status_all.sh <date> --update   # also update manifests
#   bash slurm/status_all.sh <date> -v         # verbose: full status.py output per task
#
# Run from the preselection/ directory.

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <date> [--update] [-v]"
    echo "  date: submission date to filter on, e.g. 20260416"
    exit 1
fi

DATE="$1"
shift
UPDATE_FLAG=""
VERBOSE=false
for arg in "$@"; do
    case "$arg" in
        --update) UPDATE_FLAG="--update" ;;
        -v|--verbose) VERBOSE=true ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
JOBS_DIR="$SCRIPT_DIR/jobs"

# Find all tasks matching the date
TASKS=()
for dir in "$JOBS_DIR"/*/; do
    task_name="$(basename "$dir")"
    if [[ "$task_name" == *"$DATE"* ]]; then
        TASKS+=("$task_name")
    fi
done

if [ ${#TASKS[@]} -eq 0 ]; then
    echo "No tasks found matching date '$DATE'"
    exit 1
fi

# Colors
GREEN="\033[92m"
RED="\033[91m"
YELLOW="\033[93m"
BLUE="\033[94m"
GRAY="\033[90m"
BOLD="\033[1m"
RESET="\033[0m"

GRAND_TOTAL=0
GRAND_COMPLETED=0
GRAND_RUNNING=0
GRAND_QUEUED=0
GRAND_FAILED=0
GRAND_OTHER=0

echo ""
printf "${BOLD}%-45s %6s %6s %6s %6s %6s  %s${RESET}\n" "Task" "Total" "Done" "Run" "Queue" "Fail" ""
printf "%-45s %6s %6s %6s %6s %6s\n" "---------------------------------------------" "------" "------" "------" "------" "------"

for task in "${TASKS[@]}"; do
    # Run status.py and capture output
    output=$(python3 slurm/status.py -t "$task" $UPDATE_FLAG 2>&1)

    if $VERBOSE; then
        echo "============================================================"
        echo "$output"
        echo ""
    fi

    # Strip ANSI color codes before parsing
    clean=$(echo "$output" | sed 's/\x1b\[[0-9;]*m//g')

    # Parse counts from status.py output (lines like "completed        74     100.0%")
    n_total=$(echo "$clean" | grep -oP '^Total\s+\K\d+' || echo 0)
    n_completed=$(echo "$clean" | grep -oP 'completed\s+\K\d+' || echo 0)
    n_running=$(echo "$clean" | grep -oP 'running\s+\K\d+' || echo 0)
    n_queued=$(echo "$clean" | grep -oP 'queued\s+\K\d+' || echo 0)
    n_submitted=$(echo "$clean" | grep -oP 'submitted\s+\K\d+' || echo 0)
    n_failed=$(echo "$clean" | grep -oP 'failed\s+\K\d+' || echo 0)
    n_queued=$((n_queued + n_submitted))

    GRAND_TOTAL=$((GRAND_TOTAL + n_total))
    GRAND_COMPLETED=$((GRAND_COMPLETED + n_completed))
    GRAND_RUNNING=$((GRAND_RUNNING + n_running))
    GRAND_QUEUED=$((GRAND_QUEUED + n_queued))
    GRAND_FAILED=$((GRAND_FAILED + n_failed))

    # Short task name: strip "merged_" prefix and the date+channel suffix
    short_name=$(echo "$task" | sed "s/^merged_//;s/_${DATE}[^ ]*//")

    # Determine row color
    if [ "$n_total" -gt 0 ] && [ "$n_completed" -eq "$n_total" ]; then
        ROW_COLOR="$GREEN"
    elif [ "$n_failed" -gt 0 ] && [ "$n_running" -eq 0 ] && [ "$n_queued" -eq 0 ]; then
        ROW_COLOR="$RED"
    elif [ "$n_running" -gt 0 ] || [ "$n_queued" -gt 0 ]; then
        ROW_COLOR="$BLUE"
    else
        ROW_COLOR="$GRAY"
    fi

    # Progress bar (20 chars)
    bar=""
    if [ "$n_total" -gt 0 ]; then
        pct=$((100 * n_completed / n_total))
        filled=$((20 * n_completed / n_total))
        empty=$((20 - filled))
        bar+=$'\033[92m'
        for ((i=0; i<filled; i++)); do bar+="█"; done
        bar+=$'\033[90m'
        for ((i=0; i<empty; i++)); do bar+="░"; done
        bar+=$'\033[0m'
        bar+=" ${pct}%"
    fi

    printf "${ROW_COLOR}%-45s${RESET} %6d %6d %6d %6d %6d  %s\n" \
        "$short_name" "$n_total" "$n_completed" "$n_running" "$n_queued" "$n_failed" "$bar"
done

# Grand totals
echo ""
printf "%-45s %6s %6s %6s %6s %6s\n" "---------------------------------------------" "------" "------" "------" "------" "------"
printf "${BOLD}%-45s %6d %6d %6d %6d %6d${RESET}\n" \
    "TOTAL (${#TASKS[@]} tasks)" "$GRAND_TOTAL" "$GRAND_COMPLETED" "$GRAND_RUNNING" "$GRAND_QUEUED" "$GRAND_FAILED"
echo ""
