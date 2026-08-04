#!/usr/bin/env bash
# Decay a run_ppn.sh baseline for a named nova_cases/ case, writing to
# <decay_output_base>/baseline/.
#
# Usage: decay_ppn.sh <nova_case> <run_output_base> <decay_output_base> <decay_time_seconds> [decay_ppn.jl args...]
# Example: SensitivityStudy/decay_ppn.sh ne_nova_1.15_12_X_weiss_mixed \
#              SensitivityStudy/runs SensitivityStudy/decay_runs 7200
#   -> reads SensitivityStudy/runs/baseline/final_state.csv
#   -> writes SensitivityStudy/decay_runs/baseline/decayed_state_7200s.csv
# Extra args (e.g. --option 3) are forwarded to decay_ppn.jl.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

if [ $# -lt 4 ]; then
    echo "Usage: $0 <nova_case> <run_output_base> <decay_output_base> <decay_time_seconds> [decay_ppn.jl args...]" >&2
    exit 1
fi

CASE="$1"
RUN_OUTPUT_BASE="$2"
DECAY_OUTPUT_BASE="$3"
DECAY_TIME="$4"
shift 4

CASE_DIR="$SCRIPT_DIR/nova_cases/$CASE"
[ -d "$CASE_DIR" ] || { echo "no nova_cases/$CASE directory" >&2; exit 1; }

exec julia --project="$PROJECT_ROOT" "$SCRIPT_DIR/decay_ppn.jl" \
    "$RUN_OUTPUT_BASE/baseline" "$DECAY_TIME" "$DECAY_OUTPUT_BASE/baseline" "$@"
