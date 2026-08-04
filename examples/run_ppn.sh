#!/usr/bin/env bash
# Run the baseline PPN calculation for a named nova_cases/ case, writing to
# <output_base>/baseline/.
#
# Usage: run_ppn.sh <nova_case> <output_base> [run_ppn.jl args...]
# Example: SensitivityStudy/run_ppn.sh ne_nova_1.15_12_X_weiss_mixed SensitivityStudy/runs
#   -> writes to SensitivityStudy/runs/baseline/
# Extra args (e.g. --option 3, --factor "label=value") are forwarded to run_ppn.jl.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <nova_case> <output_base> [run_ppn.jl args...]" >&2
    exit 1
fi

CASE="$1"
OUTPUT_BASE="$2"
shift 2

CASE_DIR="$SCRIPT_DIR/nova_cases/$CASE"
TRAJECTORY="$CASE_DIR/trajectory.input"
if [ -f "$CASE_DIR/initial_abundance_jch1.dat" ]; then
    ABUNDANCE="$CASE_DIR/initial_abundance_jch1.dat"
elif [ -f "$CASE_DIR/initial_abundance.dat" ]; then
    ABUNDANCE="$CASE_DIR/initial_abundance.dat"
else
    echo "no initial_abundance_jch1.dat or initial_abundance.dat in $CASE_DIR" >&2
    exit 1
fi
[ -f "$TRAJECTORY" ] || { echo "no trajectory.input in $CASE_DIR" >&2; exit 1; }

exec julia --project="$PROJECT_ROOT" "$SCRIPT_DIR/run_ppn.jl" \
    "$TRAJECTORY" "$ABUNDANCE" "$OUTPUT_BASE/baseline" "$@"
