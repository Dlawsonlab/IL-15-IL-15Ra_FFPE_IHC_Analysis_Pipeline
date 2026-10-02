#!/usr/bin/env bash
# Regenerate the cohort table, statistics, final figure and threshold summaries from the
# per-section Fiji outputs, in dependency order.
#
# Usage:  bash scripts/rebuild_all_outputs.sh
# Run after any cohort re-run, any change to data/decisions/, and any threshold change.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1

step() {
    local name="$1"; shift
    printf '\n=== %s ===\n' "$name"
    if "$@"; then printf '    ok\n'
    else printf '    *** FAILED: later steps may be stale ***\n'; FAILED+=("$name"); fi
}
FAILED=()

step "1/4 cohort table"               python3 scripts/pipeline/aggregate_cohort.py
step "2/4 cohort statistics"          python3 scripts/pipeline/plot_cohort_stats.py
step "3/4 final figure, Fig. 3 F-G"   python3 scripts/figures/plot_cohort_bars_pct.py
step "4/4 threshold scan summary"     python3 scripts/threshold_evidence/summarize_threshold_scan.py

if ((${#FAILED[@]})); then
    printf '\nFAILED STEPS: %s\n' "${FAILED[*]}"; exit 1
fi
printf '\nAll outputs regenerated.\n'
