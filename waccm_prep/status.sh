#!/bin/bash
# One line per scenario, from its prep.log.
# Usage: bash status.sh [name ...]
source "$(dirname "$0")/config.sh"
printf "%-22s %-10s %-15s %s\n" SCENARIO STATUS PERIOD SUMMARY
list_scenarios "$@" | while read -r name dir; do
    out=$(scen_dir "$name"); log=$out/prep.log
    per=$(cat "$out/period.dat" 2>/dev/null || echo -)
    if [ ! -f "$log" ]; then
        st="-"; sm="no log yet (queued, running, or not submitted)"
    elif grep -q "All requirements met" "$log"; then
        st=OK; sm=$(grep -E "^[0-9]+ failure" "$log" | tail -1)
    elif grep -q "failure(s)" "$log"; then
        st=FAIL; sm=$(grep -E "^[0-9]+ failure" "$log" | tail -1)
    else
        st=INCOMPLETE; sm=$(grep -iE "walltime|killed|error|Traceback" "$log" | tail -1)
    fi
    printf "%-22s %-10s %-15s %s\n" "$name" "$st" "$per" "$sm"
done
echo
echo "Logs: $OUT_ROOT/WACCM_data_<name>/prep.log"
