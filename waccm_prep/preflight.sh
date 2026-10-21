#!/bin/bash
# Fast check on a login node; writes nothing.
# For each scenario: input directory readable, all required variables
# present, one contiguous period, grid 192 x 70.
# Usage: bash preflight.sh [name ...]
source "$(dirname "$0")/config.sh"
activate_env

printf "%-22s %-6s %s\n" SCENARIO STATUS "PERIOD / PROBLEM"
list_scenarios "$@" | while read -r name dir; do
    out=$(python "$PREP_HOME/prep_waccm_zm.py" --input-dir "$dir" \
                 --out-dir "$(scen_dir "$name")" --dry-run 2>&1)
    if [ $? -eq 0 ]; then
        printf "%-22s %-6s %s\n" "$name" OK \
            "$(echo "$out" | sed -n 's/^Simulation period from file names: //p')"
    else
        printf "%-22s %-6s %s\n" "$name" FAIL \
            "$(echo "$out" | grep -E "DRY RUN|not found|differ|contiguous|Error" | head -2 | tr '\n' ' ')"
    fi
done
