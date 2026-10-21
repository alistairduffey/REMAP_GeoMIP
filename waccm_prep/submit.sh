#!/bin/bash
# Submit one prep job per scenario (they run in parallel as the queue allows).
# Usage: bash submit.sh [name ...]               default: all in scenarios.txt
#        OVERWRITE=1 bash submit.sh [name ...]   replace existing files
# A job that hits the walltime can simply be resubmitted: finished files
# are skipped unless OVERWRITE=1.
set -e
source "$(dirname "$0")/config.sh"
list_scenarios "$@" | while read -r name dir; do
    out=$(scen_dir "$name")
    mkdir -p "$out"
    id=$(qsub -A "$ACCOUNT" -q "$QUEUE" -l walltime="$WALLTIME" -l "$RESOURCES" \
              -N "prep_$name" -o "$out/prep.log" \
              -v NAME="$name",INDIR="$dir",PREP_HOME="$PREP_HOME",OVERWRITE="$OVERWRITE" \
              "$PREP_HOME/prep_one.pbs")
    echo "$name  $id  -> $out"
done
