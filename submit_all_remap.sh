#!/bin/bash
# Run REMAP (run_parallel.pbs) for every prepped scenario, then rename the
# outputs (finalize_outputs.py) in a job that waits for all of them.
#
# Run from anywhere, on a login node:
#   bash submit_all_remap.sh                 # all scenarios in waccm_prep/scenarios.txt
#   bash submit_all_remap.sh NAME1 NAME2     # only these
# Settings (environment): NCPUS (8), MEM (16GB), WALLTIME (02:00:00), OPT (-O2)
#
# Outputs: out/<name>/ (merged), then final/geomip_WACCM_<name>_to_NorESM2-LM_*.nc
set -e
# read settings before config.sh, which defines its own WALLTIME for the prep jobs
NCPUS=${NCPUS:-8}; MEM=${MEM:-16GB}; R_WALLTIME=${WALLTIME:-02:00:00}; OPT=${OPT:--O2}
RUN=$(cd "$(dirname "$0")" && pwd)
cd "$RUN"
source waccm_prep/config.sh          # ACCOUNT, QUEUE, list_scenarios, scen_dir
mkdir -p logs out

# ---- checks that apply to every scenario --------------------------------------------
[ -f lookup_table.dat ] || { echo "lookup_table.dat missing"; exit 1; }
for f in write_WACCM_sun.f write_WACCM_earth.f; do
    grep -q "^            wtt=wts/100d0" "$f" || {
        echo "$f does not contain the wt% fix (wtt=wts/100d0) - stopping"; exit 1; }
done

ids=(); names=()
while read -r name dir; do
    src=$(scen_dir "$name")
    if [ ! -f "$src/period.dat" ] || ! grep -qs "All requirements met" "$src/prep.log"; then
        echo "skip  $name: prep not complete/validated ($src)"; continue
    fi
    P=$(head -1 "$src/period.dat" | tr -d '[:space:]')
    cd_=$RUN/chunks/$(basename "$(realpath "$src")")_$P
    if [ -d "$cd_" ] && [ "$src/period.dat" -nt "$cd_/SOURCE" ]; then
        echo "skip  $name: $cd_ was split before the inputs were last prepped."
        echo "      Remove it first:  rm -rf $cd_"; continue
    fi
    if ls out/"$name"/ccmi_WACCM_to_ECHAM6_*.nc >/dev/null 2>&1 || \
       ls final/geomip_WACCM_"$name"_to_*.nc >/dev/null 2>&1; then
        echo "skip  $name: outputs already exist in out/$name or final/"; continue
    fi
    id=$(qsub -A "$ACCOUNT" -q "$QUEUE" -l walltime="$R_WALLTIME" \
              -l select=1:ncpus="$NCPUS":mem="$MEM" \
              -N "remap_$name" -o "$RUN/logs/remap_$name.log" \
              -v SRC="$src",OUT="$RUN/out/$name",NPROC="$NCPUS",OPT="$OPT" \
              "$RUN/run_parallel.pbs")
    echo "sub   $name ($P)  $id"
    ids+=("$id"); names+=("$name")
done < <(list_scenarios "$@")

[ ${#ids[@]} -gt 0 ] || { echo "nothing submitted"; exit 0; }

dep=$(IFS=:; echo "${ids[*]}")
fid=$(qsub -A "$ACCOUNT" -q "$QUEUE" -l walltime=03:00:00 -l select=1:ncpus=1:mem=8GB \
           -N remap_finalize -o "$RUN/logs/finalize.log" -W depend=afterany:"$dep" \
           -v NAMES="$(IFS=:; echo "${names[*]}")" "$RUN/finalize.pbs")
echo "final $fid  (runs after all ${#ids[@]} jobs; renames whatever succeeded)"
echo "logs: $RUN/logs/"

