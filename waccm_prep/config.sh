# Settings for waccm_prep. Edit here; the scripts source this file.

# Where the scenario directories are written: $OUT_ROOT/WACCM_data_<name>/
OUT_ROOT=/glade/derecho/scratch/duffey/REMAP_NorESM/scenarios

# PBS
ACCOUNT=wyom0228
QUEUE=main
WALLTIME=02:00:00
RESOURCES=select=1:ncpus=4:mem=32GB

# Python environment with xarray, dask, netCDF4
CONDA_ENV=npl

# ---- no need to edit below -------------------------------------------------
PREP_HOME=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SCEN_FILE=$PREP_HOME/scenarios.txt

scen_dir() { echo "$OUT_ROOT/WACCM_data_$1"; }

# "name input_dir" pairs from scenarios.txt; with arguments, only those names
list_scenarios() {
    awk '!/^[[:space:]]*#/ && NF>=2 {print $1, $2}' "$SCEN_FILE" |
    if [ $# -gt 0 ]; then grep -E "^($(IFS='|'; echo "$*")) "; else cat; fi
}

activate_env() {
    module load conda >/dev/null 2>&1
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
}
