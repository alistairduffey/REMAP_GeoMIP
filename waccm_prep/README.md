# waccm_prep

Turns CESM2-WACCM monthly time-series output into zonal-mean input
directories for the REMAP Mie tool, one directory per scenario.

## Output

```
$OUT_ROOT/WACCM_data_<name>/
    h0zm.PS.<period>.nc
    h0zm.num_a{1,2,3}.<period>.nc
    h0zm.dgnumwet{1,2,3}.<period>.nc
    T.<period>.nc   RELHUM.<period>.nc
    so4_a{1,2,3}.<period>.nc   SAD_AERO.<period>.nc
    period.dat      one line, YYYYMM-YYYYMM
    source.txt      input directory and files used
    prep.log        job log, including the validation report
```

`<period>` is taken from the input file names (e.g. `203501-208412`).

## Files

| File | Purpose |
|---|---|
| `config.sh` | Output location, PBS account/queue/resources, conda env |
| `scenarios.txt` | One line per scenario: `name  input_dir` |
| `prep_waccm_zm.py` | Does the work for one scenario (also usable on its own) |
| `preflight.sh` | Fast login-node check of every scenario; writes nothing |
| `submit.sh` | Submits one PBS job per scenario |
| `prep_one.pbs` | The job script `submit.sh` submits |
| `status.sh` | One-line status per scenario |

## Workflow

Run from anywhere; the scripts find their own directory.

```bash
bash waccm_prep/preflight.sh          # all rows should be OK
bash waccm_prep/submit.sh             # or: submit.sh NAME1 NAME2
bash waccm_prep/status.sh             # check progress / results
```

Resubmit failed or timed-out scenarios by name. Finished files are
skipped; use `OVERWRITE=1 bash waccm_prep/submit.sh NAME` to rewrite them.

Single scenario by hand (compute node):

```bash
python prep_waccm_zm.py --input-dir <CESM month_1 dir> --out-dir <dir>
python prep_waccm_zm.py --validate-only --out-dir <dir>
```

## Using an output directory with the current REMAP code

Until REMAP takes the input directory as an argument, it reads
`WACCM_data/` and `period.dat` from its run directory. Link a scenario in:

```bash
cd /glade/derecho/scratch/duffey/REMAP_NorESM
ln -sfn scenarios/WACCM_data_<name> WACCM_data
ln -sfn WACCM_data/period.dat period.dat
```

## Expected warnings

T (up to ~840 K) and RELHUM (up to ~10^5 %) exceed the plausible-range
checks in the mesosphere/lower thermosphere near the WACCM model top.
These are physical and harmless for the optics.
