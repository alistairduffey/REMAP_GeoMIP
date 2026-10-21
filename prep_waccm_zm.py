#!/usr/bin/env python
"""
Prepare zonal-mean WACCM inputs for the REMAP Mie tool (write_WACCM_*.f).

Reads CESM monthly time-series files (one variable per file, e.g.
<case>.cam.h0.so4_a1.203501-208412.nc), averages over longitude, and writes
the files the Fortran reads into <run-dir>/WACCM_data/, named with the real
simulation period:

    WACCM_data/<name>.<YYYYMM-YYYYMM>.nc      (see VARS below)

It also writes <run-dir>/period.dat, which the (modified) Fortran reads to
set the number of months and build the same file names. Then it validates
every output file the way the Fortran will read it.

What the Fortran requires:
  - variable shape: (time, lev, lat) for 3D fields, (time, lat) for PS
  - sizes:          lat = 192, lev = 70 (hardcoded); time = number of
                    months in period.dat (checked again by the Fortran)
  - PS file also holds lev, lat, time (with units/calendar), hyam, hybm;
    pressure is computed as hyam*1e5 + hybm*PS, i.e. P0 = 1e5 Pa

The period is taken from the input file names. If a variable is split over
several files they must be contiguous, and every variable must cover the
same period.

Usage
  # process + validate (run from anywhere; --run-dir is the REMAP directory)
  python prep_waccm_zm.py --input-dir /glade/campaign/.../month_1 --run-dir .

  # validate an existing run directory only (uses its period.dat)
  python prep_waccm_zm.py --validate-only --run-dir .

Requires xarray, dask, netCDF4 (e.g. NCAR's `npl` conda environment).
Each 3D input is several GB; run on a compute node, not a login node.
"""

import argparse
import datetime
import glob
import os
import re
import sys

import numpy as np
import xarray as xr
import netCDF4

# ----------------------------------------------------------------------------
# What the Fortran expects
# ----------------------------------------------------------------------------
NLAT, NLEV = 192, 70  # hardcoded in the Fortran
PERIOD_RE = re.compile(r"^(\d{4})(\d{2})-(\d{4})(\d{2})$")

# CESM variable -> (output file stem, expected units, required by Fortran?,
#                   plausible (min, max) range for a warning check)
VARS = {
    "PS":        ("h0zm.PS",        "Pa",      True,  (3.0e4, 1.2e5)),
    "num_a1":    ("h0zm.num_a1",    "1/kg",    True,  (0.0, 1.0e14)),
    "num_a2":    ("h0zm.num_a2",    "1/kg",    True,  (0.0, 1.0e14)),
    "num_a3":    ("h0zm.num_a3",    "1/kg",    True,  (0.0, 1.0e14)),
    "dgnumwet1": ("h0zm.dgnumwet1", "m",       True,  (0.0, 1.0e-4)),
    "dgnumwet2": ("h0zm.dgnumwet2", "m",       True,  (0.0, 1.0e-4)),
    "dgnumwet3": ("h0zm.dgnumwet3", "m",       True,  (0.0, 1.0e-4)),
    "T":         ("T",              "K",       True,  (100.0, 400.0)),
    "RELHUM":    ("RELHUM",         "percent", True,  (0.0, 200.0)),
    # read by write_WACCM_sun/earth (program stops if missing) but not used
    "so4_a1":    ("so4_a1",         "kg/kg",   True,  (0.0, 1.0)),
    "so4_a2":    ("so4_a2",         "kg/kg",   True,  (0.0, 1.0)),
    "so4_a3":    ("so4_a3",         "kg/kg",   True,  (0.0, 1.0)),
    # not read by the Fortran; useful for checking the SAD from write_WACCM_p
    "SAD_AERO":  ("SAD_AERO",       "cm2/cm3", False, (0.0, 1.0e-4)),
}
PS_EXTRA = ["hyam", "hybm"]  # read_ps also reads lev, lat, time (coords)


def out_path(outdir, var, period):
    return os.path.join(outdir, f"{VARS[var][0]}.{period}.nc")


def n_months(period):
    y1, m1, y2, m2 = map(int, PERIOD_RE.match(period).groups())
    return (y2 - y1) * 12 + (m2 - m1) + 1


def next_month(yyyymm):
    y, m = int(yyyymm[:4]), int(yyyymm[4:])
    return f"{y + m // 12:04d}{m % 12 + 1:02d}"


# ----------------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------------
def find_inputs(indir, var):
    """Files <anything>.cam.h0.<var>.YYYYMM-YYYYMM.nc, matching var exactly.

    Returns (files, period) with files sorted in time. Raises if the files
    are not contiguous."""
    pat = re.compile(rf"\.cam\.h0\.{re.escape(var)}\.(\d{{6}})-(\d{{6}})\.nc$")
    hits = []
    for f in glob.glob(os.path.join(indir, "*.nc")):
        m = pat.search(os.path.basename(f))
        if m:
            hits.append((m.group(1), m.group(2), f))
    if not hits:
        return [], None
    hits.sort()
    for (s1, e1, f1), (s2, e2, f2) in zip(hits, hits[1:]):
        if s2 != next_month(e1):
            raise ValueError(f"{var}: files not contiguous: "
                             f"{os.path.basename(f1)} ends {e1}, "
                             f"{os.path.basename(f2)} starts {s2}")
    return [h[2] for h in hits], f"{hits[0][0]}-{hits[-1][1]}"


def zonal_mean(files, var):
    keep = [var, "hyam", "hybm", "P0"]

    def pre(ds):
        return ds[[v for v in keep if v in ds.variables]]

    # decode_times=False: copy CESM's raw time values/units through untouched
    ds = xr.open_mfdataset(files, combine="by_coords", preprocess=pre,
                           decode_times=False, chunks={"time": 12},
                           data_vars="minimal", coords="minimal",
                           compat="override")
    da = ds[var]
    if "lon" in da.dims:
        da = da.mean("lon", keep_attrs=True)
    for d in ("zlon", "lon"):
        if d in da.dims and da.sizes[d] == 1:
            da = da.squeeze(d, drop=True)

    order = ("time", "lat") if var == "PS" else ("time", "lev", "lat")
    if set(da.dims) != set(order):
        raise ValueError(f"{var}: unexpected dims {da.dims}, need {order}")
    da = da.transpose(*order).astype("float32")
    da.attrs["cell_methods"] = "time: mean lon: mean"

    out = xr.Dataset({var: da})
    for v in ("hyam", "hybm", "P0"):
        if v in ds.variables:
            out[v] = ds[v]
    for c in ("lat", "lev", "time"):
        if c in ds.coords:
            out[c].attrs = ds[c].attrs
    out.attrs = dict(ds.attrs)
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    out.attrs["history"] = (f"{stamp}: zonal mean by prep_waccm_zm.py from "
                            f"{', '.join(os.path.basename(f) for f in files)}\n"
                            + str(ds.attrs.get("history", "")))
    return out


def write(ds, path):
    enc = {v: {"_FillValue": None} for v in ds.variables}
    tmp = path + ".tmp"
    # netCDF-3 64-bit offset: readable by any netCDF library the old Fortran
    # might be linked against; 972x70x192 float32 is ~52 MB, well within limits
    ds.to_netcdf(tmp, format="NETCDF3_64BIT", encoding=enc)
    os.replace(tmp, path)


# ----------------------------------------------------------------------------
# Validation: open each file raw, as the Fortran sees it
# ----------------------------------------------------------------------------
class Report:
    def __init__(self):
        self.fail, self.warn = [], []

    def F(self, msg):
        self.fail.append(msg)
        print(f"  FAIL  {msg}")

    def W(self, msg):
        self.warn.append(msg)
        print(f"  warn  {msg}")

    def ok(self, msg):
        print(f"  ok    {msg}")


def validate(outdir, period, rep):
    ntime = n_months(period)
    print(f"Period {period}: {ntime} months")
    ref = {}  # coordinate values that must agree across files
    for var, (stem, units, required, (lo, hi)) in VARS.items():
        path = out_path(outdir, var, period)
        print(f"\n{os.path.basename(path)}")
        if not os.path.exists(path):
            (rep.F if required else rep.W)(f"{var}: file missing")
            continue

        with netCDF4.Dataset(path) as nc:
            nc.set_auto_mask(False)
            if var not in nc.variables:
                rep.F(f"{var}: variable '{var}' not in file")
                continue
            v = nc.variables[var]

            # dimension names, order and sizes
            want = ("time", "lat") if var == "PS" else ("time", "lev", "lat")
            if v.dimensions != want:
                rep.F(f"{var}: dims {v.dimensions}, Fortran needs {want}")
            sizes = {"time": ntime, "lev": NLEV, "lat": NLAT}
            for d in v.dimensions:
                n = len(nc.dimensions[d])
                if d in sizes and n != sizes[d]:
                    msg = f"{var}: {d} = {n}, expected {sizes[d]}"
                    if not required:
                        rep.W(msg + " (file not read by the Fortran)")
                    else:
                        rep.F(msg)

            # units
            u = str(getattr(v, "units", "")).strip()
            if u != units:
                # so4_a* and SAD_AERO are never used in the calculation
                used = required and not var.startswith("so4")
                (rep.F if used else rep.W)(
                    f"{var}: units '{u}', expected '{units}'")

            # values
            data = v[:]
            nbad = int(np.size(data) - np.isfinite(data).sum())
            if nbad:
                rep.F(f"{var}: {nbad} non-finite values (NaN/Inf) - the "
                      "Fortran will propagate these")
            fill = getattr(v, "_FillValue", None)
            if fill is not None and np.any(data == fill):
                rep.F(f"{var}: contains _FillValue {fill}")
            fin = data[np.isfinite(data)]
            if fin.size:
                dmin, dmax = float(fin.min()), float(fin.max())
                if dmin < lo or dmax > hi:
                    rep.W(f"{var}: range [{dmin:.3g}, {dmax:.3g}] outside "
                          f"plausible [{lo:.3g}, {hi:.3g}]")
                else:
                    rep.ok(f"{var}: range [{dmin:.3g}, {dmax:.3g}]")
            if var.startswith("dgnumwet"):
                nz = int((data <= 0).sum())
                if nz:
                    rep.W(f"{var}: {nz} values <= 0 (log(0) in lookup; "
                          "can give NaN omega/g in output)")

            # coordinates must agree across all files
            for c in ("lat", "lev", "time"):
                if c in nc.variables:
                    vals = nc.variables[c][:]
                    if c not in ref:
                        ref[c] = (vals, var)
                    elif (vals.shape != ref[c][0].shape
                          or not np.allclose(vals, ref[c][0])):
                        rep.F(f"{var}: '{c}' differs from {ref[c][1]}")

            # PS file: everything read_ps() needs
            if var == "PS":
                for need in ["lev", "lat", "time"] + PS_EXTRA:
                    if need not in nc.variables:
                        rep.F(f"PS file missing '{need}' (read by read_ps)")
                if "P0" in nc.variables:
                    p0 = float(nc.variables["P0"][...])
                    if not np.isclose(p0, 1.0e5):
                        rep.F(f"P0 = {p0}; Fortran assumes 1e5 Pa")
                if {"hyam", "hybm", "lev"} <= set(nc.variables):
                    a = nc.variables["hyam"][:]
                    b = nc.variables["hybm"][:]
                    lev = nc.variables["lev"][:]
                    if not np.allclose(1000.0 * (a + b), lev, rtol=1e-3):
                        rep.F("lev != 1000*(hyam+hybm): hybrid coefficients "
                              "inconsistent with Fortran's pressure formula")
                    else:
                        rep.ok("hybrid coefficients consistent with lev")
                if "time" in nc.variables:
                    t = nc.variables["time"]
                    tv = t[:]
                    if np.any(np.diff(tv) <= 0):
                        rep.F("time not strictly increasing")
                    tu = getattr(t, "units", None)
                    tc = getattr(t, "calendar", None)
                    if not tu:
                        rep.F("time has no units attribute (the Fortran "
                              "copies it to the output and stops without it)")
                    if not tc:
                        rep.W("time has no calendar attribute; output time "
                              "will have none")
                    first = last = "?"
                    try:
                        first, last = netCDF4.num2date(
                            [tv[0], tv[-1]], tu, tc or "standard")
                    except Exception:
                        pass
                    print(f"  info  time: {len(tv)} steps, {tv[0]:g} .. "
                          f"{tv[-1]:g} '{tu}' ({tc}) = {first} .. {last}. "
                          "These values, units and calendar are copied to "
                          "the output unchanged.")


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--input-dir", help="CESM month_1 time-series directory")
    ap.add_argument("--run-dir", default=".",
                    help="REMAP run directory (gets WACCM_data/ and "
                         "period.dat); default: current directory")
    ap.add_argument("--validate-only", action="store_true",
                    help="skip processing; check existing output files "
                         "against <run-dir>/period.dat")
    ap.add_argument("--overwrite", action="store_true",
                    help="replace existing output files")
    args = ap.parse_args()

    outdir = os.path.join(args.run_dir, "WACCM_data")
    period_file = os.path.join(args.run_dir, "period.dat")
    os.makedirs(outdir, exist_ok=True)

    if not args.validate_only:
        if not args.input_dir:
            ap.error("--input-dir is required unless --validate-only")

        # find inputs and establish one common period
        inputs, periods = {}, {}
        for var, (stem, units, required, _) in VARS.items():
            files, per = find_inputs(args.input_dir, var)
            if not files:
                print(f"{'MISSING' if required else 'skip'}: no input "
                      f"files for {var}")
                continue
            inputs[var], periods[var] = files, per
        req_periods = {periods[v] for v in periods if VARS[v][2]}
        if len(req_periods) != 1:
            print("Required variables cover different periods (or none "
                  "were found):")
            for v, per in periods.items():
                print(f"  {v}: {per}")
            sys.exit(1)
        period = req_periods.pop()
        print(f"Simulation period from file names: {period} "
              f"({n_months(period)} months)\n")

        for var, files in inputs.items():
            if periods[var] != period:
                print(f"skip: {var} covers {periods[var]}, not {period}")
                continue
            dest = out_path(outdir, var, period)
            if os.path.exists(dest) and not args.overwrite:
                print(f"exists, skipping (use --overwrite): {dest}")
                continue
            print(f"{var}: {len(files)} file(s) -> {dest}")
            write(zonal_mean(files, var), dest)

        with open(period_file, "w") as f:
            f.write(period + "\n")
        print(f"\nwrote {period_file}: {period}")
    else:
        if not os.path.exists(period_file):
            sys.exit(f"{period_file} not found")
        period = open(period_file).readline().strip()

    if not PERIOD_RE.match(period):
        sys.exit(f"bad period '{period}' (need YYYYMM-YYYYMM)")

    print("\n=========== VALIDATION ===========")
    rep = Report()
    validate(outdir, period, rep)

    others = sorted(f for f in os.listdir(outdir)
                    if f.endswith(".nc") and f".{period}." not in f)
    if others:
        print(f"\nnote: {len(others)} other .nc file(s) in {outdir} "
              "(different period; ignored by the Fortran)")

    print("\n=========== SUMMARY ===========")
    print(f"{len(rep.fail)} failure(s), {len(rep.warn)} warning(s)")
    for m in rep.fail:
        print(f"  FAIL  {m}")
    if rep.fail:
        print("\nDo NOT run the Fortran until the failures are resolved.")
        sys.exit(1)
    print(f"All requirements met for period {period}.")


if __name__ == "__main__":
    main()
