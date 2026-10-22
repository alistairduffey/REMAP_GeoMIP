#!/usr/bin/env python
"""
Stratospheric aerosol optical depth (SAOD) from the REMAP outputs, compared
with WACCM's own diagnostic.

  python compare_saod.py --remap-dir <dir with ccmi_WACCM_to_ECHAM6_*.nc> \
                         --waccm-dir <CESM month_1 time-series dir> \
                         [--waccm-in <dir with the prepped T.<period>.nc>] \
                         [--band 4] [--png saod.png]

REMAP SAOD = sum over layers above the tropopause of ext_sun[band] * dz.
  ext_sun  from ccmi_WACCM_to_ECHAM6_sun.nc   [km-1]
  pressure from ccmi_WACCM_to_ECHAM6_p.nc     [Pa]
  dz       hydrostatic, from pressure and the prepped WACCM temperature
  tropopause: WACCM's TROP_P (zonal mean) if in --waccm-dir, else 100 hPa

WACCM SAOD: AODVISstdn (stratospheric AOD at 550 nm, day and night) if
present, else the first AOD*st* variable found; zonal mean, then global
mean. Band 4 (0.442-0.625 um) is the RRTMG band CAM uses for 550 nm.
"""
import argparse
import glob
import os
import sys

import numpy as np
import xarray as xr

R_D, G = 287.04, 9.80616


def open_cesm_var(d, var):
    files = sorted(glob.glob(os.path.join(d, f"*.cam.h0.{var}.*.nc")))
    if not files:
        return None
    ds = xr.open_mfdataset(files, combine="nested", concat_dim="time",
                           decode_times=False, data_vars="minimal",
                           coords="minimal", compat="override")
    ds = ds.drop_duplicates("time")
    da = ds[var]
    return da.mean("lon") if "lon" in da.dims else da


def find_waccm_saod(d):
    cands = ["AODVISstdn", "AODVISst"]
    found = sorted({os.path.basename(f).split(".cam.h0.")[1].split(".")[0]
                    for f in glob.glob(os.path.join(d, "*.cam.h0.AOD*st*.nc"))})
    for v in cands + found:
        da = open_cesm_var(d, v)
        if da is not None:
            return v, da
    return None, None


def gmean(x, lat):
    w = np.cos(np.deg2rad(lat))
    ok = np.isfinite(x)
    return np.nansum(np.where(ok, x, 0) * w, -1) / np.sum(ok * w, -1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--remap-dir", required=True)
    ap.add_argument("--waccm-dir", required=True)
    ap.add_argument("--waccm-in", help="dir with the prepped T.<period>.nc "
                    "(default: <remap-dir>/WACCM_data)")
    ap.add_argument("--band", type=int, default=4, help="sun band, 1-based")
    ap.add_argument("--png", default="saod_comparison.png")
    a = ap.parse_args()

    sun = xr.open_dataset(os.path.join(a.remap_dir, "ccmi_WACCM_to_ECHAM6_sun.nc"),
                          decode_times=False)
    pf = xr.open_dataset(os.path.join(a.remap_dir, "ccmi_WACCM_to_ECHAM6_p.nc"),
                         decode_times=False)
    tdir = a.waccm_in or os.path.join(a.remap_dir, "WACCM_data")
    tfile = glob.glob(os.path.join(tdir, "T.*.nc"))
    if not tfile:
        sys.exit(f"no T.*.nc in {tdir}; pass --waccm-in")
    T = xr.open_dataset(tfile[0], decode_times=False)["T"].values   # time, lev, lat

    ext = sun["ext_sun"].values[..., a.band - 1]                     # month, alt, lat
    p = pf["pressure"].values                                        # Pa
    lat = sun["latitude"].values
    t = sun["time"].values
    wl = (float(sun["wl1_sun"][a.band - 1]), float(sun["wl2_sun"][a.band - 1])) \
        if "wl1_sun" in sun else ("?", "?")
    if T.shape != p.shape:
        sys.exit(f"T shape {T.shape} != pressure shape {p.shape}")

    # layer thickness [km]: hydrostatic, d(ln p) between neighbouring levels
    dz = R_D * T / G * np.abs(np.gradient(np.log(p), axis=1)) / 1000.0

    # tropopause
    tp = open_cesm_var(a.waccm_dir, "TROP_P")
    if tp is not None:
        tp = tp.sel(time=t, method="nearest").values                 # month, lat
        trop_src = "WACCM TROP_P (zonal mean)"
        tp = np.where(np.isfinite(tp), tp, 100e2)
    else:
        tp = np.full(ext[:, 0, :].shape, 100e2)
        trop_src = "fixed 100 hPa (TROP_P not found)"
    strat = p < tp[:, None, :]

    saod_remap = gmean((ext * dz * strat).sum(1), lat)                # month

    vname, w = find_waccm_saod(a.waccm_dir)
    if w is None:
        print("No WACCM stratospheric AOD variable (AOD*st*) found in "
              f"{a.waccm_dir}; showing REMAP SAOD only.")
        saod_w = None
    else:
        w = w.sel(time=t, method="nearest")
        saod_w = gmean(w.values, w["lat"].values)

    print(f"band {a.band}: {wl[0]:.3f}-{wl[1]:.3f} um; tropopause: {trop_src}")
    print(f"WACCM variable: {vname}\n")
    print(f"{'year':>6} {'REMAP':>8} {'WACCM':>8} {'ratio':>7}")
    years = len(t) // 12
    import re
    m = re.search(r"since\s+(\d{4})", sun["time"].attrs.get("units", ""))
    y0 = int(m.group(1)) + int(t[0] // 365) if m else 0
    for y in range(years):
        sl = slice(12 * y, 12 * y + 12)
        r = saod_remap[sl].mean()
        if saod_w is None:
            print(f"{y0 + y:6d} {r:8.4f}")
        else:
            ww = saod_w[sl].mean()
            print(f"{y0 + y:6d} {r:8.4f} {ww:8.4f} {r / ww if ww else np.nan:7.3f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        x = np.arange(len(t)) / 12
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(x, saod_remap, label=f"REMAP ext_sun band {a.band}, above tropopause")
        if saod_w is not None:
            ax.plot(x, saod_w, label=f"WACCM {vname}")
        ax.set_xlabel("years from start of record")
        ax.set_ylabel("global-mean SAOD")
        ax.legend()
        ax.set_title(f"Tropopause: {trop_src}", fontsize=9)
        fig.tight_layout()
        fig.savefig(a.png, dpi=120)
        print(f"\nplot: {a.png}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
