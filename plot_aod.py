#!/usr/bin/env python
"""
Approximate 550 nm stratospheric AOD from the finalized REMAP outputs.

  python plot_aod.py                         # every scenario in final/
  python plot_aod.py --names G6-1p5K-SAI INJANN0N --trop none --png aod.png

AOD = sum over levels of ext_sun[band] * dz, with
  ext_sun   final/geomip_WACCM_<name>_to_NorESM2-LM_sun.nc, band 4
            (0.442-0.625 um, the RRTMG band containing 550 nm)  [km-1]
  pressure  final/geomip_WACCM_<name>_to_NorESM2-LM_p.nc        [Pa]
  dz        hydrostatic, R*T/g * d(ln p), with T from the prepped input
            scenarios/WACCM_data_<name>/T.<period>.nc (or an isothermal
            210 K if that file is missing - noted on the plot)

Vertical range (--trop):
  approx  (default) levels above p_tp(lat) = 300 - 215*cos^2(lat) hPa,
          a simple latitude-dependent tropopause (~85 hPa at the equator,
          300 hPa at the poles)
  <hPa>   levels above a fixed pressure, e.g. --trop 100
  none    whole column. Tropospheric values are not meaningful here (all
          modes are treated as sulfuric acid), so this is for comparison only.

Panels: (a) global-mean (cos-lat weighted) AOD, monthly (thin) and annual
mean (thick); (b) zonal-mean AOD averaged over each scenario's final 10
years. Scenarios have different periods, so each uses its own last decade.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_D, G = 287.04, 9.80616
# Display order: feedback-controlled runs, then single-point injections N to S.
ORDER = ["G6-1p5K-SAI", "G6-1p5K-HiLLA", "G6-1p5K-SAI-end_002",
         "INJMAM60N", "INJANN45N", "INJANN30N", "INJANN15N", "INJANN0N",
         "INJANN15S", "INJANN30S", "INJANN45S", "INJSON60S"]
FIXED = {"G6-1p5K-SAI": "black", "G6-1p5K-HiLLA": "tab:green",
         "G6-1p5K-SAI-end_002": "tab:purple"}


def discover(final, model):
    """All scenarios with a _sun file in final/, in display order."""
    found = {os.path.basename(f)[len("geomip_WACCM_"):-len(f"_to_{model}_sun.nc")]
             for f in glob.glob(os.path.join(final, f"geomip_WACCM_*_to_{model}_sun.nc"))}
    return [n for n in ORDER if n in found] + sorted(found - set(ORDER))


def style(name):
    """Colour: fixed for the feedback runs; single-point runs coloured by
    injection latitude (blue = south, red = north)."""
    if name in FIXED:
        return FIXED[name], 2.2
    m = re.search(r"INJ[A-Z]{3}(\d+)([NS])", name)
    if m:
        lat = int(m.group(1)) * (1 if m.group(2) == "N" else -1)
        if lat == 0:
            return "0.4", 1.4     # the colormap's midpoint is too pale to see
        return plt.get_cmap("coolwarm")((lat + 60) / 120), 1.4
    return "0.5", 1.4


def load(name, final, scen_root, band, model):
    fs = os.path.join(final, f"geomip_WACCM_{name}_to_{model}_sun.nc")
    fp = os.path.join(final, f"geomip_WACCM_{name}_to_{model}_p.nc")
    for f in (fs, fp):
        if not os.path.exists(f):
            raise FileNotFoundError(f)
    with netCDF4.Dataset(fs) as d:
        ext = d["ext_sun"][:, :, :, band - 1].astype("f8")   # month, alt, lat
        lat = d["latitude"][:].astype("f8")
        t = d["time"]
        tv, units, cal = t[:], t.units, getattr(t, "calendar", "standard")
        wl = (float(d["wl1_sun"][band - 1]), float(d["wl2_sun"][band - 1])) \
            if "wl1_sun" in d.variables else None
    with netCDF4.Dataset(fp) as d:
        p = d["pressure"][:].astype("f8")                    # Pa

    tfiles = glob.glob(os.path.join(scen_root, f"WACCM_data_{name}", "T.*.nc"))
    if tfiles:
        with netCDF4.Dataset(tfiles[0]) as d:
            T = d["T"][:].astype("f8")                       # time, lev, lat
        if T.shape != p.shape:
            raise ValueError(f"{name}: T {T.shape} vs pressure {p.shape}")
        tsrc = "prepped WACCM T"
    else:
        T = np.full(p.shape, 210.0)
        tsrc = "isothermal 210 K"

    # decimal years; CESM stamps monthly means at the end of the month, so
    # step back half a month to centre each value in its month
    dates = netCDF4.num2date(tv - 15.0, units, cal)
    years = np.array([dt.year + (dt.dayofyr - 0.5) / 365.0 for dt in dates])
    return dict(ext=np.ma.filled(ext, np.nan), p=np.ma.filled(p, np.nan),
                T=np.ma.filled(T, np.nan), lat=lat, years=years,
                tsrc=tsrc, wl=wl)


def aod(d, trop):
    dz = R_D * d["T"] / G * np.abs(np.gradient(np.log(d["p"]), axis=1)) / 1e3
    if trop == "none":
        mask = np.ones_like(d["p"], bool)
    elif trop == "approx":
        ptp = (300.0 - 215.0 * np.cos(np.deg2rad(d["lat"])) ** 2) * 100.0  # Pa
        mask = d["p"] < ptp[None, None, :]
    else:
        mask = d["p"] < float(trop) * 100.0
    return np.nansum(d["ext"] * dz * mask, axis=1)                 # month, lat


def gmean(x, lat):
    w = np.cos(np.deg2rad(lat))
    return (x * w).sum(-1) / w.sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--names", nargs="+",
                    help="scenarios to plot (default: every scenario in final/)")
    ap.add_argument("--final-dir", default="final")
    ap.add_argument("--scen-root", default="scenarios")
    ap.add_argument("--model", default="NorESM2-LM")
    ap.add_argument("--band", type=int, default=4)
    ap.add_argument("--trop", default="approx",
                    help="approx | <hPa> | none  (default approx)")
    ap.add_argument("--years", type=int, default=10,
                    help="length of the final-mean period (default 10)")
    ap.add_argument("--png", default="aod_check.png")
    a = ap.parse_args()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8),
                                   gridspec_kw={"width_ratios": [1.6, 1]})
    notes = set()
    print(f"{'scenario':22s} {'period':>11s} {'final-decade gm AOD':>20s}")
    names = a.names or discover(a.final_dir, a.model)
    if not names:
        sys.exit(f"no geomip_WACCM_*_to_{a.model}_sun.nc files in {a.final_dir}/")
    for name in names:
        try:
            d = load(name, a.final_dir, a.scen_root, a.band, a.model)
        except (FileNotFoundError, ValueError) as e:
            print(f"{name:22s} skipped: {e}")
            continue
        c, lw = style(name)
        tau = aod(d, a.trop)                                       # month, lat
        gm = gmean(tau, d["lat"])
        yrs = d["years"]

        ax1.plot(yrs, gm, color=c, lw=0.5, alpha=0.25)
        nyear = len(gm) // 12
        ann = gm[:12 * nyear].reshape(nyear, 12).mean(1)
        ax1.plot(yrs[:12 * nyear].reshape(nyear, 12).mean(1), ann,
                 color=c, lw=lw, label=name)

        n = min(12 * a.years, len(gm))
        if n < 12 * a.years:
            notes.add(f"{name}: only {n} months available for the final mean")
        zm = tau[-n:].mean(0)
        ax2.plot(d["lat"], zm, color=c, lw=lw)
        if d["tsrc"] != "prepped WACCM T":
            notes.add(f"{name}: dz from {d['tsrc']}")
        print(f"{name:22s} {int(yrs[0])}-{int(yrs[-1]):>4d} "
              f"{gm[-n:].mean():20.4f}")

    wl = d["wl"] if "d" in locals() and d.get("wl") else None
    band_txt = (f"band {a.band}, {wl[0]:.3f}-{wl[1]:.3f} µm" if wl
                else f"band {a.band}")
    trop_txt = {"approx": "above 300-215cos²φ hPa", "none": "whole column"}.get(
        a.trop, f"above {a.trop} hPa")
    ax1.set(xlabel="year", ylabel="global-mean AOD",
            title=f"(a) Global-mean AOD ({band_txt}, {trop_txt})")
    ax2.set(xlabel="latitude", ylabel="AOD", xlim=(-90, 90),
            title=f"(b) Zonal-mean AOD, final {a.years} years of each scenario")
    ax2.set_xticks(range(-90, 91, 30))
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)
        ax.set_ylim(bottom=0)
    ax1.legend(fontsize=7, ncol=2, loc="upper left", frameon=False)
    if notes:
        fig.text(0.01, 0.005, "; ".join(sorted(notes)), fontsize=7)
    fig.tight_layout(rect=(0, 0.03 if notes else 0, 1, 1))
    fig.savefig(a.png, dpi=130)
    print(f"\nplot: {a.png}")
    for n_ in sorted(notes):
        print(f"note: {n_}")


if __name__ == "__main__":
    sys.exit(main())

