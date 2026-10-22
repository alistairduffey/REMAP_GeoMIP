#!/usr/bin/env python
"""
Final step: give the merged REMAP outputs their release names and metadata.

For each scenario <name>, writes
    <run-dir>/final/geomip_WACCM_<name>_to_NorESM2-LM_{p,sun,earth}.nc
from
    <run-dir>/out/<name>/ccmi_WACCM_to_ECHAM6_{p,sun,earth}.nc
with an accurate title and global attributes recording the scenario,
period, input directory and processing. Data are copied unchanged.

How (and why this way): netCDF-3 keeps its metadata in a header directly
in front of the data, so lengthening it in place forces the library to
shift the whole data section - slow, and a file interrupted mid-shift is
corrupted. Instead, each file is copied once into a new file whose header
is written with the final metadata from the start, streamed 12 months at a
time, checked against the original, and only then renamed into final/.
The originals in out/<name>/ are deleted only after all three files of the
scenario have been written and checked. An interrupted run leaves the
originals intact and at most a *.tmp file, which the next run removes.

  python finalize_outputs.py NAME [NAME ...]
  python finalize_outputs.py --all          # every scenario with outputs in out/
"""
import argparse
import datetime
import glob
import os
import sys

import numpy as np
import netCDF4

KINDS = {
    "sun": "Short-wave (solar) optical properties of stratospheric aerosol: "
           "extinction, single-scattering albedo, asymmetry factor",
    "earth": "Long-wave (terrestrial) optical properties of stratospheric "
             "aerosol: extinction, single-scattering albedo, asymmetry factor",
    "p": "Pressure, surface area density and mean radius of stratospheric "
         "aerosol",
}
STEP = 12  # months per copy slice


def read_source(scen_dir):
    info = {}
    path = os.path.join(scen_dir, "source.txt")
    if os.path.exists(path):
        for line in open(path):
            if line.startswith("input directory:"):
                info["input_directory"] = line.split(":", 1)[1].strip()
            elif line.startswith("period:"):
                info["period"] = line.split(":", 1)[1].strip()
    pfile = os.path.join(scen_dir, "period.dat")
    if "period" not in info and os.path.exists(pfile):
        info["period"] = open(pfile).readline().strip()
    return info


def copy_with_attrs(src, dst, new_attrs):
    """Copy src to dst (same netCDF format) with global attributes updated."""
    with netCDF4.Dataset(src) as s, \
         netCDF4.Dataset(dst, "w", format=s.data_model) as d:
        s.set_auto_maskandscale(False)
        d.set_fill_off()

        # --- define mode: everything that goes in the header ---------------
        attrs = {k: s.getncattr(k) for k in s.ncattrs()}
        attrs.update(new_attrs)
        d.setncatts(attrs)
        for name, dim in s.dimensions.items():
            d.createDimension(name, None if dim.isunlimited() else len(dim))
        for name, v in s.variables.items():
            fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
            nv = d.createVariable(name, v.dtype, v.dimensions, fill_value=fill)
            nv.setncatts({k: v.getncattr(k) for k in v.ncattrs()
                          if k != "_FillValue"})

        # --- data: streamed in slices along the first dimension ------------
        for name, v in s.variables.items():
            nv = d.variables[name]
            nv.set_auto_maskandscale(False)
            if v.ndim == 0:
                nv[...] = v[...]
            elif v.ndim == 1 or v.shape[0] <= STEP:
                nv[:] = v[:]
            else:
                for i in range(0, v.shape[0], STEP):
                    nv[i:i + STEP] = v[i:i + STEP]


def check_copy(src, dst):
    """Shapes equal for every variable; first and last slices bit-identical."""
    with netCDF4.Dataset(src) as s, netCDF4.Dataset(dst) as d:
        s.set_auto_maskandscale(False)
        d.set_auto_maskandscale(False)
        if set(s.variables) != set(d.variables):
            return "variable lists differ"
        for name, v in s.variables.items():
            w = d.variables[name]
            if v.shape != w.shape or v.dtype != w.dtype:
                return f"{name}: shape/dtype differ"
            if v.ndim == 0:
                same = np.array_equal(v[...], w[...])
            else:
                same = (np.array_equal(v[0], w[0]) and
                        np.array_equal(v[-1], w[-1]))
            if not same:
                return f"{name}: data differ"
    return None


def finalize(name, out_root, final, scen_root, model, overwrite, stamp):
    src = {k: os.path.join(out_root, name, f"ccmi_WACCM_to_ECHAM6_{k}.nc")
           for k in KINDS}
    dst = {k: os.path.join(final, f"geomip_WACCM_{name}_to_{model}_{k}.nc")
           for k in KINDS}
    tmp = {k: dst[k] + ".tmp" for k in KINDS}
    for t in tmp.values():                      # leftovers from an interrupted run
        if os.path.exists(t):
            os.remove(t)

    missing = [k for k in KINDS if not os.path.exists(src[k])]
    if missing:
        if all(os.path.exists(dst[k]) for k in KINDS):
            return None, "already finalized"
        return None, f"missing {', '.join(missing)} in out/{name}/"
    clash = [k for k in KINDS if os.path.exists(dst[k])]
    if clash and not overwrite:
        return None, f"final/ already has {', '.join(clash)} (use --overwrite)"

    info = read_source(os.path.join(scen_root, f"WACCM_data_{name}"))
    period = info.get("period", "unknown")
    try:
        for k in KINDS:
            with netCDF4.Dataset(src[k]) as s:
                old_title = getattr(s, "title", "")
                old_hist = getattr(s, "history", "")
            new = {
                "title": (f"{KINDS[k]}, from CESM2-WACCM6 scenario {name}, on "
                          f"{model} radiation bands, monthly mean "
                          f"{period[:6]} - {period[7:]}"),
                "source_model": "CESM2-WACCM6",
                "target_model": model,
                "scenario": name,
                "period": period,
                "processing": ("REMAP Mie tool (ETH; Luo et al.), zonal-mean "
                               "WACCM inputs, H2SO4 wt% unit fix applied"),
                "history": (f"{stamp}: renamed and title/attributes set by "
                            f"finalize_outputs.py (previous title: "
                            f"'{old_title}')\n" + old_hist),
            }
            if "input_directory" in info:
                new["source_input_directory"] = info["input_directory"]
            print(f"  {name} {k}: copying ...", flush=True)
            copy_with_attrs(src[k], tmp[k], new)
            err = check_copy(src[k], tmp[k])
            if err:
                raise RuntimeError(f"{k}: copy check failed: {err}")
    except Exception as e:
        for t in tmp.values():
            if os.path.exists(t):
                os.remove(t)
        return None, f"FAILED ({e}); originals in out/{name}/ untouched"

    # all three written and checked: publish, then remove the originals
    for k in KINDS:
        os.replace(tmp[k], dst[k])
    for k in KINDS:
        os.remove(src[k])
    try:
        os.rmdir(os.path.join(out_root, name))
    except OSError:
        pass
    with open(os.path.join(final, "MANIFEST.txt"), "a") as f:
        for k in KINDS:
            f.write(f"{stamp}  {os.path.basename(dst[k])}  period={period}  "
                    f"input={info.get('input_directory', '?')}\n")
    return f"{name} ({period})", None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("names", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--run-dir", default=".")
    ap.add_argument("--scen-root", help="prep output root "
                    "(default <run-dir>/scenarios)")
    ap.add_argument("--model", default="NorESM2-LM")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    out_root = os.path.join(a.run_dir, "out")
    final = os.path.join(a.run_dir, "final")
    scen_root = a.scen_root or os.path.join(a.run_dir, "scenarios")
    names = a.names
    if a.all:
        names = sorted(os.path.basename(d.rstrip("/"))
                       for d in glob.glob(os.path.join(out_root, "*/")))
    if not names:
        ap.error("give scenario names or --all")
    os.makedirs(final, exist_ok=True)

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    done, problems = [], []
    for name in names:
        ok, why = finalize(name, out_root, final, scen_root, a.model,
                           a.overwrite, stamp)
        if ok:
            done.append(ok)
        else:
            problems.append(f"{name}: {why}")
        print(f"{'ok  ' if ok else 'skip'}  {ok or name + ': ' + why}", flush=True)

    print(f"\nFinalized {len(done)} scenario(s) into {final}/")
    for p in problems:
        print(f"  skip  {p}")
    sys.exit(1 if any("FAILED" in p for p in problems) else 0)


if __name__ == "__main__":
    main()

