#!/usr/bin/env python
"""
Final step: give the merged REMAP outputs their release names and metadata.

For each scenario <name>, moves
    <run-dir>/out/<name>/ccmi_WACCM_to_ECHAM6_{p,sun,earth}.nc
to
    <run-dir>/final/geomip_WACCM_<name>_to_NorESM2-LM_{p,sun,earth}.nc
and, in place, replaces the inherited title (which says "Surface area
density" on every file) and records the scenario, input directory, period
and processing in global attributes. Data are not touched.

A scenario is only moved if all three files are present. Existing files in
final/ are not overwritten unless --overwrite is given.

  python finalize_outputs.py NAME [NAME ...]
  python finalize_outputs.py --all          # every scenario with outputs in out/
"""
import argparse
import datetime
import glob
import os
import sys

import netCDF4

KINDS = {
    "sun": "Short-wave (solar) optical properties of stratospheric aerosol: "
           "extinction, single-scattering albedo, asymmetry factor",
    "earth": "Long-wave (terrestrial) optical properties of stratospheric "
             "aerosol: extinction, single-scattering albedo, asymmetry factor",
    "p": "Pressure, surface area density and mean radius of stratospheric "
         "aerosol",
}


def read_source(scen_dir):
    info = {}
    path = os.path.join(scen_dir, "source.txt")
    if os.path.exists(path):
        for line in open(path):
            if line.startswith("input directory:"):
                info["input_directory"] = line.split(":", 1)[1].strip()
            elif line.startswith("period:"):
                info["period"] = line.split(":", 1)[1].strip()
    if "period" not in info and os.path.exists(os.path.join(scen_dir, "period.dat")):
        info["period"] = open(os.path.join(scen_dir, "period.dat")).readline().strip()
    return info


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
        src = {k: os.path.join(out_root, name, f"ccmi_WACCM_to_ECHAM6_{k}.nc")
               for k in KINDS}
        dst = {k: os.path.join(final, f"geomip_WACCM_{name}_to_{a.model}_{k}.nc")
               for k in KINDS}
        missing = [k for k in KINDS if not os.path.exists(src[k])]
        if missing:
            already = all(os.path.exists(dst[k]) for k in KINDS)
            problems.append(f"{name}: " + ("already finalized" if already else
                            f"missing {', '.join(missing)} in out/{name}/"))
            continue
        clash = [k for k in KINDS if os.path.exists(dst[k])]
        if clash and not a.overwrite:
            problems.append(f"{name}: final/ already has {', '.join(clash)} "
                            "(use --overwrite)")
            continue

        info = read_source(os.path.join(scen_root, f"WACCM_data_{name}"))
        period = info.get("period", "unknown")
        for k in KINDS:
            os.replace(src[k], dst[k])
            with netCDF4.Dataset(dst[k], "r+") as nc:
                old_title = getattr(nc, "title", "")
                nc.title = (f"{KINDS[k]}, from CESM2-WACCM6 scenario {name}, "
                            f"on {a.model} radiation bands, monthly mean "
                            f"{period[:6]} - {period[7:]}")
                nc.source_model = "CESM2-WACCM6"
                nc.target_model = a.model
                nc.scenario = name
                nc.period = period
                if "input_directory" in info:
                    nc.source_input_directory = info["input_directory"]
                nc.processing = ("REMAP Mie tool (ETH; Luo et al.), zonal-mean "
                                 "WACCM inputs, H2SO4 wt% unit fix applied")
                nc.history = (f"{stamp}: renamed and title/attributes set by "
                              f"finalize_outputs.py (previous title: "
                              f"'{old_title}')\n" + getattr(nc, "history", ""))
        try:
            os.rmdir(os.path.join(out_root, name))
        except OSError:
            pass
        with open(os.path.join(final, "MANIFEST.txt"), "a") as f:
            for k in KINDS:
                f.write(f"{stamp}  {os.path.basename(dst[k])}  period={period}  "
                        f"input={info.get('input_directory', '?')}\n")
        done.append(f"{name} ({period})")

    print(f"Finalized {len(done)} scenario(s) into {final}/")
    for d in done:
        print(f"  ok    {d}")
    for p in problems:
        print(f"  skip  {p}")
    sys.exit(1 if problems and not done else 0)


if __name__ == "__main__":
    main()
