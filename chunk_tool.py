#!/usr/bin/env python
"""
Split a REMAP input set into time chunks, and merge the chunked outputs.

The REMAP writers treat every month independently, so a long period can be
processed as many short runs in parallel and the outputs concatenated in
time. This script does the bookkeeping; the Fortran is unchanged.

split:  python chunk_tool.py split --src <dir with *.<period>.nc> \
            --period-file <period.dat> --work <workdir> [--months 12]
        Creates <workdir>/chunk_<YYYYMM-YYYYMM>/ for each chunk, holding
        WACCM_data/<name>.<chunk>.nc (time subsets), period.dat, and
        symlinks to bands.dat, lookup_table.dat and data_base/ from the
        current directory.

merge:  python chunk_tool.py merge --work <workdir> --out <dir> [--kinds sun earth]
        Concatenates <workdir>/chunk_*/ccmi_WACCM_to_ECHAM6_<kind>.nc along
        time into <out>/ccmi_WACCM_to_ECHAM6_<kind>.nc and checks the result.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import xarray as xr

PERIOD_RE = re.compile(r"^(\d{4})(\d{2})-(\d{4})(\d{2})$")
SHARED = ["bands.dat", "lookup_table.dat", "data_base"]


def months(period):
    y1, m1, y2, m2 = map(int, PERIOD_RE.match(period).groups())
    out, y, m = [], y1, m1
    while (y, m) <= (y2, m2):
        out.append(f"{y:04d}{m:02d}")
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


def write_nc(ds, path):
    enc = {v: {"_FillValue": None} for v in ds.variables}
    tmp = path + ".tmp"
    ds.to_netcdf(tmp, format="NETCDF3_64BIT", encoding=enc)
    os.replace(tmp, path)


def split(args):
    period = open(args.period_file).readline().strip()
    if not PERIOD_RE.match(period):
        sys.exit(f"bad period '{period}' in {args.period_file}")
    ms = months(period)
    if args.months < 12:
        sys.exit("--months must be >= 12 (the Fortran requires 12 or more)")
    files = sorted(glob.glob(os.path.join(args.src, f"*.{period}.nc")))
    if not files:
        sys.exit(f"no *.{period}.nc files in {args.src}")
    for f in SHARED:
        if not os.path.exists(f):
            sys.exit(f"{f} not found in {os.getcwd()}")

    # chunk boundaries; a short remainder is merged into the last chunk so
    # every chunk has at least 12 months
    starts = list(range(0, len(ms), args.months))
    if len(starts) > 1 and len(ms) - starts[-1] < 12:
        starts.pop()
    bounds = [(s, (starts[i + 1] if i + 1 < len(starts) else len(ms)))
              for i, s in enumerate(starts)]

    datasets = {f: xr.open_dataset(f, decode_times=False) for f in files}
    for ds in datasets.values():
        if ds.sizes["time"] != len(ms):
            sys.exit(f"time length {ds.sizes['time']} != {len(ms)} months")

    os.makedirs(args.work, exist_ok=True)
    for i0, i1 in bounds:
        cp = f"{ms[i0]}-{ms[i1 - 1]}"
        d = os.path.join(args.work, f"chunk_{cp}")
        os.makedirs(os.path.join(d, "WACCM_data"), exist_ok=True)
        for f, ds in datasets.items():
            name = os.path.basename(f).replace(f".{period}.nc", f".{cp}.nc")
            write_nc(ds.isel(time=slice(i0, i1)).load(),
                     os.path.join(d, "WACCM_data", name))
        with open(os.path.join(d, "period.dat"), "w") as fh:
            fh.write(cp + "\n")
        for f in SHARED:
            link = os.path.join(d, f)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.abspath(f), link)
        print(f"{d}: {i1 - i0} months")
    print(f"{len(bounds)} chunks for {period} ({len(ms)} months)")


def merge(args):
    chunks = sorted(glob.glob(os.path.join(args.work, "chunk_*")))
    if not chunks:
        sys.exit(f"no chunk_* directories in {args.work}")
    periods = [os.path.basename(c)[len("chunk_"):] for c in chunks]
    full = f"{periods[0][:6]}-{periods[-1][7:]}"
    nexp = len(months(full))
    # chunks must tile the full period with no gaps or overlaps
    got = sum(len(months(p)) for p in periods)
    if got != nexp:
        sys.exit(f"chunks cover {got} months but {full} has {nexp}")

    os.makedirs(args.out, exist_ok=True)
    ok = True
    for kind in args.kinds:
        fname = f"ccmi_WACCM_to_ECHAM6_{kind}.nc"
        parts = [os.path.join(c, fname) for c in chunks]
        missing = [p for p in parts if not os.path.exists(p)]
        if missing:
            print(f"{kind}: {len(missing)} chunk output(s) missing, e.g. "
                  f"{missing[0]}")
            ok = False
            continue
        dss = [xr.open_dataset(p, decode_times=False) for p in parts]
        tdim = dss[0]["time"].dims[0]  # the writers call it 'month'
        merged = xr.concat(dss, dim=tdim, data_vars="minimal",
                           coords="minimal", compat="override",
                           combine_attrs="override")

        t = merged["time"].values
        nbad = sum(int(np.isnan(merged[v].values).sum())
                   for v in merged.data_vars if tdim in merged[v].dims)
        problems = []
        if merged.sizes[tdim] != nexp:
            problems.append(f"{merged.sizes[tdim]} months, expected {nexp}")
        if np.any(np.diff(t) <= 0):
            problems.append("time not strictly increasing")
        if nbad:
            problems.append(f"{nbad} NaN values")
        if problems:
            print(f"{kind}: " + "; ".join(problems))
            ok = False

        if "title" in merged.attrs:
            merged.attrs["title"] = re.sub(
                r"\d{6} - \d{6}", f"{full[:6]} - {full[7:]}",
                merged.attrs["title"])
        merged.attrs["history"] = (f"concatenated from {len(parts)} chunks "
                                   f"by chunk_tool.py\n"
                                   + merged.attrs.get("history", ""))
        out = os.path.join(args.out, fname)
        write_nc(merged.load(), out)
        print(f"{kind}: {len(parts)} chunks -> {out} ({merged.sizes[tdim]} "
              f"months{'' if not problems else ', WITH PROBLEMS'})")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("--src", required=True)
    s.add_argument("--period-file", required=True)
    s.add_argument("--work", required=True)
    s.add_argument("--months", type=int, default=12)
    m = sub.add_parser("merge")
    m.add_argument("--work", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--kinds", nargs="+", default=["sun", "earth"])
    args = ap.parse_args()
    split(args) if args.cmd == "split" else merge(args)


if __name__ == "__main__":
    main()
