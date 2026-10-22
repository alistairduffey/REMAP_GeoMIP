import xarray as xr, numpy as np, glob
T = "composition_test_run"
o = lambda v: xr.open_dataset(f"{T}/{v}/ccmi_WACCM_to_ECHAM6_sun.nc", decode_times=False)
b, f = o("base"), o("forced")
rh = xr.open_dataset(glob.glob(f"{T}/base/WACCM_data/RELHUM.*.nc")[0], decode_times=False).RELHUM

d = (b.ext_sun != f.ext_sun).any("solar_bands")            # (month, altitude, latitude)
lev = b.altitude.values
n_per_lev = d.sum(["month", "latitude"]).values
print("levels with differences [hPa]:", [(round(float(l), 3), int(n)) for l, n in zip(lev, n_per_lev) if n])

r = rh.values.transpose(0, 1, 2)                           # (time, lev, lat)
print("RELHUM at differing points: min %.1f%%, median %.1f%%" %
      (np.nanmin(r[d.values]), np.nanmedian(r[d.values])))
