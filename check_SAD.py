import xarray as xr, numpy as np, glob
p = xr.open_dataset("ccmi_WACCM_to_ECHAM6_p.nc", decode_times=False)
w_ = xr.open_dataset(glob.glob("WACCM_data/SAD_AERO.*.nc")[0], decode_times=False)["SAD_AERO"]
print("units  REMAP:", p.SAD.attrs.get("units"), "  WACCM:", w_.attrs.get("units"))
r, w = p.SAD.values, w_.values
pr, lat = p.pressure.values, p.latitude.values
cw = np.cos(np.deg2rad(lat))
strat = (pr < 100e2) & (pr > 1e2)                 # 1-100 hPa
for y in range(0, 50, 5):
    sl = slice(12*y, 12*y + 12)
    a = (r[sl] * strat[sl] * cw).sum(); b = (w[sl] * strat[sl] * cw).sum()
    print(2035 + y, f"REMAP/WACCM SAD ratio (before unit conversion): {a/b:.4g}")
