#!/usr/bin/env python3
"""Modulo-month climatology of NPP and Rh from the concatenated monthly file.

Reimplements the former pyferret-based compute_clim.sh. PyFerret is broken
on Orion (a NumPy 1.x/2.x ABI mismatch makes `import pyferret` abort), so the
climatology -- which is just the mean of each calendar month across every
year in the record -- is computed here with xarray instead.

Input:
    $MONTHLY_1X1_DIR/MiCASA_<VER>_flux_x360_y180_monthly.nc
        (concatenated multi-year monthly file; NPP/Rh as (time,lat,lon))

Output (consumed by diurnalize-ERA5.r's climatology-fallback branch):
    $MONTHLY_1X1_DIR/NPPclim.nc   variable NPPCLIM (MONTH_IRREG,lat,lon)
    $MONTHLY_1X1_DIR/Rhclim.nc    variable RHCLIM  (MONTH_IRREG,lat,lon)

The variable names, dimension order (month first), and gC m-2 s-1 units
match the old PyFerret output, so diurnalize-ERA5.r picks them up unchanged.
"""
import os
import sys
import numpy as np
import xarray as xr

work        = os.environ.get("WORK_DIR", os.getcwd())
monthly_dir = os.path.join(work, os.environ.get("MONTHLY_1X1_DIR", "monthly_1x1"))
version     = os.environ.get("MICASA_VERSION", "v1")
src         = os.path.join(monthly_dir,
                           f"MiCASA_{version}_flux_x360_y180_monthly.nc")

if not os.path.exists(src):
    sys.exit(f"compute_clim: concatenated monthly file not found: {src}\n"
             f"             run cat_monthly.sh first.")

ds = xr.open_dataset(src)
tt = ds["time"]
yrs = tt.dt.year.values
trange = f"{int(yrs.min())}-{int(yrs.max())} ({tt.sizes['time']} months)"
print(f"compute_clim: source {os.path.basename(src)}  [{trange}]")

FILL = -1.0e34

for flux, climvar in (("NPP", "NPPCLIM"), ("Rh", "RHCLIM")):
    if flux not in ds:
        sys.exit(f"compute_clim: variable {flux} missing from {src}")
    # Modulo-month mean: average every occurrence of each calendar month.
    # Equivalent to PyFerret's <var>[GT=MONTH_IRREG@MOD].
    clim = (ds[flux].groupby("time.month").mean("time")
                    .rename({"month": "MONTH_IRREG"})
                    .rename(climvar)
                    .astype("float64"))
    clim.attrs = {
        "units": "gC m-2 s-1",
        "long_name": f"{flux} modulo-month climatology",
        "climatology_time_range": trange,
        "missing_value": FILL,
    }
    out = os.path.join(monthly_dir, f"{flux}clim.nc")
    clim.to_dataset().to_netcdf(
        out, encoding={climvar: {"_FillValue": FILL}})
    finite = np.isfinite(clim.values)
    print(f"compute_clim: wrote {out}  {climvar}{tuple(clim.shape)}  "
          f"mean={np.nanmean(clim.values[finite]):.4e} gC m-2 s-1")

ds.close()
print("compute_clim: done")
