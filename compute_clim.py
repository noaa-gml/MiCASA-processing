#!/usr/bin/env python3
"""Modulo-month climatology of NPP and Rh from the concatenated monthly file.

Reimplements the former pyferret-based compute_clim.sh. PyFerret is broken
on Orion (a NumPy 1.x/2.x ABI mismatch makes `import pyferret` abort), so the
climatology -- the mean of each calendar month across every year in the
record -- is computed here.

Input:
    $MONTHLY_1X1_DIR/MiCASA_<VER>_flux_x360_y180_monthly.nc
        (concatenated multi-year monthly file; NPP/Rh as (time,lat,lon))

Output (consumed by diurnalize-ERA5.r's climatology-fallback branch):
    $MONTHLY_1X1_DIR/NPPclim.nc   variable NPPCLIM (MONTH_IRREG,lat,lon)
    $MONTHLY_1X1_DIR/Rhclim.nc    variable RHCLIM  (MONTH_IRREG,lat,lon)

The variable names, dimension order (month first), and gC m-2 s-1 units
match the old PyFerret output, so diurnalize-ERA5.r picks them up unchanged.

`modulo_month_mean` is a pure-NumPy function (no xarray / netCDF) so it can
be unit-tested standalone -- see tests/test_compute_clim.py. It now lives in
lib/climatology.py, shared with make_climatology_series.py, and is re-exported
here so existing callers and tests are unaffected.
"""
import os
import sys

import numpy as np

# lib/ holds provenance.py (CF/ACDD global attributes) and climatology.py (the
# modulo-month core). Resolve it relative to this file so the import works from
# any working directory.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from provenance import provenance_attrs
from climatology import modulo_month_mean          # noqa: F401  (re-exported)

FILL = -1.0e34


def main():
    import xarray as xr

    work        = os.environ.get("WORK_DIR", os.getcwd())
    monthly_dir = os.path.join(work, os.environ.get("MONTHLY_1X1_DIR", "monthly_1x1"))
    version     = os.environ.get("MICASA_VERSION", "v1")
    src         = os.path.join(monthly_dir,
                               f"MiCASA_{version}_flux_x360_y180_monthly.nc")

    if not os.path.exists(src):
        sys.exit(f"compute_clim: concatenated monthly file not found: {src}\n"
                 f"             run cat_monthly.sh first.")

    ds  = xr.open_dataset(src)
    yrs = ds["time"].dt.year.values
    months = ds["time"].dt.month.values
    trange = f"{int(yrs.min())}-{int(yrs.max())} ({ds.sizes['time']} months)"
    print(f"compute_clim: source {os.path.basename(src)}  [{trange}]")

    for flux, climvar in (("NPP", "NPPCLIM"), ("Rh", "RHCLIM")):
        if flux not in ds:
            sys.exit(f"compute_clim: variable {flux} missing from {src}")
        clim = modulo_month_mean(ds[flux].values, months)
        da = xr.DataArray(
            clim,
            dims=("MONTH_IRREG", "latitude", "longitude"),
            coords={"MONTH_IRREG": np.arange(1, 13),
                    "latitude":  ds["latitude"],
                    "longitude": ds["longitude"]},
            name=climvar,
            attrs={"units": "gC m-2 s-1",
                   "long_name": f"{flux} modulo-month climatology",
                   "climatology_time_range": trange,
                   "missing_value": FILL})
        out = os.path.join(monthly_dir, f"{flux}clim.nc")
        dset = da.to_dataset()
        # CF/ACDD provenance global attributes (lib/provenance.py): producing
        # software + git commit, input file + SHA-256, timestamp, host.
        dset.attrs.update(provenance_attrs(
            step="compute_clim.py", work_dir=work,
            title=f"MiCASA {version} modulo-month climatology of {flux}",
            summary=(f"Mean of each calendar month of MiCASA {flux} across "
                     f"the monthly record {trange}, on a global 1-degree "
                     f"grid. Consumed by diurnalize-ERA5.r as the climatology "
                     f"fallback for months with no published monthly file."),
            inputs={"monthly_flux_series": src},
            extra={"micasa_version": version,
                   "climatology_method": "modulo-month mean (mean per calendar month)",
                   "climatology_time_range": trange}))
        dset.to_netcdf(out, encoding={climvar: {"_FillValue": FILL}})
        finite = np.isfinite(clim)
        print(f"compute_clim: wrote {out}  {climvar}{clim.shape}  "
              f"mean={np.nanmean(clim[finite]):.4e} gC m-2 s-1")
    ds.close()
    print("compute_clim: done")


if __name__ == "__main__":
    main()
