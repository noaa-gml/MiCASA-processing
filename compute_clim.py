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
be unit-tested standalone -- see tests/test_compute_clim.py.
"""
import os
import sys
import warnings

import numpy as np

FILL = -1.0e34


def modulo_month_mean(values, months):
    """Mean of each calendar month across the time axis.

    values : ndarray, axis 0 is time.
    months : 1-D int array (calendar month 1..12), length == values.shape[0].
    returns: ndarray, axis 0 length 12 (Jan..Dec), other axes preserved.
             NaN / fill values are skipped per cell; a month with no data
             (or an all-missing cell) yields NaN.

    Equivalent to PyFerret's <var>[GT=MONTH_IRREG@MOD] and to
    xarray's `da.groupby("time.month").mean("time")`.
    """
    values = np.asarray(values, dtype="float64")
    months = np.asarray(months)
    if months.shape[0] != values.shape[0]:
        raise ValueError("months length must match values' time axis")
    out = np.full((12,) + values.shape[1:], np.nan, dtype="float64")
    for m in range(1, 13):
        sel = months == m
        if np.any(sel):
            # An all-NaN cell (e.g. ocean) yields NaN via nanmean -- that is
            # the intended result, so silence the "Mean of empty slice"
            # RuntimeWarning rather than letting it spam the log.
            with warnings.catch_warnings(), np.errstate(invalid="ignore"):
                warnings.simplefilter("ignore", RuntimeWarning)
                out[m - 1] = np.nanmean(values[sel], axis=0)
    return out


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
        da.to_dataset().to_netcdf(out, encoding={climvar: {"_FillValue": FILL}})
        finite = np.isfinite(clim)
        print(f"compute_clim: wrote {out}  {climvar}{clim.shape}  "
              f"mean={np.nanmean(clim[finite]):.4e} gC m-2 s-1")
    ds.close()
    print("compute_clim: done")


if __name__ == "__main__":
    main()
