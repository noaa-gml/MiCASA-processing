#!/usr/bin/env python3
"""Sanity-check the concatenated monthly file: print global-mean fluxes.

Reimplements the former NCO-based check_bounds.sh. The old script ran
`ncwa -O src avg` to average every dimension, then scaled the result --
but `ncwa` over the concatenated multi-year record hits an NCO chunking
bug (NC_EINVAL via nco_def_var_deflate), so the check was disabled with
a `|| true` in cat_monthly.sh and never actually ran.

This is a *crude* sanity print, not a science product: it scales the
plain (unweighted, ocean-included) mean of each tracer by the whole
Earth's area -- intentionally rough, matching the old `ncwa` behaviour.
It exists to catch gross regressions (e.g. a units flip), not to report
a precise carbon budget. The faithful per-tracer formula is preserved.
"""
import os
import sys
import numpy as np
import xarray as xr

EARTH_RADIUS    = 6378137.0                  # m
EARTH_AREA      = 4.0 * np.pi * EARTH_RADIUS ** 2
SECONDS_IN_YEAR = 3600.0 * 24.0 * 365.25

work        = os.environ.get("WORK_DIR", os.getcwd())
monthly_dir = os.path.join(work, os.environ.get("MONTHLY_1X1_DIR", "monthly_1x1"))
version     = os.environ.get("MICASA_VERSION", "v1")
src         = os.path.join(monthly_dir,
                           f"MiCASA_{version}_flux_x360_y180_monthly.nc")

if not os.path.exists(src):
    sys.exit(f"check_bounds: concatenated monthly file not found: {src}\n"
             f"             run cat_monthly.sh's ncrcat step first.")

ds = xr.open_dataset(src)
print(f"check_bounds: {os.path.basename(src)} -- crude global-mean fluxes")
for tracer in ("NPP", "Rh", "FIRE", "FUEL"):
    if tracer not in ds:
        print(f"  {tracer}: (variable absent)")
        continue
    # Unweighted mean over every dimension (time, lat, lon) -- matches the
    # old `ncwa -O` with no -a/-w options. NaN/missing are skipped.
    gmean = float(ds[tracer].mean())
    tgy   = gmean * EARTH_AREA * SECONDS_IN_YEAR / 1e15
    print(f"  {tracer}: {tgy:.6g} TgC / year")
ds.close()
