#!/bin/bash
# Build mod-month climatologies of NPP and Rh from the concatenated monthly file.
# Outputs: monthly_1x1/NPPclim.nc, monthly_1x1/Rhclim.nc

set -e

. "$(dirname "$0")/config.sh"

ferret=/apps/other/pyferret-7.6.0/ferret.bash
src="MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly.nc"

cd "${MONTHLY_1X1_DIR}"
# NPP/Rh: ferret path. As of 2026 pyferret's numpy ABI is broken on Orion
# and these calls actually fail; in our v2 tree NPPclim.nc and Rhclim.nc
# are symlinks to the v1 build (which had a working ferret at the time).
# Left here for documentation / future restoration.
for flux in NPP Rh; do
    ${ferret} -server -nojnl -memsize 1000 <<EOF || echo "WARN: ferret ${flux}clim failed (see above); relying on existing ${flux}clim.nc"
use ${src}
let/units="gC m-2 s-1" ${flux}clim=${flux}[d=1,gt=month_irreg@mod]
save/clobber/file=${flux}clim.nc ${flux}clim
quit
EOF
done

# ATMC: skip ferret entirely, use a small Python helper. ATMC was added
# 2026-04-29 (proposal #2) so we never had a working ferret-built
# ATMCclim.nc to fall back to. mod-month mean over the cat'd monthly file.
/work2/noaa/co2/miniconda3/envs/tm5/bin/python -u - <<'PYEOF2'
import netCDF4 as nc, numpy as np, cftime
src = "MiCASA_v1_flux_x360_y180_monthly.nc"
ds = nc.Dataset(src)
atmc = ds.variables["ATMC"][:]
months = np.array([d.month for d in cftime.num2pydate(
    ds.variables["time"][:], ds.variables["time"].units)])
clim = np.zeros((12, atmc.shape[1], atmc.shape[2]), dtype=np.float32)
for m in range(1, 13):
    clim[m-1] = np.nanmean(atmc[months == m], axis=0)
out = nc.Dataset("ATMCclim.nc", "w", format="NETCDF4")
out.createDimension("longitude", atmc.shape[2])
out.createDimension("latitude",  atmc.shape[1])
out.createDimension("time", 12)
out.createVariable("longitude", "f4", ("longitude",))[:] = ds.variables["longitude"][:]
out.createVariable("latitude",  "f4", ("latitude",))[:]  = ds.variables["latitude"][:]
out.createVariable("time",      "f8", ("time",))[:]      = list(range(1, 13))
v = out.createVariable("ATMCCLIM", "f4", ("time","latitude","longitude"),
                       zlib=True, complevel=4)
v[:] = clim
v.units = "gC m-2 s-1"
v.long_name = "Atmospheric correction climatology (mod-month mean)"
out.history = "Computed by compute_clim.sh / Python (pyferret broken in 2026)"
out.close(); ds.close()
print("wrote ATMCclim.nc")
PYEOF2
