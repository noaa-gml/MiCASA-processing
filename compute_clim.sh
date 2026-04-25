#!/bin/bash
# Build mod-month climatologies of NPP and Rh from the concatenated monthly file.
# Outputs: monthly_1x1/NPPclim.nc, monthly_1x1/Rhclim.nc

set -e

. "$(dirname "$0")/config.sh"

ferret=/apps/other/pyferret-7.6.0/ferret.bash
src="MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly.nc"

cd "${MONTHLY_1X1_DIR}"
for flux in NPP Rh; do
    ${ferret} -server -nojnl -memsize 1000 <<EOF
use ${src}
let/units="gC m-2 s-1" ${flux}clim=${flux}[d=1,gt=month_irreg@mod]
save/clobber/file=${flux}clim.nc ${flux}clim
quit
EOF
done
