#!/bin/bash

alias ferret='/apps/other/pyferret-7.6.0/ferret.bash -nojnl -memsize 1000'

cd monthly_1x1
for flux in NPP Rh; do
    ferret -server -nojnl <<EOF
use MiCASA_v1_flux_x360_y180_monthly.nc
!use climatological_axes
let/units="gC m-2 s-1" ${flux}clim=${flux}[d=1,gt=month_irreg@mod]
save/clobber/file=${flux}clim.nc ${flux}clim
quit
EOF
done


