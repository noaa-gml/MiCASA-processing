#!/bin/sh
cd monthly_1x1

ls MiCASA_v1_flux_x360_y180_monthly_2*.nc|sort|ncrcat -h -O -o MiCASA_v1_flux_x360_y180_monthly.nc

bash ../check_bounds.sh
