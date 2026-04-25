#!/bin/bash

# silly to dump changes in two headers

monthly_ref="/home/apera/co2/GFED-CASA/2024/MiCASA_v1/reference_2024_monthly.nc4"
new_monthly="./portal.nccs.nasa.gov/monthly/2024/MiCASA_v1_flux_x3600_y1800_monthly_202401.nc4"

diff <(ncdump -h ${monthly_ref}) <(ncdump -h ${new_monthly})
echo "\n\n"


daily_ref="/home/apera/co2/GFED-CASA/2024/MiCASA_v1/reference_2024_daily.nc4"
new_daily="./portal.nccs.nasa.gov/daily/2024/01/MiCASA_v1_flux_x3600_y1800_daily_20240101.nc4"

diff <(ncdump -h ${daily_ref}) <(ncdump -h ${new_daily})

