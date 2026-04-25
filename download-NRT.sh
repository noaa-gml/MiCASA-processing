#!/bin/sh

# server doesn't provide timestamps, so no -N :/
# only get 2025 data

wget --recursive --no-parent --no-clobber --cut-dirs=6 https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/vNRT/netcdf/daily/2025/
wget --recursive --no-parent --no-clobber --cut-dirs=6 https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/vNRT/netcdf/monthly/2025/

#wget --recursive --no-parent https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/netcdf/monthly/

