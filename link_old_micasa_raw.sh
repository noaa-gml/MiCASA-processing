#/bin/bash

for i in {2001..2023}; do
    ln -s ../../../../2024/MiCASA_v1/from_weir/portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/netcdf/daily/$i portal.nccs.nasa.gov/daily/$i
    ln -s ../../../../2024/MiCASA_v1/from_weir/portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/netcdf/monthly/$i portal.nccs.nasa.gov/monthly/$i
done