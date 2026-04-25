#/bin/bash


EARTH_RADIUS=6378137 # m
PI=3.14159265359

EARTH_AREA=$(bc <<< "4 * $PI * $EARTH_RADIUS ^ 2")
SECONDS_IN_YEAR=$(bc <<< "3600 * 24 * 365.25")

# units should all be gC m^-2 s^-1

ncwa -O ./monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc ./monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_avg.nc

# find NPP = line, find number, 
NPP_gms=$(ncdump ./monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_avg.nc | grep "NPP =" | grep -Eo "[0-9]+\.?[0-9]*[Ee][+-]?[0-9]+")
Rh_gms=$(ncdump ./monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_avg.nc | grep "Rh =" | grep -Eo "[0-9]+\.?[0-9]*[Ee][+-]?[0-9]+")
FIRE_gms=$(ncdump ./monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_avg.nc | grep "FIRE =" | grep -Eo "[0-9]+\.?[0-9]*[Ee][+-]?[0-9]+")
FUEL_gms=$(ncdump ./monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_avg.nc | grep "FUEL =" | grep -Eo "[0-9]+\.?[0-9]*[Ee][+-]?[0-9]+")

# awk understands exponetial notation
NPP_Tgy=$(awk "BEGIN { print $NPP_gms*$EARTH_AREA*$SECONDS_IN_YEAR/1e15 }" )
Rh_Tgy=$(awk "BEGIN { print $Rh_gms*$EARTH_AREA*$SECONDS_IN_YEAR/1e15 }" )
FIRE_Tgy=$(awk "BEGIN { print $FIRE_gms*$EARTH_AREA*$SECONDS_IN_YEAR/1e15 }" )
FUEL_Tgy=$(awk "BEGIN { print $FUEL_gms*$EARTH_AREA*$SECONDS_IN_YEAR/1e15 }" )

echo "NPP: $NPP_Tgy TgC / year"
echo "Rh: $Rh_Tgy TgC / year"
echo "FIRE: $FIRE_Tgy TgC / year"
echo "FUEL: $FUEL_Tgy TgC / year"

