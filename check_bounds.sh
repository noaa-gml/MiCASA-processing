#!/bin/bash
# Sanity-check the concatenated monthly file by computing global-mean fluxes
# (in TgC/yr) for NPP, Rh, FIRE, FUEL.

set -e

. "$(dirname "$0")/config.sh"

EARTH_RADIUS=6378137 # m
PI=3.14159265359

EARTH_AREA=$(bc <<< "4 * $PI * $EARTH_RADIUS ^ 2")
SECONDS_IN_YEAR=$(bc <<< "3600 * 24 * 365.25")

src="${MONTHLY_1X1_DIR}/MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly.nc"
avg="${MONTHLY_1X1_DIR}/MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly_avg.nc"

# units should all be gC m^-2 s^-1
ncwa -O "${src}" "${avg}"

for tracer in NPP Rh FIRE FUEL; do
    gms=$(ncdump "${avg}" | grep "${tracer} =" | grep -Eo "[0-9]+\.?[0-9]*[Ee][+-]?[0-9]+")
    Tgy=$(awk "BEGIN { print $gms*$EARTH_AREA*$SECONDS_IN_YEAR/1e15 }")
    echo "${tracer}: ${Tgy} TgC / year"
done
