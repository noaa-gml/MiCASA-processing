#!/bin/bash
# Symlink ingested vNRT 1° daily files as v1-named files for the same year.
# Used during the NRT phase: ingest_byyear.r writes MiCASA_vNRT_*.nc, but
# downstream consumers (CarbonTracker, etc.) want MiCASA_v1_*.nc.
# When the final v1 release lands, run with MICASA_VERSION=v1 (no-op) or
# delete the symlinks and re-ingest.

set -e

. "$(dirname "$0")/config.sh"

cd "${DAILY_1X1_DIR}"

year="${MICASA_YEAR}"
for mon in {1..12}; do
    monf=$(printf '%02d' "$mon")
    ndays=$(cal "$mon" "$year" | awk 'NF {DAYS = $NF}; END {print DAYS}')
    for day in $(seq 1 "$ndays"); do
        dayf=$(printf '%02d' "$day")
        target="MiCASA_v1_flux_x360_y180_daily_${year}${monf}${dayf}.nc"
        if [ -e "$target" ]; then
            echo "Skipping existing $target"
        else
            src="MiCASA_vNRT_flux_x360_y180_daily_${year}${monf}${dayf}.nc"
            if [ -e "$src" ]; then
                ln -s "$src" "$target"
                echo "Linked $target"
            else
                echo "WARNING source does not exist $src"
            fi
        fi
    done
done
