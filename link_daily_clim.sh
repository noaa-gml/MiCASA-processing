#!/bin/bash
# For each year in $MICASA_CLIM_YEARS (space-separated), symlink missing daily
# files to the day-of-year climatology produced by compute_daily_clim.sh.
# Used when "real" data are unavailable (early/late edges of the record).
#
# Default: 2000 (no ERA5 before this) and $MICASA_YEAR (current/NRT year).

set -e

. "$(dirname "$0")/config.sh"

: "${MICASA_CLIM_YEARS:=2000 ${MICASA_YEAR}}"

prefix="MiCASA_${MICASA_VERSION}_flux_x360_y180_daily"

cd "${DAILY_1X1_DIR}"
for yr in ${MICASA_CLIM_YEARS}; do
    for mon in $(seq 1 12); do
        monf=$(printf '%02d' "$mon")
        ndays=$(cal "$mon" "$yr" | awk 'NF {DAYS = $NF}; END {print DAYS}')
        for day in $(seq 1 "$ndays"); do
            dayf=$(printf '%02d' "$day")
            target="${prefix}_${yr}${monf}${dayf}.nc"
            clim="${prefix}_0000${monf}${dayf}.nc"
            if [ ! -e "$target" ]; then
                ln -s "$clim" "$target"
            else
                echo "File exists $target"
            fi
        done
    done
done
