#!/bin/bash
# For each year in $MICASA_CLIM_YEARS (space-separated), symlink missing
# daily files to the day-of-year climatology produced by
# compute_daily_clim.sh. Used when "real" data are unavailable —
# i.e. the year is outside the ERA5 record (2000 and earlier) or has not
# yet been fully published (the present calendar year).
#
# Default: MICASA_CLIM_YEARS = "2000 <current calendar year>".
#
# Why "current calendar year" rather than $MICASA_YEAR:
# Climatology gap-filling is driven by *what's missing on disk right now*
# (no ERA5 before 2000; the present year is still being published), not by
# which year we're processing. If you backfill an earlier year via
# `./run_year.sh 2024`, that year's data is already complete — no clim
# needed for it. Set MICASA_CLIM_YEARS explicitly to override (e.g. for
# upstream outages that left holes in a specific historical year).

set -e

. "$(dirname "$0")/config.sh"

: "${MICASA_CLIM_YEARS:=2000 $(date +%Y)}"

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
