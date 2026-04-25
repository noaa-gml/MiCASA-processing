#!/bin/bash
# Build a day-of-year climatology from all per-year daily 1° files.
# Each output MiCASA_<ver>_flux_x360_y180_daily_0000MMDD.nc averages every
# year's MMDD across the full record.

set -e

. "$(dirname "$0")/config.sh"

prefix="MiCASA_${MICASA_VERSION}_flux_x360_y180_daily"
month_names=("" "Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")

for mon in $(seq 1 12); do
    monf=$(printf '%02d' "$mon")

    # 2004 is a leap year — chosen intentionally so 29 Feb is included.
    ndays=$(cal "$mon" 2004 | awk 'NF {DAYS = $NF}; END {print DAYS}')

    echo -n "${month_names[$mon]} ["
    for day in $(seq 1 "$ndays"); do
        dayf=$(printf '%02d' "$day")

        # Match only real data years (2xxx) — must NOT include the 0000-prefixed
        # climatology output from a previous run, or it'd contaminate the average.
        # (Bug found 17 July 2025.)
        fls=$(ls "${DAILY_1X1_DIR}/${prefix}_2???${monf}${dayf}.nc")

        echo "${fls[*]}" | ncea -O -o "${DAILY_1X1_DIR}/${prefix}_0000${monf}${dayf}.nc"
        echo -n "*"
    done
    echo "]"
done
