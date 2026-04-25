#!/bin/bash
# Compare ncdump headers of newly-downloaded raw files against a reference
# from the previous year's tree. Used to catch silent metadata/units changes
# at the upstream provider (we got bitten in 2018 by a kg→g flip).

set -e

. "$(dirname "$0")/config.sh"

prev_year=$((MICASA_YEAR - 1))
prev_dir="${BASE_DIR}/${prev_year}/MiCASA_v1"

monthly_ref="${prev_dir}/reference_${prev_year}_monthly.nc4"
new_monthly="${RAW_SRC_DIR}/monthly/${MICASA_YEAR}/MiCASA_${MICASA_VERSION}_flux_x3600_y1800_monthly_${MICASA_YEAR}01.nc4"

daily_ref="${prev_dir}/reference_${prev_year}_daily.nc4"
new_daily="${RAW_SRC_DIR}/daily/${MICASA_YEAR}/01/MiCASA_${MICASA_VERSION}_flux_x3600_y1800_daily_${MICASA_YEAR}0101.nc4"

for pair in "monthly:${monthly_ref}:${new_monthly}" "daily:${daily_ref}:${new_daily}"; do
    label="${pair%%:*}"; rest="${pair#*:}"
    ref="${rest%%:*}"; new="${rest#*:}"
    if [ ! -e "${ref}" ]; then
        echo "WARNING: ${label} reference \"${ref}\" missing — skipping diff"
        continue
    fi
    if [ ! -e "${new}" ]; then
        echo "WARNING: ${label} new file \"${new}\" missing — skipping diff"
        continue
    fi
    echo "=== ${label}: ${ref} vs ${new} ==="
    diff <(ncdump -h "${ref}") <(ncdump -h "${new}") || true
    echo
done
