#!/bin/bash
# Build a day-of-year climatology from all per-year daily 1° files.
# Each output MiCASA_<ver>_flux_x360_y180_daily_0000MMDD.nc averages every
# year's MMDD across the full record.

set -e

. "$(dirname "$0")/config.sh"

prefix="MiCASA_${MICASA_VERSION}_flux_x360_y180_daily"
month_names=("" "Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")

# nullglob: an unmatched pattern expands to an empty array, not the literal
# pattern string. Lets us detect "no input files for this day" cleanly --
# the old `fls=$(ls <pattern>)` form aborted the whole script under `set -e`
# with a cryptic `ls: cannot access` whenever daily_1x1/ was mid-populate.
shopt -s nullglob

for mon in $(seq 1 12); do
    monf=$(printf '%02d' "$mon")

    # 2004 is a leap year -- chosen intentionally so 29 Feb is included.
    ndays=$(cal "$mon" 2004 | awk 'NF {DAYS = $NF}; END {print DAYS}')

    echo -n "${month_names[$mon]} ["
    for day in $(seq 1 "$ndays"); do
        dayf=$(printf '%02d' "$day")

        # Match only real data years (2xxx) -- the 0000-prefixed climatology
        # output from a previous run must NOT be averaged back in.
        # (Bug found 17 July 2025.)
        files=( "${DAILY_1X1_DIR}/${prefix}_2???${monf}${dayf}.nc" )

        if [ "${#files[@]}" -eq 0 ]; then
            echo
            echo "ERROR: compute_daily_clim: no input files match" >&2
            echo "       ${DAILY_1X1_DIR}/${prefix}_2???${monf}${dayf}.nc" >&2
            echo "       daily_1x1/ is incomplete -- run the ingest and" >&2
            echo "       link_old_micasa_finals.sh steps before this one." >&2
            exit 1
        fi

        out="${DAILY_1X1_DIR}/${prefix}_0000${monf}${dayf}.nc"
        # ncea reads the newline-separated input list from stdin.
        printf '%s\n' "${files[@]}" | ncea -O -o "$out" \
            || { echo; echo "ERROR: ncea failed building $out" >&2; exit 1; }
        echo -n "*"
    done
    echo "]"
done
