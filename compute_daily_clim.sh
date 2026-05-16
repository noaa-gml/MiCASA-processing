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

        # Match real data years (2xxx) -- the literal 0000-prefixed
        # climatology output from a previous run is excluded by the glob.
        # NOTE: the pattern is intentionally UNQUOTED so the shell glob-
        # expands it (nullglob -> empty array if nothing matches). All
        # path components are space-free, so word-splitting is safe.
        files=( ${DAILY_1X1_DIR}/${prefix}_2???${monf}${dayf}.nc )

        # ...but link_daily_clim.sh fills missing days (trailing NRT days,
        # 2025-12-22..31, all of 2026 past the real record, etc.) with
        # symlinks to the 0000MMDD climatology. Those carry a 2xxx-dated
        # NAME, so the 2??? glob matches them -- averaging last run's
        # climatology back in would contaminate the result. Keep only
        # inputs that resolve to a real (non-0000) file.
        real_inputs=()
        for f in "${files[@]}"; do
            case "$(basename "$(readlink -f "$f")")" in
                *_daily_0000*) : ;;                 # clim-fill -- skip
                *)             real_inputs+=( "$f" ) ;;
            esac
        done

        if [ "${#real_inputs[@]}" -eq 0 ]; then
            echo
            echo "ERROR: compute_daily_clim: no real input files for ${monf}${dayf}" >&2
            echo "       (${prefix}_2???${monf}${dayf}.nc matched nothing, or only" >&2
            echo "       clim-fill symlinks). daily_1x1/ is incomplete -- run the" >&2
            echo "       ingest and link_old_micasa_finals.sh steps before this." >&2
            exit 1
        fi

        out="${DAILY_1X1_DIR}/${prefix}_0000${monf}${dayf}.nc"
        # ncea reads the newline-separated input list from stdin.
        printf '%s\n' "${real_inputs[@]}" | ncea -O -o "$out" \
            || { echo; echo "ERROR: ncea failed building $out" >&2; exit 1; }
        echo -n "*"
    done
    echo "]"
done
