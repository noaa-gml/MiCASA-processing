#!/bin/bash
# Compare ncdump headers of newly-downloaded raw files against a reference
# from the previous year's tree. Used to catch silent metadata/units
# changes at the upstream provider (we got bitten in 2018 by a kg→g flip).
#
# Reference files live at:
#     ${BASE_DIR}/<YEAR>/MiCASA_v1/reference_<YEAR>_{daily,monthly}.nc4
#
# Self-healing: when the diff against last year's reference passes (i.e.
# no metadata changes between last year and this year), this script
# automatically copies this year's first daily + first monthly file as
# the *next* year's reference. That way the chain of references
# bootstraps itself once the initial 2024 reference is in place.
# Set MICASA_NO_BLESS_REFERENCE=1 to skip auto-blessing.
#
# If the previous-year reference is missing, the script prints clear
# instructions and exits 0 (non-fatal — first-time setup of a new tree).

set -e

. "$(dirname "$0")/config.sh"

prev_year=$((MICASA_YEAR - 1))
prev_dir="${BASE_DIR}/${prev_year}/MiCASA_v1"

monthly_ref="${prev_dir}/reference_${prev_year}_monthly.nc4"
new_monthly="${RAW_SRC_DIR}/monthly/${MICASA_YEAR}/MiCASA_${MICASA_VERSION}_flux_x3600_y1800_monthly_${MICASA_YEAR}01.nc4"

daily_ref="${prev_dir}/reference_${prev_year}_daily.nc4"
new_daily="${RAW_SRC_DIR}/daily/${MICASA_YEAR}/01/MiCASA_${MICASA_VERSION}_flux_x3600_y1800_daily_${MICASA_YEAR}0101.nc4"

# Where this year's "blessed" reference will live for next year's run.
next_monthly_ref="${WORK_DIR}/reference_${MICASA_YEAR}_monthly.nc4"
next_daily_ref="${WORK_DIR}/reference_${MICASA_YEAR}_daily.nc4"

any_diff=0
missing_ref=0

run_one() {
    local label="$1"
    local ref="$2"
    local new="$3"
    if [ ! -e "${ref}" ]; then
        cat <<EOF
WARNING: ${label} reference "${ref}" is missing.

  This means the previous-year tree never had a blessed reference file
  installed. Either:
    a) bootstrap by hand:
         cp "${new}" "${ref}"
       (only safe if you trust the new file's metadata!)
    b) run check_unchanged.sh once on a year for which a reference does
       exist, then it will auto-bless forward into ${MICASA_YEAR}.

EOF
        missing_ref=1
        return 0
    fi
    if [ ! -e "${new}" ]; then
        echo "WARNING: ${label} new file \"${new}\" missing — skipping diff"
        missing_ref=1
        return 0
    fi
    echo "=== ${label}: ${ref} vs ${new} ==="
    if diff <(ncdump -h "${ref}") <(ncdump -h "${new}"); then
        echo "  (no header differences)"
    else
        any_diff=1
    fi
    echo
}

run_one monthly "${monthly_ref}" "${new_monthly}"
run_one daily   "${daily_ref}"   "${new_daily}"

# Auto-bless this year's files as next year's reference, but only if we
# actually completed a clean diff (no diffs found, no missing refs).
if [ "$any_diff" -eq 0 ] && [ "$missing_ref" -eq 0 ] \
   && [ "${MICASA_NO_BLESS_REFERENCE:-0}" != "1" ]; then
    if [ -e "${new_monthly}" ] && [ ! -e "${next_monthly_ref}" ]; then
        echo "Blessing ${new_monthly} as ${next_monthly_ref}"
        cp -p "${new_monthly}" "${next_monthly_ref}"
    fi
    if [ -e "${new_daily}" ] && [ ! -e "${next_daily_ref}" ]; then
        echo "Blessing ${new_daily} as ${next_daily_ref}"
        cp -p "${new_daily}" "${next_daily_ref}"
    fi
fi

if [ "$any_diff" -ne 0 ]; then
    echo
    echo "ERROR: header changes detected — STOP and investigate before"
    echo "       running ingest. Examples of what to look for: units"
    echo "       (kg vs g), variable renames, missing-value flips."
    exit 1
fi
