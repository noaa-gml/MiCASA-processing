#!/bin/bash
# Symlink the previous year's raw NCCS portal mirror into this year's tree, so
# we don't need to re-download the entire historical record.
# Range: [$MICASA_YEAR_START, $MICASA_YEAR-1].
#
# Layout note. The directory structure under the previous year has changed
# over the dataset's history:
#   - 2024 tree (legacy): from_weir/portal.nccs.nasa.gov/datashare/gmao/
#                         geos_carb/MiCASA/<version>/netcdf/{daily,monthly}
#   - 2025 tree (current): portal.nccs.nasa.gov/{daily,monthly}
# This script auto-detects which layout the previous year uses and links
# accordingly. New layouts can be added to the `layout_candidates` array.

set -e

. "$(dirname "$0")/config.sh"

prev_year=$((MICASA_YEAR - 1))
prev_root="${BASE_DIR}/${prev_year}/MiCASA_v1"

if [ ! -d "${prev_root}" ]; then
    echo "WARNING: previous-year tree '${prev_root}' does not exist — nothing to link."
    exit 0
fi

# Candidates, in preference order. Each is the directory whose `daily/<YYYY>/`
# and `monthly/<YYYY>/` subdirs we want to link from.
layout_candidates=(
    "${prev_root}/portal.nccs.nasa.gov"
    "${prev_root}/from_weir/portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/${MICASA_VERSION}/netcdf"
)

src_abs=""
for cand in "${layout_candidates[@]}"; do
    if [ -d "${cand}/daily" ] && [ -d "${cand}/monthly" ]; then
        src_abs="${cand}"
        break
    fi
done

if [ -z "${src_abs}" ]; then
    echo "ERROR: no recognized raw layout found under ${prev_root}"
    echo "       checked: ${layout_candidates[*]}"
    exit 1
fi

echo "Linking previous-year raw data from: ${src_abs}"

mkdir -p "${RAW_SRC_DIR}/daily" "${RAW_SRC_DIR}/monthly"

# Use absolute paths for the link targets — they survive moves of the
# WORK_DIR and don't depend on relative-depth assumptions.
for year in $(seq "${MICASA_YEAR_START}" "${prev_year}"); do
    for sub in daily monthly; do
        target="${src_abs}/${sub}/${year}"
        link="${RAW_SRC_DIR}/${sub}/${year}"
        if [ -d "${target}" ]; then
            ln -sfn "${target}" "${link}"
        else
            echo "  (skip) no ${sub}/${year} in source"
        fi
    done
done
