#!/bin/bash
# Symlink the previous year's raw NCCS portal mirror into this year's tree, so
# we don't need to re-download the entire historical record.
# Range: [$MICASA_YEAR_START, $MICASA_YEAR-1].

set -e

. "$(dirname "$0")/config.sh"

prev_year=$((MICASA_YEAR - 1))
src_rel="../../../../${prev_year}/MiCASA_v1/from_weir/portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/${MICASA_VERSION}/netcdf"

mkdir -p "${RAW_SRC_DIR}/daily" "${RAW_SRC_DIR}/monthly"

for year in $(seq "${MICASA_YEAR_START}" "${prev_year}"); do
    ln -sfn "${src_rel}/daily/${year}"   "${RAW_SRC_DIR}/daily/${year}"
    ln -sfn "${src_rel}/monthly/${year}" "${RAW_SRC_DIR}/monthly/${year}"
done
