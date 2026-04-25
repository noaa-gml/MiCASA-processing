#!/bin/sh
# Concatenate all per-year-month 1° monthly files into a single time-series file.

set -e

. "$(dirname "$0")/config.sh"

prefix="MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly"

cd "${MONTHLY_1X1_DIR}"
ls ${prefix}_2*.nc | sort | ncrcat -h -O -o "${prefix}.nc"

cd ..
bash ./check_bounds.sh
