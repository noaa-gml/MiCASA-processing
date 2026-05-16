#!/bin/sh
# Concatenate all per-year-month 1° monthly files into a single time-series file.

set -e

. "$(dirname "$0")/config.sh"

prefix="MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly"

cd "${MONTHLY_1X1_DIR}"
ls ${prefix}_2*.nc | sort | ncrcat -h -O -o "${prefix}.nc"

cd ..
# check_bounds is a post-hoc sanity print -- the cat has already
# succeeded by this point, so don't let a check_bounds hiccup abort
# cat_monthly under set -e.
bash ./check_bounds.sh || echo "WARN: check_bounds sanity print failed, continuing"
