#!/bin/sh
# Concatenate all per-year-month 1° monthly files into a single time-series file.

set -e

. "$(dirname "$0")/config.sh"

prefix="MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly"

cd "${MONTHLY_1X1_DIR}"
ls ${prefix}_2*.nc | sort | ncrcat -h -O -o "${prefix}.nc"

cd ..
# check_bounds is a sanity print; ncwa hits an NCO chunking bug on the
# concatenated record. Don't let that fail the cat (mirrors the
# || true workaround already used in produce_2025_2026.sh).
bash ./check_bounds.sh || echo "WARN: check_bounds failed (known ncwa chunking issue), continuing"
