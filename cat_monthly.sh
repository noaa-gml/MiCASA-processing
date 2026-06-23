#!/bin/sh
# Concatenate all per-year-month 1° monthly files into a single time-series file.

set -e

. "$(dirname "$0")/config.sh"

prefix="MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly"

cd "${MONTHLY_1X1_DIR}"
ls ${prefix}_2*.nc | sort | ncrcat -h -O -o "${prefix}.nc"

# Helpers live in the checkout, referenced by absolute path: the cwd is now the
# monthly output dir, which may be a redirected absolute location (an env-
# overridden MONTHLY_1X1_DIR), so a relative "./stamp_provenance.py" would not
# resolve. (Was: "cd .." + "./helper", which only worked when the output dir
# sat directly under the checkout.)
WORK_DIR="${WORK_DIR:-$(cd "$(dirname "$0")" && pwd)}"

# Stamp CF/ACDD provenance onto the concatenated file (lib/provenance.py via
# stamp_provenance.py): producing software, git commit, host, timestamp.
# Non-fatal -- the concatenation itself has already succeeded.
${PYTHON:-python3} "${WORK_DIR}/stamp_provenance.py" \
    "${MONTHLY_1X1_DIR}/${prefix}.nc" \
    --step cat_monthly.sh \
    --title "MiCASA ${MICASA_VERSION} monthly land carbon flux time series" \
    --summary "Concatenation of the per-month 1-degree MiCASA monthly files (NPP, Rh, FIRE, FUEL) into a single multi-year time series." \
    || echo "WARN: stamp_provenance on ${prefix}.nc failed, continuing"

# check_bounds is a post-hoc sanity print -- the cat has already
# succeeded by this point, so don't let a check_bounds hiccup abort
# cat_monthly under set -e.
bash "${WORK_DIR}/check_bounds.sh" || echo "WARN: check_bounds sanity print failed, continuing"
