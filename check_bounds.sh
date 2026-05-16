#!/bin/bash
# Sanity-check the concatenated monthly file by printing crude global-mean
# fluxes (NPP, Rh, FIRE, FUEL). Called by cat_monthly.sh.
#
# The logic lives in check_bounds.py -- this is a thin wrapper. It was an
# NCO `ncwa` script until 2026-05; `ncwa` over the concatenated record hits
# an NCO chunking bug (NC_EINVAL), which is why cat_monthly.sh used to wrap
# this call in `|| true`. See check_bounds.py for detail.

set -e

. "$(dirname "$0")/config.sh"

${PYTHON:-python3} "$(dirname "$0")/check_bounds.py"
