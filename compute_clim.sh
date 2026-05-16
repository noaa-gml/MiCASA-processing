#!/bin/bash
# Build mod-month climatologies of NPP and Rh from the concatenated monthly
# file. Outputs: monthly_1x1/NPPclim.nc, monthly_1x1/Rhclim.nc
#
# The climatology logic lives in compute_clim.py -- this is a thin wrapper
# that sources config.sh and invokes it. (It was a PyFerret script until
# 2026-05; PyFerret is broken on Orion, see compute_clim.py for detail.)

set -e

. "$(dirname "$0")/config.sh"

${PYTHON:-python3} "$(dirname "$0")/compute_clim.py"
