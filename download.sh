#!/bin/sh
# Download raw MiCASA daily + monthly files for $MICASA_YEAR / $MICASA_VERSION
# from the NCCS portal. Replaces the old download.sh + download-NRT.sh pair.
#
# Usage:
#     ./download.sh                                  # uses config.sh defaults
#     MICASA_YEAR=2026 ./download.sh                 # different year, same version
#     MICASA_VERSION=vNRT ./download.sh              # near-real-time stream
#
# Server doesn't provide timestamps so -N (timestamping) isn't useful;
# --no-clobber gets new files only.

set -e

. "$(dirname "$0")/config.sh"

base="${PORTAL_URL_BASE}/${MICASA_VERSION}/netcdf"

echo "Downloading MiCASA ${MICASA_VERSION} for ${MICASA_YEAR} from ${base}"

wget --recursive --no-parent --no-clobber --cut-dirs=6 "${base}/daily/${MICASA_YEAR}/"
wget --recursive --no-parent --no-clobber --cut-dirs=6 "${base}/monthly/${MICASA_YEAR}/"
