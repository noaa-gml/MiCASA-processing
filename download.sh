#!/bin/sh
# Download raw MiCASA daily + monthly files for $MICASA_YEAR / $MICASA_VERSION
# from the NCCS portal. Replaces the old download.sh + download-NRT.sh pair.
#
# Usage:
#     ./download.sh                                 # uses config.sh defaults
#     MICASA_YEAR=2026 ./download.sh                # different year, same version
#     MICASA_VERSION=vNRT ./download.sh             # near-real-time stream
#     MICASA_VERSION=both ./download.sh             # v1 first, then vNRT
#                                                   # (hybrid v1+NRT stream)
#
# Hybrid v1 + vNRT stream
# -----------------------
# MiCASA publishes two parallel streams from the same upstream pipeline:
#   * v1   — final/authoritative, lags publication of the source data
#   * vNRT — near-real-time, available within days, may be revised once
#            v1 lands and supersedes it
#
# Operational policy: use v1 wherever it exists, fall back to vNRT for
# the trailing window where v1 hasn't been published yet. The two streams
# carry version-tagged basenames (MiCASA_v1_* vs MiCASA_vNRT_*) so they
# co-exist in the same local directory tree:
#
#     portal.nccs.nasa.gov/daily/<YYYY>/<MM>/
#         MiCASA_v1_flux_x3600_y1800_daily_<YYYYMMDD>.nc4   ← preferred
#         MiCASA_vNRT_flux_x3600_y1800_daily_<YYYYMMDD>.nc4 ← fallback
#         …_sha256.txt files (also version-tagged, no collision)
#
# Then ingest_byyear.r picks the v1 file when present, else vNRT
# (controlled via $MICASA_VERSION at ingest time, or via downstream
# `link_vNRT_to_v1.sh` to expose vNRT as v1 for CarbonTracker).
#
# Server doesn't provide timestamps so wget -N is useless;
# --no-clobber gets new files only.

set -e

. "$(dirname "$0")/config.sh"

download_one() {
    local version="$1"
    local base="${PORTAL_URL_BASE}/${version}/netcdf"
    echo "Downloading MiCASA ${version} for ${MICASA_YEAR} from ${base}"
    # --cut-dirs=6 strips the leading /datashare/gmao/geos_carb/MiCASA/<ver>/netcdf/
    # path, so v1 and vNRT both land under the same portal.nccs.nasa.gov/
    # mirror — version is preserved in the filename, not the directory.
    wget --recursive --no-parent --no-clobber --cut-dirs=6 \
         "${base}/daily/${MICASA_YEAR}/"
    wget --recursive --no-parent --no-clobber --cut-dirs=6 \
         "${base}/monthly/${MICASA_YEAR}/"
}

case "${MICASA_VERSION}" in
    both)
        # v1 first, then vNRT (so v1 takes precedence — vNRT only fills
        # in any remaining gaps).
        download_one v1
        download_one vNRT
        ;;
    v1|vNRT)
        download_one "${MICASA_VERSION}"
        ;;
    *)
        echo "ERROR: MICASA_VERSION must be one of {v1, vNRT, both}; got '${MICASA_VERSION}'"
        exit 2
        ;;
esac
