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
    #
    # Purge cached directory listings for $MICASA_YEAR before recursing.
    # wget --no-clobber would otherwise honor a stale index.html and miss
    # any new files NASA published since the last run -- exactly what
    # silently stalls NRT updates. The .nc4 data files themselves stay
    # cached (--no-clobber still skips them); only the index pages refresh.
    find portal.nccs.nasa.gov/daily/${MICASA_YEAR}   -name index.html -delete 2>/dev/null
    find portal.nccs.nasa.gov/monthly/${MICASA_YEAR} -name index.html -delete 2>/dev/null
    #
    # rc=8 (server 4xx) is the expected response when this version hasn't
    # published $MICASA_YEAR yet (e.g. v1 has no 2025/ until late 2026).
    # Treat as a non-fatal "not yet published" so the caller can keep going.
    local rc
    set +e
    wget --recursive --no-parent --no-clobber --cut-dirs=6 \
         "${base}/daily/${MICASA_YEAR}/"
    rc=$?
    if [ $rc -eq 8 ]; then
        echo "  (note: ${version} has not yet published ${MICASA_YEAR}; skipping)"
        set -e
        return 0
    elif [ $rc -ne 0 ]; then
        set -e
        return $rc
    fi
    wget --recursive --no-parent --no-clobber --cut-dirs=6 \
         "${base}/monthly/${MICASA_YEAR}/"
    rc=$?
    set -e
    [ $rc -eq 8 ] && return 0   # monthly may also lag; non-fatal
    return $rc
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

# Verify SHA-256 hashes of downloaded .nc4 files against the per-directory
# aggregate _sha256.txt manifests NCCS publishes alongside. Catches
# partial / corrupt / re-published files immediately rather than at
# ingest time (where a bad raw file would silently produce bad 1°
# aggregates).
#
# For tight NRT loops where you've already verified once in the same hour
# and just want a fast existence check, set MICASA_SKIP_HASH_CHECK=1.
if [ -z "${MICASA_SKIP_HASH_CHECK:-}" ]; then
    echo
    echo "Verifying SHA-256 hashes for ${MICASA_YEAR}..."
    if [ -x "$(dirname "$0")/check_hashes.py" ]; then
        # check_hashes.py reads MICASA_YEAR (set by config.sh) and
        # restricts itself to that year via the y_one branch.
        "$(dirname "$0")/check_hashes.py" || {
            echo "ERROR: SHA-256 verification failed for ${MICASA_YEAR}; see above."
            exit 1
        }
    else
        echo "WARN: check_hashes.py not executable; skipping hash verification."
    fi
fi
