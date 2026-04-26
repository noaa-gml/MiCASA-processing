#!/bin/bash
# Drive the full MiCASA_v1 pipeline for a single year.
#
# Usage:
#     ./run_year.sh YEAR [VERSION] [--skip-download] [--skip-ingest]
#                                  [--skip-aggregate] [--skip-piqs]
#                                  [--skip-diurnalize] [--skip-daysplit]
#                                  [--dry-run]
#
# Example:
#     ./run_year.sh 2026                    # full v1 pipeline for 2026
#     ./run_year.sh 2026 vNRT               # near-real-time stream
#     ./run_year.sh 2026 v1 --skip-download # data already on disk
#
# Stages (matching README flowchart):
#   1. Download           — wget MiCASA daily + monthly for YEAR
#   2. Ingest             — 0.1° → 1° aggregation (daily + monthly)
#   3. Aggregate / clim   — concat, NPP/Rh climatologies, daily climatology
#   4. PIQS fit           — fit smooth seasonal cycles
#   5. Diurnalize         — apply ERA5 hourly meteo to get hourly NEE
#   6. Day-split          — split hourly monthly files to daily NEE files
#
# SBATCH stages (ingest_byyear, diurnalize-ERA5) are submitted with
# --wait so the driver blocks until they complete.

set -e
set -o pipefail

# ---- Args -------------------------------------------------------------------

if [ $# -lt 1 ]; then
    sed -n '2,18p' "$0"
    exit 1
fi

export MICASA_YEAR="$1"
shift

if [ $# -ge 1 ] && [[ "$1" != --* ]]; then
    export MICASA_VERSION="$1"
    shift
fi

skip_download=0; skip_ingest=0; skip_aggregate=0
skip_piqs=0;     skip_diurnalize=0; skip_daysplit=0
dry_run=0

for arg in "$@"; do
    case "$arg" in
        --skip-download)   skip_download=1   ;;
        --skip-ingest)     skip_ingest=1     ;;
        --skip-aggregate)  skip_aggregate=1  ;;
        --skip-piqs)       skip_piqs=1       ;;
        --skip-diurnalize) skip_diurnalize=1 ;;
        --skip-daysplit)   skip_daysplit=1   ;;
        --dry-run)         dry_run=1         ;;
        *) echo "Unknown flag: $arg"; exit 2 ;;
    esac
done

# ---- Config -----------------------------------------------------------------

. "$(dirname "$0")/config.sh"

# Single-year mode for the SBATCH fan-outs.
export MICASA_YEAR_START="${MICASA_YEAR}"
export MICASA_YEAR_END="${MICASA_YEAR}"

cd "${WORK_DIR}"
mkdir -p "${JOBS_DIR}"

# ---- Helpers ----------------------------------------------------------------

run() {
    echo
    echo "==> $*"
    if [ "$dry_run" -eq 1 ]; then return 0; fi
    "$@"
}

# Submit an Rscript via SBATCH and block until it completes. Extra exports
# go in the second arg (comma-separated). When empty, we pass --export=ALL
# alone — sbatch rejects a trailing comma in --export="ALL,".
sbatch_wait() {
    local script="$1"
    local exports="$2"
    local jobname="$3"
    local export_arg="ALL"
    [ -n "$exports" ] && export_arg="ALL,${exports}"
    echo
    echo "==> sbatch --wait $script  (export: ${export_arg})"
    if [ "$dry_run" -eq 1 ]; then return 0; fi
    sbatch --wait \
           -J "${jobname}" \
           --mail-user="${MAIL_USER}" \
           --export="${export_arg}" \
           "${script}"
}

# ---- Banner -----------------------------------------------------------------

echo "========================================================================"
echo "MiCASA pipeline run"
echo "  YEAR        ${MICASA_YEAR}"
echo "  VERSION     ${MICASA_VERSION}"
echo "  WORK_DIR    ${WORK_DIR}"
echo "  MAIL_USER   ${MAIL_USER}"
[ "$dry_run" -eq 1 ] && echo "  *** DRY RUN — nothing will execute ***"
echo "========================================================================"

# ---- Stage 1: Download ------------------------------------------------------

if [ "$skip_download" -eq 0 ]; then
    run sh download.sh
    run Rscript check_daily_downloads.r
    run python check_hashes.py
    run sh check_unchanged.sh
else
    echo "==> [skip] download stage"
fi

# ---- Stage 2: Ingest 0.1° → 1° ---------------------------------------------

if [ "$skip_ingest" -eq 0 ]; then
    sbatch_wait ingest_byyear.r  "INGEST_YEAR=${MICASA_YEAR}" "ingest-${MICASA_YEAR}"
    sbatch_wait ingest_monthly.r "" "ingestmonthly-${MICASA_YEAR}"
else
    echo "==> [skip] ingest stage"
fi

# ---- Stage 3: Aggregate / climatology --------------------------------------

if [ "$skip_aggregate" -eq 0 ]; then
    run sh cat_monthly.sh
    run sh compute_clim.sh
    run sh compute_daily_clim.sh
    run sh link_daily_clim.sh
else
    echo "==> [skip] aggregate stage"
fi

# ---- Stage 4: PIQS fit ------------------------------------------------------

if [ "$skip_piqs" -eq 0 ]; then
    run Rscript write_piqs.r
else
    echo "==> [skip] piqs stage"
fi

# ---- Stage 5: Diurnalize (ERA5 hourly) -------------------------------------

if [ "$skip_diurnalize" -eq 0 ]; then
    sbatch_wait diurnalize-ERA5.r "diurn_year=${MICASA_YEAR}" "d-${MICASA_YEAR}-MiCASA"
else
    echo "==> [skip] diurnalize stage"
fi

# ---- Stage 6: Day-split -----------------------------------------------------

if [ "$skip_daysplit" -eq 0 ]; then
    run sh daysplitter.sh
else
    echo "==> [skip] daysplit stage"
fi

# ---- Stage 7 (optional): vNRT → v1 link ------------------------------------

if [ "${MICASA_VERSION}" = "vNRT" ]; then
    echo
    echo "Tip: when you want CarbonTracker to consume vNRT outputs as v1,"
    echo "     run:  ./link_vNRT_to_v1.sh"
fi

echo
echo "==> Pipeline finished for ${MICASA_YEAR} (${MICASA_VERSION})."
