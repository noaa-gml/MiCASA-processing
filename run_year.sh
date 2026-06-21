#!/bin/bash
# Drive the full MiCASA_v1 pipeline for a single year.
#
# Usage:
#     ./run_year.sh YEAR [VERSION] [--fitter NAME] [--skip-download]
#                                  [--skip-ingest] [--skip-aggregate]
#                                  [--skip-piqs] [--skip-diurnalize]
#                                  [--skip-daysplit] [--dry-run]
#
#   --fitter NAME  sub-monthly smoother for stage 4 (default: pchip).
#                  pchip (V2 default) | piqs (legacy V1; global solve, auto-sets
#                  MICASA_PIQS_PAD_RIGHT=2) | ppm | linmm | mss | atpk. All write
#                  fit.piqs.rda, which the diurnalize stage reads.
#
# Example:
#     ./run_year.sh 2026                      # full v2 pipeline (PCHIP) for 2026
#     ./run_year.sh 2026 vNRT                 # near-real-time stream
#     ./run_year.sh 2026 v1 --skip-download   # data already on disk
#     ./run_year.sh 2026 v1 --fitter piqs     # PIQS-fitted, V1-style product
#
# Stages (matching README flowchart):
#   1. Download           — wget MiCASA daily + monthly for YEAR
#   2. Ingest             — 0.1° → 1° aggregation (daily + monthly)
#   3. Aggregate / clim   — concat, NPP/Rh climatologies, daily climatology
#   4. Fitter             — fit smooth seasonal cycles (--fitter; default pchip)
#   5. Diurnalize         — apply ERA5 hourly meteo to get hourly NEE
#   6. Day-split          — split hourly monthly files to daily NEE files
#
# SBATCH stages (ingest_byyear, diurnalize-ERA5) are submitted with
# --wait so the driver blocks until they complete.

set -e
set -o pipefail

# ---- Args -------------------------------------------------------------------

if [ $# -lt 1 ]; then
    sed -n '2,27p' "$0"
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
fitter=pchip          # production default; --fitter swaps the stage-4 smoother
expect_fitter=0       # set when the previous token was a bare "--fitter"

for arg in "$@"; do
    if [ "$expect_fitter" -eq 1 ]; then fitter="$arg"; expect_fitter=0; continue; fi
    case "$arg" in
        --skip-download)   skip_download=1   ;;
        --skip-ingest)     skip_ingest=1     ;;
        --skip-aggregate)  skip_aggregate=1  ;;
        --skip-piqs)       skip_piqs=1       ;;
        --skip-diurnalize) skip_diurnalize=1 ;;
        --skip-daysplit)   skip_daysplit=1   ;;
        --fitter)          expect_fitter=1   ;;   # "--fitter piqs"
        --fitter=*)        fitter="${arg#*=}" ;;  # "--fitter=piqs"
        --dry-run)         dry_run=1         ;;
        *) echo "Unknown flag: $arg"; exit 2 ;;
    esac
done
if [ "$expect_fitter" -eq 1 ]; then echo "--fitter needs a value"; exit 2; fi

# Map --fitter to its writer. All writers emit fit.piqs.rda (recording the
# fitter in piqsfit.meta), which the diurnalize stage reads by default.
case "$fitter" in
    pchip) fitter_script=write_pchip.r ;;
    piqs)  fitter_script=write_piqs.r
           # PIQS is a global solve over the whole record (proposal #1, #17):
           # pad the trailing edge for NRT stability. NOTE: a PIQS refit
           # rewrites every historical month — re-diurnalize the affected tail.
           export MICASA_PIQS_PAD_RIGHT="${MICASA_PIQS_PAD_RIGHT:-2}" ;;
    ppm)   fitter_script=write_ppm.r   ;;
    linmm) fitter_script=write_linmm.r ;;
    mss)   fitter_script=write_mss.r   ;;
    atpk)  fitter_script=write_atpk.r  ;;
    *) echo "Unknown --fitter: '$fitter' (pchip|piqs|ppm|linmm|mss|atpk)"; exit 2 ;;
esac

# ---- Config -----------------------------------------------------------------

. "$(dirname "$0")/config.sh"
. "$(dirname "$0")/lib/manifest.sh"
trap 'manifest_record run_year.sh fail - "aborted (line $LINENO)"' ERR

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
    # Manifest step name: the script, not the interpreter wrapping it.
    case "$1" in
        sh|bash|Rscript|python|python3) _rs_step=$(basename "$2") ;;
        *)                              _rs_step=$(basename "$1") ;;
    esac
    _rs_t0=$(date +%s)
    if "$@"; then
        manifest_record "$_rs_step" ok "$(($(date +%s) - _rs_t0))" "run_year: $*"
    else
        _rs_rc=$?
        manifest_record "$_rs_step" fail "$(($(date +%s) - _rs_t0))" \
            "run_year: $* (exit $_rs_rc)"
        return $_rs_rc
    fi
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
    local t0; t0=$(date +%s)
    if sbatch --wait \
              -J "${jobname}" \
              --mail-user="${MAIL_USER}" \
              --export="${export_arg}" \
              "${script}"; then
        manifest_record "$(basename "$script")" ok "$(($(date +%s) - t0))" \
            "run_year sbatch: $jobname"
    else
        local rc=$?
        manifest_record "$(basename "$script")" fail "$(($(date +%s) - t0))" \
            "run_year sbatch: $jobname (exit $rc)"
        return $rc
    fi
}

# ---- Banner -----------------------------------------------------------------

echo "========================================================================"
echo "MiCASA pipeline run"
echo "  YEAR        ${MICASA_YEAR}"
echo "  VERSION     ${MICASA_VERSION}"
echo "  FITTER      ${fitter}  (stage 4: ${fitter_script})"
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

# ---- Stage 4: fitter (sub-monthly smoother; --fitter, default pchip) --------

if [ "$skip_piqs" -eq 0 ]; then
    run Rscript "$fitter_script"
else
    echo "==> [skip] fitter stage ($fitter)"
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

manifest_record run_year.sh ok - "year ${MICASA_YEAR} version ${MICASA_VERSION}"
echo
echo "==> Pipeline finished for ${MICASA_YEAR} (${MICASA_VERSION})."
