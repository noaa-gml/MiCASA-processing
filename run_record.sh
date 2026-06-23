#!/bin/bash
# Drive the full MiCASA pipeline over a MULTI-YEAR range into the configured
# output dirs. Generalizes the year-pinned produce_2025_2026.sh and the manual
# full-record orchestration used to build /work2/.../2026/MiCASA.0.
#
# Usage:
#     ./run_record.sh YEAR_START YEAR_END [options]
#
# Options:
#   --v1-through YEAR   v1 raw is available through YEAR; years after it are
#                       ingested as vNRT and then symlinked to v1 names so the
#                       fit/diurnalize/daysplit stages see one continuous v1
#                       stream. Default: YEAR_END (treat the whole range as v1).
#   --fitter NAME       stage-4 smoother: pchip (default) | piqs | ppm | linmm |
#                       mss | atpk. (Same mapping as run_year.sh.)
#   --download          run the download stage first (loops years × version).
#                       OFF by default — assumes the raw 0.1° tree is present.
#                       NOTE: only Orion login/dtn nodes have internet.
#   --skip-ingest       skip stage 2/3 (dailies + monthly ingest)
#   --skip-link         skip the vNRT→v1 symlinking (stage 4)
#   --skip-aggregate    skip stage 5 (cat_monthly / clim / link_daily_clim)
#   --skip-fit          skip stage 6 (fitter)
#   --skip-diurnalize   skip stage 7
#   --skip-daysplit     skip stage 8
#   --dry-run           print the plan; execute nothing
#
# Examples:
#   ./run_record.sh 2001 2024                       # full v1 historical record
#   ./run_record.sh 2001 2026 --v1-through 2024     # v1 2001-24 + vNRT 2025-26
#   ./run_record.sh 2020 2024 --skip-ingest         # re-fit/re-diurnalize only
#
# Honors the standard environment (config.sh): WORK_DIR, MAIL_USER, BASE_DIR,
# and the output-layout overrides (DAILY_1X1_DIR, MONTHLY_1X1_DIR, ERA5_DIR,
# JOBS_DIR, ...) so the whole product can be redirected to an absolute path.
# compute_clim needs a Python with xarray — set PYTHON to point at it
# (e.g. PYTHON=/work2/noaa/co2/miniconda3/envs/sci/bin/python3).
#
# Stages (matching run_year.sh, but multi-year + version-split):
#   1. Download         (opt-in) wget raw daily+monthly for each year/version
#   2. Ingest dailies   0.1°→1° daily, one SBATCH per year (per version group)
#   3. Ingest monthly   0.1°→1° monthly, one SBATCH per version group (range)
#   4. Link vNRT→v1     symlink vNRT-named daily+monthly files as v1 (NRT years)
#   5. Aggregate / clim cat_monthly, NPP/Rh clim, daily clim, daily-clim link
#   6. Fit              fit smooth seasonal cycles (--fitter; default pchip)
#   7. Diurnalize       ERA5 hourly meteo → hourly NEE, one SBATCH per year
#   8. Day-split        split hourly monthly files to daily NEE files
#   9. Provenance       write PROVENANCE.txt (streams, config, git, time)
#
# The fan-out stages (2, 7) submit one `sbatch --wait` per year in the
# background and block until ALL years finish — a failure in any year aborts
# the driver. Run this under nohup/screen on a login node; it blocks for the
# duration of the whole record.

set -euo pipefail

# ---- Args -------------------------------------------------------------------

if [ $# -lt 2 ]; then
    sed -n '2,46p' "$0"
    exit 1
fi

YEAR_START="$1"; shift
YEAR_END="$1";   shift

v1_through="$YEAR_END"   # default: whole range is v1
fitter=pchip
do_download=0
skip_ingest=0; skip_link=0; skip_aggregate=0
skip_fit=0;    skip_diurnalize=0; skip_daysplit=0
dry_run=0
expect_v1through=0; expect_fitter=0

for arg in "$@"; do
    if [ "$expect_v1through" -eq 1 ]; then v1_through="$arg"; expect_v1through=0; continue; fi
    if [ "$expect_fitter"    -eq 1 ]; then fitter="$arg";     expect_fitter=0;    continue; fi
    case "$arg" in
        --v1-through)      expect_v1through=1 ;;
        --v1-through=*)    v1_through="${arg#*=}" ;;
        --fitter)          expect_fitter=1 ;;
        --fitter=*)        fitter="${arg#*=}" ;;
        --download)        do_download=1 ;;
        --skip-ingest)     skip_ingest=1 ;;
        --skip-link)       skip_link=1 ;;
        --skip-aggregate)  skip_aggregate=1 ;;
        --skip-fit)        skip_fit=1 ;;
        --skip-diurnalize) skip_diurnalize=1 ;;
        --skip-daysplit)   skip_daysplit=1 ;;
        --dry-run)         dry_run=1 ;;
        *) echo "Unknown flag: $arg"; exit 2 ;;
    esac
done
[ "$expect_v1through" -eq 1 ] && { echo "--v1-through needs a value"; exit 2; }
[ "$expect_fitter"    -eq 1 ] && { echo "--fitter needs a value";     exit 2; }

# Map --fitter to its writer (same as run_year.sh). All writers emit
# fit.piqs.rda, which the diurnalize stage reads by default.
case "$fitter" in
    pchip) fitter_script=write_pchip.r ;;
    piqs)  fitter_script=write_piqs.r
           export MICASA_PIQS_PAD_RIGHT="${MICASA_PIQS_PAD_RIGHT:-2}" ;;
    ppm)   fitter_script=write_ppm.r   ;;
    linmm) fitter_script=write_linmm.r ;;
    mss)   fitter_script=write_mss.r   ;;
    atpk)  fitter_script=write_atpk.r  ;;
    *) echo "Unknown --fitter: '$fitter' (pchip|piqs|ppm|linmm|mss|atpk)"; exit 2 ;;
esac

# ---- Version split ----------------------------------------------------------
# v1 group:   [YEAR_START .. min(v1_through, YEAR_END)]
# vNRT group: [v1_through+1 .. YEAR_END]   (empty when v1_through >= YEAR_END)
years_v1=""; years_vnrt=""
for y in $(seq "$YEAR_START" "$YEAR_END"); do
    if [ "$y" -le "$v1_through" ]; then years_v1="${years_v1} $y"
    else                               years_vnrt="${years_vnrt} $y"; fi
done
years_v1="${years_v1# }"; years_vnrt="${years_vnrt# }"
years_all="$(seq "$YEAR_START" "$YEAR_END" | tr '\n' ' ')"; years_all="${years_all% }"
# Everything after the vNRT→v1 link is a single continuous v1 stream.
POST_VERSION=v1

# ---- Config -----------------------------------------------------------------

. "$(dirname "$0")/config.sh"
. "$(dirname "$0")/lib/manifest.sh"
trap 'manifest_record run_record.sh fail - "aborted (line $LINENO)"' ERR

cd "${WORK_DIR}"
mkdir -p "${JOBS_DIR}"

LOG="${JOBS_DIR}/run_record.$(date -u +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

# ---- Helpers ----------------------------------------------------------------

# Inline step (mirrors run_year.sh): run a command, record it in the manifest.
run() {
    echo; echo "==> $*"
    if [ "$dry_run" -eq 1 ]; then return 0; fi
    case "$1" in
        sh|bash|Rscript|python|python3) _rs_step=$(basename "$2") ;;
        *)                              _rs_step=$(basename "$1") ;;
    esac
    local t0; t0=$(date +%s)
    if "$@"; then
        manifest_record "$_rs_step" ok "$(($(date +%s) - t0))" "run_record: $*"
    else
        local rc=$?
        manifest_record "$_rs_step" fail "$(($(date +%s) - t0))" "run_record: $* (exit $rc)"
        return $rc
    fi
}

# Submit one SBATCH script and block until it completes (mirrors run_year.sh).
# $2 = comma-separated extra exports (empty → --export=ALL alone).
sbatch_wait() {
    local script="$1" exports="$2" jobname="$3"
    local export_arg="ALL"; [ -n "$exports" ] && export_arg="ALL,${exports}"
    echo; echo "==> sbatch --wait $script  (export: ${export_arg})"
    if [ "$dry_run" -eq 1 ]; then return 0; fi
    local t0; t0=$(date +%s)
    if sbatch --wait -J "${jobname}" --mail-user="${MAIL_USER}" \
              --export="${export_arg}" "${script}"; then
        manifest_record "$(basename "$script")" ok "$(($(date +%s) - t0))" "run_record sbatch: $jobname"
    else
        local rc=$?
        manifest_record "$(basename "$script")" fail "$(($(date +%s) - t0))" "run_record sbatch: $jobname (exit $rc)"
        return $rc
    fi
}

# Fan a per-year SBATCH script over a list of years and block until ALL finish.
#   fanout_wait SCRIPT YEAR_ENV "y1 y2 .." EXTRA_EXPORTS LABEL
# Each year is one `sbatch --wait` in the background; a failure in any year
# makes the call return non-zero (→ ERR trap aborts the driver).
fanout_wait() {
    local script="$1" yearenv="$2" years="$3" extra="$4" label="$5"
    echo; echo "==> fan-out $script over years: ${years}  (${yearenv}=<y>${extra:+,${extra}})"
    if [ "$dry_run" -eq 1 ]; then return 0; fi
    local t0; t0=$(date +%s)
    local pids=() y exp
    for y in $years; do
        exp="ALL,WORK_DIR=${WORK_DIR},${yearenv}=${y}"
        [ -n "$extra" ] && exp="${exp},${extra}"
        sbatch --wait -J "${label}-${y}" --mail-user="${MAIL_USER}" \
               --export="${exp}" "${script}" &
        pids+=("$!")
    done
    local fail=0 p
    for p in "${pids[@]}"; do wait "$p" || fail=1; done
    if [ "$fail" -eq 0 ]; then
        manifest_record "$(basename "$script")" ok "$(($(date +%s) - t0))" "run_record fanout: ${label} (${years})"
    else
        manifest_record "$(basename "$script")" fail "$(($(date +%s) - t0))" "run_record fanout: ${label} (a year failed)"
        return 1
    fi
}

# Symlink the vNRT-named monthly files of a year as v1-named (the daily side is
# handled by link_vNRT_to_v1.sh; the monthly side has no dedicated script).
link_vnrt_monthlies() {
    local y="$1"
    ( cd "${MONTHLY_1X1_DIR}" || exit 1
      local v target
      for v in MiCASA_vNRT_flux_x360_y180_monthly_"${y}"??.nc; do
          [ -e "$v" ] || continue
          target="${v/vNRT/v1}"
          if [ -e "$target" ] || [ -L "$target" ]; then rm -f "$target"; fi
          ln -s "$v" "$target"
          echo "  linked $target -> $v"
      done )
}

# ---- Banner -----------------------------------------------------------------

echo "========================================================================"
echo "MiCASA full-record run"
echo "  RANGE        ${YEAR_START}-${YEAR_END}"
echo "  v1 years     ${years_v1:-<none>}"
echo "  vNRT years   ${years_vnrt:-<none>}  (→ symlinked to v1)"
echo "  FITTER       ${fitter}  (stage 6: ${fitter_script})"
echo "  WORK_DIR     ${WORK_DIR}"
echo "  MONTHLY_DIR  ${MONTHLY_1X1_DIR}"
echo "  ERA5_DIR     ${ERA5_DIR}"
echo "  MAIL_USER    ${MAIL_USER}"
echo "  LOG          ${LOG}"
[ "$dry_run" -eq 1 ] && echo "  *** DRY RUN — nothing will execute ***"
echo "========================================================================"

# ---- Stage 1: Download (opt-in) --------------------------------------------

if [ "$do_download" -eq 1 ]; then
    for y in $years_v1;   do MICASA_YEAR="$y" MICASA_VERSION=v1   run sh download.sh; done
    for y in $years_vnrt; do MICASA_YEAR="$y" MICASA_VERSION=vNRT run sh download.sh; done
else
    echo; echo "==> [skip] download stage (assuming raw present; use --download to fetch)"
fi

# ---- Stage 2: Ingest dailies (0.1° → 1°, per-year fan-out) -----------------

if [ "$skip_ingest" -eq 0 ]; then
    [ -n "$years_v1" ]   && fanout_wait ingest_byyear.r INGEST_YEAR "$years_v1"   "MICASA_VERSION=v1"   "ingest-v1"
    [ -n "$years_vnrt" ] && fanout_wait ingest_byyear.r INGEST_YEAR "$years_vnrt" "MICASA_VERSION=vNRT" "ingest-vNRT"

    # ---- Stage 3: Ingest monthly (per version group, range loop) -----------
    [ -n "$years_v1" ]   && sbatch_wait ingest_monthly.r \
        "WORK_DIR=${WORK_DIR},MICASA_VERSION=v1,MICASA_YEAR_START=${years_v1%% *},MICASA_YEAR_END=${years_v1##* }" \
        "ingestmonthly-v1"
    [ -n "$years_vnrt" ] && sbatch_wait ingest_monthly.r \
        "WORK_DIR=${WORK_DIR},MICASA_VERSION=vNRT,MICASA_YEAR_START=${years_vnrt%% *},MICASA_YEAR_END=${years_vnrt##* }" \
        "ingestmonthly-vNRT"
else
    echo; echo "==> [skip] ingest stages (2,3)"
fi

# ---- Stage 4: Link vNRT → v1 (NRT years only) ------------------------------

if [ "$skip_link" -eq 0 ] && [ -n "$years_vnrt" ]; then
    for y in $years_vnrt; do
        echo; echo "==> link vNRT→v1 for ${y} (dailies + monthlies)"
        if [ "$dry_run" -eq 0 ]; then
            MICASA_YEAR="$y" sh link_vNRT_to_v1.sh || true   # partial final year is OK
            link_vnrt_monthlies "$y"
        fi
    done
elif [ -z "$years_vnrt" ]; then
    echo; echo "==> [skip] vNRT→v1 link (no vNRT years)"
else
    echo; echo "==> [skip] vNRT→v1 link (--skip-link)"
fi

# ---- Everything past here is a single continuous v1 stream ------------------
export MICASA_VERSION="$POST_VERSION"
export MICASA_YEAR_START="$YEAR_START"
export MICASA_YEAR_END="$YEAR_END"
export MICASA_MONTH_START=1
export MICASA_MONTH_END=12

# ---- Stage 5: Aggregate / climatology --------------------------------------

if [ "$skip_aggregate" -eq 0 ]; then
    run sh cat_monthly.sh
    run sh compute_clim.sh
    run sh compute_daily_clim.sh
    run sh link_daily_clim.sh
else
    echo; echo "==> [skip] aggregate stage (5)"
fi

# ---- Stage 6: Fit (sub-monthly smoother) -----------------------------------

if [ "$skip_fit" -eq 0 ]; then
    run Rscript "$fitter_script"
else
    echo; echo "==> [skip] fit stage (6, $fitter)"
fi

# ---- Stage 7: Diurnalize (ERA5 hourly, per-year fan-out) -------------------

if [ "$skip_diurnalize" -eq 0 ]; then
    fanout_wait diurnalize-ERA5.r diurn_year "$years_all" \
        "MICASA_VERSION=${POST_VERSION},MICASA_STRICT_PIQS=1" "d-MiCASA"
else
    echo; echo "==> [skip] diurnalize stage (7)"
fi

# ---- Stage 8: Day-split -----------------------------------------------------

if [ "$skip_daysplit" -eq 0 ]; then
    run sh daysplitter.sh
else
    echo; echo "==> [skip] daysplit stage (8)"
fi

# ---- Stage 9: Provenance stamp ---------------------------------------------
# Best-effort PROVENANCE.txt into the output dir — always exits 0.
run sh write_provenance.sh

manifest_record run_record.sh ok - "record ${YEAR_START}-${YEAR_END} (v1 ${years_v1:-none} / vNRT ${years_vnrt:-none})"
echo; echo "==> Full-record run finished for ${YEAR_START}-${YEAR_END}."
