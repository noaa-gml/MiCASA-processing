#!/bin/bash
# Build a CLIMATOLOGY PRIOR: MiCASA fluxes whose monthly magnitude is a
# multi-year climatology but whose sub-daily structure comes from each target
# day's REAL meteorology.
#
# Usage:
#     ./run_climatology_prior.sh --outdir DIR [options]
#
#   --outdir DIR       output tree (created; must not be a production tree)
#   --source FILE      concatenated monthly file to average
#                      (default: $MONTHLY_1X1_DIR/MiCASA_<ver>_flux_x360_y180_monthly.nc)
#   --baseline A-B     years averaged into the climatology   (default 2001-2020)
#   --span A:B         synthetic series span, YYYY-MM:YYYY-MM (default 2020-01:2026-12)
#   --years A-B        years actually delivered              (default 2021-2025)
#   --version V        v1 | vNRT                             (default from config.sh)
#   --fitter NAME      sub-monthly smoother                  (default pchip)
#   --account ACCT     SLURM account                         (default co2)
#   --reference DIR    ERA5 dir of the unmodified product, for the verify stage
#   --skip-series --skip-fit --skip-diurnalize --skip-daysplit --skip-verify
#   --dry-run
#
# Example — the CT2026 request (Jacobson, 2026-08-10):
#     ./run_climatology_prior.sh \
#         --outdir /work2/noaa/co2/ash/micasa-clim \
#         --source /work2/noaa/co2/GFED-CASA/2026/MiCASA.0/monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc \
#         --baseline 2001-2020 --years 2021-2025 \
#         --reference /work2/noaa/co2/GFED-CASA/2026/MiCASA.0/ERA5
#
# Stages:
#   1. Series      — make_climatology_series.py: per-month + concatenated
#                    climatological monthly files
#   2. Fit         — the sub-monthly smoother, REBUILT on that series
#   3. Diurnalize  — climatology_diurnalize.sbatch, one array task per year
#   4. Day-split   — climatology_daysplit.sbatch: stamp provenance, split to
#                    per-day NEE
#   5. Verify      — tests/verify_climatology_prior.py
#
# ** Stage 2 is not optional and not a formality. ** diurnalize's hourly flux is
# f(t) = driver(t)*mean/mean_driver - mean + qmod(t), whose monthly mean is
# mean(qmod) — taken entirely from the fit. Reusing a production fit.piqs.rda
# would reinstate the real interannual signal through qmod while producing a
# complete, healthy-looking set of files. See docs/CLIMATOLOGY_PRIOR.md.
#
# This driver blocks on SLURM stages (sbatch --wait), so for a multi-year build
# run it from a batch job or a session that can wait, not a login shell you are
# about to close.

set -e
set -o pipefail

# ---- Args -------------------------------------------------------------------

if [ $# -lt 1 ]; then
    sed -n '2,40p' "$0"
    exit 1
fi

outdir=""; source_file=""; reference=""
baseline=""; span=""; years=""
fitter=pchip; account=co2
skip_series=0; skip_fit=0; skip_diurnalize=0; skip_daysplit=0; skip_verify=0
dry_run=0

while [ $# -gt 0 ]; do
    case "$1" in
        --outdir)     outdir="$2";      shift 2 ;;
        --outdir=*)   outdir="${1#*=}"; shift   ;;
        --source)     source_file="$2"; shift 2 ;;
        --source=*)   source_file="${1#*=}"; shift ;;
        --reference)  reference="$2";   shift 2 ;;
        --reference=*) reference="${1#*=}"; shift ;;
        --baseline)   baseline="$2";    shift 2 ;;
        --baseline=*) baseline="${1#*=}"; shift ;;
        --span)       span="$2";        shift 2 ;;
        --span=*)     span="${1#*=}";   shift   ;;
        --years)      years="$2";       shift 2 ;;
        --years=*)    years="${1#*=}";  shift   ;;
        --version)    export MICASA_VERSION="$2"; shift 2 ;;
        --version=*)  export MICASA_VERSION="${1#*=}"; shift ;;
        --fitter)     fitter="$2";      shift 2 ;;
        --fitter=*)   fitter="${1#*=}"; shift   ;;
        --account)    account="$2";     shift 2 ;;
        --account=*)  account="${1#*=}"; shift  ;;
        --skip-series)     skip_series=1;     shift ;;
        --skip-fit)        skip_fit=1;        shift ;;
        --skip-diurnalize) skip_diurnalize=1; shift ;;
        --skip-daysplit)   skip_daysplit=1;   shift ;;
        --skip-verify)     skip_verify=1;     shift ;;
        --dry-run)         dry_run=1;         shift ;;
        *) echo "Unknown flag: $1"; exit 2 ;;
    esac
done

[ -n "$outdir" ] || { echo "--outdir is required"; exit 2; }

case "$fitter" in
    pchip) fitter_script=write_pchip.r ;;
    piqs)  fitter_script=write_piqs.r
           export MICASA_PIQS_PAD_RIGHT="${MICASA_PIQS_PAD_RIGHT:-2}" ;;
    ppm)   fitter_script=write_ppm.r   ;;
    linmm) fitter_script=write_linmm.r ;;
    mss)   fitter_script=write_mss.r   ;;
    atpk)  fitter_script=write_atpk.r  ;;
    *) echo "Unknown --fitter: '$fitter'"; exit 2 ;;
esac

# ---- Config -----------------------------------------------------------------

. "$(dirname "$0")/config.sh"
. "$(dirname "$0")/lib/manifest.sh"

if [ -n "$baseline" ]; then
    export MICASA_CLIM_BASELINE_START="${baseline%%-*}"
    export MICASA_CLIM_BASELINE_END="${baseline##*-}"
fi
if [ -n "$span" ]; then
    export MICASA_CLIM_SPAN_START="${span%%:*}"
    export MICASA_CLIM_SPAN_END="${span##*:}"
fi
if [ -n "$years" ]; then
    export MICASA_CLIM_YEAR_START="${years%%-*}"
    export MICASA_CLIM_YEAR_END="${years##*-}"
fi

# Default source: the concatenation in whatever monthly dir config.sh resolved
# to (i.e. the tree being climatologized) — captured BEFORE we repoint
# MONTHLY_1X1_DIR at the output tree.
if [ -z "$source_file" ]; then
    source_file="${MICASA_CLIM_SOURCE:-${MONTHLY_1X1_DIR}/MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly.nc}"
fi
export MICASA_CLIM_SOURCE="$source_file"

# Everything the pipeline writes goes into the climatology tree.
export MONTHLY_1X1_DIR="${outdir}/monthly_1x1"
export ERA5_DIR="${outdir}/ERA5"
export JOBS_DIR="${outdir}/jobs"
export MICASA_FIT_RDA="${outdir}/fit.piqs.rda"      # absolute is fine (see diurnalize-ERA5.r)

# Production diurnalization settings, stated rather than inherited, so the
# climatology is the ONLY difference from the ordinary product.
export MICASA_RESP_DRIVER="${MICASA_RESP_DRIVER:-airtemp}"
export MICASA_RESP_TEMPFUN="${MICASA_RESP_TEMPFUN:-q10}"
export MICASA_POLAR_CLIP="${MICASA_POLAR_CLIP:-conserve}"
export MICASA_STRICT_PIQS="${MICASA_STRICT_PIQS:-1}"

Y0="$MICASA_CLIM_YEAR_START"; Y1="$MICASA_CLIM_YEAR_END"
MARKER="${outdir}/.micasa-climatology-prior"

echo "========================================================================"
echo "MiCASA climatology prior"
echo "  OUTDIR      ${outdir}"
echo "  SOURCE      ${MICASA_CLIM_SOURCE}"
echo "  BASELINE    ${MICASA_CLIM_BASELINE_START}-${MICASA_CLIM_BASELINE_END}"
echo "  SPAN        ${MICASA_CLIM_SPAN_START} .. ${MICASA_CLIM_SPAN_END}"
echo "  DELIVER     ${Y0}-${Y1}"
echo "  VERSION     ${MICASA_VERSION}"
echo "  FITTER      ${fitter}  (${fitter_script})"
echo "  WORK_DIR    ${WORK_DIR}"
echo "  DIURNAL     driver=${MICASA_RESP_DRIVER} tempfun=${MICASA_RESP_TEMPFUN} clip=${MICASA_POLAR_CLIP}"
[ "$dry_run" -eq 1 ] && echo "  *** DRY RUN — nothing will execute ***"
echo "========================================================================"

# ---- Stage 0: preflight -----------------------------------------------------
# Refuse to write into a tree we did not create. The output is a synthetic
# product; dropping it on top of a real one would be unrecoverable and the
# resulting mixture unattributable.
if [ -e "$outdir" ] && [ ! -e "$MARKER" ]; then
    if [ -n "$(ls -A "$outdir" 2>/dev/null)" ]; then
        echo "FATAL: ${outdir} exists, is not empty, and carries no"
        echo "       ${MARKER##*/} marker — refusing to write a synthetic"
        echo "       product into a tree this driver did not create."
        exit 1
    fi
fi
[ -s "$MICASA_CLIM_SOURCE" ] || {
    echo "FATAL: source monthly file not found: ${MICASA_CLIM_SOURCE}"
    echo "       run cat_monthly.sh on the source tree first."
    exit 1
}

if [ "$dry_run" -eq 0 ]; then
    mkdir -p "$MONTHLY_1X1_DIR" "$ERA5_DIR" "$JOBS_DIR"
    if [ ! -e "$MARKER" ]; then
        {
            echo "MiCASA climatology prior — written by run_climatology_prior.sh"
            echo "created   $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
            echo "source    ${MICASA_CLIM_SOURCE}"
            echo "baseline  ${MICASA_CLIM_BASELINE_START}-${MICASA_CLIM_BASELINE_END}"
            echo "span      ${MICASA_CLIM_SPAN_START} .. ${MICASA_CLIM_SPAN_END}"
            echo "deliver   ${Y0}-${Y1}"
            echo "SYNTHETIC PRODUCT — not observational. See docs/CLIMATOLOGY_PRIOR.md"
        } > "$MARKER"
    fi
fi

trap 'manifest_record run_climatology_prior.sh fail - "aborted (line $LINENO)"' ERR
_t0=$(date +%s)

run() {
    echo; echo "==> $*"
    [ "$dry_run" -eq 1 ] && return 0
    "$@"
}

# Submit a per-year array and block until EVERY task has reached a terminal
# state, then require every one of them to have succeeded.
#
# ⚠ Do NOT use `sbatch --wait` here. On a job array it returns once SOME tasks
# have finished, not all of them: observed 2026-08-10 returning while two of
# five diurnalize tasks were still running, which started the day-split against
# a year whose last monthly file was still being written (ncatted failed on a
# half-written fluxes_202412.nc). The per-task artefact gates caught it, but the
# driver must not create that race in the first place.
sbatch_array_wait() {
    local script="$1" jobname="$2" logtag="$3"
    echo; echo "==> sbatch --array=${Y0}-${Y1} ${script}  (then wait for ALL tasks)"
    [ "$dry_run" -eq 1 ] && return 0

    local jid
    jid=$(sbatch --parsable \
                 --array="${Y0}-${Y1}" \
                 -J "${jobname}" \
                 -A "${account}" \
                 --mail-user="${MAIL_USER}" \
                 --mail-type=FAIL \
                 --export=ALL \
                 --output="${JOBS_DIR}/${logtag}-%a.o%A" \
                 "${WORK_DIR}/${script}")
    if [ -z "$jid" ]; then echo "FATAL: sbatch did not return a job id"; return 1; fi
    echo "    job ${jid}, ${Y0}-${Y1}"

    # Poll until nothing for this array is queued, running or completing.
    # (This runs on a compute node inside the driver's own job, not on a login
    # node, so a poll loop is appropriate here.)
    while squeue -h -j "$jid" -t PENDING,RUNNING,CONFIGURING,COMPLETING 2>/dev/null \
          | grep -q .; do
        sleep 30
    done

    # Then require every task to have succeeded. -X reports the allocation only,
    # not the .batch/.extern steps.
    local states bad n
    states=$(sacct -j "$jid" --format=JobID%20,State,ExitCode -n -X 2>/dev/null)
    echo "$states" | sed 's/^/    /'
    # ⚠ Count with `wc -l`, never `grep -c`: grep exits 1 when it counts zero,
    # and under `set -e` a command substitution in an assignment propagates that
    # status -- so `bad=$(... | grep -c .)` aborted the driver in exactly the
    # case where every task had SUCCEEDED. (Observed 2026-08-10.)
    n=$(echo "$states" | awk 'NF' | wc -l)
    bad=$(echo "$states" | awk 'NF && $2 != "COMPLETED"' | wc -l)
    local want=$((Y1 - Y0 + 1))
    if [ "$n" -ne "$want" ]; then
        echo "FATAL: ${script}: expected ${want} array tasks, sacct reports ${n}"
        return 1
    fi
    if [ "$bad" -ne 0 ]; then
        echo "FATAL: ${script}: ${bad}/${n} array tasks did not complete cleanly"
        return 1
    fi
    echo "    all ${n} tasks COMPLETED"
}

cd "${WORK_DIR}"

# ---- Stage 1: climatological monthly series ---------------------------------
if [ "$skip_series" -eq 0 ]; then
    run ${PYTHON:-python3} ./make_climatology_series.py
else
    echo "==> [skip] series stage"
fi

# ---- Stage 2: refit the sub-monthly smoother on that series -----------------
# write_pchip.r and friends save fit.piqs.rda into the CURRENT directory.
if [ "$skip_fit" -eq 0 ]; then
    echo; echo "==> (cd ${outdir} && Rscript ${fitter_script})"
    if [ "$dry_run" -eq 0 ]; then
        ( cd "${outdir}" && Rscript "${WORK_DIR}/${fitter_script}" )
        [ -s "${outdir}/fit.piqs.rda" ] || {
            echo "FATAL: ${outdir}/fit.piqs.rda absent or empty after ${fitter_script}"
            exit 1
        }
        ls -l "${outdir}/fit.piqs.rda"
    fi
else
    echo "==> [skip] fit stage"
fi

# ---- Stage 3: diurnalize ----------------------------------------------------
if [ "$skip_diurnalize" -eq 0 ]; then
    sbatch_array_wait climatology_diurnalize.sbatch \
        "climprior-diurn" "clim-diurnalize"
else
    echo "==> [skip] diurnalize stage"
fi

# ---- Stage 4: stamp + day-split ---------------------------------------------
if [ "$skip_daysplit" -eq 0 ]; then
    sbatch_array_wait climatology_daysplit.sbatch \
        "climprior-daysplit" "clim-daysplit"
else
    echo "==> [skip] daysplit stage"
fi

# ---- Stage 5: verify --------------------------------------------------------
if [ "$skip_verify" -eq 0 ]; then
    verify_args=(--outdir "${outdir}" --years "${Y0}-${Y1}"
                 --source "${MICASA_CLIM_SOURCE}"
                 --baseline "${MICASA_CLIM_BASELINE_START}-${MICASA_CLIM_BASELINE_END}")
    [ -n "$reference" ] && verify_args+=(--reference "${reference}")
    run ${PYTHON:-python3} ./tests/verify_climatology_prior.py "${verify_args[@]}"
else
    echo "==> [skip] verify stage"
fi

manifest_record run_climatology_prior.sh ok "$(($(date +%s) - _t0))" \
    "outdir=${outdir} baseline=${MICASA_CLIM_BASELINE_START}-${MICASA_CLIM_BASELINE_END} years=${Y0}-${Y1}"
trap - ERR

echo
echo "========================================================================"
echo "climatology prior complete: ${outdir}"
echo "  CarbonTracker uses it by pointing ONE rc key at the ERA5 dir:"
echo "      bio.input.dir : ${ERA5_DIR}"
echo "  fire/fuel (ct.wildfire -> fires.input.dir) are NOT part of this"
echo "  product and should keep pointing at the unmodified daily_1x1 tree."
echo "========================================================================"
