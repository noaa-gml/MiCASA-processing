#!/bin/bash
# spinup_year.sh -- generate climatology spin-up flux files for pre-record years.
#
# The real MiCASA record starts 2001-01. For an inversion that needs earlier
# boundary / spin-up years, this fills a year with the day-of-year NPP/Rh
# climatology diurnalized against that year's REAL ERA5 meteo (each output is
# flagged flux_from_climatology="yes"). ERA5 meteo is available back to 1989.
#
# Per year, in the configured output dirs (honors run.env overrides):
#   daily_1x1/ ...daily_YYYYMMDD.nc      day-of-year clim symlinks (link_daily_clim.sh)
#   ERA5/      fluxes_YYYYMM.nc          diurnalized hourly        (diurnalize-ERA5.r)
#   ERA5/      MiCASA_v1.nee.YYYYMMDD.nc per-day NEE               (daysplitter.sh)
# then refreshes PROVENANCE.txt once every year's daysplit has succeeded.
#
# Usage:
#   cd <checkout>; source <product>/run.env      # point outputs at the product
#   ./spinup_year.sh [--dry-run] YYYY [YYYY ...]
#
# Example (the MiCASA.0 product):
#   source /work2/noaa/co2/GFED-CASA/2026/MiCASA.0/run.env
#   ./spinup_year.sh 1998 1997 1996
#
# SBATCH chain per year: diurnalize -> (afterok) daysplit; one final
# write_provenance after all daysplits (afterok). Short walltimes so the jobs
# slot in around maintenance reservations (a single year takes ~20 min).
#
# Env (beyond the standard config.sh / run.env set):
#   MICASA_FIT_RDA   sub-monthly fit, RELATIVE to WORK_DIR (default fit.piqs.rda;
#                    for a redirected product set ../2026/MiCASA.0/fit.piqs.rda --
#                    diurnalize-ERA5.r does file.path(WORK_DIR, MICASA_FIT_RDA),
#                    so an ABSOLUTE path mis-concatenates).
#   SPINUP_DIURN_TIME / SPINUP_SPLIT_TIME   SBATCH walltimes (def 02:00:00/00:45:00)

set -eu

DRY=0
if [ "${1:-}" = "--dry-run" ]; then DRY=1; shift; fi
[ "$#" -ge 1 ] || { sed -n '2,33p' "$0"; exit 1; }

. "$(dirname "$0")/config.sh"

WORK_DIR="${WORK_DIR:-$(cd "$(dirname "$0")" && pwd)}"
FIT_RDA="${MICASA_FIT_RDA:-fit.piqs.rda}"
DIURN_TIME="${SPINUP_DIURN_TIME:-02:00:00}"
SPLIT_TIME="${SPINUP_SPLIT_TIME:-00:45:00}"
RECORD_START="${MICASA_RECORD_START:-2001}"
METEO_ROOT="${MICASA_ERA5_DIR:-${CARBONTRACKER:-}/METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100}"

# sbatch wrapper: in --dry-run, print the command and return a placeholder id.
_n=0
submit() {
    if [ "$DRY" = 1 ]; then
        _n=$((_n + 1)); echo "  [dry] sbatch $*" >&2; echo "DRY$_n"
    else
        sbatch --parsable "$@"
    fi
}

echo "spin-up: WORK_DIR=$WORK_DIR  ERA5_DIR=${ERA5_DIR}  FIT=$FIT_RDA${DRY:+  (dry-run=$DRY)}"
splits=""
for Y in "$@"; do
    case "$Y" in *[!0-9]*|"") echo "skip '$Y': not a 4-digit year"; continue;; esac
    if [ "$Y" -ge "$RECORD_START" ]; then
        echo "WARN $Y >= $RECORD_START: real MiCASA record exists -- spin-up is for pre-record years; skipping"
        continue
    fi
    # Soft meteo guard: only skip if we can see the meteo root and the year is absent.
    if [ -d "$METEO_ROOT" ] && [ ! -d "$METEO_ROOT/$Y" ]; then
        echo "WARN no ERA5 meteo for $Y under $METEO_ROOT -- cannot diurnalize; skipping"
        continue
    fi

    echo "=== spin-up $Y ==="
    # 1. daily_1x1 day-of-year climatology links (instant, idempotent).
    if [ "$DRY" = 1 ]; then
        echo "  [dry] MICASA_CLIM_YEARS=$Y bash link_daily_clim.sh"
    else
        MICASA_VERSION="${MICASA_VERSION:-v1}" MICASA_CLIM_YEARS="$Y" \
            bash "$WORK_DIR/link_daily_clim.sh" >/dev/null
        echo "  daily_1x1 clim links done"
    fi

    # 2. diurnalize (worker mode: diurn_year set => this year only).
    diurn=$(submit --account=co2 --time="$DIURN_TIME" --ntasks=1 --cpus-per-task=1 \
        --mem=40g --partition=orion -J "diurn-$Y-MiCASA" \
        --output="$JOBS_DIR/diurn-$Y.o%j" --mail-type=FAIL --mail-user="$MAIL_USER" \
        --chdir="$WORK_DIR" \
        --export="ALL,WORK_DIR=$WORK_DIR,diurn_year=$Y,MICASA_YEAR_START=$Y,MICASA_YEAR_END=$Y,MICASA_VERSION=v1,MICASA_FIT_RDA=$FIT_RDA,MICASA_MONTH_START=1,MICASA_MONTH_END=12" \
        --wrap='Rscript diurnalize-ERA5.r')
    echo "  diurnalize job: $diurn"

    # 3. daysplit (afterok on diurnalize).
    split=$(submit --dependency="afterok:$diurn" --account=co2 --time="$SPLIT_TIME" \
        --ntasks=1 --cpus-per-task=1 --mem=8g --partition=orion -J "daysplit-$Y-MiCASA" \
        --output="$JOBS_DIR/daysplit-$Y.o%j" --mail-type=FAIL --mail-user="$MAIL_USER" \
        --chdir="$WORK_DIR" \
        --export="ALL,WORK_DIR=$WORK_DIR,MICASA_YEAR_START=$Y,MICASA_YEAR_END=$Y,MICASA_VERSION=v1,MICASA_MONTH_START=1,MICASA_MONTH_END=12" \
        --wrap='module load nco && bash daysplitter.sh')
    echo "  daysplit job:   $split (afterok:$diurn)"
    splits="${splits:+$splits:}$split"
done

[ -n "$splits" ] || { echo "no years submitted"; exit 0; }

# Final: refresh PROVENANCE.txt once every daysplit has succeeded.
prov=$(submit --dependency="afterok:$splits" --account=co2 --time=00:05:00 \
    --ntasks=1 --cpus-per-task=1 --mem=2g --partition=orion -J "provstamp-spinup" \
    --output="$JOBS_DIR/provstamp-spinup.o%j" --mail-type=FAIL --mail-user="$MAIL_USER" \
    --chdir="$WORK_DIR" --export="ALL,WORK_DIR=$WORK_DIR" \
    --wrap='sh write_provenance.sh')
echo "provenance re-stamp: $prov (afterok:$splits)"
echo "done."
