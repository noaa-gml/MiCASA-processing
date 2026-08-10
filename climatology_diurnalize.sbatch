#!/bin/bash
#
# climatology_diurnalize.sbatch — diurnalize the climatological monthly series,
# one array task per year. Stage 3 of run_climatology_prior.sh.
#
# Usage (run from the checkout; the driver does this for you):
#   sbatch --array=2021-2025 --export=ALL climatology_diurnalize.sbatch
#
# Everything here is the stock diurnalization at production settings. The ONLY
# difference from an ordinary run is that MONTHLY_1X1_DIR and MICASA_FIT_RDA
# point at a CLIMATOLOGICAL series and its matching fit. The meteorology is the
# real, year-specific ERA5 — that is the whole point of the product.
#
#SBATCH --job-name=clim-diurnalize
#SBATCH --account=co2
#SBATCH --time=8:00:00
#SBATCH --ntasks=1 --mem=48g
#SBATCH --output=jobs/clim-diurnalize-%a.o%A

set -u

cd "$SLURM_SUBMIT_DIR"
. ./config.sh

YR=${SLURM_ARRAY_TASK_ID:?"missing SLURM_ARRAY_TASK_ID; submit with --array=YYYY-YYYY"}

export MICASA_YEAR=$YR
export MICASA_YEAR_START=$YR
export MICASA_YEAR_END=$YR
export diurn_year=$YR

# A month with real monthly data that falls outside the fit window would get
# coefficient-CLIMATOLOGY instead of its own coefficients, silently. For this
# product the fit is built to span the series, so that must never happen.
export MICASA_STRICT_PIQS="${MICASA_STRICT_PIQS:-1}"

echo "=== clim-diurnalize year=${YR} ==="
echo "host          $(hostname)"
echo "monthly dir   ${MONTHLY_1X1_DIR}"
echo "era5 out dir  ${ERA5_DIR}"
echo "fit           ${MICASA_FIT_RDA:-fit.piqs.rda}"

# ---- PRE-GATE ---------------------------------------------------------------
# Point-at-the-right-reference guard. The fit and the monthly series must both
# be the climatological ones; a run that quietly used the production fit would
# reinstate the trend and still produce every expected file.
FIT="${MICASA_FIT_RDA:-fit.piqs.rda}"
case "$FIT" in /*) FIT_ABS="$FIT" ;; *) FIT_ABS="${WORK_DIR}/${FIT}" ;; esac
[ -s "$FIT_ABS" ] || { echo "FATAL: fit not found: $FIT_ABS"; exit 1; }

MON="${MONTHLY_1X1_DIR}/MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly_${YR}07.nc"
[ -s "$MON" ] || { echo "FATAL: monthly input not found: $MON"; exit 1; }

# The series must actually be climatological. Reading the attribute the builder
# stamps is cheap and turns a silent wrong-input into a hard stop.
FLAG=$(ncdump -h "$MON" 2>/dev/null |
       sed -n 's/.*:flux_from_climatology = "\(.*\)" ;.*/\1/p' | head -1)
if [ "$FLAG" != "yes" ]; then
    echo "FATAL: $MON is not a climatological series (flux_from_climatology='$FLAG')."
    echo "       Run make_climatology_series.py first, or point MONTHLY_1X1_DIR"
    echo "       at the climatology tree."
    exit 1
fi

# ---- RUN --------------------------------------------------------------------
Rscript ./diurnalize-ERA5.r
rc=$?
echo "diurnalize-ERA5.r exit $rc"

# ---- POST-GATE: the artefact count, never the exit code ---------------------
MS=${MICASA_MONTH_START:-1}
ME=${MICASA_MONTH_END:-12}
want=$((ME - MS + 1))
got=$(ls "${ERA5_DIR}"/fluxes_${YR}??.nc 2>/dev/null | wc -l)
echo "artefacts: ${got}/${want} fluxes_${YR}MM.nc"
[ "$got" -eq "$want" ] || { echo "FATAL: expected $want monthly flux files, found $got"; exit 1; }

z=$(find "${ERA5_DIR}" -maxdepth 1 -name "fluxes_${YR}??.nc" -size 0 | wc -l)
[ "$z" -eq 0 ] || { echo "FATAL: $z zero-byte flux files for $YR"; exit 1; }

echo "=== Done year ${YR}: ${got} monthly hourly flux files ==="
