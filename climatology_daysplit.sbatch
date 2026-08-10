#!/bin/bash
#
# climatology_daysplit.sbatch — stamp the climatology provenance onto the
# monthly hourly files, then split them into the per-day NEE files that
# CarbonTracker reads. Stage 4 of run_climatology_prior.sh.
#
# Usage (run from the checkout; the driver does this for you):
#   sbatch --array=2021-2025 --export=ALL climatology_daysplit.sbatch
#
# The stamp goes on BEFORE the split so ncks propagates it into every daily
# file. diurnalize-ERA5.r writes flux_from_climatology="no" here, because it
# was handed ordinary per-month files and only sets "yes" on its own
# NPPclim/Rhclim fallback path — true of the code path it took, misleading
# about the data. The attributes below are the ones a downstream user (or a
# future us) will actually read off the file to find out what it is.
#
#SBATCH --job-name=clim-daysplit
#SBATCH --account=co2
#SBATCH --time=4:00:00
#SBATCH --ntasks=1 --cpus-per-task=2 --mem=8g
#SBATCH --output=jobs/clim-daysplit-%a.o%A

set -u

cd "$SLURM_SUBMIT_DIR"
. ./config.sh

module load nco/5.1.6 2>/dev/null || true

YR=${SLURM_ARRAY_TASK_ID:?"missing SLURM_ARRAY_TASK_ID; submit with --array=YYYY-YYYY"}
MS=${MICASA_MONTH_START:-1}
ME=${MICASA_MONTH_END:-12}

echo "=== clim-daysplit ${MICASA_VERSION} year=${YR} months=${MS}-${ME} ==="
command -v ncks    >/dev/null || { echo "FATAL: ncks not on PATH"; exit 1; }
command -v ncatted >/dev/null || { echo "FATAL: ncatted not on PATH"; exit 1; }

want=$((ME - MS + 1))
got=$(ls "${ERA5_DIR}"/fluxes_${YR}??.nc 2>/dev/null | wc -l)
[ "$got" -eq "$want" ] || { echo "FATAL: expected $want fluxes_${YR}MM.nc, found $got"; exit 1; }

# ---- Stamp ------------------------------------------------------------------
# Inherited from the climatological monthly series so the two agree; fall back
# to the configured baseline if the attribute is somehow absent.
SRC_MON="${MONTHLY_1X1_DIR}/MiCASA_${MICASA_VERSION}_flux_x360_y180_monthly_${YR}07.nc"
BASE=$(ncdump -h "$SRC_MON" 2>/dev/null |
       sed -n 's/.*:climatology_baseline = "\(.*\)" ;.*/\1/p' | head -1)
BASE=${BASE:-${MICASA_CLIM_BASELINE_START}-${MICASA_CLIM_BASELINE_END}}

VARIANT="MiCASA climatology ${BASE} (bio only)"
METHOD="modulo-month mean of MiCASA monthly NPP and Rh over ${BASE}, diurnalized with each target day's REAL ERA5 hourly meteorology (ssrd drives GPP; t2m through Q10 drives respiration). The monthly means AND the sub-monthly PCHIP coefficient fit are both climatological; only the meteorology is year-specific."
SCOPE="BIO CHANNEL ONLY (NEE = Rh - NPP). FIRE and FUEL are NOT part of this product and are NOT climatological: CarbonTracker reads them from a separate rc key (ct.wildfire -> fires.input.dir), which continues to point at the unmodified daily_1x1 tree."
NOTE="Interannual variability and trend are removed BY CONSTRUCTION: every calendar month carries the same monthly mean in every year. Day-to-day and hour-to-hour structure comes entirely from real meteorology. Not an observational product."

for f in "${ERA5_DIR}"/fluxes_${YR}??.nc; do
    ncatted -O \
        -a flux_from_climatology,global,o,c,"yes" \
        -a prior_variant,global,o,c,"$VARIANT" \
        -a climatology_baseline,global,o,c,"$BASE" \
        -a climatology_method,global,o,c,"$METHOD" \
        -a scope_note,global,o,c,"$SCOPE" \
        -a provenance_note,global,o,c,"$NOTE" \
        "$f" || { echo "FATAL: ncatted failed on $f"; exit 1; }
done
echo "stamped ${got} monthly files with baseline ${BASE}"

# ---- Split ------------------------------------------------------------------
MICASA_YEAR_START=$YR MICASA_YEAR_END=$YR \
MICASA_MONTH_START=$MS MICASA_MONTH_END=$ME \
  ./daysplitter.sh
rc=$?
echo "daysplitter.sh exit $rc"

# ---- POST-GATE: the artefact count, never the exit code ---------------------
# Expected days across the requested month range, leap-aware.
exp=$(python3 - "$YR" "$MS" "$ME" <<'EOF'
import calendar, sys
y, a, b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
print(sum(calendar.monthrange(y, m)[1] for m in range(a, b + 1)))
EOF
)
n=$(ls "${ERA5_DIR}/MiCASA_${MICASA_VERSION}.nee.${YR}"????.nc 2>/dev/null | wc -l)
echo "artefacts: ${n}/${exp} daily NEE files"
[ "$n" -eq "$exp" ] || { echo "FATAL: expected $exp daily files, found $n"; exit 1; }

z=$(find "${ERA5_DIR}" -maxdepth 1 -name "MiCASA_${MICASA_VERSION}.nee.${YR}????.nc" -size 0 | wc -l)
[ "$z" -eq 0 ] || { echo "FATAL: $z zero-byte daily files"; exit 1; }

# Every daily file must carry a full 24-hour axis. A truncated file would be
# read without complaint by a consumer that just indexes hours.
bad=0
for f in "${ERA5_DIR}/MiCASA_${MICASA_VERSION}.nee.${YR}"????.nc; do
    ns=$(ncdump -h "$f" | sed -n 's/.*time = UNLIMITED ; \/\/ (\([0-9]*\) currently).*/\1/p')
    [ "$ns" = "24" ] || { echo "  BAD slots=$ns $f"; bad=$((bad+1)); }
done
[ "$bad" -eq 0 ] || { echo "FATAL: $bad daily files without 24 hourly slots"; exit 1; }

echo "=== Done year ${YR}: ${n} daily NEE files, all with 24 slots ==="
