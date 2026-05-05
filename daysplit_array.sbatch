#!/bin/bash
#
# daysplit_array.sbatch — parallel daysplitter via SLURM array
#
# Runs daysplitter.sh in parallel, one array task per year. Splits each
# ERA5/fluxes_YYYYMM.nc into per-day MiCASA_<ver>.nee.YYYYMMDD.nc files.
#
# Usage (run from this directory):
#   sbatch --array=2001-2025          daysplit_array.sbatch   # full record
#   sbatch --array=2025               daysplit_array.sbatch   # single year
#   sbatch --array=2025-2026%4        daysplit_array.sbatch   # cap concurrency at 4
#
# Optional env overrides (passed through to daysplitter.sh):
#   MICASA_VERSION         "v1" / "v2" / "vNRT" (default from config.sh = v1)
#   MICASA_MONTH_START     restrict to a month range within each year
#   MICASA_MONTH_END
#
# Examples:
#   # Backfill v1 daily NEE for the gap 2025-06..2025-12
#   MICASA_VERSION=v1 MICASA_MONTH_START=6 MICASA_MONTH_END=12 \
#     sbatch --array=2025 daysplit_array.sbatch
#
#   # Backfill all of v2 (matches what 22987776 did)
#   MICASA_VERSION=v2 sbatch --array=2001-2025 daysplit_array.sbatch
#
# Files that already exist are overwritten (ncks -O), so re-running over
# already-split months is safe but does redo the work; use MICASA_MONTH_*
# to restrict if you want only the gap.
#
#SBATCH --job-name=daysplit
#SBATCH --account=co2
#SBATCH --time=1:00:00
#SBATCH --ntasks=1 --cpus-per-task=2 --mem=4g
#SBATCH --output=jobs/daysplit-%a.o%A
#SBATCH --error=jobs/daysplit-%a.o%A

set -e

cd "$SLURM_SUBMIT_DIR"
. ./config.sh

module load nco/5.1.6 2>/dev/null || true

YR=${SLURM_ARRAY_TASK_ID:?"missing SLURM_ARRAY_TASK_ID; submit with --array=YYYY-YYYY"}
MS=${MICASA_MONTH_START:-1}
ME=${MICASA_MONTH_END:-12}

echo "=== daysplit ${MICASA_VERSION} year=${YR} months=${MS}-${ME} ==="

MICASA_YEAR_START=$YR MICASA_YEAR_END=$YR \
MICASA_MONTH_START=$MS MICASA_MONTH_END=$ME \
  ./daysplitter.sh

n=$(ls "${ERA5_DIR}/MiCASA_${MICASA_VERSION}.nee.${YR}"*.nc 2>/dev/null | wc -l)
echo "=== Done year ${YR}; ${n} daily NEE files for MiCASA_${MICASA_VERSION} ==="
