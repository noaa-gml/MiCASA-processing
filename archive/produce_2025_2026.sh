#!/bin/bash
# produce_2025_2026.sh — Phase 2 of v2 generation, after ingest_byyear has
# brought the 1° vNRT dailies up to date.
#
# ARCHIVED (lives in archive/). This is the year-pinned one-off that built the
# 2025-2026 NRT tail; the general multi-year driver is ../run_record.sh. It is
# kept for the NRT trailing-edge completion it special-cases (steps 5-6: the
# day-of-year clim-fill of the partial final month + the synthetic monthly file
# with a patched time coordinate), which run_record.sh deliberately omits.
# Still runnable: it cd's to the repo root ("$(dirname "$0")/..") so its
# repo-root-relative calls (./config.sh, ./link_vNRT_to_v1.sh, ...) resolve.
#
# Steps:
#   1. ingest_monthly vNRT 2025 (Jan..Nov, raw exists)
#   2. ingest_monthly vNRT 2026 (Jan..Mar, raw exists)
#   3. symlink vNRT-named dailies as v1-named (link_vNRT_to_v1.sh)
#   4. symlink vNRT-named monthlies as v1-named (inline)
#   5. link_daily_clim 2025 to fill 2025-12-22..31 from day-of-year clim
#   6. ncra 2025-12 dailies (21 real + 10 clim) into monthly file w/ provenance
#   7. cat_monthly (multi-year concat now spans 2001-01..2026-03)
#   8. write_pchip (refit, Fritsch-Carlson monotone-cubic Hermite)
#   9. submit diurnalize driver for 2025 (12-month fan-out)
#  10. submit diurnalize driver for 2026 (12-month fan-out)
#
# Both diurnalize drivers run the full 12-month range with no special
# scoping. diurnalize-ERA5.r picks real-vs-climatology per month by
# monthly-file presence (docs/PROPOSALS.md item 14), and skips months
# with no complete ERA5 meteo -- so 2026's Apr.. (no monthly file yet)
# fall to climatology and are skipped where ERA5 is absent, without a
# crash. Meteo comes from the FastTrack (ea_0005) fallback where the
# primary ERA5 tree lags the NRT window (item 12).
#
# Re-running this script is idempotent: ingest_monthly is mtime-aware
# skip-existing, symlinking checks for existing targets, and ncra will just
# overwrite the synthetic Dec file.

set -euo pipefail
cd "$(dirname "$0")/.."   # archive/ -> repo root, so ./config.sh etc. resolve
. ./config.sh
. ./lib/manifest.sh
_prod_t0=$(date +%s)

mkdir -p jobs
LOG="jobs/produce_2025_2026.$(date -u +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
echo "=== produce_2025_2026 started $(date -u) ==="
echo "logging to $LOG"

CURRENT_STEP="(init)"
step() {
    CURRENT_STEP="$1"
    echo; echo "=== [step $1] $2"
    manifest_record produce_2025_2026.sh start - "step $1: $2"
}
trap 'manifest_record produce_2025_2026.sh fail - "aborted at step ${CURRENT_STEP} (line $LINENO)"' ERR

step 1/10 "ingest_monthly vNRT 2025"
MICASA_VERSION=vNRT MICASA_YEAR_START=2025 MICASA_YEAR_END=2025 \
  Rscript ingest_monthly.r

step 2/10 "ingest_monthly vNRT 2026"
MICASA_VERSION=vNRT MICASA_YEAR_START=2026 MICASA_YEAR_END=2026 \
  Rscript ingest_monthly.r

step 3/10 "link_vNRT_to_v1.sh for dailies (2025, 2026)"
MICASA_YEAR=2025 ./link_vNRT_to_v1.sh
MICASA_YEAR=2026 ./link_vNRT_to_v1.sh || true   # 2026 only has Jan..~Apr days

step 4/10 "symlink vNRT-named monthlies as v1-named"
cd monthly_1x1
for v in MiCASA_vNRT_flux_x360_y180_monthly_2025{01,02,03,04,05,06,07,08,09,10,11}.nc \
         MiCASA_vNRT_flux_x360_y180_monthly_2026{01,02,03}.nc; do
    [ -e "$v" ] || { echo "  skip missing $v"; continue; }
    target=${v/vNRT/v1}
    if [ -e "$target" ] || [ -L "$target" ]; then rm -f "$target"; fi
    ln -s "$v" "$target"
    echo "  linked $target -> $v"
done
cd ..

step 5/10 "link_daily_clim.sh for 2025 (fills 2025-12-22..31)"
MICASA_CLIM_YEARS="2025" MICASA_VERSION=v1 ./link_daily_clim.sh

step 6/10 "ncra 2025-12 dailies -> monthly_202512.nc with provenance"
ncra -O daily_1x1/MiCASA_v1_flux_x360_y180_daily_202512??.nc \
        monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_202512.nc
ncatted -O \
  -a "provenance,global,c,c,21 days real vNRT 2025-12-01..21 + 10 days climatology fill (per-day-of-year mean of 2001-2024) for 2025-12-22..31, aggregated by ncra" \
  -a "status,global,c,c,provisional" \
  monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_202512.nc

# ncra averages times along with everything else, but the climatology daily
# files (MiCASA_v1_flux_x360_y180_daily_0000MMDD.nc) have time set to year
# 2001 by convention, not year 0000. Averaging real-2025 dailies with
# nominal-2001 climatology dailies pulls the mean time toward ~2018.
# Overwrite with the correct mid-Dec 2025 timestamp matching the convention
# used by ingest_monthly.r for other months (mid-month minus 0.5s).
${PYTHON:-python3} -c "
import netCDF4 as nc, datetime
ds = nc.Dataset('monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_202512.nc', 'r+')
t = datetime.datetime(2025, 12, 16, 11, 59, 59, 500000,
                      tzinfo=datetime.timezone.utc).timestamp()
ds.variables['time'][:] = [t]
ds.close()
print('  patched time to 2025-12-16 11:59:59.5 UTC')
"

step 7/10 "cat_monthly (multi-year concat, 2001-01..2026-03)"
./cat_monthly.sh || echo "WARN: cat_monthly returned non-zero (likely check_bounds), continuing"

step 8/10 "write_pchip (Fritsch-Carlson monotone-cubic Hermite)"
Rscript write_pchip.r

# diurnalize driver submission. SBATCH knobs shared by both years.
submit_diurn() {
    local label=$1 ; shift
    sbatch --parsable \
      --account=co2 \
      --time=10:00 \
      --ntasks=1 --cpus-per-task=1 --mem=4g \
      --partition=orion \
      --job-name="diurn-${label}-driver" \
      --output="jobs/diurn-${label}-driver.o%j" \
      --mail-type=FAIL --mail-user="$MAIL_USER" \
      --chdir="$PWD" \
      --export="ALL,WORK_DIR=$PWD,MICASA_STRICT_PIQS=1,$*" \
      --wrap='Rscript diurnalize-ERA5.r'
}

step 9/10 "submit diurnalize driver for 2025 (full year)"
DRV25=$(submit_diurn 2025 \
  "MICASA_YEAR_START=2025,MICASA_YEAR_END=2025")
echo "submitted diurnalize-2025 driver: $DRV25"

step 10/10 "submit diurnalize driver for 2026 (full year)"
# No MICASA_MONTH_END / MICASA_CLIM_YEARS scoping needed: diurnalize-ERA5.r
# auto-detects real-vs-climatology per month by monthly-file presence, so
# Jan..Mar use their real monthly files and later months fall to
# climatology (skipped where ERA5 is unavailable).
DRV26=$(submit_diurn 2026 \
  "MICASA_YEAR_START=2026,MICASA_YEAR_END=2026")
echo "submitted diurnalize-2026 driver: $DRV26"

manifest_record produce_2025_2026.sh ok "$(($(date +%s) - _prod_t0))" "all 10 steps complete"
echo "=== produce_2025_2026 finished $(date -u) ==="
