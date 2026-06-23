#!/bin/sh
# write_provenance.sh -- drop a human-readable PROVENANCE.txt into the output
# directory: what the data streams mean (v1 / vNRT / FastTrack), when and where
# it was generated, the run configuration, and the git commit.
#
# Complements (does not replace) the per-file CF/ACDD netCDF global attributes
# stamped by lib/provenance.r and the run log in $JOBS_DIR/run_manifest.tsv --
# this is the top-level orientation file for someone browsing the output dir.
#
# Usage:
#     ./write_provenance.sh [OUTPUT_DIR]
#
# OUTPUT_DIR defaults to $ERA5_DIR (the shipped hourly + daily-NEE product).
# Reads the exported MICASA_* environment (as set by config.sh / run_year.sh);
# it does NOT require MAIL_USER / BASE_DIR, so it is safe to run standalone, and
# it always exits 0 so it can never abort a pipeline run.

# Repo / work dir (holds .git and lib/). Default to this script's location.
WORK_DIR="${WORK_DIR:-$(cd "$(dirname "$0")" && pwd)}"

# Citation constants -- single source of truth, valid POSIX shell.
if [ -f "$WORK_DIR/lib/provenance.conf" ]; then
    . "$WORK_DIR/lib/provenance.conf"
fi

# Resolved configuration (exported env, falling back to the config.sh defaults).
VERSION="${MICASA_VERSION:-v1}"
Y0="${MICASA_YEAR_START:-2001}"
Y1="${MICASA_YEAR_END:-${MICASA_YEAR:-2025}}"
M0="${MICASA_MONTH_START:-1}"
M1="${MICASA_MONTH_END:-12}"
ERA5_DIR="${ERA5_DIR:-ERA5}"
DAILY_1X1_DIR="${DAILY_1X1_DIR:-daily_1x1}"
MONTHLY_1X1_DIR="${MONTHLY_1X1_DIR:-monthly_1x1}"
JOBS_DIR="${JOBS_DIR:-jobs}"
RESP_DRIVER="${MICASA_RESP_DRIVER:-airtemp}"
RESP_TEMPFUN="${MICASA_RESP_TEMPFUN:-q10}"
POLAR_CLIP="${MICASA_POLAR_CLIP:-conserve}"
METEO_ROOT="${MICASA_ERA5_DIR:-${CARBONTRACKER:-<CARBONTRACKER unset>}/METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100}"

OUTDIR="${1:-$ERA5_DIR}"
mkdir -p "$OUTDIR" 2>/dev/null
OUT="$OUTDIR/PROVENANCE.txt"

# Code version (matches lib/provenance.r's git logic).
GIT_COMMIT="$(git -C "$WORK_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
GIT_DESCRIBE="$(git -C "$WORK_DIR" describe --tags --always --dirty 2>/dev/null || echo unknown)"
GIT_BRANCH="$(git -C "$WORK_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
NOW="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo unknown)"
HOST="$(hostname 2>/dev/null || echo unknown)"
WHO="$(id -un 2>/dev/null || echo unknown)"

# Best-effort: echo a few global attrs off the newest output file, so the
# stamp reflects what is actually on disk (which may predate the current env).
SAMPLE="$(ls -t "$ERA5_DIR"/fluxes_*.nc 2>/dev/null | head -1)"
SAMPLE_ATTRS=""
if [ -n "$SAMPLE" ] && command -v ncdump >/dev/null 2>&1; then
    SAMPLE_ATTRS="$(ncdump -h "$SAMPLE" 2>/dev/null \
        | grep -E ':(respiration_temperature_(driver|function)|status|source_commit|software_version|date_created|meteo_fallback_used)' \
        | sed 's/^[[:space:]]*:/    /; s/ ;[[:space:]]*$//')"
fi

cat > "$OUT" <<EOF
================================================================================
 MiCASA processed carbon-flux product -- PROVENANCE
================================================================================
Pipeline  : ${MICASA_PROV_PIPELINE:-MiCASA-processing}
Source    : ${MICASA_PROV_PIPELINE_URL:-https://github.com/noaa-gml/MiCASA-processing}
Institute : ${MICASA_PROV_INSTITUTION:-NOAA Global Monitoring Laboratory}
Upstream  : MiCASA Land Carbon Flux v1 (NASA GSFC; Weir et al.)

This directory holds hourly 1-degree GPP / respiration / NEE (fluxes_YYYYMM.nc)
and per-day NEE (MiCASA_<ver>.nee.YYYYMMDD.nc) derived from the MiCASA model for
ingestion by NOAA CarbonTracker. Per-file detail lives in the CF/ACDD netCDF
global attributes ("ncdump -h <file>") and in ${JOBS_DIR}/run_manifest.tsv.

-- Data streams --------------------------------------------------------------
 v1         Final / authoritative MiCASA stream. Lags the source data by some
            weeks; this is the production-quality product.
 vNRT (NRT) Near-real-time stream: available within days of the source data,
            but may be revised once the v1 version lands and supersedes it.
            (Exposed to CarbonTracker as v1 via link_vNRT_to_v1.sh.)
 FastTrack  ERA5 meteo fallback (the "ea_0005" tree), used only for the NRT
            trailing window where the primary ERA5 tree is not yet populated.
            Resolved per day; each output records which tree fed it in its
            "meteo_source_*" global attributes (meteo_fallback_used=TRUE if any
            day used FastTrack). A month missing some ERA5 days is written with
            status="provisional".

-- Generated -----------------------------------------------------------------
 when : $NOW
 host : $HOST
 by   : $WHO

-- Configuration -------------------------------------------------------------
 version             : $VERSION
 year range          : $Y0 .. $Y1
 month range         : $M0 .. $M1
 respiration driver  : $RESP_DRIVER   (MICASA_RESP_DRIVER; soiltemp = opt-in)
 resp. response fn   : $RESP_TEMPFUN   (MICASA_RESP_TEMPFUN; lloydtaylor = opt-in)
 polar-night clip    : $POLAR_CLIP   (MICASA_POLAR_CLIP; plain = legacy zero-clip)
 sub-monthly fitter  : recorded in fit.piqs.rda (piqsfit.meta\$fitter) and the
                       per-file attributes; default pchip (PCHIP-on-cumulative)
 ERA5 meteo (input)  : $METEO_ROOT
 output locations    : ERA5_DIR        = $ERA5_DIR
                       DAILY_1X1_DIR   = $DAILY_1X1_DIR
                       MONTHLY_1X1_DIR = $MONTHLY_1X1_DIR
                       JOBS_DIR        = $JOBS_DIR

-- Code version --------------------------------------------------------------
 git commit   : $GIT_COMMIT
 git describe : $GIT_DESCRIBE
 git branch   : $GIT_BRANCH
EOF

if [ -n "$SAMPLE_ATTRS" ]; then
    {
        echo ""
        echo "-- Attributes on the newest output file ($(basename "$SAMPLE")) --"
        echo "$SAMPLE_ATTRS"
    } >> "$OUT"
fi

cat >> "$OUT" <<EOF

-- Citation ------------------------------------------------------------------
 DOI          : ${MICASA_DOI:-PENDING}
 license      : ${MICASA_PROV_LICENSE:-CC0-1.0}
 conventions  : ${MICASA_PROV_CONVENTIONS:-CF-1.10, ACDD-1.3}
================================================================================
EOF

echo "Wrote $OUT"
exit 0
