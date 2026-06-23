#!/bin/sh
# write_provenance.sh -- drop a human-readable PROVENANCE.txt into the output
# directory: what the data are (variables / units / sign conventions / files),
# what the streams mean (v1 / vNRT / FastTrack), what is on disk, the run
# configuration, the upstream source, and the git commit.
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

# Citation / upstream constants -- single source of truth, valid POSIX shell.
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
PORTAL="${PORTAL_URL_BASE:-https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA}"
METEO_ROOT="${MICASA_ERA5_DIR:-${CARBONTRACKER:-<CARBONTRACKER unset>}/METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100}"

OUTDIR="${1:-$ERA5_DIR}"
mkdir -p "$OUTDIR" 2>/dev/null
OUT="$OUTDIR/PROVENANCE.txt"

# Code version (matches lib/provenance.r's git logic).
GIT_COMMIT="$(git -C "$WORK_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
GIT_DESCRIBE="$(git -C "$WORK_DIR" describe --tags --always --dirty 2>/dev/null || echo unknown)"
GIT_BRANCH="$(git -C "$WORK_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
case "$GIT_DESCRIBE" in
    *-dirty|*dirty) DIRTY_WARN="
 *** WARNING: built from a MODIFIED working tree (git reports -dirty); the commit
 ***          above does not fully capture the code that produced these files. ***" ;;
    *)              DIRTY_WARN="" ;;
esac

NOW="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo unknown)"
HOST="$(hostname 2>/dev/null || echo unknown)"
WHO="$(id -un 2>/dev/null || echo unknown)"

# Software stack (best-effort; one R startup, skipped cleanly if R is absent).
R_INFO="(R not on PATH; not captured)"
if command -v Rscript >/dev/null 2>&1; then
    _ri="$(Rscript -e 'cat(R.version.string); cat(" | ncdf4", tryCatch(as.character(packageVersion("ncdf4")), error = function(e) "n/a"))' 2>/dev/null)"
    [ -n "$_ri" ] && R_INFO="$_ri"
fi

# On-disk coverage (cheap: filename parse, no per-file ncdump).
fmt_ym() { echo "$1" | sed -E 's/.*fluxes_([0-9]{4})([0-9]{2})\.nc$/\1-\2/'; }
N_HOURLY="$(ls "$ERA5_DIR"/fluxes_*.nc 2>/dev/null | wc -l | tr -d ' ')"
FIRST="$(ls "$ERA5_DIR"/fluxes_*.nc 2>/dev/null | head -1)"
LAST="$(ls "$ERA5_DIR"/fluxes_*.nc 2>/dev/null | tail -1)"
[ -n "$FIRST" ] && FIRST="$(fmt_ym "$FIRST")" || FIRST="(none)"
[ -n "$LAST" ]  && LAST="$(fmt_ym "$LAST")"   || LAST="(none)"
N_DAILY="$(ls "$ERA5_DIR"/MiCASA_*.nee.*.nc 2>/dev/null | wc -l | tr -d ' ')"

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

Hourly 1-degree GPP / respiration / NEE derived from the MiCASA land model for
ingestion by NOAA CarbonTracker. Per-file detail lives in the CF/ACDD netCDF
global attributes ("ncdump -h <file>") and in ${JOBS_DIR}/run_manifest.tsv.

-- Contents & conventions ----------------------------------------------------
 Files:
   fluxes_YYYYMM.nc                hourly, all components + carried meteo
   MiCASA_<ver>.nee.YYYYMMDD.nc    per-day NEE only -- what CarbonTracker ingests
 Variables (hourly file): GPP, resp, NEE; QGPP/qresp (sub-monthly fit terms);
   carried meteo ssr (ssrd), t2m, stl1, swvl1.
 Units : mol m-2 s-1 (carbon flux) -- see each variable's "units" attribute.
 Signs : NEE = Rh - NPP.  POSITIVE = source to the atmosphere (net release);
         negative = net uptake. GPP is carried NEGATIVE (uptake); resp positive.
 Grid  : 1 deg x 1 deg global, 360 lon x 180 lat.

-- Data streams --------------------------------------------------------------
 v1         Final / authoritative MiCASA stream. Lags the source by some weeks;
            this is the production-quality product.
 vNRT (NRT) Near-real-time stream: available within days of the source, but may
            be revised once the v1 version lands and supersedes it. (Exposed to
            CarbonTracker as v1 via link_vNRT_to_v1.sh.)
 FastTrack  ERA5 meteo fallback (the "ea_0005" tree), used only for the NRT
            trailing window where the primary ERA5 tree is not yet populated.
            Resolved per day; each output records which tree fed it in its
            "meteo_source_*" attributes (meteo_fallback_used=TRUE if any day did).
            A month missing some ERA5 days is written with status="provisional".

-- Coverage (on disk, this directory) ----------------------------------------
 hourly files : $N_HOURLY  ($FIRST .. $LAST)
 daily files  : $N_DAILY  (per-day NEE)
 Climatology-filled and provisional (partial-meteo) months are flagged per file
 (status / meteo_partial attributes), not here.

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
 output locations    : ERA5_DIR=$ERA5_DIR  DAILY_1X1_DIR=$DAILY_1X1_DIR
                       MONTHLY_1X1_DIR=$MONTHLY_1X1_DIR  JOBS_DIR=$JOBS_DIR

 reproduce (per year YYYY):
   MICASA_VERSION=$VERSION MICASA_RESP_DRIVER=$RESP_DRIVER \\
   MICASA_RESP_TEMPFUN=$RESP_TEMPFUN MICASA_POLAR_CLIP=$POLAR_CLIP \\
   ./run_year.sh YYYY            # + export the output-dir vars above as needed

-- Inputs --------------------------------------------------------------------
 upstream model : ${MICASA_UPSTREAM_NAME:-MiCASA Land Carbon Flux v1 (NASA GSFC)}
 upstream DOI   : ${MICASA_UPSTREAM_DOI:-10.5067/ZBXSA1LEN453}
 download from  : $PORTAL
 ERA5 meteo     : $METEO_ROOT
 (per-file input paths + SHA-256 checksums are in the netCDF attributes.)

-- Caveats -------------------------------------------------------------------
 * Fire and fuel-wood emissions are passed through from the MiCASA daily product
   -- they are NOT fitted or diurnalized here.
 * NEE excludes the MiCASA ATMC term (NEE = Rh - NPP, not - ATMC): ATMC is an
   atmospheric-inversion correction; subtracting it would double-count the
   constraint a downstream CO2 inversion applies. See docs/METHODOLOGY.md.
 * Single deterministic realization -- no per-pixel uncertainty is provided.
 * Monthly means are preserved exactly (fit + mass-conserving polar clip); the
   sub-monthly and diurnal SHAPE is reconstructed, not native to MiCASA.

-- Code version --------------------------------------------------------------
 git commit   : $GIT_COMMIT
 git describe : $GIT_DESCRIBE
 git branch   : $GIT_BRANCH
 software     : $R_INFO$DIRTY_WARN
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
 DOI          : ${MICASA_DOI:-PENDING}    (this processed product)
 license      : ${MICASA_PROV_LICENSE:-CC0-1.0}
 conventions  : ${MICASA_PROV_CONVENTIONS:-CF-1.10, ACDD-1.3}
================================================================================
EOF

echo "Wrote $OUT"
exit 0
