README written by Ash Pera 2025-05-13 17:03:37
Updated 2026-04-26 — Tier-1 refactor + latent-bug sweep documented

##########################
# Overview
##########################

This directory holds scripts to take raw MiCASA data and process it for use in CarbonTracker
https://gml.noaa.gov/ccgg/carbontracker/documentation.php

MiCASA Land Carbon Flux
Global, daily and monthly mean 0.1 degree resolution carbon fluxes from net primary production (NPP),
heterotrophic respiration (Rh), wildfire emissions (FIRE), fuel wood burning emissions (FUEL),
net ecosystem exchange (NEE), and net biosphere exchange (NBE) derived from the MiCASA model, version 1
https://earth.gov/ghgcenter/data-catalog/micasa-carbonflux-grid-v1

##########################
# Quick start
##########################

Run the whole pipeline for one year:

    ./run_year.sh 2026                    # full v1 pipeline for 2026
    ./run_year.sh 2026 vNRT               # near-real-time stream
    ./run_year.sh 2026 v1 --skip-download # data already on disk
    ./run_year.sh 2026 --dry-run          # show stages without running

Stage skip flags: --skip-download, --skip-ingest, --skip-aggregate,
--skip-piqs, --skip-diurnalize, --skip-daysplit. SBATCH stages are
submitted with --wait so the driver blocks until completion.

##########################
# Versions: v1 vs vNRT (the hybrid stream)
##########################

MiCASA publishes two parallel streams from the same upstream pipeline:

  * v1     The final/authoritative stream. Lags the source data by some
           weeks; what you want for production-quality NEE.
  * vNRT   Near-real-time. Available within days of the source data, but
           may be revised once v1 lands and supersedes it.

Both streams use version-tagged basenames so they coexist in one tree:

    portal.nccs.nasa.gov/daily/<YYYY>/<MM>/
        MiCASA_v1_flux_x3600_y1800_daily_<YYYYMMDD>.nc4    ← preferred
        MiCASA_vNRT_flux_x3600_y1800_daily_<YYYYMMDD>.nc4  ← fallback
        …_sha256.txt                                       (also versioned)

Operational policy:  use v1 wherever it exists, fall back to vNRT for
the trailing window where v1 has not yet been published. Concretely:

  1. Download both streams for the current year:
        MICASA_VERSION=both ./download.sh
     (v1 is downloaded first; --no-clobber means vNRT only fills gaps.)

  2. Ingest with whichever version you intend to use *for the gap days*.
     v1 days that already exist on disk will be ingested as v1; vNRT
     days fill the rest. Set MICASA_VERSION=v1 to ingest both as v1
     (the file basename then determines which you read), or
     MICASA_VERSION=vNRT to keep the vNRT-tagged outputs.

  3. After ingest, expose vNRT-tagged 1° outputs as v1-tagged for any
     downstream consumer that doesn't know about vNRT:
        ./link_vNRT_to_v1.sh                  (per current year)
     This skips days where MiCASA_v1_*.nc already exists, so v1 always
     wins.

  4. When the upstream v1 record catches up, re-run ingest for those
     days with MICASA_VERSION=v1 — the v1 outputs replace the symlinks,
     and CarbonTracker silently switches over.

##########################
# Configuration (env-driven)
##########################

config.sh and config.r read the same environment variables, so any knob
can be set on the command line or in run_year.sh. Defaults match the
operational 2025+ setup.

  MICASA_YEAR          single-year focus, used by SBATCH workers and
                       by run_year.sh
  MICASA_VERSION       v1 (default), vNRT, or both (download.sh only)
  MICASA_YEAR_START    first year for multi-year stages (default 2001)
  MICASA_YEAR_END      last year  (default 2025)
  MICASA_MONTH_START   first month for diurnalize (default 1)
  MICASA_MONTH_END     last month  (default 12)
  MICASA_CLIM_YEARS    space-separated years that should use day-of-year
                       climatology instead of real ERA5 data.
                       Default: "2000 <current calendar year>"
                       — i.e. the years where ERA5 is either not yet
                       (pre-2000) or not yet fully (current NRT year)
                       available. Independent of MICASA_YEAR so
                       backfilling an earlier year doesn't accidentally
                       clim it.
  MAIL_USER            SBATCH --mail-user
  BASE_DIR             /work2/noaa/co2/GFED-CASA   (parent of YYYY trees)
  WORK_DIR             auto-detected from script path; this directory
  PORTAL_URL_BASE      NCCS download base URL
  DAILY_1X1_DIR        daily_1x1
  MONTHLY_1X1_DIR      monthly_1x1
  ERA5_DIR             ERA5
  RAW_SRC_DIR          portal.nccs.nasa.gov   (raw 0.1° mirror)
  JOBS_DIR             jobs

Runtime-only knobs:
  INGEST_YEAR          set by ingest_byyear.r driver; do not set manually
  diurn_year           set by diurnalize-ERA5.r driver; do not set manually
  RECOMPUTE_EXISTING   1 to force ingest_monthly.r to re-write existing
                       outputs (default: skip them)
  MICASA_NO_BLESS_REFERENCE
                       1 to skip auto-blessing this year's downloaded
                       file as next year's reference in check_unchanged.sh

##########################
# Flowchart
##########################

run_year.sh
   |
   |--- link_old_micasa_raw.sh (auto-detect 2024 vs 2025+ layout)
   |
   |--- download.sh ---> portal.nccs.nasa.gov/{daily,monthly}/YYYY/...
   |       check_daily_downloads.r
   |       check_hashes.py            (year range from MICASA_YEAR_*)
   |       check_unchanged.sh         (auto-blesses next year's reference)
   |
   |--- ingest_monthly.r       ingest_byyear.r
   |       (both source lib/ingest_common.r — area-weighted aggregator
   |        bug-fixed 2026-04-26; see lib/test_aggregate.r)
   |
   |--- cat_monthly.sh                compute_daily_clim.sh
   |       check_bounds.sh                    |
   |                                  link_daily_clim.sh
   |       compute_clim.sh
   |       write_piqs.r
   |
   |--- diurnalize-ERA5.r
   |
   `--- daysplitter.sh
                |
        link_vNRT_to_v1.sh  (only if MICASA_VERSION=vNRT was used)

##########################
# Programs
##########################

run_year.sh:
    Top-level driver. Sets MICASA_YEAR, sources config.sh, and calls
    each pipeline stage in order. SBATCH stages submitted with --wait.

config.sh / config.r:
    Single source of truth for env-driven knobs (see Configuration).
    Sourced by every other script.

lib/ingest_common.r:
    Shared helpers between ingest_byyear.r and ingest_monthly.r.
    Defines: archimedes(), compute.gca(), aggregate.to.1x1(),
    micasa.dim.lon(), micasa.dim.lat(), micasa.time.dim(),
    make.tracer.vars(), write.netcdf(), and the constants
    micasa.tracers and EARTH_RADIUS_M.

lib/test_aggregate.r:
    Self-contained Rscript test harness for aggregate.to.1x1, with a
    regression test against the pre-2026-04-26 buggy implementation.

download.sh:
    wget MiCASA daily + monthly files for $MICASA_YEAR / $MICASA_VERSION
    from the NCCS portal. Supports MICASA_VERSION=both for the hybrid
    v1+vNRT stream (v1 downloads first, then vNRT fills gaps via
    --no-clobber). Writes to portal.nccs.nasa.gov/{daily,monthly}/<YYYY>/
    in the work directory.

check_daily_downloads.r:
    Verify NPP, Rh, FIRE, FUEL exist for every day in
    [MICASA_YEAR_START, MICASA_YEAR_END].

check_hashes.py:
    Verify SHA-256 of each downloaded file. Year range from
    MICASA_YEAR_START / MICASA_YEAR_END (silently skips years outside).
    Handles both v1 and vNRT files in the same directory.

check_unchanged.sh:
    Diff ncdump -h headers of new vs reference (previous year's tree).
    Catches silent provider-side metadata changes (e.g. the 2018 kg→g
    units flip). On a clean diff, automatically blesses this year's
    first daily/monthly as the *next* year's reference, so the chain
    self-bootstraps after the initial 2024 reference.
    Set MICASA_NO_BLESS_REFERENCE=1 to skip the auto-bless step.

check_bounds.sh:
    Simple unweighted-area average sanity check, called by cat_monthly.sh.
    NOT used in production aggregation (that's aggregate.to.1x1).

ingest_byyear.r:
    For a given INGEST_YEAR, aggregate every day's raw 0.1° NPP/Rh/FIRE/
    FUEL to 1° via lib/ingest_common.r:aggregate.to.1x1, write to
    daily_1x1/MiCASA_<VER>_flux_x360_y180_daily_<YYYYMMDD>.nc.
    Driver mode (no INGEST_YEAR): fans out one SBATCH per year in
    [MICASA_YEAR_START, MICASA_YEAR_END].

ingest_monthly.r:
    Plain year-loop monthly aggregator. Skips outputs that already
    exist unless RECOMPUTE_EXISTING=1.

cat_monthly.sh:
    Concatenate monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly_<YYYYMM>.nc
    into a single time-stacked monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly.nc.
    Runs check_bounds.sh.

compute_clim.sh:
    Ferret-driven modulo-month average of the concatenated monthly file.
    Writes monthly_1x1/{NPP,Rh}clim.nc.

compute_daily_clim.sh:
    ncea across-year average per day-of-year, writing
    daily_1x1/MiCASA_<VER>_flux_x360_y180_daily_0000<MMDD>.nc.

link_daily_clim.sh:
    For each year in $MICASA_CLIM_YEARS (default: 2000 + current
    calendar year), symlink missing daily files to the 0000<MMDD> clim.

link_old_micasa_raw.sh:
    Auto-detect the previous year's raw layout (legacy from_weir/...
    or current portal.nccs.nasa.gov/...) and absolute-path-symlink
    daily/monthly into this year's tree, range
    [MICASA_YEAR_START, MICASA_YEAR-1]. New layouts can be added
    to its layout_candidates array.

link_old_micasa_finals.sh:
    Same idea but for the 1° outputs.

link_vNRT_to_v1.sh:
    Symlink ingested vNRT daily files as v1-named files for the same
    year. Run this once vNRT-stream ingest is complete so downstream
    consumers (CarbonTracker, etc.) read MiCASA_v1_*.nc transparently.
    Skips days where MiCASA_v1_*.nc already exists, so v1 always wins.

write_piqs.r:
    Source config.r and load monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly.nc
    via micasa.out.monthly.cat(cfg). Per grid cell, fit GPP and rtot
    with piecewise integral quadratic splines (PIQS), save to
    fit.piqs.rda.
    https://gml.noaa.gov/ccgg/carbontracker/documentation.php#tth_sEc2.2

diurnalize-ERA5.r:
    Apply ERA5 hourly meteo (ssrd, t2m, stl1, swvl1) to the PIQS-smoothed
    monthly fluxes to get hourly GPP/RESP/NEE per (year, month).
    Writes ERA5/fluxes_<YYYYMM>.nc.
    Driver mode (no diurn_year): fans out per year in
    [MICASA_YEAR_START, MICASA_YEAR_END].
    Years in MICASA_CLIM_YEARS use day-of-year climatology (NPPclim.nc,
    Rhclim.nc) instead of monthly real data.

daysplitter.sh:
    For every (year, month) in ERA5/fluxes_<YYYYMM>.nc, split into
    ERA5/MiCASA_v1.nee.<YYYYMMDD>.nc dailies, keeping only NEE.
    Range from MICASA_YEAR_START..END and MICASA_MONTH_START..END.

Deprecated / kept-for-reference:
    ingest.r                       — superseded by ingest_byyear/monthly
    download_and_check.sh          — superseded by run_year.sh stage 1
    test_gca.r                     — geometry sanity test (pre-refactor)
    create_era5_move.py            — one-time data-move script
    download-NRT.sh                — merged into download.sh

##########################
# Latent-bug sweep — 2026-04-26
##########################

Six bugs found and fixed during the post-refactor audit (commits
9ea6970 and following):

  1. lib/ingest_common.r:aggregate.to.1x1
     Latitude-area weights were being recycled column-major across a
     10×10 sub-block, applying them along the LONGITUDE axis instead
     of latitude. The inner `for (inlon in inlons)` loop was also dead
     (×10 then ÷10). Fix: build a flat length-100 weight vector that
     correctly assigns gca[inlats[k]] to every cell at lat-position k.
     Magnitude depends on field gradient within a 1° block — typically
     <0.01% for smooth fields, growing toward the poles. See
     lib/test_aggregate.r for verification + regression test.

  2. run_year.sh:sbatch_wait
     `--export="ALL,${exports}"` produced a trailing comma when called
     with empty exports (e.g. for ingest_monthly.r). Now passes
     "ALL" alone in that case.

  3. write_piqs.r
     load.ncdf() path was hardcoded to MiCASA_v1_*.nc. Now sources
     config.r and uses micasa.out.monthly.cat(cfg), so it works under
     MICASA_VERSION=vNRT too.

  4. check_hashes.py
     Directory glob 202[4-5] silently skipped any year outside
     2024-2025 — verification would pass vacuously for 2026+. Now
     reads MICASA_YEAR_START/END from env and globs the requested
     years. Also added a missing-checksum-file warning.

  5. link_old_micasa_raw.sh
     Hardcoded the legacy from_weir/portal.nccs.nasa.gov/... path,
     which only existed in the 2024 layout. Now auto-detects between
     legacy and 2025+ layouts via a `layout_candidates` array, and
     uses absolute paths so the link survives WORK_DIR moves.

  6. check_unchanged.sh
     Used to silently warn-and-continue when the previous-year
     reference was missing; new years would slip through unchecked.
     Now: clearer warning with bootstrap instructions, and on a
     successful clean diff it auto-blesses the new year's file as
     next year's reference — chain bootstraps itself once the initial
     2024 reference is in place.

Also: link_daily_clim.sh and diurnalize-ERA5.r now default
MICASA_CLIM_YEARS to "2000 <current calendar year>" instead of
"2000 $MICASA_YEAR" — climatology fallback should track *what's
missing on disk right now*, not which year you happen to be
processing.

##########################
# Data layout
##########################

portal.nccs.nasa.gov/{daily/<YYYY>/<MM>,monthly/<YYYY>}/
                MiCASA_<VER>_flux_x3600_y1800_<freq>_<...>.nc4
    Created by download.sh. Both v1 and vNRT files coexist here.

daily_1x1/MiCASA_<VER>_flux_x360_y180_daily_<YYYYMMDD>.nc
    Created by ingest_byyear.r (1° area-weighted aggregate of NPP,
    Rh, FIRE, FUEL).

daily_1x1/MiCASA_<VER>_flux_x360_y180_daily_0000<MMDD>.nc
    Day-of-year climatology, created by compute_daily_clim.sh.

monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly_<YYYYMM>.nc
    Created by ingest_monthly.r.

monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly.nc
    Concatenated multi-year monthly file (cat_monthly.sh).

monthly_1x1/{NPP,Rh}clim.nc
    Climatology, created by compute_clim.sh (Ferret).

ERA5/fluxes_<YYYYMM>.nc
    Hourly diurnalized monthly file, created by diurnalize-ERA5.r.

ERA5/MiCASA_v1.nee.<YYYYMMDD>.nc
    Daily NEE-only files, created by daysplitter.sh.

##########################
# Extra Notes
##########################

  - Climatology fallback applies to (a) years before ERA5 starts (2000
    and earlier) and (b) the current calendar year (NRT phase, ERA5
    not yet fully published). See MICASA_CLIM_YEARS above.

  - Diurnalize is described in the CT documentation. It needs monthly
    Rh and NPP, from which it generates temporally-downscaled GPP and
    total respiration.

  - Fire and (bio)fuel emissions are taken from the daily files
    provided by MiCASA. This is the only source for fairly
    high-resolution-in-time emissions for those processes.

  - MiCASA provides temporally-downscaled (non-fire and -fuel) fluxes,
    using a method similar to ours, but with different meteorology
    (NASA's MERRA2 reanalysis). Our method is a little bit better due
    to the PIQS part of the scheme, which smooths out abrupt changes
    at monthly boundaries, and our meteo comes from ERA5. That makes
    the downscaling consistent with the atmospheric transport provided
    by TM5. So, we do not use the MiCASA temporally-downscaled fluxes;
    we start with the monthlies and apply the downscaling ourselves.

https://nco.sourceforge.net/nco.pdf

ncea (netCDF Ensemble Average)
    performs gridpoint averages of variables across an arbitrary number
    (an ensemble) of input files, with each file receiving an equal
    weight in the average.
    -O, overwrite output if it exists
    Note: ncea is deprecated for nces (netCDF Ensemble Statistics)

ncks (netCDF Kitchen Sink)
    extracts (a subset of the) data from input-file, regrids it
    according to map-file if specified, then writes in netCDF format
    to output-file.

VSEM-ET (possibly unrelated, but nice picture)
    https://insightmaker.com/insight/6DkHwGgVTkedbnCviUX8bD/Clone-of-Very-Simple-Ecosystem-Model-with-Evapotranspiration-VSEM-ET


Original daily header dump:
netcdf MiCASA_v1_flux_x3600_y1800_daily_20130507 {
dimensions:
	lat = 1800 ;
	lon = 3600 ;
	time = UNLIMITED ; // (1 currently)
	nv = 2 ;
variables:
	double lat(lat) ;
		lat:units = "degrees_north" ;
		lat:long_name = "latitude" ;
	double lon(lon) ;
		lon:units = "degrees_east" ;
		lon:long_name = "longitude" ;
	double time(time) ;
		time:units = "days since 1980-01-01" ;
		time:long_name = "time" ;
		time:bounds = "time_bnds" ;
	double time_bnds(time, nv) ;
		time_bnds:units = "days since 1980-01-01" ;
		time_bnds:long_name = "time bounds" ;
	float NPP(time, lat, lon) ;
		NPP:units = "kg m-2 s-1" ;
		NPP:expressed_as = "carbon" ;
		NPP:long_name = "Net primary productivity" ;
	float Rh(time, lat, lon) ;
		Rh:units = "kg m-2 s-1" ;
		Rh:expressed_as = "carbon" ;
		Rh:long_name = "Heterotrophic respiration" ;
	float FIRE(time, lat, lon) ;
		FIRE:units = "kg m-2 s-1" ;
		FIRE:expressed_as = "carbon" ;
		FIRE:long_name = "Fire emission" ;
	float FUEL(time, lat, lon) ;
		FUEL:units = "kg m-2 s-1" ;
		FUEL:expressed_as = "carbon" ;
		FUEL:long_name = "Fuel wood emission" ;
	float ATMC(time, lat, lon) ;
		ATMC:units = "kg m-2 s-1" ;
		ATMC:expressed_as = "carbon" ;
		ATMC:long_name = "Atmospheric correction" ;
	float NEE(time, lat, lon) ;
		NEE:units = "kg m-2 s-1" ;
		NEE:expressed_as = "carbon" ;
		NEE:long_name = "Net ecosystem exchange" ;

// global attributes:
		:Conventions = "CF-1.9" ;
		:contact = "Brad Weir <brad.weir@nasa.gov>" ;
		:institution = "NASA Goddard Space Flight Center" ;
		:title = "MiCASA Daily NPP Rh ATMC NEE FIRE FUEL Fluxes 0.1 degree x 0.1 degree v1" ;
		:LongName = "MiCASA Daily NPP Rh ATMC NEE FIRE FUEL Fluxes 0.1 degree x 0.1 degree" ;
		:ShortName = "MICASA_FLUX_D" ;
		:VersionID = "1" ;
		:GranuleID = "MiCASA_v1_flux_x3600_y1800_daily_20130507.nc4" ;
		:Format = "netCDF" ;
		:ProcessingLevel = "4" ;
		:IdentifierProductDOIAuthority = "https://doi.org/" ;
		:IdentifierProductDOI = "10.5067/ZBXSA1LEN453" ;
		:ReadMeURL = "https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/MiCASA_README.pdf" ;
		:RangeBeginningDate = "2013-05-07" ;
		:RangeBeginningTime = "00:00:00.000000" ;
		:RangeEndingDate = "2013-05-07" ;
		:RangeEndingTime = "23:59:59.999999" ;
		:NorthernmostLatiude = "90.0" ;
		:WesternmostLongitude = "-180.0" ;
		:SouthernmostLatitude = "-90.0" ;
		:EasternmostLongitude = "180.0" ;
		:comment = "Positive NPP indicates uptake by vegetation. Positive Rh indicates emission to the atmosphere. NEE = Rh - NPP - ATMC, and NBE = NEE + FIRE + FUEL. ATMC adjusts net exchange to account for missing processes and better match long-term atmospheric budgets." ;
		:ProductionDateTime = "2024-09-23T01:58:21Z" ;
}

##########################
# Performance: compression-level tuning — 2026-04-26
##########################

`diurnalize-ERA5.r` writes 12 ~660 MB files per year (9 hourly vars at
1°). Default deflate level was 9; bench results on a real
fluxes_202401.nc (lib/bench_compression_diurnal.r):

  level  time/file   size_MB   per-year writes
    9      108 s      632       1298 s   (= reference)
    6       72 s      633        870 s   (-33%)
    4       65 s      634        786 s   (-39%, +0.3% size)
    3       60 s      646        715 s   (-45%, +2.2% size)
    1       55 s      646        654 s   (-50%, +2.2% size)

Chose level 4: nearly identical file size to level 9 (+0.3%) for ~9 min
saved per year on the diurnalize stage. Levels 1-3 buy another ~2 min
but cost +14 MB/file = +170 MB/year, not worth it for archived output.

Ingest paths (`lib/ingest_common.r`, `ingest.r`) left at level 9
because per-file output is only ~164 KB and the prior bench
(lib/bench_compression.r) showed ~9 s/year savings — not worth the
file-size cost for users who pull the daily 1° aggregates.

##########################
# Performance: ingest_byyear skip-existing + read-only-needed — 2026-04-26
##########################

Two changes to `ingest_byyear.r` (and a smaller one to `ingest_monthly.r`):

  1. Skip-existing — `RECOMPUTE_EXISTING=1` to override (default off).
     A daily NRT cycle that adds 1 new day previously deleted and
     rebuilt all 365 daily 1° outputs. Now: re-run skips finished
     days, processes only what's missing.

  2. Read only the 4 needed tracers (NPP, Rh, FIRE, FUEL) instead of
     the full 6-var raw file (which also has ATMC and NEE). Done by
     passing `vars = micasa.tracers` to `load.ncdf()`.

Measured impact on ingest_byyear 2024 (full year, 366 days):

  | run                                  |  wall-time |
  | ------------------------------------ | ---------- |
  | baseline (vectorized aggregator)     |     610 s  |
  | + read-only-needed (RECOMPUTE=1)     |     504 s  | -17%
  | + skip-existing (cached re-run)      |       4 s  | -99%

Output is bit-identical (`ncdiff` on 4 sample days × 4 tracers: max
|Δ| = 0). Only the `:history` attribute timestamp differs on rewrite,
as expected.

The vectorized aggregator (commit ce1bccc) was the big win that
collapsed ingest_byyear from 3.6 hr to ~10 min/year. These two
changes shave another ~17% of the throughput case and ~99% of the
NRT-rerun case.

Verified by:
  - lib/test_ingest_bitident.r  — read-path bit-identity
  - lib/profile_ingest_day.r    — per-step cost breakdown
  - lib/test_aggregate.r        — aggregator regression test (earlier)
