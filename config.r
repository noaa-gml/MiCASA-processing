## MiCASA pipeline configuration for R scripts.
## Source this from any R script in this repo:
##     source(file.path(Sys.getenv("WORK_DIR", getwd()), "config.r"))
##
## All knobs come from environment variables (set by config.sh / SBATCH
## --export=ALL). Defaults match config.sh so an R script can be run
## standalone for testing.

micasa.config <- function() {
  list(
    year         = as.integer(Sys.getenv("MICASA_YEAR",        "2025")),
    version      = Sys.getenv("MICASA_VERSION",                "v1"),
    year.start   = as.integer(Sys.getenv("MICASA_YEAR_START",  "2001")),
    year.end     = as.integer(Sys.getenv("MICASA_YEAR_END",    "2025")),
    month.start  = as.integer(Sys.getenv("MICASA_MONTH_START", "1")),
    month.end    = as.integer(Sys.getenv("MICASA_MONTH_END",   "12")),
    mail.user    = Sys.getenv("MAIL_USER",                     "ashley.pera@noaa.gov"),
    base.dir     = Sys.getenv("BASE_DIR",                      "/work2/noaa/co2/GFED-CASA"),
    work.dir     = Sys.getenv("WORK_DIR",                      getwd()),
    portal.url   = Sys.getenv("PORTAL_URL_BASE",
                              "https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA"),
    daily.1x1    = Sys.getenv("DAILY_1X1_DIR",   "daily_1x1"),
    monthly.1x1  = Sys.getenv("MONTHLY_1X1_DIR", "monthly_1x1"),
    era5.dir     = Sys.getenv("ERA5_DIR",        "ERA5"),
    raw.dir      = Sys.getenv("RAW_SRC_DIR",     "portal.nccs.nasa.gov"),
    jobs.dir     = Sys.getenv("JOBS_DIR",        "jobs")
  )
}

## ---- Filename helpers ------------------------------------------------------
## All filenames embed the version (v1 or vNRT) so callers don't need to
## special-case it.

# Raw 0.1° daily file as downloaded from the NCCS portal.
micasa.raw.daily <- function(cfg, year, month, day) {
  sprintf("%s/daily/%d/%02d/MiCASA_%s_flux_x3600_y1800_daily_%d%02d%02d.nc4",
          cfg$raw.dir, year, month, cfg$version, year, month, day)
}

# Raw 0.1° monthly file as downloaded from the NCCS portal.
micasa.raw.monthly <- function(cfg, year, month) {
  sprintf("%s/monthly/%d/MiCASA_%s_flux_x3600_y1800_monthly_%d%02d.nc4",
          cfg$raw.dir, year, cfg$version, year, month)
}

# Aggregated 1° daily file written by ingest_byyear.r.
micasa.out.daily <- function(cfg, year, month, day) {
  sprintf("%s/MiCASA_%s_flux_x360_y180_daily_%d%02d%02d.nc",
          cfg$daily.1x1, cfg$version, year, month, day)
}

# Aggregated 1° monthly file written by ingest_monthly.r.
micasa.out.monthly <- function(cfg, year, month) {
  sprintf("%s/MiCASA_%s_flux_x360_y180_monthly_%d%02d.nc",
          cfg$monthly.1x1, cfg$version, year, month)
}

# Concatenated multi-year monthly file (output of cat_monthly.sh).
micasa.out.monthly.cat <- function(cfg) {
  sprintf("%s/MiCASA_%s_flux_x360_y180_monthly.nc",
          cfg$monthly.1x1, cfg$version)
}
