#!/usr/bin/env Rscript

#SBATCH --account co2
#SBATCH --time 8:00:00
#SBATCH --ntasks 1
#SBATCH --mem 20g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --mail-type=FAIL

## Ingest raw 0.1° MiCASA daily files → aggregated 1° daily files.
##
## Year-parallel: when launched without $INGEST_YEAR set, self-submits one
## SBATCH job per year in [$MICASA_YEAR_START, $MICASA_YEAR_END]. When run
## inside an SBATCH job (with $INGEST_YEAR set), processes that single year.
##
## Tracers ingested: NPP, Rh, FIRE, FUEL (FIRE and FUEL come from the daily
## stream; NPP and Rh from the daily stream are also retained for QC).
##
## Skip-existing is mtime-aware: a day is re-ingested if the source
## .nc4 is newer than the existing 1° output (NASA can republish files).
## Set RECOMPUTE_EXISTING=1 to force re-ingest unconditionally.
##
## Reads only the 4 tracers we need from the raw file (saves ~22%
## per-day read time).

ct.setup()
script.name <- "ingest_byyear.r"

work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
source(file.path(work.dir, "lib", "ingest_common.r"))
cfg <- micasa.config()

year.env           <- Sys.getenv("INGEST_YEAR")
recompute.existing <- nchar(Sys.getenv("RECOMPUTE_EXISTING")) > 0

if (nchar(year.env) == 0) {
  ## Driver mode — fan out one SBATCH per year. Inherit env so the workers
  ## see the same MICASA_YEAR/VERSION/etc. as the driver.
  for (year in cfg$year.start:cfg$year.end) {
    cmd <- sprintf("sbatch -J ingest-%d --mail-user=%s --export=ALL,INGEST_YEAR=%d ingest_byyear.r",
                   year, cfg$mail.user, year)
    cat(cmd, "\n")
    system(cmd)
  }
  quit(save = "no")
}

## ---- Worker mode -----------------------------------------------------------

year <- as.integer(year.env)

setwd(work.dir)
if (!dir.exists(cfg$daily.1x1)) {
  cat(sprintf("Creating output dir \"%s\"\n", cfg$daily.1x1))
  dir.create(cfg$daily.1x1, recursive = TRUE, showWarnings = TRUE)
}

## Pre-compute lat-cell-area vector once (depends only on the input grid).
## Use any existing raw daily file to read the lat coordinate.
probe.year  <- cfg$year.start
probe       <- load.ncdf(micasa.raw.daily(cfg, probe.year, 1, 1))
gca         <- compute.gca(probe$lat)

dim.lon <- micasa.dim.lon()
dim.lat <- micasa.dim.lat()

dpm <- days.in.month(year)
pb  <- progress.bar.start(message = sprintf("%d: %d days", year, sum(dpm)),
                          nx = sum(dpm))

iday <- 0
for (month in 1:12) {
  for (day in 1:dpm[month]) {
    iday <- iday + 1

    this.date <- ISOdatetime(year, month, day, 0, 0, 0, tz = "UTC") + 86400 / 2
    dim.time  <- micasa.time.dim(this.date)

    srcnm <- micasa.raw.daily(cfg, year, month, day)
    ncout <- micasa.out.daily(cfg, year, month, day)

    if (!recompute.existing && out.is.fresh(ncout, srcnm)) {
      cat(sprintf("skipping (fresh) \"%s\"\n", ncout))
      pb <- progress.bar.print(pb, iday)
      next
    }
    if (file.exists(ncout)) {
      cat(sprintf("re-ingesting (source newer or RECOMPUTE_EXISTING=1) \"%s\"\n",
                  ncout))
      file.remove(ncout)
    }

    ## Read only the 4 tracers we aggregate (raw file also has ATMC, NEE).
    ncin <- load.ncdf(srcnm, vars = micasa.tracers, quiet = TRUE)
    vars <- make.tracer.vars(ncin, dim.lon, dim.lat, dim.time)

    vals <- list()
    for (nm in micasa.tracers) {
      ## Source units kg C m-2 s-1 → output gC m-2 s-1.
      vals[[nm]] <- aggregate.to.1x1(ncin[[nm]], gca) * 1e3
    }

    write.netcdf(ncout, vars, vals, srcnm, script.name)
    pb <- progress.bar.print(pb, iday)
  }
}
progress.bar.end(pb)
