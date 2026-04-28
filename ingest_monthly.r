#!/usr/bin/env Rscript

#SBATCH --account co2
#SBATCH --time 8:00:00
#SBATCH --ntasks 1
#SBATCH --mem 20g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --mail-type=FAIL

## Ingest raw 0.1° MiCASA monthly files → aggregated 1° monthly files.
##
## Loops over [$MICASA_YEAR_START, $MICASA_YEAR_END]. Skip-existing is
## mtime-aware: a month is re-ingested if the source .nc4 is newer than
## the existing 1° output. Set RECOMPUTE_EXISTING=1 to force re-ingest.

ct.setup()
script.name <- "ingest_monthly.r"

work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
source(file.path(work.dir, "lib", "ingest_common.r"))
cfg <- micasa.config()

recompute.existing <- nchar(Sys.getenv("RECOMPUTE_EXISTING")) > 0

setwd(work.dir)
if (!dir.exists(cfg$monthly.1x1)) {
  cat(sprintf("Creating output dir \"%s\"\n", cfg$monthly.1x1))
  dir.create(cfg$monthly.1x1, recursive = TRUE, showWarnings = TRUE)
}

## Pre-compute lat-cell-area vector once.
probe.year <- cfg$year.start
probe      <- load.ncdf(micasa.raw.monthly(cfg, probe.year, 1))
gca        <- compute.gca(probe$lat)

dim.lon <- micasa.dim.lon()
dim.lat <- micasa.dim.lat()

for (year in cfg$year.start:cfg$year.end) {
  for (month in 1:12) {

    this.date <- seq.midmon(year, month, year, month)
    dim.time  <- micasa.time.dim(this.date)

    srcnm <- micasa.raw.monthly(cfg, year, month)
    ncout <- micasa.out.monthly(cfg, year, month)

    cat(sprintf("Processing %s...", basename(srcnm)))

    if (!file.exists(srcnm)) {
      cat(sprintf("source missing, skipping: %s\n", srcnm))
      next
    }
    if (!recompute.existing && out.is.fresh(ncout, srcnm)) {
      cat(sprintf("skipping (fresh) \"%s\"\n", ncout))
      next
    }
    if (file.exists(ncout)) {
      cat(sprintf("re-ingesting (source newer or RECOMPUTE_EXISTING=1) \"%s\"\n",
                  ncout))
      file.remove(ncout)
    }

    ncin <- load.ncdf(srcnm, vars = micasa.tracers, quiet = TRUE)
    vars <- make.tracer.vars(ncin, dim.lon, dim.lat, dim.time)

    vals <- list()
    for (nm in micasa.tracers) {
      ## Source units kg C m-2 s-1 → output gC m-2 s-1.
      vals[[nm]] <- aggregate.to.1x1(ncin[[nm]], gca) * 1e3
    }

    write.netcdf(ncout, vars, vals, srcnm, script.name)
    cat("\n")
  }
}
