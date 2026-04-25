#!/usr/bin/env Rscript

## Verify that every expected daily MiCASA file is present in the local
## NCCS portal mirror for [$MICASA_YEAR_START, $MICASA_YEAR_END].

ct.setup()

work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
cfg <- micasa.config()

setwd(work.dir)

for (year in cfg$year.start:cfg$year.end) {
  dpm <- days.in.month(year)
  pb  <- progress.bar.start(message = sprintf("%d: %d days", year, sum(dpm)),
                            nx = sum(dpm))
  iday <- 0
  for (month in 1:12) {
    for (day in 1:dpm[month]) {
      iday <- iday + 1
      ## Some years are .nc, others .nc4 — accept either.
      src.nc4 <- micasa.raw.daily(cfg, year, month, day)
      src.nc  <- sub("\\.nc4$", ".nc", src.nc4)
      if (!file.exists(src.nc) && !file.exists(src.nc4)) {
        cat(sprintf("no file %s\n", src.nc4))
      }
      pb <- progress.bar.print(pb, iday)
    }
  }
  progress.bar.end(pb)
}
