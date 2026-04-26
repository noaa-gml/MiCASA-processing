#!/usr/bin/env Rscript
## Compare write time + file size across deflate compression levels
## using a real ingested 1° daily file as input data.
##
## Usage: Rscript lib/bench_compression.r

suppressPackageStartupMessages({
  ct.setup()
})
work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
source(file.path(work.dir, "lib", "ingest_common.r"))

# Pick a representative existing daily file as the source data.
src <- "daily_1x1/MiCASA_v1_flux_x360_y180_daily_20240101.nc"
if (!file.exists(src)) stop(sprintf("Need %s to exist for benchmark", src))

din <- load.ncdf(src)
nlon <- length(din$longitude); nlat <- length(din$latitude)
cat(sprintf("Source: %s  (%d×%d, %d tracers)\n", src, nlon, nlat,
            length(intersect(names(din), micasa.tracers))))

dim.lon <- micasa.dim.lon()
dim.lat <- micasa.dim.lat()
dim.time <- micasa.time.dim(ISOdatetime(2024, 1, 1, 12, 0, 0, tz = "UTC"))

# Replicate make.tracer.vars but with a configurable compression level.
make.vars.lvl <- function(lvl) {
  vars <- list()
  for (nm in micasa.tracers) {
    vars[[nm]] <- ncvar_def(name = nm, units = "gC m^-2 s^-1",
                            dim = list(dim.lon, dim.lat, dim.time),
                            missval = -1e34, compression = lvl,
                            longname = nm, prec = "float")
  }
  vars
}

## Time-resolved write: do N writes per level, drop the first (warmup),
## report mean of the remaining.
bench_one <- function(lvl, n = 10) {
  out_path <- sprintf("/tmp/micasa_bench_lvl%d.nc", lvl)
  times <- numeric(n)
  for (i in seq_len(n)) {
    if (file.exists(out_path)) file.remove(out_path)
    vars <- make.vars.lvl(lvl)
    t0 <- proc.time()[3]
    ncf <- nc_create(out_path, vars = vars)
    for (nm in micasa.tracers) ncvar_put(ncf, vars[[nm]], din[[nm]])
    nc_close(ncf)
    times[i] <- proc.time()[3] - t0
  }
  fsize <- file.info(out_path)$size
  list(lvl = lvl,
       mean_t = mean(times[-1]),
       sd_t   = sd(times[-1]),
       size_b = fsize)
}

cat("\n  level   mean ± sd (s)        bytes        ratio vs lvl9\n")
results <- list()
for (lvl in c(1, 2, 4, 6, 9)) {
  r <- bench_one(lvl)
  results[[as.character(lvl)]] <- r
}
ref9 <- results[["9"]]
for (lvl in names(results)) {
  r <- results[[lvl]]
  cat(sprintf("  %5s   %.4f ± %.4f s   %9d   size %.2fx  time %.2fx\n",
              lvl, r$mean_t, r$sd_t, r$size_b,
              r$size_b / ref9$size_b,
              r$mean_t / ref9$mean_t))
}

# Project the wall-time delta onto a full ingest_byyear year:
# 365 daily writes per year (each writes 4 tracers in one nc_create call,
# which is what bench_one measures).
cat(sprintf("\nProjected per-year write cost (365 days):\n"))
for (lvl in names(results)) {
  r <- results[[lvl]]
  cat(sprintf("  level %s: %.0f s = %.1f min\n",
              lvl, 365 * r$mean_t, 365 * r$mean_t / 60))
}
