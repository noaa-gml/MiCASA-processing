#!/usr/bin/env Rscript
## bench_compression_diurnal.r
##
## Benchmark netCDF compression levels on a payload representative of
## diurnalize-ERA5.r output by reading a real production file
## and rewriting at each level.
##
## Reading real flux/met fields is critical -- synthetic Gaussian noise
## fills the full dynamic range and defeats compression at every level.

library(ncdf4)

## Pass the path to the input fluxes file as the first argument; e.g.:
##     Rscript lib/bench_compression_diurnal.r ERA5/fluxes_202401.nc
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript lib/bench_compression_diurnal.r <fluxes_YYYYMM.nc>")
}
src <- args[1]
stopifnot(file.exists(src))
cat(sprintf("source: %s  (%.1f MB)\n", src, file.info(src)$size / 1024^2))

bench_dir <- tempfile("bench_diurnal_real_")
dir.create(bench_dir)
cat(sprintf("bench dir: %s\n", bench_dir))

cat("reading source ... ")
t0 <- proc.time()[3]
ncf <- nc_open(src)

dim_names <- names(ncf$dim)
lon_name <- intersect(c("lon", "longitude"), dim_names)[1]
lat_name <- intersect(c("lat", "latitude"),  dim_names)[1]
stopifnot(!is.na(lon_name), !is.na(lat_name))

lon <- ncvar_get(ncf, lon_name)
lat <- ncvar_get(ncf, lat_name)
tim <- ncvar_get(ncf, "time")
timeunits <- ncatt_get(ncf, "time", "units")$value
# variables = anything in $var that is not a dim coord
varnames <- setdiff(names(ncf$var), c(lon_name, lat_name, "time"))
data <- list()
for (nm in varnames) data[[nm]] <- ncvar_get(ncf, nm)
nc_close(ncf)
cat(sprintf("%.1fs   nlon=%d nlat=%d nt=%d   vars=%s\n",
            proc.time()[3] - t0, length(lon), length(lat), length(tim),
            paste(varnames, collapse=",")))

write_one <- function(comp) {
  fn <- sprintf("%s/comp%d.nc", bench_dir, comp)
  if (file.exists(fn)) file.remove(fn)

  lon.dim  <- ncdim_def("lon",  "degrees_east",  vals = lon)
  lat.dim  <- ncdim_def("lat",  "degrees_north", vals = lat)
  date.dim <- ncdim_def("time", timeunits, vals = tim, unlim = TRUE)

  vars <- list()
  for (nm in varnames) {
    a <- data[[nm]]
    dims <- if (length(dim(a)) <= 1) list(date.dim) else
            list(lon.dim, lat.dim, date.dim)
    vars[[nm]] <- ncvar_def(name = nm, units = "x", dim = dims,
                            missval = -1e34, compression = comp, longname = nm)
  }

  t0 <- proc.time()[3]
  ncf <- nc_create(fn, vars = vars)
  for (nm in varnames) ncvar_put(ncf, vars[[nm]], vals = data[[nm]])
  nc_close(ncf)
  elapsed <- proc.time()[3] - t0
  size_mb <- file.info(fn)$size / 1024^2
  list(comp = comp, time_s = elapsed, size_mb = size_mb)
}

# warm cache
invisible(write_one(1))

results <- list()
for (comp in c(1, 2, 3, 4, 6, 9)) {
  r <- write_one(comp)
  cat(sprintf("  level %d: %6.2fs   %7.1f MB\n", r$comp, r$time_s, r$size_mb))
  results[[as.character(comp)]] <- r
}

cat("\n=== Per-month (this file) ===\n")
ref <- results[["9"]]
orig_mb <- file.info(src)$size / 1024^2
cat(sprintf("%-7s %8s %8s %10s %10s %10s\n",
            "level", "time_s", "size_MB", "speedup", "size_x", "vs_orig"))
for (r in results) {
  cat(sprintf("%-7d %8.2f %8.1f %9.2fx %9.2fx %9.2fx\n",
              r$comp, r$time_s, r$size_mb,
              ref$time_s / r$time_s, r$size_mb / ref$size_mb,
              r$size_mb / orig_mb))
}

cat("\n=== Per-year (12 months) projection ===\n")
for (r in results) {
  cat(sprintf("level %d  write-only: %5.0fs   total: %5.1f GB\n",
              r$comp, 12 * r$time_s, 12 * r$size_mb / 1024))
}

cat(sprintf("\nbench files left in: %s\n", bench_dir))
