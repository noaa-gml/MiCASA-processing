## Shared helpers for ingest_byyear.r and ingest_monthly.r.
##
## These were duplicated byte-for-byte in both scripts before. Source via:
##     source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "ingest_common.r"))
##
## Assumes ct.setup() has been called (provides ncdim_def, ncvar_def, etc.).

## ---- Constants -------------------------------------------------------------

# MiCASA tracers we ingest from the raw 0.1° files.
micasa.tracers <- c("NPP", "Rh", "FIRE", "FUEL")

# Earth mean radius (m), per the MiCASA dataset documentation.
EARTH_RADIUS_M <- 6371007.2

# Output 1° grid axes — same for daily and monthly.
micasa.dim.lon <- function() {
  ncdim_def("longitude", "degrees_east", vals = seq(-179.5, 179.5, 1))
}
micasa.dim.lat <- function() {
  ncdim_def("latitude",  "degrees_north", vals = seq(-89.5, 89.5, 1))
}

# UNIX epoch as POSIX, matching how time is encoded in the 1° outputs.
micasa.epoch <- function() ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC")
micasa.timeunits          <- "seconds"
micasa.timeunits.difftime <- "secs"

## ---- Geometry / aggregation ------------------------------------------------

# Area (m^2) of a single grid cell with corners at `lons` and `lats` (radians).
archimedes <- function(lons, lats) {
  if (length(lons) != 2) stop("Lons vector length not 2")
  if (length(lats) != 2) stop("Lats vector length not 2")
  if (any(abs(range(lons)) >  pi))     stop("abs(lons) vector exceeds pi")
  if (any(abs(range(lats)) > (pi / 2))) stop("abs(lats) vector exceeds pi/2")
  (sin(lats[2]) - sin(lats[1])) * (lons[2] - lons[1]) * EARTH_RADIUS_M^2
}

# Compute the 1800-element latitude-cell-area vector for the 0.1° MiCASA grid.
# `lats` is a vector of cell-center latitudes in degrees (length 1800).
compute.gca <- function(lats) {
  gca <- rep(NA_real_, length(lats))
  lon.rad <- c(-0.05, 0.05) * (pi / 180)  # 0.1° wide cell at the equator
  for (ilat in seq_along(lats)) {
    lat.rad <- (pi / 180) * (lats[ilat] + c(-0.05, 0.05))
    gca[ilat] <- archimedes(lon.rad, lat.rad)
  }
  gca
}

# Aggregate a 3600x1800 0.1° field to 360x180 1° using cell-area weights.
# `gca` is the 1800-element latitude-area vector from compute.gca().
#
# Vectorized 2026-04-26. Decomposition:
#   * The 0.1° grid factors as (10 lon-in × 360 lon-out) × (10 lat-in × 180 lat-out).
#   * Each output cell averages 100 input cells with weight gca[lat] (constant in lon).
#   * The unnormalized sum factors: collapse the 10 lon-in's first (uniform weight),
#     then weight each row by gca, then collapse 10 lat-in's into the 180 lat-out's.
#   * NA handling: build a 0/1 mask, run the same pipeline on the mask, divide.
#     This matches weighted.mean(..., na.rm = TRUE), which renormalizes by the
#     remaining weights; an all-NA output cell becomes NaN.
#
# History: previously a triple-loop in R. Pre-2026-04-26 versions also had a
# numerical bug — see lib/test_aggregate.r regression test — where lat-area
# weights were recycled along the lon axis instead of the lat axis. The
# vectorized form here is ~4.6× faster than the bug-fixed scalar version on
# 3600×1800 fields (random + 1%% NA, single thread, Orion login node), and
# matches it to machine precision (max |err| ~2e-16).
aggregate.to.1x1 <- function(fld, gca) {
  mask        <- !is.na(fld)
  fld_clean   <- fld
  fld_clean[!mask] <- 0

  # 1) Sum over the 10 lon-in cells per lon-block (uniform weight, since
  #    lon spacing is constant): 3600×1800 → 360×1800.
  s_lon <- matrix(colSums(array(fld_clean, dim = c(10, 360, 1800))), 360, 1800)
  m_lon <- matrix(colSums(array(mask + 0,  dim = c(10, 360, 1800))), 360, 1800)

  # 2) Apply lat-area weight along the lat axis. Doing this AFTER step 1
  #    (not before) means we sweep 360×1800 = 0.65M elements instead of
  #    3600×1800 = 6.5M — a ~10× reduction in scalar multiplies.
  s_lon_w <- sweep(s_lon, 2, gca, "*")
  m_lon_w <- sweep(m_lon, 2, gca, "*")

  # 3) Sum over the 10 lat-in cells per lat-block: 360×1800 → 360×180.
  #    Reshape (360, 10, 180) and unroll the 10-element sum (faster than apply).
  arr_n <- array(s_lon_w, dim = c(360, 10, 180))
  arr_d <- array(m_lon_w, dim = c(360, 10, 180))
  num   <- arr_n[, 1, ]; for (k in 2:10) num   <- num   + arr_n[, k, ]
  denom <- arr_d[, 1, ]; for (k in 2:10) denom <- denom + arr_d[, k, ]

  out <- num / denom
  out[denom == 0] <- NaN  # all-NA block → NaN (matches weighted.mean na.rm=TRUE)
  out
}

## ---- Time dimension --------------------------------------------------------

# Build a netCDF unlimited time dimension whose value is `date` (a POSIXct),
# encoded as "seconds since 1970-01-01 00:00:00 UTC".
micasa.time.dim <- function(date) {
  time.vals <- as.numeric(difftime(date, micasa.epoch(),
                                   units = micasa.timeunits.difftime))
  ncdim_def("time",
            sprintf("%s since %s",
                    micasa.timeunits,
                    format(micasa.epoch(), format = "%Y-%m-%d %H:%M:%S UTC")),
            vals = time.vals, unlim = TRUE)
}

## ---- netCDF write ----------------------------------------------------------

# Build the per-tracer ncvar_def list for a daily/monthly output file.
# `ncin` is the loaded raw input (used only for long_name passthrough).
make.tracer.vars <- function(ncin, dim.lon, dim.lat, dim.time) {
  vars <- list()
  for (nm in micasa.tracers) {
    vars[[nm]] <- ncvar_def(name = nm, units = "gC m^-2 s^-1",
                            dim = list(dim.lon, dim.lat, dim.time),
                            missval = -1e34, compression = 9,
                            longname = attributes(ncin[[nm]])$long_name,
                            prec = "float")
  }
  vars
}

# Write `vals` (a list keyed by tracer name) to `ncout` with provenance attrs.
# `srcnm`       — path of the raw input file (recorded in :Source)
# `script.name` — the calling script's Time-stamp (recorded in :history)
write.netcdf <- function(ncout, vars, vals, srcnm, script.name) {
  if (file.exists(ncout)) file.remove(ncout)
  ncf <- nc_create(ncout, vars = vars)
  ncatt_put(ncf, 0, "history",
            attval = sprintf("Created on %s\nby script '%s'",
                             format(Sys.time(), "%a %b %d %Y %H:%M:%S %Z"),
                             script.name),
            prec = "text")
  ncatt_put(ncf, 0, "Source", attval = srcnm, prec = "text")
  for (nm in names(vars)) ncvar_put(ncf, vars[[nm]], vals[[nm]])
  nc_close(ncf)
}
