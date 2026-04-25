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
aggregate.to.1x1 <- function(fld, gca) {
  retval <- matrix(0, 360, 180)
  for (jlat in 1:180) {
    inlats <- 1:10 + 10 * (jlat - 1)
    for (ilon in 1:360) {
      inlons <- 1:10 + 10 * (ilon - 1)
      retval[ilon, jlat] <- 0
      for (inlon in inlons) {
        retval[ilon, jlat] <- retval[ilon, jlat] +
          weighted.mean(fld[inlons, inlats], weights = gca[inlats], na.rm = TRUE)
      }
      retval[ilon, jlat] <- retval[ilon, jlat] / 10
    }
  }
  retval
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
