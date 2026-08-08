## Shared helpers for ingest_byyear.r and ingest_monthly.r.
##
## These were duplicated byte-for-byte in both scripts before. Source via:
##     source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "ingest_common.r"))
##
## Assumes ct.setup() has been called (provides ncdim_def, ncvar_def, etc.).

## ---- Constants -------------------------------------------------------------

# MiCASA fields we ingest from the raw 0.1° files, in TWO classes.
#
# micasa.tracers -- the carbon flux components. Physical fluxes; the
# consumer combines them (NEE = Rh - NPP, with fire/fuel added separately).
#
# micasa.diagnostics -- CARRIED FOR REFERENCE, NEVER ADDED TO A FLUX.
# ATMC ("Atmospheric correction") is the LoFI empirical sink of Weir et
# al. 2021a (ACP, doi:10.5194/acp-21-9609-2021), rescaled ANNUALLY so the
# global biospheric total matches the observed atmospheric CO2 growth
# rate. MiCASA's own NEE is (Rh - NPP - ATMC); ours is deliberately
# (Rh - NPP) -- ATMC is NOT subtracted -- because these fluxes are the
# PRIOR to a global inversion that itself assimilates atmospheric CO2.
# Pre-correcting would smuggle observational information from the same
# data class into the prior (data leakage), after which the inversion can
# no longer independently constrain the long-term sink. Subtracting it
# was tried 2026-04-29 (proposal #2) and reverted the same day; full
# reasoning in docs/V1_TO_V2_JUSTIFICATION.md section 5.2.
#
# It IS ingested and shipped as its own variable so the choice is
# REVERSIBLE by the consumer without re-reading the 0.1° archive: a
# forward/diagnostic user outside an inversion (no double-dipping) can
# form MiCASA's NEE as (Rh - NPP - ATMC) themselves. Nothing in this
# pipeline ever adds or subtracts it. Same units as the tracers
# (kg m-2 s-1 carbon upstream -> gC m^-2 s^-1 here), so it rides the
# generic aggregation and x1e3 conversion unchanged.
micasa.tracers     <- c("NPP", "Rh", "FIRE", "FUEL")
micasa.diagnostics <- c("ATMC")
micasa.ingest.vars <- c(micasa.tracers, micasa.diagnostics)

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

## ---- Skip-existing freshness check -----------------------------------------

# Make-style staleness check: returns TRUE iff `ncout` exists AND was last
# modified after `srcnm`. Use to gate skip-existing in ingest hot loops:
#
#     if (!recompute.existing && out.is.fresh(ncout, srcnm)) next
#
# Why mtime, not just file.exists: NASA republishes source files (especially
# vNRT). A pure file.exists skip would silently keep stale aggregates.
#
# *** WARNING -- this check is NOT sufficient after a MICASA_REFRESH=1 run. ***
# Plain wget (--no-clobber) sets local mtime to the DOWNLOAD time, so a
# re-download did make mtime(srcnm) > mtime(ncout). But `wget --timestamping`
# (the MICASA_REFRESH path added after the 2025-06-17 republication incident)
# sets mtime to UPSTREAM's Last-Modified, which is typically OLDER than the
# 1° output it must replace -- so this returns "fresh" and the re-ingest is
# silently skipped. Measured 2026-08-07: refreshed raw 2025-06-17 vs output
# 2026-06-23. ALWAYS pass RECOMPUTE_EXISTING=1 when re-ingesting after a
# refresh; see CHANGELOG 2026-08-07.
out.is.fresh <- function(ncout, srcnm) {
  file.exists(ncout) && file.mtime(ncout) > file.mtime(srcnm)
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
  for (nm in micasa.ingest.vars) {
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
  # Diagnostics carry an explicit "do not add this" note, so the distinction
  # survives into the file and does not live only in this repo's comments.
  for (nm in intersect(micasa.diagnostics, names(vars))) {
    ncatt_put(ncf, nm, "usage", prec = "text", attval = paste(
      "DIAGNOSTIC -- carried for reference only. NOT subtracted from, or added",
      "to, any flux in this product. MiCASA's own published NEE is",
      "(Rh - NPP - ATMC); this product's NEE is (Rh - NPP), because these fluxes",
      "are the prior to an inversion that assimilates the same atmospheric CO2",
      "the ATMC correction was tuned to. Subtract it yourself ONLY for",
      "forward/diagnostic use outside an inversion. See",
      "docs/V1_TO_V2_JUSTIFICATION.md 5.2."))
  }
  nc_close(ncf)
}
