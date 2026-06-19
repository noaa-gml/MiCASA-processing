#!/usr/bin/env Rscript

#SBATCH --account co2
#SBATCH --time 8:00:00
#SBATCH --ntasks 1
#SBATCH --mem 40g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --mail-type=FAIL

## Diurnalize monthly NPP/Rh into hourly GPP/RESP/NEE on ERA5 meteo, writing
## ERA5/fluxes_YYYYMM.nc per (year, month) in the configured ranges.
##
## Year-parallel: when launched without $diurn_year, fans out one SBATCH job
## per year in [$MICASA_YEAR_START, $MICASA_YEAR_END]. Months processed are
## [$MICASA_MONTH_START, $MICASA_MONTH_END] (default 1..12).
##
## Climatology fallback: per (year, month), the monthly NPP/Rh file is
## used if present; otherwise the Rh/NPP day-of-year climatology
## (NPPclim.nc / Rhclim.nc) stands in. The choice is made per month by
## file existence -- no year list to maintain -- so a partially-published
## year is handled automatically (real early months, climatology for the
## rest).

script.name <- "diurnalize-ERA5.r"

work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
cfg <- micasa.config()

product.name <- sprintf("MiCASA_%s_flux_x360_y180_monthly", cfg$version)

## Hourly 1° ERA5 from the TM5 meteo tree. The primary tree is tried
## first; the FastTrack tree (ea_0005) is a fallback for the NRT
## trailing window where the primary has not yet been populated. A day
## is resolved wholly to the first tree that holds all of varnms for it,
## so provenance stays clean; the per-day attribution is written to the
## output (meteo_source_* global attributes). Either root can be
## overridden via MICASA_ERA5_DIR / MICASA_ERA5_DIR_FALLBACK.
.ct <- Sys.getenv("CARBONTRACKER", "")
era5dirs <- c(
  primary   = sprintf("%s/METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100",     .ct),
  fasttrack = sprintf("%s/METEO/tm5-nc/ec/ea_0005/h06h18tr1/sfc/glb100x100", .ct))
if (nchar(Sys.getenv("MICASA_ERA5_DIR")) > 0)
  era5dirs["primary"]   <- Sys.getenv("MICASA_ERA5_DIR")
if (nchar(Sys.getenv("MICASA_ERA5_DIR_FALLBACK")) > 0)
  era5dirs["fasttrack"] <- Sys.getenv("MICASA_ERA5_DIR_FALLBACK")
era5template <- "YYYY/MM/VVV_YYYYMMDD_00p01.nc"
cat(sprintf("ERA5 meteo search order:\n%s\n",
            paste(sprintf("  %-9s %s", names(era5dirs), era5dirs),
                  collapse = "\n")))

## ERA5 meteo path-resolution helpers (era5.relpath, resolve.era5.source,
## encode.day.runs) live in lib/era5_meteo.r so they can be unit-tested
## standalone (tests/test_era5_meteo.r).
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "era5_meteo.r"))

## Provenance helpers (provenance.attrs, nc.write.provenance): the CF/ACDD
## global-attribute set -- producing software, git commit, input checksums,
## host, timestamp -- stamped onto every output. See lib/provenance.r.
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "provenance.r"))

## Run-manifest helper (manifest.record): the worker records start / ok /
## fail to jobs/run_manifest.tsv, which verify_v2 reads. See lib/manifest.r.
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "manifest.r"))

## Diurnalize flux transform (diurnal.flux, polar.night.clip): the pure
## arithmetic core, extracted so its invariants can be unit-tested
## (tests/test_diurnal.r). See lib/diurnal.r.
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "diurnal.r"))

yr.env <- Sys.getenv("diurn_year")

if (nchar(yr.env) == 0) {
  ## Driver mode — fan out per year. Inherit env so workers see same config.
  for (yr in cfg$year.start:cfg$year.end) {
    cmd <- sprintf("sbatch -J d-%d-MiCASA --mail-user=%s --export=ALL,diurn_year=%d diurnalize-ERA5.r",
                   yr, cfg$mail.user, yr)
    cat(cmd, "\n")
    system(cmd)
  }
  quit(save = "no")
}

## ---- Worker mode -----------------------------------------------------------

yr <- as.integer(yr.env)

## Run manifest: record this worker's start, and arm an error handler that
## records `fail` on any uncaught error. The matching `ok` (and disarm) is
## written after the month loop.
diurn.t0 <- Sys.time()
manifest.record("diurnalize-ERA5.r", "start", detail = sprintf("year=%d", yr))
options(error = function() {
  manifest.record("diurnalize-ERA5.r", "fail",
                  elapsed = as.integer(difftime(Sys.time(), diurn.t0,
                                                units = "secs")),
                  detail  = sprintf("year=%d uncaught error", yr))
})

setwd(work.dir)
in.dir  <- cfg$monthly.1x1
out.dir <- cfg$era5.dir
## A/B test mode: shadow output dir + single-month restriction. Default
## (env vars unset): no-op, full production loop. Set
## MICASA_DIURN_OUT_DIR=ERA5_atmc_test/ and MICASA_DIURN_ONLY_MONTH=7 to
## diurnalize one month into a shadow dir for diff-against-canonical
## comparisons (used 2026-04 to validate the ATMC integration before
## fanning out the full re-diurnalize).
.diurn.out.override <- Sys.getenv("MICASA_DIURN_OUT_DIR", unset = "")
if (nchar(.diurn.out.override) > 0) {
  out.dir <- .diurn.out.override
  cat(sprintf("** Test mode: out.dir overridden to %s **\n", out.dir))
}
.diurn.only.mon <- as.integer(Sys.getenv("MICASA_DIURN_ONLY_MONTH", unset = "0"))

ct.setup()

## Load PIQS coefficients to smooth month-month variability. The fit file
## defaults to fit.piqs.rda; MICASA_FIT_RDA overrides it (e.g. fit.linmm.rda)
## to A/B an alternative fitter without clobbering the production fit.
fit.rda.path <- file.path(work.dir, Sys.getenv("MICASA_FIT_RDA", "fit.piqs.rda"))
load(fit.rda.path)
piqsfit.time <- epoch.seconds.to.POSIX(piqsfit.time)
piqsfit.lts  <- as.POSIXlt(piqsfit.time)

## Optional prior-uncertainty output: when the fit carries per-piece kriging
## variance ($var, written by write_atpk.r) AND MICASA_WRITE_FLUX_SD is set,
## diurnalize emits an NEE_sd field (1-sigma prior uncertainty on the smoothed
## sub-monthly NEE flux, constant within a month). Default off; no effect on the
## PCHIP/PPM/PIQS path (those fits have no $var).
write.flux.sd <- (!is.null(piqsfit.gpp$var)) && nchar(Sys.getenv("MICASA_WRITE_FLUX_SD")) > 0
if (write.flux.sd) cat("MICASA_WRITE_FLUX_SD: emitting NEE_sd prior-uncertainty (ATP kriging variance)\n")

## Provenance of the coefficient fit: path + SHA-256 (hashed once here, not
## once per month -- fit.piqs.rda is ~190 MB) and the fitter that wrote it.
fit.rda.sha256 <- prov.file.sha256(fit.rda.path)
fit.method     <- if (exists("piqsfit.meta") && !is.null(piqsfit.meta$fitter))
                     piqsfit.meta$fitter else "unknown"

## ---- Fit-window banner (proposal #2 in README.ash) ------------------------
## Always report the loaded fit's coverage and any padding metadata, and warn
## if the active diurnalization year extends past the right edge -- the
## climatology-fallback branch below would silently substitute a smoother
## coefficient field for those months. Set MICASA_STRICT_PIQS=1 to escalate
## this from a warning to a hard error (recommended for NRT operations).
cat(sprintf("PIQS fit window: %s -- %s (%d months)\n",
            format(min(piqsfit.time), "%Y-%m"),
            format(max(piqsfit.time), "%Y-%m"),
            length(piqsfit.time)))
if (exists("piqsfit.meta")) {
  cat(sprintf("PIQS fit padding: left=%d, right=%d (written %s)\n",
              piqsfit.meta$pad.left, piqsfit.meta$pad.right,
              piqsfit.meta$written.at))
} else {
  cat("PIQS fit padding: unknown (no piqsfit.meta in .rda; assume 0/0)\n")
}
cat(sprintf("Active diurnalization year: %d\n", yr))

strict.piqs <- as.integer(Sys.getenv("MICASA_STRICT_PIQS", unset = "0")) == 1
## The fit-edge check is per month, in the loop below: a month with real
## monthly data but past the fit window gets climatological sub-monthly
## coefficients (MICASA_STRICT_PIQS=1 escalates that to a hard error).
## --------------------------------------------------------------------------

lon.dim <- ncdim_def("longitude", "degrees_east", vals = seq(-179.5, 179.5, 1))
lat.dim <- ncdim_def("latitude",  "degrees_north", vals = seq(-89.5, 89.5, 1))

epoch <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC")
timeunits <- "days"

dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

mon.range <- cfg$month.start:cfg$month.end
if (.diurn.only.mon > 0) {
  mon.range <- .diurn.only.mon
  cat(sprintf("** Test mode: restricted to month %d **\n", .diurn.only.mon))
}
n.written <- 0L
for (mon in mon.range) {

  cat(sprintf("%d/%02d\n", yr, mon))
  monstr <- sprintf("%d%02d", yr, mon)
  current.time <- ISOdatetime(yr, mon, 1, 0, 0, 0, tz = "UTC")
  ncname.out <- sprintf("%s/fluxes_%s.nc", out.dir, monstr)

  ## Per-month clim decision: use the real monthly file when it exists,
  ## else fall back to the NPP/Rh day-of-year climatology. Decided by file
  ## existence, so a partially-published year needs no special handling.
  real.monthly <- sprintf("%s/%s_%s.nc", in.dir, product.name, monstr)
  use.clim <- !file.exists(real.monthly)
  in.fit   <- current.time >= min(piqsfit.time) &&
              current.time <= max(piqsfit.time)
  if (!use.clim && !in.fit) {
    msg <- sprintf("%d-%02d has real monthly data but lies past the PIQS fit window [%s..%s]; its sub-monthly shape will come from coefficient-climatology. Re-run write_pchip.r with the latest monthly data.",
                   yr, mon, format(min(piqsfit.time), "%Y-%m"),
                   format(max(piqsfit.time), "%Y-%m"))
    if (strict.piqs) stop(msg, " Or unset MICASA_STRICT_PIQS to allow.")
    else warning(msg, " Set MICASA_STRICT_PIQS=1 to make this an error.",
                 immediate. = TRUE)
  }

  if (use.clim) {
    fname <- sprintf("%s/NPPclim.nc", in.dir)
    if (!file.exists(fname)) stop(sprintf("%d-%02d:  %s does not exist.", yr, mon, fname))
    foo <- load.ncdf(fname)
    ## Sign flip: negative = sink. Units: gC m-2 s-1 → mol m-2 s-1.
    gpp.clim <- -2 * foo$NPPCLIM / 12

    fname <- sprintf("%s/Rhclim.nc", in.dir)
    if (!file.exists(fname)) stop(sprintf("%d-%02d:  %s does not exist.", yr, mon, fname))
    foo <- load.ncdf(fname)
    rh.clim   <- foo$RHCLIM / 12
    rtot.clim <- rh.clim - 0.5 * gpp.clim   # autotrophic = NPP = -GPP/2

    rtot.mn <- rtot.clim[, , mon]
    gpp.mn  <- gpp.clim[, , mon]
    rm(foo)
  } else {
    fname <- sprintf("%s/%s_%s.nc", in.dir, product.name, monstr)
    foo <- load.ncdf(fname)
    gpp.mn  <- -2 * foo$NPP / 12
    rh.mn   <- foo$Rh / 12
    rtot.mn <- rh.mn - 0.5 * gpp.mn
    rm(foo)
  }
  cat(sprintf("Finished reading %s...\n", fname))

  ## ---- Read ERA5 meteo for this month ----
  ##
  ## ERA5 sometimes has gaps near year boundaries while ECMWF processes the
  ## final day of the year (e.g. 2025-12-31 was missing on Orion as of
  ## 2026-04-28). Detect available days first, drop missing days from the
  ## hourly time axis, and mark the output as provisional. Older behaviour
  ## (crash on first nc_open) is preserved when every day is present.
  varnms <- c("t2m", "ssrd", "stl1", "swvl1")
  dpm.full <- days.in.month(yr)[mon]
  ## Per-day meteo source: "primary", "fasttrack", or NA (day missing
  ## from every tree). day.source is indexed by day number.
  day.source <- vapply(1:dpm.full,
                       function(d) resolve.era5.source(era5dirs, era5template, yr, mon, d, varnms),
                       character(1))
  available.days <- which(!is.na(day.source))
  missing.days   <- which(is.na(day.source))
  if (length(available.days) == 0) {
    cat(sprintf("WARN: %d/%02d has no complete ERA5 days; skipping month.\n", yr, mon))
    next
  }
  partial.month <- length(missing.days) > 0
  if (partial.month) {
    cat(sprintf("WARN: %d/%02d: only %d/%d days have complete ERA5 meteo (missing day(s): %s); writing partial month.\n",
                yr, mon, length(available.days), dpm.full,
                paste(missing.days, collapse=",")))
  }
  cat(sprintf("ERA5 meteo source: %s\n",
              encode.day.runs(available.days, day.source)))
  dpm <- length(available.days)
  mets <- list()
  times <- rep(NA, dpm * 24)
  for (k in seq_along(available.days)) {
    day <- available.days[k]
    src <- day.source[day]
    for (varnm in varnms) {
      ncname.in <- file.path(unname(era5dirs[[src]]),
                             era5.relpath(era5template, yr, mon, day, varnm))
      foo <- load.ncdf(ncname.in)
      if (is.null(mets[[varnm]])) {
        mets[[varnm]] <- array(NA, dim = c(360, 180, 24 * dpm))
      }
      k0 <- 1 + (k - 1) * 24
      k1 <- 23 + k0
      mets[[varnm]][, , k0:k1] <- foo[[varnm]]
      if (varnm == varnms[1]) {
        ## ERA5 time axes were faulty 2022-04-30..2024-06-29; rebuild explicitly.
        times[k0:k1] <- seq(ISOdatetime(yr, mon, day,  0, 0, 0),
                            ISOdatetime(yr, mon, day, 23, 0, 0), by = "1 hour")
      }
    }
  }
  times <- epoch.seconds.to.POSIX(times) + 1800   # shift to mid-hour
  nslots <- length(times)

  q10    <- 1.5 ^ ((mets$t2m - 273.15) / 10.0)
  q10.mn <- apply(q10,         c(1, 2), mean)
  ssr.mn <- apply(mets$ssrd,   c(1, 2), mean)

  ## ERA5 SSRD is identically 0 above ~70°N — guard against div-by-zero.
  lx <- which(ssr.mn == 0)
  if (length(lx) > 0) ssr.mn[lx] <- 1e-16

  cat("GPP mean summary:\n");  print(summary(as.vector(gpp.mn)))
  cat("RTOT mean summary:\n"); print(summary(as.vector(rtot.mn)))
  factor.gpp  <- gpp.mn  / ssr.mn
  factor.resp <- rtot.mn / q10.mn
  lx <- which(is.na(factor.gpp))
  if (length(lx) > 0) factor.gpp[lx] <- 0

  gpp   <- array(NA, dim = dim(mets$ssrd))
  resp  <- array(NA, dim = dim(mets$ssrd))
  nee   <- array(NA, dim = dim(mets$ssrd))
  qgpp  <- array(NA, dim = dim(mets$ssrd))
  qresp <- array(NA, dim = dim(mets$ssrd))

  ## Subtract monthly mean and insert smoothed PIQS fit.
  if ((current.time >= min(piqsfit.time)) & (current.time <= max(piqsfit.time))) {
    imon <- which(piqsfit.time == current.time)
    gpp.a  <- piqsfit.gpp$a[, , imon];  gpp.b  <- piqsfit.gpp$b[, , imon];  gpp.c  <- piqsfit.gpp$c[, , imon]
    resp.a <- piqsfit.resp$a[, , imon]; resp.b <- piqsfit.resp$b[, , imon]; resp.c <- piqsfit.resp$c[, , imon]
    nee.sd.mn <- if (write.flux.sd) sqrt(piqsfit.gpp$var[, , imon] + piqsfit.resp$var[, , imon]) / 12 else NULL
  } else {
    monseq <- which((piqsfit.lts$mon + 1) == mon)
    gpp.a  <- apply(piqsfit.gpp$a[, , monseq],  c(1, 2), mean)
    gpp.b  <- apply(piqsfit.gpp$b[, , monseq],  c(1, 2), mean)
    gpp.c  <- apply(piqsfit.gpp$c[, , monseq],  c(1, 2), mean)
    resp.a <- apply(piqsfit.resp$a[, , monseq], c(1, 2), mean)
    resp.b <- apply(piqsfit.resp$b[, , monseq], c(1, 2), mean)
    resp.c <- apply(piqsfit.resp$c[, , monseq], c(1, 2), mean)
    nee.sd.mn <- if (write.flux.sd)
      sqrt(apply(piqsfit.gpp$var[, , monseq], c(1,2), mean) +
           apply(piqsfit.resp$var[, , monseq], c(1,2), mean)) / 12 else NULL
  }
  ## Sign-flip diagnostics (proposal #4 in README.ash). GPP is negative-for-
  ## uptake here (gpp = -2*NPP), so a positive qmod.gpp at sub-monthly
  ## resolution is unphysical and indicates the PIQS quadratic overshot above
  ## zero. Respiration is positive-typical, so the symmetric concern is
  ## qmod.resp going negative.
  land.gpp  <- is.finite(gpp.mn)  & gpp.mn  < 0
  land.resp <- is.finite(rtot.mn) & rtot.mn > 0
  any.flip.gpp  <- array(FALSE, dim = c(360, 180))
  any.flip.resp <- array(FALSE, dim = c(360, 180))
  flip.hours.gpp  <- 0L
  flip.hours.resp <- 0L

  for (islot in 1:nslots) {
    dt <- as.numeric(times[islot]) - as.numeric(times[1])
    qmod.gpp  <- (gpp.a  * dt^2 + gpp.b  * dt + gpp.c)  / 12
    qmod.resp <- (resp.a * dt^2 + resp.b * dt + resp.c) / 12

    flip.gpp.now  <- land.gpp  & is.finite(qmod.gpp)  & qmod.gpp  > 0
    flip.resp.now <- land.resp & is.finite(qmod.resp) & qmod.resp < 0
    any.flip.gpp    <- any.flip.gpp  | flip.gpp.now
    any.flip.resp   <- any.flip.resp | flip.resp.now
    flip.hours.gpp  <- flip.hours.gpp  + sum(flip.gpp.now)
    flip.hours.resp <- flip.hours.resp + sum(flip.resp.now)

    ## Hourly flux: the meteo driver scales the monthly mean, and the flat
    ## monthly mean is swapped for the fitted sub-monthly shape (qmod).
    ## See lib/diurnal.r :: diurnal.flux.
    gpp[ , , islot] <- diurnal.flux(mets$ssrd[, , islot], gpp.mn,  ssr.mn, qmod.gpp)
    resp[, , islot] <- diurnal.flux(q10[      , , islot], rtot.mn, q10.mn, qmod.resp)

    ## Polar-night clip (Check 12.2). No incoming shortwave => no
    ## photosynthesis. Without this, the sub-monthly quadratic leaks a
    ## small residual into hours where ssrd is identically 0 (~2.6% of
    ## cells in fluxes_202512.nc, max |GPP|=9.4e-9 mol m-2 s-1). resp /
    ## qgpp / qresp are unaffected (Rh has no reason to vanish in dark).
    gpp[ , , islot] <- polar.night.clip(gpp[, , islot], mets$ssrd[, , islot])
    nee[ , , islot] <- gpp[, , islot] + resp[, , islot]
    qgpp[ , , islot] <- qmod.gpp
    qresp[, , islot] <- qmod.resp
  }

  ## One-line per-month diagnostic. Denominators are land cells with the
  ## expected sign on the monthly mean.
  n.land.gpp  <- sum(land.gpp)
  n.land.resp <- sum(land.resp)
  cat(sprintf("PIQS sign-flip [GPP > 0]:  %d / %d cells (%.2f%%), %d / %d cell-hours (%.3f%%)\n",
              sum(any.flip.gpp), n.land.gpp,
              100 * sum(any.flip.gpp) / max(n.land.gpp, 1L),
              flip.hours.gpp, n.land.gpp * nslots,
              100 * flip.hours.gpp / max(n.land.gpp * nslots, 1L)))
  cat(sprintf("PIQS sign-flip [resp < 0]: %d / %d cells (%.2f%%), %d / %d cell-hours (%.3f%%)\n",
              sum(any.flip.resp), n.land.resp,
              100 * sum(any.flip.resp) / max(n.land.resp, 1L),
              flip.hours.resp, n.land.resp * nslots,
              100 * flip.hours.resp / max(n.land.resp * nslots, 1L)))

  lx <- which(is.na(nee))
  if (length(lx) > 0) nee[lx] <- 0

  cat(sprintf("  %s\n", ncname.out))
  date.vals    <- as.numeric(difftime(times, epoch, units = timeunits))
  decimal.date <- POSIX.to.decimal(times)
  date.dim <- ncdim_def("time",
                        sprintf("%s since %s", timeunits,
                                format(epoch, format = "%Y-%m-%d %H:%M:%S UTC")),
                        vals = date.vals, unlim = TRUE)

  ## NOTE: deflate level chosen via lib/bench_compression_diurnal.r on a real
  ## fluxes_YYYYMM.nc. Level 4 vs level 9 saves ~40%% wall-clock on writes
  ## (108s -> 65s per file, ~9 min/year) at +0.3%% file size. Lower levels
  ## (1-3) save another ~10s/file but cost ~+2%% file size.
  ncvar <- function(name, units, longname, dims = list(lon.dim, lat.dim, date.dim)) {
    ncvar_def(name = name, units = units, dim = dims,
              missval = -1e34, compression = 4, longname = longname)
  }

  vars <- list()
  vars$dd    <- ncvar("decimal_date", "years",        "decimal_date", dims = list(date.dim))
  vars$gpp   <- ncvar("GPP",   "mol m-2 s-1", "gross_primary_production, twice the modeled NPP, positive is source to atm (contrary to conventional definition)")
  vars$resp  <- ncvar("resp",  "mol m-2 s-1", "ecosystem_respiration, as sum of Rhetero and Rauto, positive is source to atm")
  vars$nee   <- ncvar("NEE",   "mol m-2 s-1", "NEE=GPP+RESP, positive is source to atm, as is each component")
  vars$qgpp  <- ncvar("QGPP",  "mol m-2 s-1", "gross_primary_production model")
  vars$qresp <- ncvar("qresp", "mol m-2 s-1", "ecosystem_respiration model")
  vars$ssr   <- ncvar("ssr",   "W/m2",        "ERA5 surface shortwave radiation downward")
  vars$t2m   <- ncvar("t2m",   "K",           "ERA5 2-meter air temperature")
  vars$stl1  <- ncvar("stl1",  "K",           "ERA5 soil level 1 temperature (0-7 cm)")
  vars$swvl1 <- ncvar("swvl1", "m3/m3",       "ERA5 soil level 1 volumetric moisture content (0-7 cm)")
  if (write.flux.sd)
    vars$nee_sd <- ncvar("NEE_sd", "mol m-2 s-1", "prior 1-sigma uncertainty on the smoothed sub-monthly NEE flux (ATP kriging variance, sqrt(var_GPP+var_RESP); constant within a month)")

  if (file.exists(ncname.out)) {
    cat(sprintf("Removing existing output file \"%s\"\n", ncname.out))
    file.remove(ncname.out)
  }
  ncf <- nc_create(ncname.out, vars = vars)

  ## CF/ACDD provenance: producing software (git commit/version), input
  ## files (path + SHA-256), the fitter, host and timestamp -- and the CF
  ## `history` line. The fit file goes through `extra` so its ~190 MB are
  ## hashed once (before the loop), not once per month.
  flux.inputs <- if (use.clim)
    list(npp_climatology = sprintf("%s/NPPclim.nc", in.dir),
         rh_climatology  = sprintf("%s/Rhclim.nc",  in.dir))
  else
    list(monthly_flux = real.monthly)
  nc.write.provenance(
    ncf, step = script.name, work.dir = work.dir,
    title = sprintf("MiCASA %s hourly land carbon flux (GPP, ecosystem respiration, NEE), %d-%02d",
                    cfg$version, yr, mon),
    summary = paste("Hourly gross primary production, ecosystem respiration",
                    "and net ecosystem exchange on a global 1-degree grid,",
                    "from monthly MiCASA NPP/Rh diurnalized with ERA5 hourly",
                    "meteorology. Sign convention: positive = flux to atmosphere."),
    inputs = flux.inputs,
    extra = list(micasa_version         = cfg$version,
                 flux_fit_method        = fit.method,
                 flux_from_climatology  = if (use.clim) "yes" else "no",
                 input_pchip_fit        = fit.rda.path,
                 input_pchip_fit_sha256 = if (is.na(fit.rda.sha256))
                                            "unavailable" else fit.rda.sha256))

  ## Meteo provenance. Each day is read wholly from one tree; record the
  ## tree paths and the run-length per-day attribution. meteo_source_directory
  ## is kept (pointing at the tree that supplied the most days) for any
  ## downstream reader that still expects the old single-path attribute.
  for (src in names(era5dirs))
    ncatt_put(ncf, 0, sprintf("meteo_source_%s", src),
              attval = unname(era5dirs[[src]]), prec = "text")
  ncatt_put(ncf, 0, "meteo_source_by_day",
            attval = encode.day.runs(available.days, day.source), prec = "text")
  .src.counts <- table(day.source[available.days])
  .dominant   <- names(.src.counts)[which.max(.src.counts)]
  ncatt_put(ncf, 0, "meteo_source_directory",
            attval = unname(era5dirs[[.dominant]]), prec = "text")
  ## "yes" if any day came from a non-primary (fallback) tree.
  ncatt_put(ncf, 0, "meteo_fallback_used",
            attval = if (any(day.source[available.days] != names(era5dirs)[1]))
                       "yes" else "no",
            prec = "text")
  if (partial.month) {
    ncatt_put(ncf, 0, "status", attval = "provisional", prec = "text")
    ncatt_put(ncf, 0, "meteo_partial",
              attval = sprintf("only %d/%d days have ERA5 meteo; missing day(s) %s excluded from output",
                               dpm, dpm.full, paste(missing.days, collapse=",")),
              prec = "text")
  }

  ncvar_put(ncf, vars$dd,    vals = decimal.date)
  ncvar_put(ncf, vars$gpp,   vals = gpp)
  ncvar_put(ncf, vars$resp,  vals = resp)
  ncvar_put(ncf, vars$nee,   vals = nee)
  ncvar_put(ncf, vars$qgpp,  vals = qgpp)
  ncvar_put(ncf, vars$qresp, vals = qresp)
  ncvar_put(ncf, vars$ssr,   vals = mets$ssrd)
  ncvar_put(ncf, vars$t2m,   vals = mets$t2m)
  ncvar_put(ncf, vars$stl1,  vals = mets$stl1)
  ncvar_put(ncf, vars$swvl1, vals = mets$swvl1)
  if (write.flux.sd)                                    # per-month sd broadcast to all hour slots
    ncvar_put(ncf, vars$nee_sd, vals = array(nee.sd.mn, dim = dim(nee)))
  nc_close(ncf)
  n.written <- n.written + 1L
  cat("\n")
}

## Worker finished cleanly: record `ok` and disarm the error handler.
options(error = NULL)
manifest.record("diurnalize-ERA5.r", "ok",
                elapsed = as.integer(difftime(Sys.time(), diurn.t0,
                                              units = "secs")),
                detail  = sprintf("year=%d months_written=%d", yr, n.written))
