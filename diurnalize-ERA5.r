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
## Climatology fallback: years in $MICASA_CLIM_YEARS (space-separated env var,
## default "2000 $MICASA_YEAR") use Rh/NPP climatology instead of monthly files.

script.name <- "diurnalize-ERA5.r"

work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
cfg <- micasa.config()

product.name <- sprintf("MiCASA_%s_flux_x360_y180_monthly", cfg$version)

clim.yrs <- as.integer(strsplit(
  Sys.getenv("MICASA_CLIM_YEARS", sprintf("2000 %d", cfg$year)), "\\s+")[[1]])

## Hourly 1° ERA5 from the TM5 meteo tree.
era5dir <- sprintf("%s/METEO/tm5-nc/ec/ea_0005/h06h18tr1/sfc/glb100x100",
                   Sys.getenv("CARBONTRACKER", ""))
era5template <- "YYYY/MM/VVV_YYYYMMDD_00p01.nc"
metstr <- "ERA5"

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

setwd(work.dir)
in.dir  <- cfg$monthly.1x1
out.dir <- cfg$era5.dir

ct.setup()

## Load PIQS coefficients to smooth month-month variability.
load(file.path(work.dir, "fit.piqs.rda"))
piqsfit.time <- epoch.seconds.to.POSIX(piqsfit.time)
piqsfit.lts  <- as.POSIXlt(piqsfit.time)

lon.dim <- ncdim_def("longitude", "degrees_east", vals = seq(-179.5, 179.5, 1))
lat.dim <- ncdim_def("latitude",  "degrees_north", vals = seq(-89.5, 89.5, 1))

epoch <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC")
timeunits <- "days"

dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

for (mon in cfg$month.start:cfg$month.end) {

  cat(sprintf("%d/%02d\n", yr, mon))
  monstr <- sprintf("%d%02d", yr, mon)
  current.time <- ISOdatetime(yr, mon, 1, 0, 0, 0, tz = "UTC")
  ncname.out <- sprintf("%s/fluxes_%s.nc", out.dir, monstr)

  if (yr %in% clim.yrs) {
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
  varnms <- c("t2m", "ssrd", "stl1", "swvl1")
  mets <- list()
  dpm <- days.in.month(yr)[mon]
  times <- rep(NA, dpm * 24)
  for (day in 1:dpm) {
    for (varnm in varnms) {
      e5nm <- gsub("YYYY", sprintf("%d",   yr),    era5template)
      e5nm <- gsub("MM",   sprintf("%02d", mon),   e5nm)
      e5nm <- gsub("DD",   sprintf("%02d", day),   e5nm)
      e5nm <- gsub("VVV",  varnm,                  e5nm)
      ncname.in <- sprintf("%s/%s", era5dir, e5nm)
      foo <- load.ncdf(ncname.in)
      if (is.null(mets[[varnm]])) {
        mets[[varnm]] <- array(NA, dim = c(360, 180, 24 * dpm))
      }
      k0 <- 1 + (day - 1) * 24
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
  } else {
    monseq <- which((piqsfit.lts$mon + 1) == mon)
    gpp.a  <- apply(piqsfit.gpp$a[, , monseq],  c(1, 2), mean)
    gpp.b  <- apply(piqsfit.gpp$b[, , monseq],  c(1, 2), mean)
    gpp.c  <- apply(piqsfit.gpp$c[, , monseq],  c(1, 2), mean)
    resp.a <- apply(piqsfit.resp$a[, , monseq], c(1, 2), mean)
    resp.b <- apply(piqsfit.resp$b[, , monseq], c(1, 2), mean)
    resp.c <- apply(piqsfit.resp$c[, , monseq], c(1, 2), mean)
  }
  for (islot in 1:nslots) {
    dt <- as.numeric(times[islot]) - as.numeric(times[1])
    qmod.gpp  <- (gpp.a  * dt^2 + gpp.b  * dt + gpp.c)  / 12
    qmod.resp <- (resp.a * dt^2 + resp.b * dt + resp.c) / 12

    gpp[ , , islot] <- mets$ssrd[, , islot] * gpp.mn  / ssr.mn
    resp[, , islot] <- q10[      , , islot] * rtot.mn / q10.mn

    gpp[ , , islot] <- gpp[ , , islot] - gpp.mn  + qmod.gpp
    resp[, , islot] <- resp[, , islot] - rtot.mn + qmod.resp
    nee[ , , islot] <- gpp[, , islot] + resp[, , islot]
    qgpp[ , , islot] <- qmod.gpp
    qresp[, , islot] <- qmod.resp
  }

  lx <- which(is.na(nee))
  if (length(lx) > 0) nee[lx] <- 0

  cat(sprintf("  %s\n", ncname.out))
  date.vals    <- as.numeric(difftime(times, epoch, units = timeunits))
  decimal.date <- POSIX.to.decimal(times)
  date.dim <- ncdim_def("time",
                        sprintf("%s since %s", timeunits,
                                format(epoch, format = "%Y-%m-%d %H:%M:%S UTC")),
                        vals = date.vals, unlim = TRUE)

  ncvar <- function(name, units, longname, dims = list(lon.dim, lat.dim, date.dim)) {
    ncvar_def(name = name, units = units, dim = dims,
              missval = -1e34, compression = 9, longname = longname)
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

  if (file.exists(ncname.out)) {
    cat(sprintf("Removing existing output file \"%s\"\n", ncname.out))
    file.remove(ncname.out)
  }
  ncf <- nc_create(ncname.out, vars = vars)
  ncatt_put(ncf, 0, "history",
            attval = sprintf("Created on %s\nby script '%s'",
                             format(Sys.time(), "%a %b %d %Y %H:%M:%S %Z"),
                             script.name),
            prec = "text")
  ncatt_put(ncf, 0, "meteo_source_directory", attval = era5dir, prec = "text")

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
  nc_close(ncf)
  cat("\n")
}
