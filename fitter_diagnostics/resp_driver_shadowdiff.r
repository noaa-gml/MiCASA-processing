#!/usr/bin/env Rscript
#SBATCH --account co2
#SBATCH --time 0:30:00
#SBATCH --ntasks 1
#SBATCH --mem 24g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
## Shadow-diff: air-temp vs soil-temp respiration driver (prototype #1,
## docs/DIURNALIZATION_ALTERNATIVES.md). Quantifies the rectifier-relevant
## change: (1) mass conservation (monthly means unchanged), (2) respiration &
## NEE diurnal AMPLITUDE, (3) diurnal PHASE (hour of max). Reads the two
## single-month shadow outputs and writes a text report + a global-land mean
## diurnal-cycle figure. Pure ncdf4 + base R.
suppressMessages(library(ncdf4))

fair  <- "ERA5_resp_airtemp/fluxes_202007.nc"
fsoil <- "ERA5_resp_soiltemp/fluxes_202007.nc"
stopifnot(file.exists(fair), file.exists(fsoil))

rd <- function(f, v) { nc <- nc_open(f); on.exit(nc_close(nc)); ncvar_get(nc, v) }
nc <- nc_open(fair)
lon <- ncvar_get(nc, "longitude"); lat <- ncvar_get(nc, "latitude")
tsec <- ncvar_get(nc, "time")
tunit <- ncatt_get(nc, "time", "units")$value
nc_close(nc)
## epoch from "<unit> since YYYY-MM-DD HH:MM:SS UTC"
ep <- as.POSIXct(sub(".*since ", "", tunit), tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
uw <- sub(" since.*", "", tunit)
scale <- switch(uw, seconds=1, second=1, sec=1, hours=3600, hour=3600, days=86400, day=86400,
                stop(sprintf("unknown time unit '%s'", uw)))
tt <- ep + tsec*scale
hod <- as.integer(format(tt, "%H", tz = "UTC"))      # hour-of-day (UTC), mid-hour
nt  <- length(tsec)
cat(sprintf("file: %d timesteps, hours %d..%d, %d unique\n", nt, min(hod), max(hod), length(unique(hod))))

## area weight (cos lat), broadcast to [lon,lat]
wlat <- cos(lat * pi/180); wlat[wlat < 0] <- 0
W <- matrix(rep(wlat, each = length(lon)), length(lon), length(lat))

## mean diurnal cycle per cell: average each variable over timesteps sharing
## an hour-of-day. Returns [lon,lat,24].
diurnal.cycle <- function(x3) {
  d <- array(0, c(length(lon), length(lat), 24))
  for (h in 0:23) {
    sl <- which(hod == h)
    d[ , , h + 1] <- apply(x3[ , , sl, drop = FALSE], c(1, 2), mean)
  }
  d
}

report <- c(); CYC <- list()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); report <<- c(report, s) }

say("=== Respiration-driver shadow-diff: July 2020 (air vs soil temp) ===")
for (vv in c("resp", "NEE", "GPP")) {
  xa <- rd(fair,  vv); xs <- rd(fsoil, vv)
  ## (1) mass conservation: monthly mean per cell
  ma <- apply(xa, c(1,2), mean); ms <- apply(xs, c(1,2), mean)
  mass.max <- max(abs(ma - ms))
  mass.rel <- mass.max / max(abs(ma), 1e-30)
  say("")
  say("[%s] monthly-mean |air-soil| max = %.3e (rel %.2e)  <- expect ~0 (pure redistribution)",
      vv, mass.max, mass.rel)
  if (vv == "GPP") { say("  (GPP driver unchanged -> identical by construction)"); next }

  ## land mask: cells with appreciable |monthly mean| of this flux
  thr  <- 1e-9
  land <- is.finite(ma) & abs(ma) > thr
  ## (2)+(3) diurnal amplitude & phase, per cell
  da <- diurnal.cycle(xa); ds <- diurnal.cycle(xs)
  amp.a <- apply(da, c(1,2), function(z) max(z) - min(z))
  amp.s <- apply(ds, c(1,2), function(z) max(z) - min(z))
  ph.a  <- apply(da, c(1,2), which.max) - 1L     # hour of max (0..23)
  ph.s  <- apply(ds, c(1,2), which.max) - 1L
  ## amplitude ratio soil/air and phase shift, on land
  ratio <- amp.s[land] / pmax(amp.a[land], 1e-30)
  ## circular phase difference in hours, wrapped to (-12,12]
  dph <- ((ph.s[land] - ph.a[land] + 12) %% 24) - 12
  wl  <- W[land]
  wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o])/sum(w); x[o][which.max(cw >= p)] }
  say("  diurnal AMPLITUDE ratio soil/air (area-wtd):  median %.3f  [p25 %.3f, p75 %.3f]",
      wq(ratio, wl, .5), wq(ratio, wl, .25), wq(ratio, wl, .75))
  say("  diurnal PHASE shift soil-air (hours, area-wtd): median %+.2f  [p25 %+.2f, p75 %+.2f]",
      wq(dph, wl, .5), wq(dph, wl, .25), wq(dph, wl, .75))

  ## latitude-band breakdown
  bands <- list("boreal 50-70N"=c(50,70), "NH temp 25-50N"=c(25,50),
                "tropics 25S-25N"=c(-25,25), "SH temp 25-50S"=c(-50,-25))
  for (bn in names(bands)) {
    br <- bands[[bn]]; latsel <- lat >= br[1] & lat <= br[2]
    bm <- land & matrix(rep(latsel, each=length(lon)), length(lon), length(lat))
    if (sum(bm) < 5) next
    r <- amp.s[bm]/pmax(amp.a[bm],1e-30); dp <- ((ph.s[bm]-ph.a[bm]+12)%%24)-12; w2 <- W[bm]
    say("    %-16s amp ratio %.3f   phase %+.2f h   (n=%d)", bn, wq(r,w2,.5), wq(dp,w2,.5), sum(bm))
  }

  ## global-land mean diurnal cycle (area-weighted over land)
  gca <- sapply(1:24, function(h) sum(da[,,h][land]*wl)/sum(wl))
  gcs <- sapply(1:24, function(h) sum(ds[,,h][land]*wl)/sum(wl))
  say("  global-land mean diurnal range: air %.3e  soil %.3e  (soil/air %.3f)",
      max(gca)-min(gca), max(gcs)-min(gcs), (max(gcs)-min(gcs))/(max(gca)-min(gca)))
  say("  global-land peak hour (UTC):    air %d  soil %d", which.max(gca)-1L, which.max(gcs)-1L)
  CYC[[vv]] <- data.frame(hour=0:23, air=gca, soil=gcs)  # stash for CSV/figure
}

writeLines(report, "fitter_diagnostics/resp_driver_shadowdiff.txt")

## dump global-land mean diurnal cycles (resp + NEE, air vs soil) for plotting
cyc <- data.frame(hour=0:23,
  resp_air=CYC$resp$air, resp_soil=CYC$resp$soil,
  nee_air =CYC$NEE$air,  nee_soil =CYC$NEE$soil)
write.csv(cyc, "fitter_diagnostics/resp_driver_diurnal.csv", row.names=FALSE)
cat("\nWrote resp_driver_shadowdiff.txt and resp_driver_diurnal.csv\n")
