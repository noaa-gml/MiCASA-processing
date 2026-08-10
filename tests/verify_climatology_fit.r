#!/usr/bin/env Rscript
## verify_climatology_fit.r -- assertions about the climatology prior's
## sub-monthly coefficient fit. Invoked by tests/verify_climatology_prior.py;
## runnable standalone.
##
## Usage:
##   Rscript tests/verify_climatology_fit.r FIT.rda SERIES.nc YEAR0 YEAR1
##
## The fit is the component that actually sets the delivered monthly mean --
## diurnalize's hourly flux has monthly mean mean(qmod), taken entirely from
## here -- so it is the one place a silent reinstatement of the real
## interannual signal could hide behind a complete set of output files.
##
##   WINDOW  every delivered month has its own coefficients, so diurnalize
##           never falls back to the coefficient-climatology branch.
##   MEAN    the analytic mean of the fitted quadratic over each month
##           reproduces that month's climatological mean:
##           mean = a*h^2/3 + b*h/2 + c.
##   FLAT    the same calendar month carries identical coefficients in every
##           delivered year, i.e. no interannual structure survives.
##
## FLAT is measured RELATIVE to the coefficient scale. Bit-exact equality is
## not achievable: cumsum() accumulates from a different absolute offset in
## each year, so the last bits differ by pure round-off. A climatological fit
## sits near 1e-13 relative; a fit carrying real structure sits near 1. The
## separation is reported so a marginal pass is visible rather than being
## reported as a bare PASS.
##
## Emits CHECK|<id>|<PASS|FAIL>|<detail> lines and exits non-zero on failure.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("usage: verify_climatology_fit.r FIT.rda SERIES.nc YEAR0 YEAR1\n")
  quit(status = 2)
}
fit.path <- args[1]
mon.path <- args[2]
y0 <- as.integer(args[3])
y1 <- as.integer(args[4])

FLAT.TOL <- 1e-9
MEAN.TOL <- 1e-6

ct.setup()

fail <- 0L
emit <- function(id, ok, detail) {
  cat(sprintf("CHECK|%s|%s|%s\n", id, if (ok) "PASS" else "FAIL", detail))
  if (!ok) fail <<- fail + 1L
}

if (!file.exists(fit.path)) { emit("fit.present", FALSE,
    sprintf("fit not found: %s", fit.path)); quit(status = 1) }
if (!file.exists(mon.path)) { emit("series.present", FALSE,
    sprintf("series not found: %s", mon.path)); quit(status = 1) }

load(fit.path)
ft <- as.POSIXct(piqsfit.time, origin = "1970-01-01", tz = "UTC")
cat(sprintf("INFO|fit months=%d span=%s..%s fitter=%s\n",
            length(ft), format(min(ft), "%Y-%m"), format(max(ft), "%Y-%m"),
            if (exists("piqsfit.meta")) piqsfit.meta$fitter else "unknown"))

## ---- WINDOW ---------------------------------------------------------------
want <- seq(ISOdatetime(y0, 1, 1, 0, 0, 0, tz = "UTC"),
            ISOdatetime(y1, 12, 1, 0, 0, 0, tz = "UTC"), by = "1 month")
nfound <- sum(want %in% ft)
emit("fit.window", nfound == length(want),
     sprintf("%d/%d delivered months have their own coefficients",
             nfound, length(want)))

## ---- MEAN-PRESERVATION ----------------------------------------------------
din  <- load.ncdf(mon.path)
gpp  <- -2 * din$NPP                 # gC m-2 s-1, as the fitters build it
rtot <- din$Rh + din$NPP
dt.mon <- as.POSIXct(din$time, tz = "UTC")

mstart <- as.POSIXct(format(ft, "%Y-%m-01"), tz = "UTC")
nxt <- seq(mstart[length(mstart)], by = "1 month", length.out = 2)[2]
h.all <- as.numeric(difftime(c(mstart[-1], nxt), mstart, units = "secs"))

land <- which(apply(abs(gpp), c(1, 2), max) > 1e-15, arr.ind = TRUE)
set.seed(1)
samp <- land[sample(nrow(land), min(400L, nrow(land))), , drop = FALSE]

worst.gpp <- 0; worst.resp <- 0
ks <- which(ft %in% want)
for (k in ks) {
  h  <- h.all[k]
  im <- which(format(dt.mon, "%Y-%m") == format(ft[k], "%Y-%m"))[1]
  if (is.na(im)) next
  for (r in seq_len(nrow(samp))) {
    i <- samp[r, 1]; j <- samp[r, 2]
    mg <- piqsfit.gpp$a[i, j, k]  * h^2 / 3 + piqsfit.gpp$b[i, j, k]  * h / 2 +
          piqsfit.gpp$c[i, j, k]
    mr <- piqsfit.resp$a[i, j, k] * h^2 / 3 + piqsfit.resp$b[i, j, k] * h / 2 +
          piqsfit.resp$c[i, j, k]
    worst.gpp  <- max(worst.gpp,
                      abs(mg - gpp[i, j, im])  / max(abs(gpp[i, j, im]),  1e-12))
    worst.resp <- max(worst.resp,
                      abs(mr - rtot[i, j, im]) / max(abs(rtot[i, j, im]), 1e-12))
  }
}
emit("fit.mean_preserving", worst.gpp < MEAN.TOL && worst.resp < MEAN.TOL,
     sprintf("worst relative deviation gpp=%.3e resp=%.3e over %d cells x %d months (tol %.0e)",
             worst.gpp, worst.resp, nrow(samp), length(ks), MEAN.TOL))

## ---- FLATNESS -------------------------------------------------------------
maxdev <- 0; worst.m <- NA
for (m in 1:12) {
  kk <- which(format(ft, "%m") == sprintf("%02d", m) &
              as.integer(format(ft, "%Y")) >= y0 &
              as.integer(format(ft, "%Y")) <= y1)
  if (length(kk) < 2) next
  hs <- h.all[kk]
  kk <- kk[hs == hs[1]]           # equal-length instances only (Feb: non-leap)
  if (length(kk) < 2) next
  for (k2 in kk[-1]) {
    for (nm in c("a", "b", "c")) {
      ref <- piqsfit.gpp[[nm]][, , kk[1]]
      cur <- piqsfit.gpp[[nm]][, , k2]
      sc  <- max(abs(ref))
      if (sc > 0) {
        d <- max(abs(cur - ref)) / sc
        if (d > maxdev) { maxdev <- d; worst.m <- m }
      }
    }
  }
}
emit("fit.flat", maxdev < FLAT.TOL,
     sprintf("worst relative deviation %.3e (month %s), tol %.0e, margin %.3g x",
             maxdev, worst.m, FLAT.TOL,
             if (maxdev > 0) FLAT.TOL / maxdev else Inf))

cat(sprintf("SUMMARY|%d\n", fail))
quit(status = if (fail == 0L) 0 else 1)
