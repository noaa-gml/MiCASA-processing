#!/usr/bin/env Rscript
##
## Selectable alternative to write_pchip.r / write_piqs.r using a
## minmod-limited integral-preserving piecewise-LINEAR reconstruction
## (finite-volume / MUSCL). Storage layout matches the others (a, b, c per
## piece per cell, here with a == 0), so diurnalize-ERA5.r consumes the
## output transparently via MICASA_FIT_RDA.
##
## Output path: $MICASA_FIT_OUT (default fit.piqs.rda). Set it to e.g.
## fit.linmm.rda to A/B against the production PCHIP fit without clobbering.
##
## Vectorized: the whole 360x180 grid is fit in one linmm.fit.grid() call
## rather than a per-cell loop (grid==cell bit-for-bit, tests/test_linmm_fit.r).
##
## See lib/linmm_fit.r for the core and docs/PROPOSALS.md (17) for rationale:
## integral-preserving + provably no overshoot, at the cost of a flux
## discontinuity at each month boundary.
ct.setup()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "config.r"))
cfg <- micasa.config()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "linmm_fit.r"))
din <- load.ncdf(micasa.out.monthly.cat(cfg))
gpp  <- -2 * din$NPP
rtot <- din$Rh + din$NPP
nmon <- length(din$time)
plt.start <- as.POSIXlt(din$time[1]); y0<-plt.start$year+1900; m0<-plt.start$mon+1
x.time <- as.numeric(seq(ISOdatetime(y0,m0,1,0,0,0,tz="UTC"), by="1 month", length.out=nmon+1))

N <- 360 * 180
COEF_ZERO_THRESHOLD <- 1e-15
NPPm <- matrix(din$NPP, N, nmon); Rhm <- matrix(din$Rh, N, nmon)
skip <- (apply(abs(NPPm), 1, max) < COEF_ZERO_THRESHOLD) &
        (apply(abs(Rhm),  1, max) < COEF_ZERO_THRESHOLD)

fitone <- function(field) {
  f <- linmm.fit.grid(x.time, matrix(field, N, nmon))
  f$a[skip, ] <- 0; f$b[skip, ] <- 0; f$c[skip, ] <- 0
  list(a = array(f$a, c(360,180,nmon)),
       b = array(f$b, c(360,180,nmon)),
       c = array(f$c, c(360,180,nmon)))
}
piqsfit.gpp  <- fitone(gpp)
piqsfit.resp <- fitone(rtot)

piqsfit.time <- x.time[1:(length(x.time)-1)]
piqsfit.meta <- list(fitter      = "linmm",
                     pad.left    = 0L,
                     pad.right   = 0L,
                     fit.range   = range(x.time),
                     saved.range = range(piqsfit.time),
                     written.at  = format(Sys.time(), tz = "UTC", usetz = TRUE),
                     notes       = "minmod-limited integral-preserving piecewise-linear (MUSCL); flux a==0; no overshoot; discontinuous at month edges")
out <- Sys.getenv("MICASA_FIT_OUT", "fit.piqs.rda")
save(file = out, piqsfit.gpp, piqsfit.resp, piqsfit.time, piqsfit.meta)
cat(sprintf("Wrote %s (minmod-linear MUSCL, vectorized)\n", out))
