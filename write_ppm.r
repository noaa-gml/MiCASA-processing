#!/usr/bin/env Rscript
## PPM (Colella & Woodward 1984) limited piecewise-parabolic integral-preserving
## fitter driver. Storage matches write_pchip.r (a,b,c per piece). Output path
## $MICASA_FIT_OUT (default fit.piqs.rda). See lib/ppm_fit.r + PROPOSALS (17).
##
## Vectorized: the whole 360x180 grid is fit in one ppm.fit.grid() call rather
## than a per-cell loop (~6x faster; grid==cell to FP, tests/test_ppm_fit.r).
ct.setup()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "config.r"))
cfg <- micasa.config()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "ppm_fit.r"))
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
  f <- ppm.fit.grid(x.time, matrix(field, N, nmon))
  f$a[skip, ] <- 0; f$b[skip, ] <- 0; f$c[skip, ] <- 0   # match loop's all-zero skip
  list(a = array(f$a, c(360,180,nmon)),
       b = array(f$b, c(360,180,nmon)),
       c = array(f$c, c(360,180,nmon)))
}
piqsfit.gpp  <- fitone(gpp)
piqsfit.resp <- fitone(rtot)

piqsfit.time <- x.time[1:(length(x.time)-1)]
piqsfit.meta <- list(fitter="ppm", pad.left=0L, pad.right=0L, fit.range=range(x.time),
  saved.range=range(piqsfit.time), written.at=format(Sys.time(),tz="UTC",usetz=TRUE),
  notes="PPM (Colella & Woodward 1984) limited piecewise-parabolic; integral-preserving; no overshoot; mostly C0")
out <- Sys.getenv("MICASA_FIT_OUT","fit.piqs.rda")
save(file=out, piqsfit.gpp, piqsfit.resp, piqsfit.time, piqsfit.meta)
cat(sprintf("Wrote %s (PPM limited-parabolic, vectorized)\n", out))
